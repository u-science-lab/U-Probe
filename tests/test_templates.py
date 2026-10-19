import copy
import random
import re

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from uprobe.core.api import UProbeAPI
from uprobe.core.attributes import add_attributes
from uprobe.core.cli import cli, _extract_probe_targets, _validate_and_normalize_protocol
from uprobe.core.gene_expression import validate_gene_expression
from uprobe.core.hybridization import FIELDS as HYBRIDIZATION_FIELDS, TM_DEFAULTS, effective_conditions
from uprobe.core.layout import validate_layout_references
from uprobe.core.offtarget import mapped_genes_settings
from uprobe.core.process import post_process
from uprobe.core.sampling import normalize_sampling
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import render_text, validate_files
from uprobe.templates import list_templates, load_template, read_template, template_names

NAMES = template_names()
REQUIRED = {
    'smfish_stellaris', 'merfish', 'seqfish_plus', 'starmap', 'padlock_iss',
    'mipseq', 'hcr_v3', 'oligopaint', 'double_hyb_rca',
}
META_KEYS = ['name', 'technology', 'description', 'citations', 'parameter_sources',
             'hybridization', 'replace']
# Attribute types that need no genome, aligner index or external tool.
OFFLINE_TYPES = {'gc_content', 'annealing_temperature', 'fold_score', 'self_match', 'length'}
SEQ = ''.join(random.Random(0).choice('ACGT') for _ in range(400))


def test_library_has_every_technology():
    assert REQUIRED <= set(NAMES)


@pytest.mark.parametrize('name', NAMES)
def test_meta_documents_sources(name):
    template = load_template(name)
    meta = template['meta']
    assert meta['name'] == name
    assert all(meta.get(key) for key in META_KEYS)
    # The primary citation is a paper; protocols and code may be cited by URL.
    assert 'doi' in meta['citations'][0]
    for citation in meta['citations']:
        assert citation['short'] and citation['title']
        assert citation.get('doi') or citation.get('url')
        if 'doi' in citation:
            assert re.fullmatch(r'10\.\d{4,9}/\S+', citation['doi'])
    # The source protocol's conditions apply to the workflow; meta records where they come from.
    assert {'temperature_c', 'na_mM', 'formamide_pct'} <= set(template['hybridization'])
    assert meta['hybridization']['source']
    for row in meta['parameter_sources']:
        assert row['parameter'] and row['source']


@pytest.mark.parametrize('name', NAMES)
def test_template_passes_uprobe_validate(name, tmp_path):
    path = tmp_path / f'{name}.yaml'
    path.write_text(read_template(name), encoding='utf-8')
    findings = validate_files(path)
    assert not [f for f in findings if f.severity == 'error'], render_text(findings)


@pytest.mark.parametrize('name', [name for name in NAMES if load_template(name)['meta'].get('tm_window')])
def test_tm_window_applied_only_when_verified(name):
    protocol = load_template(name)
    window = protocol['meta']['tm_window']
    verified = isinstance(window['conditions'], dict) and 'TODO' not in str(window)
    assert window['applied'] is verified
    attributes = protocol['attributes']
    if not verified:
        assert not any('hybridization' in attr for attr in attributes.values())
        return
    name = window['attribute']
    design = {key: value for key, value in window['conditions'].items() if key in HYBRIDIZATION_FIELDS}
    assert attributes[name]['hybridization'] == design
    assert effective_conditions(protocol)[name] == {**TM_DEFAULTS, **design}
    low, high = window['value']
    assert protocol['post_process']['filters'][name]['condition'] == f'{name} >= {low} & {name} <= {high}'


def test_double_hyb_rca_tm_filters_keep_legacy_conditions():
    # The 35-45 C filters were set under the legacy conditions, not the lifted experimental ones.
    conditions = effective_conditions(load_template('double_hyb_rca'))
    assert conditions['tm1'] == conditions['tm2'] == conditions['tm3'] == TM_DEFAULTS
    assert conditions['circle_fold_score'] == {'temperature_c': 37}


def offline_api(protocol, seq=SEQ, encoding=None):
    """Swap the template targets for one direct sequence; no genome needed."""
    protocol = copy.deepcopy(protocol)
    api = object.__new__(UProbeAPI)
    api.protocol = protocol
    api.genome = None
    protocol['targets'] = [{'TEST': seq}]
    protocol['encoding'] = {'TEST': encoding or next(iter(protocol['encoding'].values()))}
    return api


def design(name, seq, encoding):
    api = offline_api(_validate_and_normalize_protocol(load_template(name)), seq, encoding)
    targets = api.generate_target_seqs()
    return targets, api.construct_probes(targets)


@pytest.mark.parametrize('name', NAMES)
def test_template_validates_and_builds_probes(name, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    protocol = _validate_and_normalize_protocol(load_template(name))
    validate_layout_references(protocol)
    api = offline_api(protocol)
    targets = api.generate_target_seqs()
    assert not targets.empty
    probes = api.construct_probes(targets)
    expected = {target for target, _ in _extract_probe_targets(protocol['probes'])}
    assert set(probes.columns) == expected
    assert (probes != '').all().all() and len(probes) == len(targets)

    combined = pd.concat([targets.reset_index(drop=True), probes], axis=1)
    attributes = protocol['attributes']
    for attr in attributes.values():
        target = attr['target']
        assert target in combined or target.replace(':', '.') in combined, target
    pp = protocol['post_process']
    for spec in pp.get('filters', {}).values():
        if 'condition' in spec:
            columns = set(re.findall(r'([A-Za-z_]\w*)\s*[<>=]', spec['condition']))
            assert columns <= set(attributes), spec
        else:
            assert spec['target'] in combined, spec
    sorts = pp.get('sorts', {})
    assert set(sorts.get('is_ascending', []) + sorts.get('is_descending', [])) <= set(attributes)
    assert set(protocol['summary'].get('attributes', [])) <= set(attributes)

    offline = {k: v for k, v in attributes.items() if v['type'] in OFFLINE_TYPES}
    result = add_attributes(combined, {**protocol, 'attributes': offline}, {})
    assert result[list(offline)].notna().all().all()
    # Post-processing runs on the offline attributes (genome-only filters skip missing columns).
    offline_pp = {k: v for k, v in pp.items() if k not in ('avoid_otp', 'off_target_pool')}
    processed = post_process(result, {**protocol, 'post_process': offline_pp})
    if 'equal_space' in pp:
        assert len(processed) <= pp['equal_space']['number_desired']


@pytest.mark.parametrize('name', NAMES)
def test_off_target_search_sees_the_template_sequences(name):
    # Bowtie2's local mode aligns nothing shorter than 22 nt and misses many
    # sites of 25-30-mers, so windows up to 30 nt are searched with Bowtie 1.
    protocol = load_template(name)
    shortest = normalize_sampling(protocol['extracts']['target_region']['length'])[0]
    for attr in protocol['attributes'].values():
        if attr['type'] != 'mapped_genes':
            continue
        # target_region is the RNA sequence; the probes bind where it occurs.
        assert attr['target'] == 'target_region' and attr['strand'] == 'sense'
        if shortest <= 30:
            assert attr['aligner'] == 'bowtie' and attr['max_mismatches'] == 2
        else:
            assert attr['aligner'] == 'bowtie2'


def test_smfish_allows_no_off_target_gene_and_documents_the_pool_option(tmp_path):
    # No off-target gene per probe by default; the pool (one off-target gene
    # per probe, 2 probes per gene) is a commented-out option.
    text = read_template('smfish_stellaris')
    pp = load_template('smfish_stellaris')['post_process']
    assert pp['filters']['target_mapped_genes']['condition'] == 'target_mapped_genes <= 1'
    assert 'off_target_pool' not in pp and '# gene_expression:' in text
    option = '  # off_target_pool:\n  #   attribute: target_mapped_genes\n  #   max_probes_per_gene: 2\n'
    assert text.count(option) == 1
    enabled = text.replace(option, option.replace('# ', '')).replace(
        'condition: target_mapped_genes <= 1', 'condition: target_mapped_genes <= 2')
    path = tmp_path / 'smfish_pool.yaml'
    path.write_text(enabled, encoding='utf-8')
    assert yaml.safe_load(enabled)['post_process']['off_target_pool'] == {
        'attribute': 'target_mapped_genes', 'max_probes_per_gene': 2}
    findings = validate_files(path)
    assert not findings, render_text(findings)


def commented_block(text, first_line):
    """The YAML of a commented-out top-level block: `# <first_line>` and the commented lines after it."""
    lines = text[text.index(f'\n# {first_line}\n') + 1:].splitlines()
    end = next((i for i, line in enumerate(lines) if not line.startswith('#')), len(lines))
    return yaml.safe_load('\n'.join(line[2:] for line in lines[:end]))


def test_smfish_gene_expression_example_uses_the_recommended_settings():
    example = commented_block(read_template('smfish_stellaris'), 'gene_expression:')['gene_expression']
    settings = validate_gene_expression(example)
    assert (settings['min_value'], settings['min_fraction_of_target']) == (1, 0.05)


def test_smfish_tiles_the_reference_mrna_and_exempts_copies(tmp_path):
    # The whole mRNA of the MANE Select (else Ensembl canonical) transcript; copies of the
    # target are exempt by default, which needs no expression table.
    protocol = load_template('smfish_stellaris')
    region = protocol['extracts']['target_region']
    assert (region['source'], region['transcripts']) == ('mRNA', 'mane')
    attr = protocol['attributes']['target_mapped_genes']
    assert mapped_genes_settings(attr)['exempt_copies'] == {'kmer': 30, 'min_shared': 0.5}
    # The gene_expression example, uncommented, is valid as it stands.
    example = commented_block(read_template('smfish_stellaris'), 'gene_expression:')['gene_expression']
    example['table'] = str(tmp_path / 'expression.tsv')
    (tmp_path / 'expression.tsv').write_text('gene_id\tTPM\nENSG00000111640\t2000\n')
    path = tmp_path / 'smfish_expression.yaml'
    path.write_text(yaml.safe_dump({**protocol, 'gene_expression': example}), encoding='utf-8')
    findings = validate_files(path)
    assert not findings, render_text(findings)


RNA_TEMPLATES = [name for name in NAMES if load_template(name)['extracts']['target_region']['source'] != 'genome']


@pytest.mark.parametrize('name', RNA_TEMPLATES)
def test_rna_templates_choose_source_and_transcripts(name):
    # Every RNA template names what it extracts (not the ambiguous exon) and designs on
    # the reference transcript; meta.parameter_sources gives the evidence.
    protocol = load_template(name)
    region = protocol['extracts']['target_region']
    assert region['source'] in ('mRNA', 'CDS') and region['transcripts'] == 'mane'
    rows = [row for row in protocol['meta']['parameter_sources'] if 'transcripts' in row['parameter']]
    assert len(rows) == 1 and rows[0]['source']


@pytest.mark.parametrize('name, count', [('smfish_stellaris', 48), ('merfish', 92), ('seqfish_plus', 24),
                                         ('hcr_v3', 20), ('starmap', 4), ('padlock_iss', 5)])
def test_published_probe_counts_per_gene(name, count):
    protocol = _validate_and_normalize_protocol(load_template(name))
    assert protocol['post_process']['equal_space'] == {'number_desired': count}
    sources = {row['parameter']: row for row in protocol['meta']['parameter_sources']}
    assert sources.get('post_process.equal_space.number_desired', {'value': count})['value'] == count


# Published target windows and oligos: the templates must rebuild them exactly.
def test_padlock_iss_reproduces_hybiss_pvalb():
    # Gyllborg et al. 2020, Supplementary Table 1 (mouse Pvalb, BpID_0015).
    targets, probes = design('padlock_iss', 'CCCATCTCAGTCCTTCTCGCGGCCCTCCTG', {'ID': 'TAGACGTAGTGAGCATGACT'})
    assert probes['padlock'].tolist() == [
        'CTCGCGGCCCTCCTGTAGACGTAGTGAGCATGACTTCTACGAGTTTGCAGTCACGCCCATCTCAGTCCTT']
    # The published target passes the applied Tm window under its design conditions.
    protocol = load_template('padlock_iss')
    combined = pd.concat([targets.reset_index(drop=True), probes], axis=1)
    tm = add_attributes(combined, {**protocol, 'attributes': {'target_tm': protocol['attributes']['target_tm']}}, {})
    assert 65 <= tm['target_tm'][0] <= 75


def test_mipseq_reproduces_published_pair():
    # Wu et al. 2024 Gad2-2 probes on mouse Gad2 mRNA (NM_008078.2).
    _, probes = design('mipseq', 'AGCGTGTGATGGAGAAAGGCCCACTCTGGCATTTCTACAA', {'B1': 'GCA', 'B2': 'GCA'})
    assert probes['padlock'].tolist() == ['TCCATCACACGCTCCGCATGCGTCTATTTGTGCATAGTGGAGCCTGAGTGGGCCTTTC']
    assert probes['rca_primer'].tolist() == ['TTGTAGAAATGCCTCAGGCTCCACTA']


def test_starmap_reproduces_published_pair():
    # Zeng et al. 2023 ALDOC_01/ALDOC_11 (64-gene set) on mouse Aldoc mRNA (NM_009657.4).
    _, probes = design('starmap', 'CCAGGATAAGGGCATTCTCGTAGGCATCAAGGTTGACAAGGGT', {'code': 'AAA'})
    pairs = set(zip(probes['padlock'], probes['primer']))
    assert ('AAAATAACGAGAATGCCCTTATCCTGGAATTATTACAAACATACACTAAAGATA',
            'ACCCTTGTCAACCTTGATGCTATTTTTATCTT') in pairs


def test_hcr_matches_insitu_probe_generator():
    # maker37cb.py reverse complements the mRNA (W) and prints
    # up + upspc + W[27:52] and W[0:25] + dnspc + dn for each 52-nt window.
    encoding = load_template('hcr_v3')['encoding']['GAPDH']
    targets, probes = design('hcr_v3', SEQ[:120], encoding)
    for window, odd, even in zip(targets['target_region'], probes['probe_odd'], probes['probe_even']):
        w = reverse_complement(window)
        assert odd == encoding['init_odd'] + encoding['spacer_odd'] + w[27:52]
        assert even == w[0:25] + encoding['spacer_even'] + encoding['init_even']


def test_double_hyb_rca_primer_binds_rna_at_5_prime_and_circle_at_3_prime():
    encoding = {'barcode1': 'A' * 17, 'barcode2': 'ACGTTGCAACGTTGCAA'}
    targets, probes = design('double_hyb_rca', SEQ[:40], encoding)
    window = targets['target_region'][0]
    assert probes['amp_probe'][0] == reverse_complement(window[-13:]) + 'N' + reverse_complement(encoding['barcode2'])
    assert probes['circle_probe'][0].startswith(reverse_complement(window[:13]))
    assert probes['circle_probe'][0].endswith(reverse_complement(window[13:26]))


def test_list_command():
    result = CliRunner().invoke(cli, ['templates', 'list'])
    assert result.exit_code == 0, result.output
    lines = result.output.splitlines()
    assert lines[0].split()[:3] == ['NAME', 'TECHNOLOGY', 'CITATION']
    for row in list_templates():
        assert any(line.startswith(row['name'] + ' ') and row['citation'] in line for line in lines)
        assert 'doi:10.' in row['citation']


def test_show_command(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, ['templates', 'show', 'hcr_v3'])
    assert result.exit_code == 0
    assert result.output == read_template('hcr_v3')

    out = tmp_path / 'protocol.yaml'
    result = runner.invoke(cli, ['templates', 'show', 'merfish', '-o', str(out)])
    assert result.exit_code == 0
    assert yaml.safe_load(out.read_text()) == load_template('merfish')
    result = runner.invoke(cli, ['templates', 'show', 'merfish', '-o', str(out)])
    assert result.exit_code != 0 and 'exists' in result.output
    assert runner.invoke(cli, ['templates', 'show', 'merfish', '-o', str(out), '--force']).exit_code == 0

    result = runner.invoke(cli, ['templates', 'show', 'no_such_template'])
    assert result.exit_code == 2 and 'hcr_v3' in result.output
