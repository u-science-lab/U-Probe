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
OFFLINE_TYPES = {'gc_content', 'annealing_temperature', 'fold_score', 'self_match', 'length', 'difference'}
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


def _parts_with_literal_n(parts, prefix=''):
    """Paths of probes and parts whose own template text holds a literal N."""
    for name, part in (parts or {}).items():
        if 'N' in re.sub(r'\{[^}]*\}', '', part.get('template') or ''):
            yield prefix + name
        yield from _parts_with_literal_n(part.get('parts'), f'{prefix}{name}.')


@pytest.mark.parametrize('name', NAMES)
def test_meta_replace_names_every_literal_n(name):
    # A synthesized N is a mix of all four bases: meta.replace must tell users to replace
    # each one (double_hyb_rca's 1-nt spacers were missing, so its probes kept an N).
    template = load_template(name)
    replace = ' '.join(template['meta']['replace'])
    for path in _parts_with_literal_n(template['probes']):
        assert path in replace or path.rsplit('.', 1)[-1] in replace, path


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


@pytest.mark.parametrize('name', NAMES)
def test_templates_rank_lower_scores_first(name, tmp_path):
    # fold_score is -MFE: higher = more structure. Fewer self-matches, k-mer
    # counts and mapped genes are better too.
    protocol = load_template(name)
    sorts = protocol['post_process'].get('sorts') or {}
    lower_is_better = {key for key, attr in protocol['attributes'].items()
                       if attr['type'] in ('fold_score', 'self_match', 'kmer_count', 'mapped_genes')}
    assert not lower_is_better & set(sorts.get('is_descending') or [])
    path = tmp_path / f'{name}.yaml'
    path.write_text(read_template(name), encoding='utf-8')
    assert not [f for f in validate_files(path) if 'sorted descending' in f.message]


def test_double_hyb_rca_tm_filters_keep_legacy_conditions():
    # The 35-45 C filters were set under the legacy conditions, not the lifted experimental ones;
    # without Mg2+ their 0.6 mM dNTP changed no Tm, so the template's 0 keeps the window.
    conditions = effective_conditions(load_template('double_hyb_rca'))
    assert conditions['tm1'] == conditions['tm2'] == conditions['tm3'] == {**TM_DEFAULTS, 'dntp_mM': 0}
    assert conditions['circle_fold_score'] == {'temperature_c': 37}


@pytest.mark.parametrize('name', NAMES)
def test_templates_compute_tm_without_dntps(name):
    # FISH buffers hold no dNTPs; unset, the legacy PCR default of 0.6 mM applied to every Tm
    # and, with padlock_iss's 10 mM Mg2+, lowered the free Mg2+ of its arm Tm.
    protocol = load_template(name)
    assert protocol['hybridization']['dntp_mM'] == 0
    conditions = effective_conditions(protocol)
    for attr, spec in protocol['attributes'].items():
        if spec['type'] == 'annealing_temperature':
            assert conditions[attr]['dntp_mM'] == 0, attr


@pytest.mark.parametrize('name', NAMES)
def test_tm_filters_filter_as_under_the_legacy_dntp(name):
    # dNTPs only bind Mg2+: every filtered Tm has no Mg2+ or a formula without dNTPs, so the
    # published windows (MERFISH, HybISS, Oligopaint, double_hyb_rca) filter as before.
    from uprobe.core.attributes._attributes import cal_temp, cal_temp_multi_padlock_design
    from uprobe.core.hybridization import multi_padlock_design_arguments, tm_arguments
    protocol = load_template(name)
    conditions = effective_conditions(protocol)
    filtered = set(protocol['post_process'].get('filters') or {})
    rng = random.Random(5)
    for attr, spec in protocol['attributes'].items():
        if spec['type'] != 'annealing_temperature' or attr not in filtered:
            continue
        legacy = {**conditions[attr], 'dntp_mM': 0.6}
        for length in (13, 15, 20, 25, 30, 41):
            seq = ''.join(rng.choice('ACGT') for _ in range(length))
            if spec.get('tm_method') == 'multi_padlock_design':
                now = cal_temp_multi_padlock_design(seq, **multi_padlock_design_arguments(conditions[attr]))
                before = cal_temp_multi_padlock_design(seq, **multi_padlock_design_arguments(legacy))
            else:
                now, before = cal_temp(seq, **tm_arguments(conditions[attr])), cal_temp(seq, **tm_arguments(legacy))
            assert now == pytest.approx(before), (attr, seq)


@pytest.mark.parametrize('name, floors, weak', [
    # The weakest probe each template kept on the e2e genes (GAPDH, TFRC, EEF2; mouse Gad2).
    ('smfish_stellaris', {'probe_tm': 37}, {'probe_tm': 'ATAATTATATTCTATTAAAA'}),
    ('hcr_v3', {'odd_tm': 37, 'even_tm': 37}, {}),
    ('mipseq', {'arm_5_tm': 28, 'arm_3_tm': 28, 'primer_tm': 28},
     {'arm_5_tm': 'TTAATTGCATATT', 'primer_tm': 'TGATTAAAATATC'}),
])
def test_tm_floors_from_the_published_probe_sets(name, floors, weak):
    # Without a GC or Tm filter these templates kept probes that cannot bind at the hybridization
    # temperature (TFRC 5% GC, Tm 33.7 C at 37 C; a MiP-seq arm of 23.6 C and primer of 20.4 C).
    # The floors and the share of published probes they reject are in meta.parameter_sources.
    from uprobe.core.attributes._attributes import cal_temp
    from uprobe.core.hybridization import tm_arguments
    protocol = load_template(name)
    filters = protocol['post_process']['filters']
    conditions = effective_conditions(protocol)
    for attr, floor in floors.items():
        assert filters[attr]['condition'] == f'{attr} >= {floor}'
        assert any(f'filters.{attr}' in row['parameter'] for row in protocol['meta']['parameter_sources'])
    for attr, seq in weak.items():
        assert cal_temp(seq, **tm_arguments(conditions[attr])) < floors[attr]


@pytest.mark.parametrize('name', ['smfish_stellaris', 'merfish', 'seqfish_plus'])
def test_many_probe_templates_exempt_copies_of_the_target(name):
    # A false spot needs many probes of the set on one RNA, so a copy of the target (a processed
    # pseudogene, a retrocopy) is not counted: without it mouse Gapdh kept 2 of 24 seqFISH+ probes,
    # with it 24 (Gm13394 and 43 other copies, NIH/3T3 expression).
    protocol = load_template(name)
    assert protocol['attributes']['target_mapped_genes']['exempt_copies'] == {'kmer': 30, 'min_shared': 0.5}
    rows = [row for row in protocol['meta']['parameter_sources'] if 'exempt_copies' in row['parameter']]
    assert rows and rows[0]['value'] == '{kmer: 30, min_shared: 0.5}'


def test_smfish_tm_uses_the_stellaris_concentration_per_oligo():
    # Stellaris hybridizes the 48-oligo set at 125 nM in total, about 2.6 nM per oligo.
    assert load_template('smfish_stellaris')['hybridization']['probe_nM'] == 2.6


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
        if attr['type'] == 'difference':  # of two attributes computed before it
            assert set(attr['of']) <= set(attributes), attr
            continue
        target = attr['target']
        assert target in combined or target.replace(':', '.') in combined, target
    pp = protocol['post_process']
    # A mapped_sites attribute adds <attribute>_num, the number of sites.
    columns_of_attributes = set(attributes) | {f'{name}_num' for name, attr in attributes.items()
                                               if attr['type'] == 'mapped_sites'}
    for spec in pp.get('filters', {}).values():
        if 'condition' in spec:
            columns = set(re.findall(r'([A-Za-z_]\w*)\s*[<>=]', spec['condition']))
            # Attributes, or a position column of the extracted windows (direct sequences have none,
            # and a condition on a missing column is skipped).
            assert columns <= columns_of_attributes | {'target_start', 'target_end', 'start', 'end'}, spec
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
    region = protocol['extracts']['target_region']
    for attr in protocol['attributes'].values():
        if attr['type'] != 'mapped_genes':
            continue
        # target_region, or a layout part of it, is the RNA sequence; the probes bind where it occurs.
        assert attr['strand'] == 'sense'
        if attr['target'] == 'target_region':
            shortest = normalize_sampling(region['length'])[0]
        else:
            part = attr['target'].removeprefix('target_parts.')
            assert attr['target'] == f'target_parts.{part}' and part in region['layout']['parts']
            shortest = normalize_sampling(region['layout']['parts'][part]['length'])[0]
        if shortest <= 30:
            assert attr['aligner'] == 'bowtie'
            # A whole window within 2 substitutions; a part may set its own limit (meta.parameter_sources).
            assert attr['max_mismatches'] == 2 or attr['target'] != 'target_region'
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


def test_smfish_meta_gives_raj_probe_lengths():
    # Raj 2008's probe list: the d2EGFP (GFP) set is 48 17-mers and the Actb set 48 20-mers;
    # only Supplementary Fig. 4 labels the GFP probes "20 bp".
    meta = load_template('smfish_stellaris')['meta']
    description = ' '.join(meta['description'].split())
    assert 'd2EGFP test set is 48 17-mers' in description and '20-mers in the GFP test' not in description
    row = next(r for r in meta['parameter_sources'] if r['parameter'] == 'extracts.target_region.length')
    assert 'GFP) set 48 17-mers' in ' '.join(row['source'].split()) and '20 bp' not in row['quote']


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
    # The published target passes the applied Tm window under the formula and conditions of the
    # HybISS design code (multi_padlock_design lib/screenseq.py gives 70.50 C; primer3 67.2 C).
    protocol = load_template('padlock_iss')
    combined = pd.concat([targets.reset_index(drop=True), probes], axis=1)
    tm = add_attributes(combined, {**protocol, 'attributes': {'target_tm': protocol['attributes']['target_tm']}}, {})
    assert tm['target_tm'][0] == pytest.approx(70.50057357343485)
    assert protocol['post_process']['filters']['target_tm']['condition'] == 'target_tm >= 65 & target_tm <= 75'
    # The junction search covers the 20 nt around the ligation junction (mRNA nt 15|16).
    junction = targets['target_parts.junction'][0]
    assert junction == 'CCCATCTCAGTCCTTCTCGCGGCCCTCCTG'[5:25]
    assert probes['padlock'][0].endswith(junction[:10]) and probes['padlock'][0].startswith(junction[10:])


def test_mipseq_reproduces_published_pair():
    # Wu et al. 2024, 217-gene set (sheet "Fig 3"): Fig3_463_C1qb on the window at nt 24-63 of
    # mouse C1qb (ENSMUST00000031132). Its primer is 25 nt, as 247 of the 432 published ones are.
    _, probes = design('mipseq', 'GTTCTTCCTGCCTCTAGGGACCCAGACTTCCGCTTTCTGA', {'B1': 'AGA', 'B2': 'TTT'})
    assert probes['padlock'].tolist() == ['AGGCAGGAAGAACCCAGATGCGTCTATTTGTTTTTAGTGGAGCCTTCTGGGTCCCTAG']
    assert probes['rca_primer'].tolist() == ['TCAGAAAGCGGAACAGGCTCCACTA']


def test_mipseq_keeps_two_sites_near_the_5_prime_end():
    # Wu 2023 Methods: 2 sites per gene, in an exon close to the 5' end (the 1,000-nt cap is U-Probe's).
    protocol = _validate_and_normalize_protocol(load_template('mipseq'))
    pp = protocol['post_process']
    assert pp['equal_space'] == {'number_desired': 2}
    assert pp['filters']['near_5p']['condition'] == 'target_start <= 1000'
    # A direct sequence has no target_start, so number the windows as an mRNA target's are numbered.
    api = offline_api(protocol, SEQ * 4)
    windows = api.generate_target_seqs()
    windows['target_start'], windows['target_end'] = windows['start'], windows['end']
    assert windows['target_start'].max() > 1000  # the filter has something to remove
    kept = post_process(windows, {'post_process': {k: pp[k] for k in ('filters', 'remove_overlap', 'equal_space')}})
    assert len(kept) == 2 and kept['target_start'].max() <= 1000


def test_starmap_reproduces_published_pair():
    # Zeng et al. 2023 ALDOC_01/ALDOC_11 (64-gene set) on mouse Aldoc mRNA (NM_009657.4);
    # its two sites are 1 nt apart, as every published pair's are.
    targets, probes = design('starmap', 'CCAGGATAAGGGCATTCTCGTAGGCATCAAGGTTGACAAGGGT', {'code': 'AAA'})
    pairs = set(zip(probes['padlock'], probes['primer']))
    assert ('AAAATAACGAGAATGCCCTTATCCTGGAATTATTACAAACATACACTAAAGATA',
            'ACCCTTGTCAACCTTGATGCTATTTTTATCTT') in pairs


def test_starmap_forbids_abutting_sites():
    # Wang 2018 allows a 0-nt gap, but none of the 4,178 located published pairs has one.
    protocol = load_template('starmap')
    parts = protocol['extracts']['target_region']['layout']['parts']
    assert parts['gap']['length'] == [1, 2]
    targets, _ = design('starmap', SEQ[:60], {'code': 'AAA'})
    assert set(targets['target_parts.gap.length']) == {1, 2}
    # The designed 60 °C of the halves is recorded, not applied (the paper's conditions are unknown).
    window = protocol['meta']['tm_window']
    assert (window['applied'], window['value']) == (False, 60)


def test_validate_checks_difference_attributes_and_prefer():
    from uprobe.core.validate import validate_protocol
    protocol = load_template('starmap')
    assert not [f for f in validate_protocol(protocol) if f.severity == 'error']
    protocol['attributes']['tm_balance']['of'] = ['padlock_tm', 'target_mapped_genes']  # defined after it
    protocol['post_process']['remove_overlap']['prefer'] = 'tm_balanse'
    errors = {f.path: f.message for f in validate_protocol(protocol) if f.severity == 'error'}
    assert "'target_mapped_genes' is not an attribute defined before this one" in errors['attributes.tm_balance.of']
    assert 'prefer must name an attribute' in errors['post_process.remove_overlap.prefer']


def test_starmap_chooses_the_layout_of_each_site_by_tm_balance(tmp_path, monkeypatch):
    # The published halves are 19-25 nt each, Tm-matched (median difference 0.9 C) around 58 C under
    # the template's conditions; remove_overlap kept the shortest layout, 19 + 1 + 20 nt, everywhere.
    monkeypatch.chdir(tmp_path)
    protocol = _validate_and_normalize_protocol(load_template('starmap'))
    filters = protocol['post_process']['filters']
    assert filters['padlock_tm']['condition'] == 'padlock_tm >= 50 & padlock_tm <= 65'
    assert filters['primer_tm']['condition'] == 'primer_tm >= 50 & primer_tm <= 65'
    assert protocol['post_process']['remove_overlap']['prefer'] == 'tm_balance'
    protocol['attributes'] = {k: v for k, v in protocol['attributes'].items() if v['type'] in OFFLINE_TYPES | {'difference'}}
    protocol['post_process']['filters'] = {k: filters[k] for k in ('padlock_tm', 'primer_tm')}
    protocol['post_process'].pop('equal_space')
    rng = random.Random(3)  # SEQ is poly-T: one Random(0) per base
    api = offline_api(protocol, ''.join(rng.choice('ACGT') for _ in range(400)), {'code': 'AAA'})
    targets = api.generate_target_seqs()
    combined = pd.concat([targets.reset_index(drop=True), api.construct_probes(targets).reset_index(drop=True)], axis=1)
    final = post_process(add_attributes(combined, protocol, None), protocol)
    assert len(final) >= 4
    assert final['target_region'].str.len().nunique() > 1   # not one fixed layout
    assert final['tm_balance'].max() < final['tm_balance'].quantile(0.5) + 5
    assert final['tm_balance'].median() < 1.5


def test_seqfish_plus_reproduces_lohoff_probes():
    # Lohoff et al. 2022, Supplementary Table 1A: the first two mouse Abcc4 probes (30-nt target
    # regions; seqFISH+ readouts and primer). The table lists the readout sites as they appear in
    # the probes, so the readout probes are their reverse complements.
    sites = ['TAATCAAAGGCCGCA', 'GGGCCAGCGCAATAA', 'GGCTAGCAGTCCGTC', 'TATCATACGACGTGG']
    primary = ['CGAGGGCTGTCCATTGGAGGTGTTCATAAC', 'TGTTTCGGTGAGGGCAGCGGCTTCGCCCTT']
    published = [
        'GTCTTATCATGTGCCTTTCCTTAATCAAAGGCCGCATGGGCCAGCGCAATAATCGAGGGCTGTCCATTGGAGGTGTTCATAACTGGCTAGCAGTCCGTCT'
        'TATCATACGACGTGGCCCTATAGTGAGTCGTATTA',
        'GTCTTATCATGTGCCTTTCCTGGGCCAGCGCAATAATGGCTAGCAGTCCGTCTTGTTTCGGTGAGGGCAGCGGCTTCGCCCTTTTATCATACGACGTGGT'
        'TAATCAAAGGCCGCACCCTATAGTGAGTCGTATTA',
    ]
    protocol = load_template('seqfish_plus')
    protocol['extracts']['target_region'].update(length=30, step=30)
    api = offline_api(_validate_and_normalize_protocol(protocol), ''.join(reverse_complement(p) for p in primary),
                      {'readouts': [reverse_complement(s) for s in sites]})
    probes = api.construct_probes(api.generate_target_seqs())
    assert probes['primary_probe'].tolist() == published  # probe_index 0 and 1: slots 12|34, then 23|41


def test_hcr_matches_insitu_probe_generator():
    # maker37cb.py reverse complements the mRNA (W) and prints
    # up + upspc + W[27:52] and W[0:25] + dnspc + dn for each 52-nt window.
    encoding = load_template('hcr_v3')['encoding']['GAPDH']
    targets, probes = design('hcr_v3', SEQ[:120], encoding)
    for window, odd, even in zip(targets['target_region'], probes['probe_odd'], probes['probe_even']):
        w = reverse_complement(window)
        assert odd == encoding['init_odd'] + encoding['spacer_odd'] + w[27:52]
        assert even == w[0:25] + encoding['spacer_even'] + encoding['init_even']


def test_hcr_pairs_are_two_nt_apart_as_in_choi_2018():
    # The published Choi 2018 sets tile pairs on a 54-nt grid: 52-nt windows, 2 nt apart.
    protocol = _validate_and_normalize_protocol(load_template('hcr_v3'))
    assert protocol['post_process']['remove_overlap'] == {'location_interval': 2}
    api = offline_api(protocol, SEQ[:270])
    windows = api.generate_target_seqs()
    kept = post_process(windows, {'post_process': {'remove_overlap': protocol['post_process']['remove_overlap']}})
    assert kept['start'].tolist() == [1, 55, 109, 163, 217]


def test_double_hyb_rca_keeps_units_without_off_target_genes():
    # One unit makes one RCA product, so no off-target gene is allowed on the window or the arm site;
    # fold_score only sorts, because a threshold would depend on the barcodes.
    protocol = load_template('double_hyb_rca')
    filters = protocol['post_process']['filters']
    assert filters['n_mapped_genes']['condition'] == 'n_mapped_genes <= 1'
    assert filters['arm_mapped_genes']['condition'] == 'arm_mapped_genes <= 1'
    assert protocol['attributes']['arm_mapped_genes']['target'] == 'target_parts.arm_site'
    assert protocol['extracts']['target_region']['layout']['parts']['arm_site']['length'] == 26
    assert not [name for name in filters if 'fold' in name]
    assert protocol['post_process']['sorts']['is_ascending'][:2] == ['circle_fold_score', 'amp_fold_score']


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
