import copy
import inspect
import json
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from uprobe.core.cli import _validate_and_normalize_protocol, cli
from uprobe.core.gen.expression import check_probe_expressions
from uprobe.core.validate import ATTRIBUTE_TYPES, render_text, validate_files, validate_protocol

DATA = Path(__file__).parent / 'data'
DEFAULT_PROTOCOL = Path(__file__).parents[1] / 'uprobe' / 'core' / 'agent' / 'templates' / 'DEFAULT_PROTOCOL.yaml'
DNA_TARGET = 'NC_000017.11:40304334-40362489'


def load(name):
    return yaml.safe_load((DATA / name).read_text(encoding='utf-8'))


def errors(findings):
    return [f for f in findings if f.severity == 'error']


def find(findings, path, severity='error'):
    matches = [f for f in findings if f.path == path and f.severity == severity]
    assert matches, f'no {severity} at {path}:\n{render_text(findings)}'
    return matches[0]


@pytest.fixture
def dna():
    return load('DNA_format.yaml')


@pytest.fixture
def rna():
    """RNA_format.yaml without its ignored min_mapq option."""
    protocol = load('RNA_format.yaml')
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    return protocol


@pytest.fixture
def genomes(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>1\n' + 'ACGT' * 100 + '\n')
    gtf = tmp_path / 'ref.gtf'
    gtf.write_text('1\ttest\texon\t1\t200\t.\t+\t.\tgene_id "ID1"; gene_name "Cryl1"; transcript_id "T1";\n'
                   '1\ttest\texon\t201\t400\t.\t-\t.\tgene_id "ID2"; gene_name "g42115"; transcript_id "T2";\n')
    return {'hg38': {'fasta': str(fasta), 'gtf': str(gtf), 'align_index': ['bowtie2'], 'jellyfish': True}}


# -- example protocols --------------------------------------------------

def test_dna_format_is_valid(dna, genomes):
    assert validate_protocol(dna) == []
    assert errors(validate_protocol(dna, genomes)) == []


def test_fixed_rna_format_is_valid(rna):
    assert validate_protocol(rna) == []


@pytest.mark.parametrize('name', ['DNA_format.yaml', 'RNA_format.yaml', 'double_hyb_rca.yaml'])
def test_test_data_protocols_are_valid(name):
    assert errors(validate_files(DATA / name)) == []


def test_rna_format_yields_candidates_for_every_target(tmp_path):
    from uprobe.core.api import UProbeAPI
    api = UProbeAPI(load('RNA_format.yaml'), {}, tmp_path, require_genome=False)
    targets = api.generate_target_seqs()
    assert set(targets['target']) == {'xxx', 'ccc'}
    probes = api.construct_probes(targets)
    assert len(probes) == len(targets) and (probes != '').all().all()


def test_legacy_names_get_fixes(tmp_path):
    # Protocols written before the renames gene_id -> target and n_mapped_genes ->
    # mapped_genes (as the test data were) get the current name as the fix.
    path = tmp_path / 'legacy.yaml'
    path.write_text((DATA / 'double_hyb_rca.yaml').read_text(encoding='utf-8')
                    .replace('encoding[target]', 'encoding[gene_id]')
                    .replace('type: mapped_genes', 'type: n_mapped_genes'))
    rca = validate_files(path)
    barcode1 = 'probes.circle_probe.parts.part2.parts.barcode1.expr'
    assert {f.path for f in errors(rca)} == {
        barcode1, 'probes.circle_probe.parts.part2.parts.barcode2.expr', 'attributes.n_mapped_genes.type'}
    assert "encoding[target]['barcode1']" in find(rca, barcode1).suggestion
    assert find(rca, barcode1).line == 35
    assert "'mapped_genes'" in find(rca, 'attributes.n_mapped_genes.type').suggestion


def test_default_protocol_template():
    findings = validate_files(DEFAULT_PROTOCOL)
    # Only the required placeholders of the template are errors; `uprobe run`
    # generates a blank `attributes:` and skips an all-blank post_process step.
    assert {f.path for f in errors(findings)} == {
        'genome', 'probes', 'extracts.target_region.source', 'extracts.target_region.length'}
    assert 'skips it' in find(findings, 'post_process.remove_overlap', 'warning').message
    protocol = yaml.safe_load(DEFAULT_PROTOCOL.read_text(encoding='utf-8'))
    protocol.update(name='filled', genome='hg38', encoding={'xxx1': {'BC1': 'ACGTACGTAC', 'BC2': 'TTGCAAGGTC'}})
    protocol['extracts']['target_region'].update(source='exon', length=40, step=20)
    protocol['probes'] = {'probe_1': {'template': '{part1}{part2}', 'parts': {
        'part1': {'expr': 'rc(target_region[0:20])'}, 'part2': {'expr': "encoding[target]['BC1']"}}}}
    protocol['summary']['report_name'] = 'rna_report'
    # As `uprobe run` normalizes it: defaults for the blank attributes and post_process.
    normalized = _validate_and_normalize_protocol(copy.deepcopy(protocol))
    assert normalized['attributes'] and normalized['post_process']['filters']
    assert 'remove_overlap' not in normalized['post_process']
    protocol['post_process']['remove_overlap']['location_interval'] = 0
    assert validate_protocol(protocol) == []


# -- YAML ---------------------------------------------------------------

def test_yaml_syntax_error(tmp_path):
    path = tmp_path / 'p.yaml'
    path.write_text("probes:\n  probe_1:\n    expr: 'encoding[target]['BC1']'\n")
    [finding] = validate_files(path)
    assert finding.severity == 'error' and finding.message.startswith('YAML syntax error')
    assert finding.line == 3 and finding.suggestion


def test_line_numbers_and_duplicate_keys(tmp_path):
    text = (DATA / 'DNA_format.yaml').read_text(encoding='utf-8').replace('source: genome', 'source: genomic')
    path = tmp_path / 'p.yaml'
    path.write_text(text + '\nname: again\n')
    findings = validate_files(path)
    source = find(findings, 'extracts.target_region.source')
    assert source.line == text.splitlines().index('    source: genomic') + 1
    assert source.file == str(path) and "'genome'" in source.suggestion
    assert 'duplicate' in find(findings, 'name', 'warning').message


# -- sections -----------------------------------------------------------

def test_top_level_sections(dna):
    del dna['probes']
    dna['post_proces'] = {}
    dna['hybridization'] = {'na_mM': 50}
    dna['meta'] = {'template': 'x'}
    findings = validate_protocol(dna)
    assert 'missing required section' in find(findings, 'probes').message
    assert "'post_process'" in find(findings, 'post_proces', 'warning').suggestion
    assert not [f for f in findings if f.path in ('hybridization', 'meta')]
    dna['probes'] = [{'template': '{a}'}]
    assert 'must be a mapping' in find(validate_protocol(dna), 'probes').message


@pytest.mark.parametrize('change,path,hint', [
    ({'source': 'exons'}, 'source', "did you mean 'exon'"),
    ({'length': 0}, 'length', 'length: 40'),
    ({'length': [40, 30]}, 'length', 'length: 40'),
    ({'step': 0}, 'step', 'step: 20'),
    ({'length': [60, 70]}, 'overlap', 'step'),
    ({'layout': {'template': '{a}{b}', 'parts': {'a': {'length': 30}}}}, 'layout', 'template'),
])
def test_target_region_checks(dna, change, path, hint):
    dna['extracts']['target_region'].update(change)
    finding = find(validate_protocol(dna), f'extracts.target_region.{path}')
    assert hint in finding.suggestion


def test_step_and_overlap_warns(dna):
    dna['extracts']['target_region']['step'] = 10
    assert find(validate_protocol(dna), 'extracts.target_region.overlap', 'warning')


def test_target_checks(dna, rna):
    dna['targets'].append('chr1')
    assert 'chrom:start-end' in find(validate_protocol(dna), 'targets[2]').message
    rna['targets'][0]['xxx'] = 'ACGT'
    assert 'no candidates' in find(validate_protocol(rna), 'targets[0].xxx').message


def test_probe_structure_checks(dna):
    parts = dna['probes']['probe_1']['parts']
    dna['probes']['probe_1']['template'] = '{part1}{prt2}{part3}'
    parts['part1']['expr'] = 'rc(target_regoin[0:20])'
    parts['part3']['expr'] = 'rc(target_region[0:20]'
    dna['probes']['probe_2'] = {'expr': 'rc(probe_1)'}
    dna['probes']['probe_3'] = {'template': '{a}', 'parts': [{'expr': 'target_region'}]}
    findings = validate_protocol(dna)
    assert "did you mean 'part2'" in find(findings, 'probes.probe_1.template').suggestion
    assert "did you mean 'target_region'" in find(findings, 'probes.probe_1.parts.part1.expr').suggestion
    assert 'syntax' in find(findings, 'probes.probe_1.parts.part3.expr').message
    assert 'expr: "rc(probe_1)"' in find(findings, 'probes.probe_2').suggestion
    assert find(findings, 'probes.probe_3.parts')
    assert find(findings, 'probes.probe_1.parts.part2', 'warning')  # not used in the template


def test_probe_reference_checks(dna):
    dna['probes']['probe_2'] = {'template': '{a}{b}', 'parts': {
        'a': {'expr': 'rc(probe_1.part9)'}, 'b': {'expr': "target_parts['left']"}}}
    findings = validate_protocol(dna)
    assert "one of 'probe_1.part1', 'probe_1.part2'" in find(findings, 'probes.probe_2.parts.a.expr').suggestion
    assert 'layout' in find(findings, 'probes.probe_2.parts.b.expr').suggestion


def test_dependency_cycle(dna):
    dna['probes']['probe_1']['parts']['part2']['expr'] = 'rc(probe_1)'
    finding = find(validate_protocol(dna), 'probes.probe_1.parts.part2.expr')
    assert finding.message == 'dependency cycle: probe_1.part2 -> probe_1 -> probe_1.part2'


def test_encoding_lookups(dna):
    del dna['encoding']['NC_000003.12:197795559-198295559']
    dna['encoding'][DNA_TARGET]['BC_2'] = dna['encoding'][DNA_TARGET].pop('BC2')
    findings = validate_protocol(dna)
    assert 'NC_000003.12:197795559-198295559' in find(findings, 'encoding').message
    assert "did you mean 'BC_2'" in find(findings, f"encoding['{DNA_TARGET}']").suggestion


@pytest.mark.parametrize('expr,message', [
    ('target_region[0:20].lower()', 'method or attribute call is not allowed'),
    ("'A' * 20", "operator '*' is not allowed"),
    ("__import__('os')", "call to '__import__()' is not allowed"),
    ('[b for b in target_region]', 'list comprehension is not allowed'),
    ('target_region.__class__', 'dunder name is not allowed'),
    ('target.upper', 'attribute access is not allowed'),
    ('rc(target_regoin[0:20])', "unknown name 'target_regoin'"),
])
def test_expressions_checked_like_the_evaluator(dna, expr, message):
    # The validator uses the grammar of uprobe.core.gen.expression and reports each problem once.
    dna['probes']['probe_1']['parts']['part2']['expr'] = expr
    path = 'probes.probe_1.parts.part2.expr'
    [finding] = [f for f in errors(validate_protocol(dna)) if f.path == path]
    assert message in finding.message
    assert check_probe_expressions(dna['probes'])  # `uprobe run` and the HTTP API reject it too


# -- attributes ----------------------------------------------------------

def test_attribute_checks(dna):
    attributes = dna['attributes']
    attributes['probe1_part2_tm']['type'] = 'annealing_temp'
    attributes['probe1_part2_gcContent']['target'] = 'probe_1:prt2'
    del attributes['probe1_part2_mappedSites']['aligner']
    findings = validate_protocol(dna)
    assert "'annealing_temperature'" in find(findings, 'attributes.probe1_part2_tm.type').suggestion
    assert "'probe_1.part2'" in find(findings, 'attributes.probe1_part2_gcContent.target').suggestion
    assert find(findings, 'attributes.probe1_part2_mappedSites.aligner').suggestion == 'aligner: bowtie2'


def test_aligner_requirements(dna, genomes):
    genomes['hg38'].update(jellyfish=False, align_index=['blast'])
    findings = validate_protocol(dna, genomes)
    assert 'jellyfish: true' in find(findings, 'attributes.probe1_part2_kmerCount').suggestion
    assert 'bowtie2' in find(findings, 'attributes.probe1_part2_mappedSites').suggestion
    for key in ('attributes', 'post_process', 'summary'):  # generated from the probes instead
        del dna[key]
    genomes['hg38']['align_index'] = ['bowtie2']
    finding = find(validate_protocol(dna, genomes), 'attributes')
    assert 'generated' in finding.message and 'probe_1_part1_kmer_count' in finding.message
    assert 'jellyfish: true' in finding.suggestion


def test_hybridization_checks(dna):
    attributes = dna['attributes']
    dna['hybridization'] = {'temperature_c': 37, 'na_mM': 390, 'formamid_pct': 30}
    attributes['probe1_part2_tm']['hybridization'] = {'na_mM': 300, 'probe_nM': 5}
    attributes['probe1_part2_foldScore']['hybridization'] = {'temperature_c': 'warm'}
    attributes['probe1_part2_gcContent']['hybridization'] = {'na_mM': 300}
    findings = validate_protocol(dna)
    assert 'formamid_pct' in find(findings, 'hybridization').message
    assert 'must be a number' in find(findings, 'attributes.probe1_part2_foldScore.hybridization').message
    assert 'only supported' in find(findings, 'attributes.probe1_part2_gcContent.hybridization').message
    assert not [f for f in findings if f.path.startswith('attributes.probe1_part2_tm')]
    dna['hybridization']['formamide_pct'] = dna['hybridization'].pop('formamid_pct')
    attributes['probe1_part2_foldScore']['hybridization'] = {'temperature_c': 30}
    del attributes['probe1_part2_gcContent']['hybridization']
    assert validate_protocol(dna) == []


def test_attribute_types_match_add_attributes():
    from uprobe.core.attributes import _add_attributes
    source = inspect.getsource(_add_attributes)
    for kind in ATTRIBUTE_TYPES:
        assert f'"{kind}"' in source


# -- post_process and summary -------------------------------------------

@pytest.mark.parametrize('condition,message', [
    ('probe1_part2_kmerCnt <= 100', 'unknown column'),
    ('probe1_part2_kmerCount => 100', "operator '=>'"),
    ('probe1_part2_kmerCount <= 100 | probe1_part2_tm >= 40', 'cannot read'),
    ('probe1_part2_kmerCount <= 1e5', 'read as'),
    ('probe1_part2_kmerCount >= -5', 'Invalid condition'),
])
def test_filter_conditions(dna, condition, message):
    dna['post_process']['filters']['probe1_part2_kmerCount']['condition'] = condition
    finding = find(validate_protocol(dna), 'post_process.filters.probe1_part2_kmerCount.condition')
    assert message in finding.message and finding.suggestion


def test_post_process_checks(dna):
    pp = dna['post_process']
    pp['filter'] = pp.pop('filters')
    pp['sorts']['is_ascending'].append('probe1_part2_gc')
    pp['remove_overlap'] = True
    pp['equal_space']['NC_000017.11:40304334-40362488'] = pp['equal_space'].pop(DNA_TARGET)
    findings = validate_protocol(dna)
    assert "'filters'" in find(findings, 'post_process.filter').suggestion
    assert find(findings, 'post_process.sorts.is_ascending[2]')
    assert 'location_interval' in find(findings, 'post_process.remove_overlap').suggestion
    assert DNA_TARGET in find(findings, "post_process.equal_space['NC_000017.11:40304334-40362488']",
                              'warning').suggestion


def test_avoid_otp_and_sequence_filters(dna):
    del dna['attributes']['probe1_part2_mappedSites']
    dna['post_process']['filters']['ends'] = {'type': 'sequence_pattern', 'target': 'probe_1',
                                              'exclude_patterns': ['[GC']}
    findings = validate_protocol(dna)
    assert 'mapped_sites' in find(findings, 'post_process.avoid_otp').message
    assert 'pattern' in find(findings, 'post_process.filters.ends').message


def test_rna_mode_drops_dna_only_steps(rna):
    rna['post_process']['equal_space'] = {'number_desired': 5}
    findings = validate_protocol(rna)
    assert errors(findings) == [] and find(findings, 'post_process.equal_space', 'warning')


def test_summary_attributes_exist(dna):
    dna['summary']['attributes'].append('probe1_part2_gcContnt')
    finding = find(validate_protocol(dna), 'summary.attributes[3]')
    assert finding.suggestion == "did you mean 'probe1_part2_gcContent'?"


# -- genome config -------------------------------------------------------

def test_genome_config_checks(dna, genomes):
    dna['genome'] = 'hg_38'
    assert "'hg38'" in find(validate_protocol(dna, genomes), 'genome').suggestion
    dna['genome'] = 'hg38'
    genomes['hg38']['fasta'] += '.missing'
    del genomes['hg38']['align_index']
    findings = validate_protocol(dna, genomes)
    assert find(findings, 'hg38.fasta').file == 'genomes'
    assert 'align_index' in find(findings, 'hg38').message


def test_genome_fetch_entry_keys_are_accepted(dna, genomes):
    # Entries written by `uprobe genome fetch` also carry description, species, out and source.
    expected = validate_protocol(dna, genomes)
    genomes['hg38'].update(
        description='Human GRCh38 (Ensembl release 115, primary_assembly)', species='homo_sapiens',
        out=str(Path(genomes['hg38']['fasta']).parent),
        source={'provider': 'Ensembl', 'division': 'vertebrates', 'release': 115, 'assembly': 'GRCh38',
                'sequence': 'primary_assembly', 'fasta_url': 'https://ftp.ensembl.org/pub/x.fa.gz',
                'gtf_url': 'https://ftp.ensembl.org/pub/x.gtf.gz', 'manifest': 'genome_source.json'})
    assert validate_protocol(dna, genomes) == expected


def test_check_targets_against_gtf_and_fasta(genomes):
    protocol = load('double_hyb_rca.yaml')
    protocol['genome'] = 'hg38'
    protocol['targets'] = ['Cryl1', 'cryl1', 'Nope1']
    findings = validate_protocol(protocol, genomes, check_targets=True)
    assert "did you mean 'Cryl1'" in find(findings, 'targets[1]').suggestion
    assert find(findings, 'targets[2]')
    assert not [f for f in findings if f.path == 'targets[0]']
    dna = load('DNA_format.yaml')
    dna['targets'] = ['1:1-100', '1:300-500', 'chr2:1-100']
    dna['encoding'] = {target: {'BC1': 'ACGT', 'BC2': 'ACGT'} for target in dna['targets']}
    dna['post_process'] = {}
    findings = validate_protocol(dna, genomes, check_targets=True)
    assert 'ends after' in find(findings, 'targets[1]').message
    assert 'not in the FASTA' in find(findings, 'targets[2]').message
    assert not [f for f in findings if f.path == 'targets[0]']


# -- probe construction on synthetic targets -----------------------------

def test_dry_run_reports_runtime_problems(rna):
    parts = rna['probes']['probe_1']['parts']
    parts['part1']['expr'] = 'rc(target_region[50])'
    parts['part3']['length'] = 12
    rna['probes']['probe_2']['parts']['part2']['expr'] = 'encoding[target]'
    rna['probes']['probe_3'] = {'template': '{a}', 'parts': {'a': {'expr': 'target_region[45:60]'}}}
    findings = validate_protocol(rna)
    assert 'index 50 out of range' in find(findings, 'probes.probe_1.parts.part1.expr').message
    assert 'declared length 12' in find(findings, 'probes.probe_1.parts.part3.length', 'warning').message
    assert 'not a nucleotide' in find(findings, 'probes.probe_2.parts.part2.expr', 'warning').message
    assert 'empty' in find(findings, 'probes.probe_3.parts.a.expr', 'warning').message


@pytest.mark.parametrize('template', ["open({!r}, 'w').name", "open({!r}, 'w')"])
def test_dry_run_never_runs_arbitrary_code(rna, tmp_path, template):
    marker = tmp_path / 'marker'
    rna['probes']['probe_2']['parts']['part2']['expr'] = template.format(str(marker))
    findings = validate_protocol(rna)
    assert 'not allowed' in find(findings, 'probes.probe_2.parts.part2.expr').message
    assert not marker.exists()


def test_validator_never_raises(dna, genomes):
    def paths(value, prefix=()):
        yield prefix
        children = value.items() if isinstance(value, dict) else enumerate(value) if isinstance(value, list) else []
        for key, child in children:
            yield from paths(child, prefix + (key,))

    for path in list(paths(dna))[1:]:
        for bad in (None, 'x', -1, ['x'], {'x': 1}):
            protocol = yaml.safe_load(yaml.safe_dump(dna))
            parent = protocol
            for key in path[:-1]:
                parent = parent[key]
            parent[path[-1]] = bad
            validate_protocol(protocol, genomes)


# -- command line --------------------------------------------------------

def test_cli_json_and_exit_codes(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, ['validate', '-p', str(DATA / 'DNA_format.yaml'), '--json'])
    assert result.exit_code == 0
    assert json.loads(result.output) == {'summary': {'valid': True, 'errors': 0, 'warnings': 0}, 'findings': []}
    result = runner.invoke(cli, ['validate', '-p', str(DATA / 'RNA_format.yaml'), '--json'])
    assert result.exit_code == 0
    assert json.loads(result.output)['summary'] == {'valid': True, 'errors': 0, 'warnings': 1}  # min_mapq
    broken = tmp_path / 'broken.yaml'
    broken.write_text((DATA / 'RNA_format.yaml').read_text(encoding='utf-8')
                      .replace('type: mapped_genes', 'type: n_mapped_genes'))
    result = runner.invoke(cli, ['validate', '-p', str(broken), '--json'])
    assert result.exit_code == 1
    report = json.loads(result.output)
    assert report['summary'] == {'valid': False, 'errors': 1, 'warnings': 1}
    assert set(report['findings'][0]) == {'severity', 'path', 'message', 'suggestion', 'file', 'line'}
    result = runner.invoke(cli, ['validate', '-p', str(broken)])
    assert result.exit_code == 1 and '  fix: ' in result.output and 'FAILED: 1 error(s)' in result.output
    result = runner.invoke(cli, ['validate', '-p', str(DATA / 'DNA_format.yaml'), '--check-targets'])
    assert result.exit_code == 2
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text('hg38: [')
    result = runner.invoke(cli, ['validate', '-p', str(DATA / 'DNA_format.yaml'), '-g', str(genomes), '--json'])
    assert result.exit_code == 1 and json.loads(result.output)['findings'][0]['file'] == str(genomes)
