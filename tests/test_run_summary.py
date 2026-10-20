import copy
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner
from pyfaidx import Fasta

from uprobe.core.api import InvalidTargetsError, UProbeAPI
from uprobe.core.cli import cli
from uprobe.core.gen.fun import extract_trans_seqs, generate_target_seqs
from uprobe.core.hybridization import TM_DEFAULTS
from uprobe.core.process import filter_table
from uprobe.core.sampling import iter_target_windows

ROOT = Path(__file__).parents[1]
TOY = ROOT / 'tests/examples/toy_genome'
MIXED = 'ATGCGTACGTTAGCCGATCGATGCTAGCTAGGCTAGCATCGATCGGCTAGCTAGCGATCGATCGTAGCTAGCTGACTGATCGATGCTAGCTAGCTAGG' * 2
POLY_A = 'A' * 150
FILTERS = {'target_gc': {'condition': 'target_gc >= 0.3'}}


def make_api(tmp_path, targets=None, filters=FILTERS, report_name='rna_report'):
    protocol = {
        'name': 'summary_test', 'genome': 'none', 'encoding': {},
        'targets': targets or [{'G1': MIXED}, {'G2': POLY_A}],
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 7}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': {'target_gc': {'target': 'target_region', 'type': 'gc_content'},
                       'target_tm': {'target': 'target_region', 'type': 'annealing_temperature'}},
        'post_process': {'filters': copy.deepcopy(filters), 'remove_overlap': {'location_interval': 0},
                         'sorts': {'is_ascending': ['target_tm']}},
        'summary': {'report_name': report_name, 'attributes': ['target_gc', 'target_tm']},
    }
    api = UProbeAPI(protocol, {}, tmp_path, require_genome=False)
    api.build_genome_index = lambda **kwargs: None
    return api


def read_summary(directory):
    return json.loads((Path(directory) / 'summary.json').read_text(encoding='utf-8'))


def test_summary_matches_result_tables(tmp_path):
    api = make_api(tmp_path)
    result = api.run_workflow(raw_csv=True)
    summary = read_summary(tmp_path)
    assert summary == json.loads(json.dumps(api.summary))
    assert summary['status'] == 'partial' and summary['stopped_at'] is None and summary['error'] is None
    assert summary['protocol_name'] == 'summary_test'
    assert summary['inputs']['requested_targets'] == ['G1', 'G2']
    assert summary['inputs']['valid_targets'] == ['G1', 'G2'] and summary['inputs']['invalid_targets'] == []
    # No hybridization block: the legacy (fisheye-compatible) Tm conditions are in effect.
    assert summary['inputs']['hybridization'] is None
    assert summary['inputs']['hybridization_conditions'] == {'target_tm': TM_DEFAULTS}
    assert summary['inputs']['genome_source'] is None

    raw = api._report_raw_data
    after_filters = filter_table(raw.copy(), FILTERS).groupby('target').size()
    final = result.groupby('target').size()
    for name, seq in [('G1', MIXED), ('G2', POLY_A)]:
        entry = summary['targets'][name]
        windows = len(list(iter_target_windows(seq, 40, 7)))
        assert entry['candidate_windows'] == entry['probes_constructed'] == windows
        assert entry['after_stage']['Filters'] == after_filters.get(name, 0)
        assert entry['final'] == entry['after_stage']['Sort'] == final.get(name, 0)
        counts = [windows] + list(entry['after_stage'].values())
        drops = {stage: before - after for stage, before, after in zip(entry['after_stage'], counts, counts[1:])}
        top = max(drops, key=drops.get)
        assert entry['top_removal_stage'] == (top if drops[top] > 0 else None)
    assert summary['targets']['G2']['final'] == 0 and summary['targets']['G2']['top_removal_stage'] == 'Filters'
    assert summary['targets']['G1']['top_removal_stage'] == 'Remove overlap'
    assert summary['counts'] == {'requested_targets': 2, 'valid_targets': 2, 'candidate_windows': len(raw),
                                 'probes_constructed': len(raw), 'final_probes': len(result),
                                 'targets_with_probes': 1}
    assert [s['stage'] for s in summary['stages']] == ['Filters', 'Remove overlap', 'Sort']
    assert summary['stages'][-1]['after'] == len(result)

    for name in ('target_gc', 'target_tm'):
        stats = summary['attributes'][name]
        assert stats['count'] == len(result)
        assert stats['min'] == pytest.approx(result[name].min(), abs=1e-4)
        assert stats['median'] == pytest.approx(result[name].median(), abs=1e-4)
        assert stats['max'] == pytest.approx(result[name].max(), abs=1e-4)

    files = summary['output_files']
    assert all(files[key] and (tmp_path / files[key]).is_file() for key in ('xlsx', 'raw_xlsx', 'tsv', 'html', 'summary'))
    assert files['tsv'] == files['xlsx'].replace('.xlsx', '_probes.tsv')
    table = pd.read_csv(tmp_path / files['tsv'], sep='\t')
    assert list(table.columns) == list(result.columns)
    assert table['probe_id'].tolist() == result['probe_id'].tolist()
    assert table['target_gc'].tolist() == pytest.approx(result['target_gc'].tolist())
    assert summary['flags'] == {'no_filtered_probes': False, 'raw_file': None}
    assert summary['warnings'] == [{'code': 'targets_without_probes', 'targets': ['G2'],
                                    'eliminated_by': {'G2': 'Filters'},
                                    'message': '1 target(s) have no final probes.'}]


def test_summary_lists_files_for_a_relative_output_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    make_api(Path('out')).run_workflow()
    files = read_summary(tmp_path / 'out')['output_files']
    assert all(files[key] and (tmp_path / 'out' / files[key]).is_file() for key in ('xlsx', 'tsv', 'html'))


def test_summary_records_conditions_and_genome_source(tmp_path):
    api = make_api(tmp_path)
    api.protocol['hybridization'] = {'temperature_c': 37, 'na_mM': 390, 'formamide_pct': 10}
    api.protocol['attributes']['target_tm']['hybridization'] = {'probe_nM': 5}
    source = {'provider': 'Ensembl', 'division': 'vertebrates', 'release': 115, 'assembly': 'GRCh38'}
    api.genome = {'source': source}
    api.run_workflow()
    inputs = read_summary(tmp_path)['inputs']
    assert inputs['hybridization'] == {'temperature_c': 37, 'na_mM': 390, 'formamide_pct': 10}
    assert inputs['hybridization_conditions'] == {
        'target_tm': {**TM_DEFAULTS, 'na_mM': 390, 'probe_nM': 5, 'formamide_pct': 10}}
    assert inputs['genome_source'] == source


def test_summary_when_nothing_passes_post_processing(tmp_path):
    api = make_api(tmp_path, filters={'target_gc': {'condition': 'target_gc >= 2'}}, report_name='')
    result = api.run_workflow()
    summary = read_summary(tmp_path)
    assert result.empty
    assert summary['status'] == 'no_probes'
    assert summary['flags'] == {'no_filtered_probes': True, 'raw_file': api.raw_file}
    assert (tmp_path / api.raw_file).is_file()
    assert summary['output_files'] == {'xlsx': None, 'raw_xlsx': api.raw_file, 'tsv': None, 'raw_tsv': None,
                                       'html': None, 'summary': 'summary.json'}
    assert not list(tmp_path.glob('*.tsv'))
    assert summary['counts']['final_probes'] == 0 and summary['counts']['candidate_windows'] > 0
    assert all(t['final'] == 0 and t['top_removal_stage'] == 'Filters' for t in summary['targets'].values())
    assert summary['attributes']['target_gc'] == {'count': 0, 'min': None, 'median': None, 'max': None}
    codes = {w['code']: w for w in summary['warnings']}
    assert codes['stage_removed_all']['stage'] == 'Filters'
    assert codes['stage_removed_all']['removed'] == summary['counts']['candidate_windows']
    assert codes['targets_without_probes']['eliminated_by'] == {'G1': 'Filters', 'G2': 'Filters'}


def test_summary_when_no_target_sequences(tmp_path):
    api = make_api(tmp_path, targets=[{'Short': 'ACGT' * 5}])
    assert api.run_workflow().empty
    summary = read_summary(tmp_path)
    assert summary['status'] == 'no_target_sequences'
    assert summary['stopped_at'] == 'extract_target_regions'
    assert summary['counts']['candidate_windows'] == 0 and summary['counts']['probes_constructed'] is None
    assert summary['targets']['Short']['candidate_windows'] == 0
    assert [w['code'] for w in summary['warnings']] == ['no_candidate_windows']
    assert summary['output_files']['xlsx'] is None


def test_summary_when_no_probe_can_be_built(tmp_path):
    # A blank encoding: every probe fails to build, so the run stops instead of
    # reporting candidates with empty probe sequences as designed probes.
    api = make_api(tmp_path)
    api.protocol['probes']['pad']['parts']['part1']['expr'] = "encoding[target]['BC1']"
    api.protocol['encoding'] = None
    assert api.run_workflow().empty
    summary = read_summary(tmp_path)
    assert summary['status'] == 'no_probes_constructed' and summary['stopped_at'] == 'construct_probes'
    assert summary['counts']['probes_constructed'] == 0 and summary['counts']['final_probes'] == 0
    windows = summary['counts']['candidate_windows']
    assert [w for w in summary['warnings'] if w['code'] == 'probe_construction_failed'][0]['count'] == windows
    assert summary['output_files']['tsv'] is None and not list(tmp_path.glob('*.tsv'))


def test_summary_for_invalid_targets(tmp_path):
    gtf = tmp_path / 'ref.gtf'
    gtf.write_text('1\ttest\texon\t1\t48\t.\t-\t.\tgene_name "G"; transcript_id "T";\n')
    api = make_api(tmp_path / 'out', targets=['Missing'])
    api.genome = {'gtf': str(gtf)}
    assert api.run_workflow(continue_on_invalid_targets=True).empty
    summary = read_summary(tmp_path / 'out')
    assert summary['status'] == 'no_valid_targets' and summary['stopped_at'] == 'validate_targets'
    assert summary['inputs']['valid_targets'] == [] and summary['inputs']['invalid_targets'] == ['Missing']
    assert summary['warnings'][0]['code'] == 'invalid_targets'

    api = make_api(tmp_path / 'out', targets=['Missing'])
    api.genome = {'gtf': str(gtf)}
    with pytest.raises(InvalidTargetsError, match='1 target.s. not found in the genome annotation: Missing'):
        api.run_workflow()
    summary = read_summary(tmp_path / 'out')
    assert summary['status'] == 'error' and summary['stopped_at'] == 'validate_targets'
    assert 'Missing' in summary['error'] and '--continue-invalid' in summary['error']
    assert summary['inputs']['invalid_targets'] == ['Missing']


def test_any_invalid_target_stops_the_run_unless_continue_invalid(tmp_path):
    # One gene target is not in the GTF; the direct sequence G1 is valid.
    gtf = tmp_path / 'ref.gtf'
    gtf.write_text('1\ttest\texon\t1\t48\t.\t-\t.\tgene_name "G"; transcript_id "T";\n')
    api = make_api(tmp_path / 'stop', targets=[{'G1': MIXED}, 'Missing'])
    api.genome = {'gtf': str(gtf)}
    indexed = []
    api.build_genome_index = lambda **kwargs: indexed.append(kwargs)
    with pytest.raises(InvalidTargetsError) as error:
        api.run_workflow()
    assert error.value.targets == ['Missing'] and '--continue-invalid' in str(error.value)
    assert indexed == []  # stopped before the indexes are built and anything is designed
    summary = read_summary(tmp_path / 'stop')
    assert summary['status'] == 'error' and summary['stopped_at'] == 'validate_targets'
    assert summary['inputs']['invalid_targets'] == ['Missing'] and summary['counts']['candidate_windows'] is None
    warning = summary['warnings'][0]
    assert warning['code'] == 'invalid_targets' and 'stopped the run' in warning['message']
    assert not list((tmp_path / 'stop').glob('*.tsv'))

    api = make_api(tmp_path / 'skip', targets=[{'G1': MIXED}, 'Missing'])
    api.genome = {'gtf': str(gtf)}
    result = api.run_workflow(continue_on_invalid_targets=True)
    summary = read_summary(tmp_path / 'skip')
    assert set(result['target']) == {'G1'} and summary['status'] == 'partial'
    assert summary['inputs']['invalid_targets'] == ['Missing'] and summary['targets']['G1']['final'] > 0
    assert summary['warnings'][0] == {'code': 'invalid_targets', 'targets': ['Missing'],
                                      'message': '1 target(s) not found in the genome annotation and skipped.'}


def test_cli_stops_on_an_invalid_target_without_continue_invalid(tmp_path):
    genomes = toy_genome(tmp_path, align_index=[])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'}})
    output = tmp_path / 'out'
    result = CliRunner().invoke(cli, ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(output), '--json'])
    assert result.exit_code == 1
    summary = json.loads(result.stdout)
    assert summary['status'] == 'error' and summary['stopped_at'] == 'validate_targets'
    assert 'NotAGene' in summary['error'] and '--continue-invalid' in summary['error']
    assert summary['inputs']['invalid_targets'] == ['NotAGene'] and summary['output_files']['tsv'] is None
    # The step commands stop too; validate-targets accepts the flag for a partly valid list.
    args = ['-p', str(protocol), '-g', str(genomes)]
    assert CliRunner().invoke(cli, ['validate-targets', *args]).exit_code == 1
    assert CliRunner().invoke(cli, ['validate-targets', *args, '--continue-invalid']).exit_code == 0
    assert CliRunner().invoke(cli, ['generate-targets', *args, '-o', str(tmp_path / 'g')]).exit_code == 1
    result = CliRunner().invoke(cli, ['generate-targets', *args, '-o', str(tmp_path / 'g'), '--continue-invalid'])
    assert result.exit_code == 0 and (tmp_path / 'g' / 'target_sequences.csv').is_file()


def toy_genome(tmp_path, align_index):
    """Copy the toy reference; its FASTA is a Git LFS pointer unless LFS objects were fetched."""
    fasta = tmp_path / 'toy_uprobe.fa'
    source = (TOY / 'toy_uprobe.fa').read_text()
    if source.startswith('version https://git-lfs'):
        # Regenerate with the seed from prepare_wsl.py and check the pointer's checksum.
        rng = random.Random(20261004)
        sequences = {'1': ''.join(rng.choices('ACGT', k=12000)), '2': ''.join(rng.choices('ACGT', k=8000))}
        source = ''.join(f'>{c}\n' + ''.join(s[i:i + 80] + '\n' for i in range(0, len(s), 80))
                         for c, s in sequences.items())
        oid = (TOY / 'toy_uprobe.fa').read_text().split('sha256:')[1].split()[0]
        if hashlib.sha256(source.encode()).hexdigest() != oid:
            pytest.skip('toy FASTA unavailable')
    fasta.write_text(source)
    shutil.copy(TOY / 'toy_uprobe.gtf', tmp_path / 'toy_uprobe.gtf')
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump({'toy': {'fasta': str(fasta), 'gtf': str(tmp_path / 'toy_uprobe.gtf'),
                                               'align_index': align_index}}))
    return genomes


def toy_protocol(tmp_path, attributes, targets=('ToyGeneA', 'ToyGeneB', 'ToyGeneC', 'NotAGene'), post_process=None):
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({
        'name': 'toy', 'genome': 'toy', 'encoding': {},
        'targets': list(targets),
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 60}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': attributes,
        'post_process': post_process or {'filters': {'target_gc': {'condition': 'target_gc >= 0.4 & target_gc <= 0.6'}}},
        'summary': {'report_name': 'rna_report', 'attributes': list(attributes)},
    }))
    return protocol


def check_toy_summary(summary, output):
    assert summary == read_summary(output)
    assert summary['status'] == 'partial'
    assert summary['inputs']['genome'] == 'toy' and summary['inputs']['extract_source'] == 'exon'
    assert summary['inputs']['valid_targets'] == ['ToyGeneA', 'ToyGeneB', 'ToyGeneC']
    assert summary['inputs']['invalid_targets'] == ['NotAGene']
    table = pd.read_csv(output / summary['output_files']['tsv'], sep='\t')
    assert summary['counts']['final_probes'] == len(table) > 0
    assert {name: t['final'] for name, t in summary['targets'].items()} == table.groupby('target').size().to_dict()
    assert sum(t['candidate_windows'] for t in summary['targets'].values()) == summary['counts']['candidate_windows']
    # List cells are JSON, unquoted.
    assert all(json.loads(cell) for cell in table['transcript_names'])


def test_cli_json_on_toy_genome_and_failure(tmp_path):
    genomes = toy_genome(tmp_path, align_index=[])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'},
                                       'target_tm': {'target': 'target_region', 'type': 'annealing_temperature'}})
    output = tmp_path / 'out'
    args = ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(output), '--continue-invalid', '--json']
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0, result.output
    check_toy_summary(json.loads(result.stdout), output)

    broken = tmp_path / 'broken.yaml'
    broken.write_text(protocol.read_text().replace('genome: toy', 'genome: absent'))
    result = CliRunner().invoke(cli, ['run', '-p', str(broken), '-g', str(genomes), '-o', str(output), '--json'])
    assert result.exit_code == 1
    summary = json.loads(result.stdout)
    assert summary['status'] == 'error' and "Genome 'absent' not found" in summary['error']
    assert read_summary(output) == summary


def test_cli_json_stdout_stays_clean_with_tool_output(tmp_path):
    """bowtie2-build writes to stdout; with --json it must land on stderr."""
    env = dict(os.environ, PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', ''))
    if not shutil.which('bowtie2-build', path=env['PATH']):
        pytest.skip('bowtie2 not installed')
    genomes = toy_genome(tmp_path, align_index=['bowtie2'])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'},
                                       'target_mapped_genes': {'target': 'target_region', 'type': 'mapped_genes',
                                                               'aligner': 'bowtie2'}})
    output = tmp_path / 'out'
    process = subprocess.run(
        [sys.executable, '-m', 'uprobe', 'run', '-p', str(protocol), '-g', str(genomes), '-o', str(output),
         '--continue-invalid', '--threads', '2', '--json'],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
    assert process.returncode == 0, process.stderr[-2000:]
    assert 'bowtie2' in process.stderr
    summary = json.loads(process.stdout)
    check_toy_summary(summary, output)
    assert summary['attributes']['target_mapped_genes']['count'] == summary['counts']['final_probes']


def run_toy_cli(tmp_path, post_process):
    genomes = toy_genome(tmp_path, align_index=[])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
                            targets=('ToyGeneA', 'ToyGeneB', 'ToyGeneC'), post_process=post_process)
    output = tmp_path / 'out'
    result = CliRunner().invoke(cli, ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(output), '--json'])
    return result, output


def test_rna_windows_are_numbered_along_the_spliced_rna(tmp_path):
    # Exon windows restart at 1 in every exon; target_start/target_end count along the
    # spliced RNA, 5' to 3', for ToyGeneA (+ strand) and ToyGeneB (- strand, 3 exons each).
    genomes = yaml.safe_load(toy_genome(tmp_path, align_index=[]).read_text())['toy']
    windows = generate_target_seqs('exon', ['ToyGeneA', 'ToyGeneB'], genomes['fasta'], genomes['gtf'], 40, step=97)
    extract_trans_seqs(genomes['gtf'], genomes['fasta'], str(tmp_path / 'transcripts.fa'))
    with Fasta(str(tmp_path / 'transcripts.fa')) as transcripts:
        rna = {gene: str(transcripts[f'{gene}_{gene}.t1']).upper() for gene in ('ToyGeneA', 'ToyGeneB')}
    for row in windows.itertuples():
        assert rna[row.target][row.target_start - 1:row.target_end] == row.target_region
    assert windows.groupby('target')['exon_name'].nunique().tolist() == [3, 3]
    assert (windows['start'] == 1).sum() == 6 and windows['target_start'].max() > 1800


def test_cli_keeps_probes_per_target_spread_along_the_rna(tmp_path):
    # 45 windows per 3-exon gene (15 per exon), 40 for ToyGeneC; 6 are kept, evenly spaced ranks.
    result, output = run_toy_cli(tmp_path, {'remove_overlap': {'location_interval': 0},
                                            'equal_space': {'number_desired': 6}})
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout)
    assert summary['status'] == 'ok' and summary['message'] == '18 probes designed for all 3 target(s).'
    assert {name: (t['final'], t['requested_probes']) for name, t in summary['targets'].items()} == {
        'ToyGeneA': (6, 6), 'ToyGeneB': (6, 6), 'ToyGeneC': (6, 6)}
    assert [s['stage'] for s in summary['stages']] == ['Remove overlap', 'Equal spacing']
    table = pd.read_csv(output / summary['output_files']['tsv'], sep='\t')
    for gene in ('ToyGeneA', 'ToyGeneB'):
        probes = table[table['target'] == gene]
        # Ranks 0, 8, 17, 26, 35, 44 along the target: two probes in each exon.
        assert probes['target_start'].tolist() == [1, 481, 1021, 1561, 2101, 2641]
        assert probes['exon_name'].value_counts().tolist() == [2, 2, 2]


def test_cli_reports_targets_with_fewer_probes_than_requested(tmp_path, monkeypatch):
    from unittest import mock
    from uprobe.core import cli as cli_module
    log = mock.Mock()
    monkeypatch.setattr(cli_module, 'log', log)
    result, output = run_toy_cli(tmp_path, {'remove_overlap': {'location_interval': 0},
                                            'equal_space': {'number_desired': 6, 'ToyGeneC': {'number_desired': 100}}})
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout)
    assert summary == read_summary(output)
    assert summary['status'] == 'partial'
    assert summary['message'] == '52 probes designed; 1 target(s) have fewer probes than requested.'
    assert summary['targets']['ToyGeneC']['final'] == 40 and summary['targets']['ToyGeneC']['requested_probes'] == 100
    warning = summary['warnings'][-1]
    assert warning['code'] == 'fewer_probes_than_requested' and warning['targets'] == ['ToyGeneC']
    assert warning['counts'] == {'ToyGeneC': {'requested': 100, 'final': 40}}
    logged = [str(call) for call in log.method_calls]
    assert any('ToyGeneA 6/6, ToyGeneB 6/6, ToyGeneC 40/100' in line for line in logged)
    assert any('warning' in line and 'ToyGeneC 40 of 100' in line for line in logged)


NOTHING_PASSES = {'filters': {'target_gc': {'condition': 'target_gc >= 2'}}}


def test_targets_without_probes_are_not_reported_as_short(tmp_path):
    # A run where nothing passed also warned that 'all their passing probes were kept' for its
    # targets; a target with no probe is reported once, as targets_without_probes.
    result, output = run_toy_cli(tmp_path, {**NOTHING_PASSES, 'equal_space': {'number_desired': 6}})
    assert result.exit_code == 3, result.output
    summary = read_summary(output)
    assert summary['status'] == 'no_probes'
    codes = [warning['code'] for warning in summary['warnings']]
    assert 'targets_without_probes' in codes and 'fewer_probes_than_requested' not in codes
    assert summary['targets']['ToyGeneA']['requested_probes'] == 6


@pytest.mark.parametrize('json_flag', [['--json'], []])
def test_cli_exits_with_3_when_no_probe_passes(tmp_path, json_flag):
    genomes = toy_genome(tmp_path, align_index=[])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
                            targets=('ToyGeneA',), post_process=NOTHING_PASSES)
    output = tmp_path / 'out'
    result = CliRunner().invoke(cli, ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(output), *json_flag])
    assert result.exit_code == 3, result.output
    summary = read_summary(output)
    assert summary['status'] == 'no_probes' and summary['counts']['candidate_windows'] > 0
    assert json.loads(result.stdout) == summary if json_flag else result.stdout == ''


def test_cli_exit_codes_for_construction_failure_usage_and_success(tmp_path):
    genomes = toy_genome(tmp_path, align_index=[])
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
                            targets=('ToyGeneA',))
    output = tmp_path / 'out'
    args = ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(output), '--json']
    assert CliRunner().invoke(cli, args).exit_code == 0
    # A blank encoding: no probe can be built.
    config = yaml.safe_load(protocol.read_text())
    config['probes']['pad']['parts']['part1']['expr'] = "encoding[target]['BC1']"
    config['encoding'] = None
    protocol.write_text(yaml.safe_dump(config))
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 3 and json.loads(result.stdout)['status'] == 'no_probes_constructed'
    # All targets invalid, continued with --continue-invalid.
    config = yaml.safe_load(protocol.read_text().replace('ToyGeneA', 'NotAGene'))
    protocol.write_text(yaml.safe_dump(config))
    result = CliRunner().invoke(cli, args + ['--continue-invalid'])
    assert result.exit_code == 3 and json.loads(result.stdout)['status'] == 'no_valid_targets'
    # Without --continue-invalid this is an error (1); a missing option is a usage error (2).
    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 1 and json.loads(result.stdout)['status'] == 'error'
    assert CliRunner().invoke(cli, ['run', '-g', str(genomes)]).exit_code == 2


def test_http_runner_is_unaffected_by_the_exit_code(tmp_path, monkeypatch):
    # The web server runs the workflow through UProbeAPI, not the CLI: no SystemExit.
    from uprobe.http.utils import task_control
    from uprobe.http.utils.paths import get_genomes_yaml
    from uprobe.http.utils.uprobe_runner import run_uprobe_workflow
    monkeypatch.setenv('UPROBE_SERVER_ROOT', str(tmp_path / 'server'))
    monkeypatch.setenv('PATH', os.environ['PATH'])  # restored after the runner prepends to it
    # It detaches the worker process (setsid) and reads /proc.
    monkeypatch.setattr(task_control, 'register_task_process', lambda *args, **kwargs: None)
    genomes = toy_genome(tmp_path, align_index=[])
    get_genomes_yaml().write_text(genomes.read_text())
    protocol = toy_protocol(tmp_path, {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
                            targets=('ToyGeneA',), post_process=NOTHING_PASSES)
    result = run_uprobe_workflow(protocol_yaml=protocol.read_text(), username='tester', task_id='t1',
                                 output_dir=str(tmp_path / 'task'), threads=1, raw_csv=False,
                                 continue_invalid_targets=False, log_path=str(tmp_path / 'task.log'))
    assert result['ok'] and result['no_filtered_probes'], result.get('error')
    assert read_summary(tmp_path / 'task')['status'] == 'no_probes'


def test_web_tasks_skip_invalid_targets(tmp_path, monkeypatch):
    # The web form has no --continue-invalid switch: web tasks skip invalid targets, as before.
    import asyncio
    from concurrent.futures import ThreadPoolExecutor
    from types import SimpleNamespace
    from uprobe.http.routers import task as task_router
    from uprobe.http.utils import task_queue, uprobe_runner
    monkeypatch.setenv('UPROBE_SERVER_ROOT', str(tmp_path / 'server'))
    monkeypatch.setattr(task_queue, '_task_semaphore', None)
    calls = []
    monkeypatch.setattr(uprobe_runner, 'run_uprobe_workflow',
                        lambda **kwargs: calls.append(kwargs) or {'ok': False, 'error': 'stub'})
    record = SimpleNamespace(status='pending', progress=0, yaml_content='name: x')
    monkeypatch.setattr(task_router, 'find_task_by_id', lambda username, task_id: record)
    monkeypatch.setattr(task_router, 'update_task_in_db', lambda username, task: None)
    with ThreadPoolExecutor(1) as pool:
        monkeypatch.setattr(task_router, 'get_process_pool', lambda: pool)
        asyncio.run(task_router._run_uprobe_task('tester', 't1'))
    assert [call['continue_invalid_targets'] for call in calls] == [True]
    assert record.status == 'failed' and record.error_message == 'stub'
