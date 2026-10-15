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

from uprobe.core.api import UProbeAPI
from uprobe.core.cli import cli
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
    assert summary['output_files'] == {'xlsx': None, 'raw_xlsx': api.raw_file, 'tsv': None,
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
    with pytest.raises(ValueError, match='No valid targets'):
        api.run_workflow()
    summary = read_summary(tmp_path / 'out')
    assert summary['status'] == 'error' and summary['stopped_at'] == 'validate_targets'
    assert 'No valid targets' in summary['error']
    assert summary['inputs']['invalid_targets'] == ['Missing']


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


def toy_protocol(tmp_path, attributes):
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({
        'name': 'toy', 'genome': 'toy', 'encoding': {},
        'targets': ['ToyGeneA', 'ToyGeneB', 'ToyGeneC', 'NotAGene'],
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 60}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': attributes,
        'post_process': {'filters': {'target_gc': {'condition': 'target_gc >= 0.4 & target_gc <= 0.6'}}},
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
