"""Output files are named after the protocol's `name`; the name never leaves the output directory."""
from pathlib import Path

import pytest

from uprobe.core.api import UProbeAPI, _file_stem

SEQ = ('ATGCGTACGTTAGCCGATCGATGCTAGCTAGGCTAGCATCGATCGGCTAGCTAGCGATCGATCGTAGCTAGCTGACTGATCGATGCTAGCTAGCTAGG'
       'GGATCGATTCGATCGGATTAGCTAGCATCGATCGATGCTAGCTAGGATCGATCGTAGCTAGCTAGGCTAGCTAGCATCGATGCATCGAT')


@pytest.mark.parametrize('name, stem', [
    ('summary_test', 'summary_test'),
    ('my-probes.v2', 'my-probes.v2'),
    ('../../ESCAPED', '_.._ESCAPED'),
    ('/etc/passwd', '_etc_passwd'),
    ('a\\b', 'a_b'),
    ('..', 'probes'),
    ('', 'probes'),
    (None, 'None'),
    ('探针 设计', '探针_设计'),
])
def test_file_stem(name, stem):
    assert _file_stem(name) == stem
    assert '/' not in _file_stem(name) and '\\' not in _file_stem(name)


def test_outputs_stay_in_the_output_directory(tmp_path):
    output_dir = tmp_path / 'task' / 'out'
    protocol = {
        'name': '../../ESCAPED', 'genome': 'none', 'encoding': {}, 'targets': [{'G1': SEQ}],
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 7}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
        'post_process': {'filters': {'target_gc': {'condition': 'target_gc >= 0.1'}},
                         'remove_overlap': {'location_interval': 0}},
        'summary': {'report_name': 'rna_report', 'attributes': ['target_gc']},
    }
    api = UProbeAPI(protocol, {}, output_dir, require_genome=False)
    api.build_genome_index = lambda **kwargs: None
    api.run_workflow(raw_csv=True)
    named = [p for p in tmp_path.rglob('*ESCAPED*') if p.is_file()]
    assert named, 'the run wrote no file named after the protocol'
    assert all(p.parent == output_dir for p in named), [str(p) for p in named]
    assert not (tmp_path / 'ESCAPED').exists() and not list(tmp_path.parent.glob('*ESCAPED*'))
