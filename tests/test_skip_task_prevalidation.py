
import pandas as pd
import pytest
from uprobe.core.api import UProbeAPI

@pytest.mark.parametrize('skip', [True, False])
def test_skip_target_validation_reaches_design_without_changing_targets(monkeypatch, skip):
    from uprobe.core import cli
    api = object.__new__(UProbeAPI)
    api.protocol = {'targets': ['unknown'], 'extracts': {'target_region': {'source': 'transcript'}}}
    api.genome = {}
    api.transcripts_used = []
    api._parse_targets = lambda: (['unknown'], {})
    monkeypatch.setattr(cli, '_validate_and_normalize_protocol', lambda protocol, genome: protocol)
    calls = []
    api.validate_targets = lambda **kwargs: calls.append('validate') or False
    api.build_genome_index = lambda **kwargs: calls.append('design')
    api.generate_target_seqs = lambda: pd.DataFrame()
    run = {}
    api._run_workflow(run, False, False, 1, skip_target_validation=skip)
    assert calls == (['design'] if skip else ['validate'])
    assert api.protocol['targets'] == ['unknown']
    assert run['status'] == ('no_target_sequences' if skip else 'no_valid_targets')
