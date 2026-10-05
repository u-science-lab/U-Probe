import pandas as pd
import pytest
from uprobe.core.process import filter_table, post_process
from uprobe.core.report.compact import build_report_data


def test_generic_terminal_exclusions_and_report():
    raw = pd.DataFrame({'probe_id': list('abcdefg'), 'custom.part': ['GCAT', 'ATCG', 'ATAT', 'atat', None, '', 'GGCC'], 'score': [1]*7})
    filters = {'ends': {'type': 'sequence_pattern', 'target': 'custom:part', 'exclude_patterns': ['^[GC]{2}', '[GC]{2}$']}}
    result = filter_table(raw, filters)
    assert result.probe_id.tolist() == ['c', 'd']
    data = build_report_data(result, {'post_process': {'filters': filters}}, raw)
    assert data['filterCount'] == 2
    assert data['diagnostics'][0]['rejected'] == 5
    assert 'custom:part' in data['diagnostics'][0]['condition']


@pytest.mark.parametrize('config', [
    {'target': 'missing', 'exclude_patterns': ['GC']},
    {'target': 'seq', 'exclude_patterns': ['[']},
    {'target': 'seq', 'exclude_patterns': []},
])
def test_invalid_sequence_filters_fail_explicitly(config):
    with pytest.raises(ValueError):
        filter_table(pd.DataFrame({'seq': ['ACGT']}), {'bad': {'type': 'sequence_pattern', **config}})


def test_sequence_and_numeric_filters_combine():
    raw = pd.DataFrame({'seq': ['GCAT', 'ATAT', 'ATAT'], 'score': [1, 1, 9]})
    filters = {'score': {'condition': 'score <= 4'}, 'ends': {'type': 'sequence_pattern', 'target': 'seq', 'exclude_patterns': ['^[GC]{2}']}}
    assert filter_table(raw, filters).index.tolist() == [1]


def test_post_processing_counts_reach_report():
    raw = pd.DataFrame({'probe_id': ['a', 'b', 'c'], 'seq': ['GCAT', 'ATAT', 'ATAT'], 'score': [1, 2, 1]})
    protocol = {'post_process': {'filters': {'ends': {'type': 'sequence_pattern', 'target': 'seq', 'exclude_patterns': ['^[GC]{2}']}}, 'sorts': {'is_ascending': ['score']}}}
    result = post_process(raw, protocol)
    assert result.probe_id.tolist() == ['c', 'b']
    stages = build_report_data(result, protocol, raw)['postProcessStages']
    assert stages == [
        {'stage': 'Filters', 'before': 3, 'after': 2, 'removed': 1},
        {'stage': 'Sort', 'before': 2, 'after': 2, 'removed': 0},
    ]
