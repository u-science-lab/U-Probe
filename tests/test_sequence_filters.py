import pandas as pd
import pytest
from uprobe.core.attributes import add_attributes
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


def test_gc_window_applies_to_the_exact_fraction(tmp_path, monkeypatch):
    # 19 of 30 nt G or C is 63.3% GC, outside the 43-63% MERFISH window; rounded
    # to 0.63 it passed, and 50 of 214 MERFISH probes of GAPDH, TFRC and EEF2 had it.
    monkeypatch.chdir(tmp_path)
    raw = pd.DataFrame({'probe_id': ['gc19', 'gc18', 'gc13'],
                        'target_region': ['G' * 19 + 'A' * 11, 'G' * 18 + 'A' * 12, 'G' * 13 + 'A' * 17]})
    protocol = {'attributes': {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
                'post_process': {'filters': {'target_gc': {'condition': 'target_gc >= 0.43 & target_gc <= 0.63'}}}}
    df = add_attributes(raw, protocol, genome={})
    assert df.target_gc.tolist() == [19 / 30, 18 / 30, 13 / 30]
    assert post_process(df, protocol).probe_id.tolist() == ['gc18', 'gc13']
