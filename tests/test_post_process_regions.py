import pandas as pd
import pytest

from uprobe.core.process import post_process
from uprobe.core.process.equal_space import equal_space
from uprobe.core.process.regions import region_bounds


def test_region_bounds_accepts_every_source_format():
    df = pd.DataFrame({
        'sub_region': ['101-140', '1_40', '5-44', None],
        'start': [None, None, 7, 9],
        'end': [None, None, 46, 48],
    })
    bounds = region_bounds(df)
    assert bounds['start'].tolist() == [101, 1, 7, 9]
    assert bounds['end'].tolist() == [140, 40, 46, 48]


def test_region_bounds_rejects_unparseable_rows():
    with pytest.raises(ValueError):
        region_bounds(pd.DataFrame({'sub_region': ['chr1:1-40']}))


def test_equal_space_orders_exon_sub_regions_by_start():
    # Exon/UTR extraction writes "start_end"; int("2_41") parses as 241, so the old
    # split('-') ordering sorted by the concatenated digits (1_100 -> 1100 after 241).
    df = pd.DataFrame({
        'probe_id': ['a', 'b', 'c'],
        'target': ['G'] * 3,
        'sub_region': ['1_100', '2_41', '50_89'],
    })
    result = equal_space(df, {'number_desired': 2})
    assert result['probe_id'].tolist() == ['a', 'c']


def test_remove_overlap_orders_genome_windows_numerically():
    # As strings "100-139" < "41-80", so the non-overlapping 41-80 window was dropped.
    df = pd.DataFrame({
        'probe_id': ['a', 'b', 'c'],
        'target': ['chr1:0-200'] * 3,
        'sub_region': ['1-40', '41-80', '100-139'],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a', 'b', 'c']
    assert '_start' not in result.columns


def test_remove_overlap_genome_windows_with_empty_input():
    df = pd.DataFrame(columns=['probe_id', 'target', 'sub_region'])
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result.empty and list(result.columns) == ['probe_id', 'target', 'sub_region']


def test_remove_overlap_keeps_windows_of_each_direct_target():
    # Direct-sequence targets all have transcript_names 'custom'; windows of
    # different targets at the same positions must not remove each other.
    df = pd.DataFrame({
        'probe_id': ['a1', 'a2', 'b1', 'b2'],
        'target': ['A', 'A', 'B', 'B'],
        'transcript_names': ['custom'] * 4,
        'start': [0, 20, 0, 20],
        'end': [39, 59, 39, 59],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a1', 'b1']
    df['start'], df['end'] = [0, 40, 0, 40], [39, 79, 39, 79]
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a1', 'a2', 'b1', 'b2']


def test_remove_overlap_keeps_windows_of_each_exon():
    # Exon and UTR windows are numbered within their exon or UTR, so windows of
    # different exons of one transcript at the same positions do not overlap.
    df = pd.DataFrame({
        'probe_id': ['e1a', 'e1b', 'e1c', 'e2a', 'e2b', 'u1', 'c1'],
        'target': ['G'] * 6 + ['C'],
        'exon_name': ['1_101_200_+'] * 3 + ['1_301_400_+'] * 2 + [None, '.'],
        'utr_name': [None] * 5 + ['1_501_600_+', None],
        'transcript_names': [['T1']] * 5 + [['T1']] + ['custom'],
        'start': [1, 21, 41, 1, 41, 1, 1],
        'end': [40, 60, 80, 40, 80, 40, 40],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert sorted(result['probe_id']) == ['c1', 'e1a', 'e1c', 'e2a', 'e2b', 'u1']
