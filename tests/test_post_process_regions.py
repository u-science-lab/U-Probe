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
