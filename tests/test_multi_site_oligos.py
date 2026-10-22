import pandas as pd
import pytest

from uprobe.core.process import post_process, remove_overlap
from uprobe.core.process.equal_space import equal_space

A, B, C, D = 'A' * 20, 'C' * 20, 'G' * 20, 'T' * 20


def windows(spec, target='G'):
    """spec: (sequence, start); 20-nt windows on the target's merged model."""
    return pd.DataFrame([{'probe_id': f'p{i}', 'target': target, 'target_region': seq,
                          'target_start': start, 'target_end': start + len(seq) - 1,
                          'transcript_name': 'T1'} for i, (seq, start) in enumerate(spec)])


def test_unique_oligos_are_unchanged():
    df = windows([(A, 1), (B, 15), (C, 41), (D, 61)])
    out = remove_overlap(df, {}, 0)
    assert out['probe_id'].tolist() == ['p0', 'p2', 'p3']
    assert 'n_sites' not in out.columns and 'multi_site_oligos' not in out.attrs


def test_identical_oligo_is_counted_once_and_occupies_all_its_sites():
    # A repeats at 1 and 101; B overlaps the second copy of A and is first in line after the first copy.
    df = windows([(A, 1), (B, 95), (A, 101), (C, 201)])
    out = remove_overlap(df, {}, 0)
    # Per position, B (95-114) would be kept after A (1-20); A's second site (101-120) now blocks it.
    assert out['probe_id'].tolist() == ['p0', 'p3']
    assert out['n_sites'].tolist() == [2, 1]
    assert out.loc[0, 'oligo_sites'] == '1-20;101-120'
    report = out.attrs['multi_site_oligos']
    assert report == [{'target': 'G', 'sequence': A, 'n_sites': 2, 'sites': ['1-20', '101-120']}]


def test_location_interval_applies_to_every_site():
    df = windows([(A, 1), (B, 125), (A, 101), (C, 300)])  # B starts after a 4-nt gap beyond the second A site
    assert remove_overlap(df, {}, 3)['probe_id'].tolist() == ['p0', 'p1', 'p3']
    assert remove_overlap(df, {}, 5)['probe_id'].tolist() == ['p0', 'p3']


def test_oligo_refused_when_any_site_is_blocked():
    # C (50-69) is kept first; A (1-20) would be fine alone but its second site (60-79) overlaps C.
    df = windows([(A, 1), (C, 50), (A, 60)])
    # By start A is visited first and takes both sites, so C is the one refused.
    assert remove_overlap(df, {}, 0)['target_region'].tolist() == [A]
    # C (1-20) is first; A's site at 25 is free but its other site (10-29) overlaps C: A is refused whole.
    assert remove_overlap(windows([(C, 1), (A, 25)]), {}, 0)['target_region'].tolist() == [C, A]
    assert remove_overlap(windows([(C, 1), (A, 25), (A, 10)]), {}, 0)['target_region'].tolist() == [C]


def test_overlapping_sites_of_one_oligo_do_not_block_each_other():
    # tandem period 10 < length 20: sites 1-20 and 11-30 are one oligo
    df = windows([(A, 1), (A, 11), (B, 31)])
    out = remove_overlap(df, {}, 0)
    assert out['target_region'].tolist() == [A, B]


def test_layouts_of_a_window_are_different_oligos():
    df = windows([(A, 1), (A, 1)])
    df['layout_id'] = ['left:10', 'left:12']
    assert len(remove_overlap(df, {}, 0)) == 1  # same window: overlap keeps one
    assert len(equal_space(df, {'number_desired': 5})) == 2  # distinct probes: both kept


def test_equal_space_counts_distinct_oligos():
    df = windows([(A, 1), (A, 101), (B, 201), (C, 301), (D, 401)])
    out = equal_space(df, {'number_desired': 3})
    assert len(set(out['target_region'])) == 3 and out['target_region'].tolist()[0] == A
    assert out.attrs['multi_site_oligos'][0]['n_sites'] == 2


def test_post_process_reports_multi_site_oligos():
    df = windows([(A, 1), (A, 101), (B, 201)])
    out = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0},
                                             'equal_space': {'number_desired': 5}}})
    assert len(out) == 2
    assert out.attrs['multi_site_oligos'][0]['sequence'] == A


def test_run_summary_warns_about_multi_site_oligos():
    from uprobe.core.report.run_summary import build_run_summary
    df = windows([(A, 1), (A, 101), (B, 201)])
    final = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    summary = build_run_summary({'name': 'x'}, {'final': final, 'requested_targets': ['G'], 'targets': df})
    warning = next(w for w in summary['warnings'] if w['code'] == 'multi_site_oligos')
    assert warning['count'] == 1 and warning['oligos'][0]['sites'] == ['1-20', '101-120']
