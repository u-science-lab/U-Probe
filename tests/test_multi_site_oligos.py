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
    assert out.attrs['multi_site_oligos'] == [{'target': 'G', 'sequence': A, 'n_sites': 2, 'sites': ['1-20', '101-120']}]
    assert out['n_sites'].tolist() == [2, 1, 1] and out['oligo_sites'].tolist()[1:] == ['201-220', '401-420']


E = 'ACGT' * 5
MULTI = {'target': 'G', 'sequence': A, 'n_sites': 2, 'sites': ['1-20', '101-120']}


def summary_warning(final, df):
    from uprobe.core.report.run_summary import build_run_summary
    summary = build_run_summary({'name': 'x'}, {'final': final, 'requested_targets': ['G'], 'targets': df})
    return [w for w in summary['warnings'] if w['code'] == 'multi_site_oligos']


def test_post_process_reports_multi_site_oligos_once():
    # remove_overlap reports A, and equal_space keeps it: one entry, not one per step.
    df = windows([(A, 1), (A, 101), (B, 201)])
    out = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0},
                                             'equal_space': {'number_desired': 5}}})
    assert out['target_region'].tolist() == [A, B]
    assert out.attrs['multi_site_oligos'] == [MULTI]
    [warning] = summary_warning(out, df)
    assert warning['count'] == 1 and warning['oligos'] == [MULTI]


def test_only_oligos_of_the_final_table_are_reported():
    # remove_overlap keeps A (two sites); equal_space then picks the first and last candidates, not A.
    df = windows([(B, 1), (C, 201), (A, 401), (A, 601), (D, 801), (E, 1001)])
    out = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0},
                                             'equal_space': {'number_desired': 2}}})
    assert out['target_region'].tolist() == [B, E]
    assert 'multi_site_oligos' not in out.attrs and 'n_sites' not in out.columns
    assert summary_warning(out, df) == []
    # equal_space alone reports the repeated oligos it keeps, not every repeated candidate.
    df = windows([(B, 1), (A, 101), (A, 201), (C, 301), (D, 401), (D, 501), (E, 601)])
    out = equal_space(df, {'number_desired': 3})
    assert out['target_region'].tolist() == [B, C, E]
    assert 'multi_site_oligos' not in out.attrs and 'n_sites' not in out.columns
    out = equal_space(df, {'number_desired': 5})
    assert out['target_region'].tolist() == [B, A, C, D, E]
    assert out.attrs['multi_site_oligos'] == [{'target': 'G', 'sequence': A, 'n_sites': 2, 'sites': ['101-120', '201-220']},
                                              {'target': 'G', 'sequence': D, 'n_sites': 2, 'sites': ['401-420', '501-520']}]
    assert out['n_sites'].tolist() == [2 if seq in (A, D) else 1 for seq in out['target_region']]


def test_report_lists_every_kept_multi_site_oligo_with_all_its_sites():
    # Two targets; F repeats three times on G2, A twice on G; the summary counts the kept rows with n_sites > 1.
    F = 'GATC' * 5
    df = pd.concat([windows([(A, 1), (B, 50), (A, 101), (C, 201)]),
                    windows([(F, 1), (F, 301), (F, 601), (D, 900)], target='G2')], ignore_index=True)
    out = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0},
                                             'equal_space': {'number_desired': 3}}})
    expected = [MULTI, {'target': 'G2', 'sequence': F, 'n_sites': 3, 'sites': ['1-20', '301-320', '601-620']}]
    assert out.attrs['multi_site_oligos'] == expected
    assert int((out['n_sites'] > 1).sum()) == 2
    [warning] = summary_warning(out, df)
    assert warning['count'] == 2 and warning['oligos'] == expected


def test_report_names_an_oligo_once_when_the_table_holds_it_twice():
    from uprobe.core.process.oligo_sites import report_from_table
    table = pd.DataFrame({'target': ['G', 'G'], 'target_region': [A, A], 'n_sites': [2, 2],
                          'oligo_sites': ['1-20;101-120', '101-120;301-320']})
    assert report_from_table(table) == [{'target': 'G', 'sequence': A, 'n_sites': 3,
                                         'sites': ['1-20', '101-120', '301-320']}]


def test_identical_windows_at_one_site_are_one_probe_listing_all_transcripts():
    # transcripts: all without remove_overlap: one window shared by two transcripts is one candidate.
    df = windows([(A, 1), (A, 1), (B, 101), (C, 201)])
    df['transcript_name'] = ['T1', 'T2', 'T1', 'T1']
    df['n_trans'] = 1
    out = equal_space(df, {'number_desired': 5})
    assert out['target_region'].tolist() == [A, B, C]
    assert out['transcript_name'].tolist() == ['T1, T2', 'T1', 'T1'] and out['n_trans'].tolist() == [2, 1, 1]
    assert 'n_sites' not in out.columns and 'multi_site_oligos' not in out.attrs  # one site: not a repeat
    df['transcript_names'] = [['T1'], ['T2'], ['T1'], ['T1']]
    out = equal_space(df.drop(columns='transcript_name'), {'number_desired': 5})
    assert out['transcript_names'].tolist()[0] == ['T1', 'T2']


def naive_keep(spec, gap):
    """Reference: visit (sequence, start) in start order; an oligo is kept when none of its sites is
    within gap of a kept site (all its sites then count as kept)."""
    sites = {}
    for seq, start in spec:
        sites.setdefault(seq, set()).add((start, start + len(seq) - 1))
    kept_sites, kept, decided = [], [], set()
    for i in sorted(range(len(spec)), key=lambda i: spec[i][1]):
        seq = spec[i][0]
        if seq in decided:
            continue
        decided.add(seq)
        own = sorted(sites[seq])
        merged = []
        for s, e in own:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
            else:
                merged.append((s, e))
        if all(s > ke + gap or ks > e + gap for s, e in merged for ks, ke in kept_sites):
            kept.append(i)
            kept_sites += merged
    return kept


def test_remove_overlap_matches_a_naive_reference_on_random_repeats():
    import random
    rng = random.Random(3)
    for _ in range(200):
        pool = [''.join(rng.choice('ACGT') for _ in range(20)) for _ in range(5)]
        spec = []
        for _ in range(rng.randint(2, 40)):
            seq = rng.choice(pool) if rng.random() < 0.4 else ''.join(rng.choice('ACGT') for _ in range(20))
            spec.append((seq, rng.randint(1, 600)))
        spec = list(dict.fromkeys(spec))  # one row per (sequence, start)
        gap = rng.choice([0, 2, 7])
        out = remove_overlap(windows(spec), {}, gap)
        assert out['probe_id'].tolist() == [f'p{i}' for i in sorted(naive_keep(spec, gap), key=lambda i: spec[i][1])]
