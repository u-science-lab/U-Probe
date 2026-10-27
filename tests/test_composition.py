import math
import random

import pandas as pd
import pytest

from uprobe.core.attributes import add_attributes
from uprobe.core.composition import base_fraction, composition_settings, max_homopolymer, max_window_count
from uprobe.core.process import post_process
from uprobe.core.validate import validate_protocol

# The five rules of Xu et al. 2009 (PNAS 106:2289), as Oligostan applies them to smiFISH probes.
XU2009_ATTRIBUTES = {
    'a_frac': {'target': 'seq', 'type': 'base_fraction', 'bases': 'A'},
    'a_run': {'target': 'seq', 'type': 'max_homopolymer', 'bases': 'A'},
    'c_frac': {'target': 'seq', 'type': 'base_fraction', 'bases': 'C'},
    'c_run': {'target': 'seq', 'type': 'max_homopolymer', 'bases': 'C'},
    'c_5prime': {'target': 'seq', 'type': 'max_window_count', 'bases': 'C', 'window': 6, 'positions': [1, 12]},
}
XU2009_FILTERS = {
    'rule1': {'condition': 'a_frac < 0.28'},
    'rule2': {'condition': 'a_run <= 3'},
    'rule3': {'condition': 'c_frac >= 0.22 & c_frac <= 0.28'},
    'rule4': {'condition': 'c_run <= 3'},
    'rule5': {'condition': 'c_5prime <= 3'},
}


def xu2009_by_hand(seq):
    n = len(seq)
    return (seq.count('A') / n < 0.28 and 'AAAA' not in seq and 0.22 <= seq.count('C') / n <= 0.28
            and 'CCCC' not in seq and all(seq[i:i + 6].count('C') < 4 for i in range(0, 7)))


def test_base_fraction():
    assert base_fraction('AACG', 'A') == 0.5
    assert base_fraction('AACG', 'GC') == 0.5
    assert base_fraction('aacg', 'C') == 0.25  # case-insensitive
    assert base_fraction('ACGU', 'T') == 0.25  # U counts as T
    assert base_fraction('AAAACCCC', 'A', (5, 8)) == 0
    assert base_fraction('AAAACCCC', 'C', (-2, -1)) == 1
    assert math.isnan(base_fraction(None, 'A'))
    assert math.isnan(base_fraction('', 'A'))


def test_max_homopolymer():
    assert max_homopolymer('ACCCAAG') == 3
    assert max_homopolymer('ACCCAAG', 'A') == 2
    assert max_homopolymer('ACCCAAG', 'G') == 1
    assert max_homopolymer('ACCCAAG', 'T') == 0
    assert max_homopolymer('AAAACCAA', 'A', (3, 8)) == 2  # the run is clipped at the part's start
    assert math.isnan(max_homopolymer(float('nan')))


def test_max_window_count():
    #                 positions 1-12: CACCTC GACCAC ...
    seq = 'CACCTCGACCACCCCC'
    assert max_window_count(seq, 'C', 6, (1, 12)) == 4  # CACCTC
    assert max_window_count('ACACACGTTTTTCCCC', 'C', 6, (1, 12)) == 3
    assert max_window_count('ACACACGTTTTTCCCC', 'C', 6) == 4  # beyond 12 when no positions
    assert max_window_count('ACACACGTTTTTCCCC', 'C', 6, (1, 13)) == 3  # window 8-13 holds 1 C
    assert max_window_count('ACACACGTTTTTCCCC', 'C', 6, (1, 14)) == 3  # window 9-14 holds 2
    assert max_window_count('CCC', 'C', 6) == 3  # a part shorter than the window is counted whole
    assert max_window_count('GCGCAT', 'GC', 3) == 3


@pytest.mark.parametrize('attr, message', [
    ({'type': 'base_fraction'}, 'bases is required'),
    ({'type': 'base_fraction', 'bases': 'AN'}, 'only A, C, G, T'),
    ({'type': 'max_window_count', 'bases': 'C'}, 'window must be a positive integer'),
    ({'type': 'max_window_count', 'bases': 'C', 'window': 0}, 'window must be a positive integer'),
    ({'type': 'max_window_count', 'bases': 'C', 'window': 6, 'positions': [0, 12]}, 'two non-zero integers'),
    ({'type': 'max_window_count', 'bases': 'C', 'window': 6, 'positions': [12, 1]}, 'is after end'),
    ({'type': 'max_homopolymer', 'positions': 12}, 'two non-zero integers'),
])
def test_invalid_settings(attr, message):
    with pytest.raises(ValueError, match=message):
        composition_settings(attr['type'], attr, 'attributes.x')


def test_default_homopolymer_bases_and_list_bases():
    assert composition_settings('max_homopolymer', {}, 'x')['bases'] == 'ACGT'
    assert composition_settings('base_fraction', {'bases': ['g', 'C']}, 'x')['bases'] == 'GC'


def test_xu2009_rules_match_a_direct_implementation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rng = random.Random(7)
    seqs = [''.join(rng.choice('ACGT') for _ in range(rng.randint(26, 32))) for _ in range(3000)]
    # Sequences built to sit on the rule boundaries.
    seqs += ['CACCTC' + 'GTGTGT' + 'ACAGCTAGTCAG',  # rule 5: 4 C in positions 1-6
             'GTGTGT' + 'CACCTC' + 'ACAGTGAGTGAG',  # rule 5: 4 C in positions 7-12
             'GTGTGTGTGTGC' + 'CACCTCACAGTG',  # 4 C in a 6-nt window from 13 on: allowed
             'ACGT' * 7]  # A exactly 25 %, C 25 %
    df = add_attributes(pd.DataFrame({'probe_id': range(len(seqs)), 'seq': seqs}),
                        {'attributes': XU2009_ATTRIBUTES}, genome={})
    kept = post_process(df, {'post_process': {'filters': XU2009_FILTERS}})
    expected = [i for i, seq in enumerate(seqs) if xu2009_by_hand(seq)]
    assert kept.probe_id.tolist() == expected
    assert 0 < len(expected) < len(seqs)
    n = len(seqs)
    assert n - 4 not in expected and n - 3 not in expected  # rule 5 on the first 12 nt
    assert n - 1 in expected


def test_filters_run_before_remove_overlap():
    # Three overlapping windows of one target: the first fails the A rule. Filters run first, so
    # remove_overlap packs from the passing windows instead of keeping the failing one.
    seqs = ['AAAAAAAAAA', 'CGCGCGCGCG', 'GCGCGCGCGC', 'CGCGCGCGCG']
    df = pd.DataFrame({'probe_id': list('abcd'), 'target': 'X', 'seq': seqs,
                       'target_start': [0, 2, 8, 13], 'target_end': [10, 12, 18, 23]})
    df = add_attributes(df, {'attributes': {'a_frac': XU2009_ATTRIBUTES['a_frac']}}, genome={})
    protocol = {'post_process': {'remove_overlap': {'location_interval': 0},
                                 'filters': {'rule1': XU2009_FILTERS['rule1']}}}  # YAML order does not matter
    assert post_process(df, protocol).probe_id.tolist() == ['b', 'd']


def minimal(**attrs):
    return {'name': 'comp', 'genome': 'test', 'targets': [{'x': 'ACGT' * 20}],
            'extracts': {'target_region': {'source': 'exon', 'length': 30}}, 'encoding': {},
            'probes': {'probe_1': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
            'attributes': attrs,
            'post_process': {'filters': {name: {'condition': f'{name} <= 3'} for name in attrs}}}


def test_validate_accepts_the_attributes_and_reports_bad_options():
    attrs = {name: {**attr, 'target': 'probe_1.part1'} for name, attr in XU2009_ATTRIBUTES.items()}
    assert not [f for f in validate_protocol(minimal(**attrs)) if f.severity == 'error']
    errors = [f for f in validate_protocol(minimal(c=dict(attrs['c_5prime'], window=None)))
              if f.severity == 'error']
    assert [f.path for f in errors] == ['attributes.c'] and 'window' in errors[0].message
    warnings = [f for f in validate_protocol(minimal(gc={'target': 'probe_1', 'type': 'gc_content', 'bases': 'A'}))
                if f.path == 'attributes.gc.bases']
    assert len(warnings) == 1 and warnings[0].severity == 'warning'
