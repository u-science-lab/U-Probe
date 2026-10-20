"""The attributes, filters and sorts generated for a protocol without them."""
import copy

import pandas as pd
import pytest

from uprobe.core.attributes import cal_target_fold_score
from uprobe.core.cli import _generate_default_attributes, _generate_default_post_process
from uprobe.core.process import post_process

PROBES = {'probe_1': {'template': '{part1}{part2}',
                      'parts': {'part1': {'expr': 'rc(target_region)'}, 'part2': {'expr': "'ACGTACGTAC'"}}}}


def protocol(source):
    return {'genome': 'g', 'targets': ['chr1:1-200' if source == 'genome' else 'GAPDH'], 'encoding': {},
            'extracts': {'target_region': {'source': source, 'length': 30, 'step': 1}},
            'probes': copy.deepcopy(PROBES)}


@pytest.mark.parametrize('source', ['mRNA', 'genome'])
def test_generated_sorts_rank_the_lowest_scores_first(source):
    # fold_score is -MFE (higher = more structure); fewer self-matches and
    # lower genome k-mer counts are better too, so every score sorts ascending.
    config = protocol(source)
    attributes = _generate_default_attributes(config)
    sorts = _generate_default_post_process(config, attributes)['sorts']
    scores = [name for name in attributes if name.endswith(('_fold', '_self_match', '_kmer_count'))]
    assert scores and set(scores) <= set(sorts['is_ascending'])
    assert sorts['is_descending'] == []


def test_generated_sort_puts_the_least_structured_window_first():
    hairpin = 'GGGGCCCCGAAAGGGGCCCC'
    loose = 'ATATATTTATATTAATATAT'
    assert cal_target_fold_score(hairpin) > cal_target_fold_score(loose)
    config = protocol('mRNA')
    attributes = {'target_fold': {'target': 'target_region', 'type': 'fold_score'}}
    config['post_process'] = {'sorts': _generate_default_post_process(config, attributes)['sorts']}
    table = pd.DataFrame({'target': ['G', 'G'], 'target_region': [hairpin, loose],
                          'target_fold': [cal_target_fold_score(s) for s in (hairpin, loose)]})
    assert post_process(table, config)['target_region'].tolist() == [loose, hairpin]
