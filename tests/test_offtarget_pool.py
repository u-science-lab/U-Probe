"""post_process.off_target_pool: at most k probes of a target's set per off-target gene."""
import copy
from pathlib import Path

import pandas as pd
import pytest
import yaml

from uprobe.core.cli import _validate_and_normalize_protocol
from uprobe.core.process import STEP_ORDER, post_process
from uprobe.core.process.off_target_pool import limit_off_target_pool, pool_settings, rank_order
from uprobe.core.report.run_summary import build_run_summary
from uprobe.core.validate import validate_protocol

DATA = Path(__file__).parent / 'data'


def candidates():
    # Ranked by score (higher is better): p1 > p2 > ... ; p1-p4 hit ACTG1.
    return pd.DataFrame({
        'probe_id': ['p6', 'p1', 'p3', 'p2', 'p5', 'p4', 'q1', 'q2', 'q3'],
        'target': ['A'] * 6 + ['B'] * 3,
        'score': [4, 9, 7, 8, 5, 6, 3, 2, 1],
        'hits_off_targets': ['', 'ACTG1:0', 'ACTG1:1, POTEF:2', 'ACTG1:2', 'POTEF:0', 'ACTG1:0',
                             'ACTG1:0', 'ACTG1:1', 'ACTG1:2'],
    })


SORTS = {'is_descending': ['score']}


def test_keeps_the_best_ranked_probes_per_gene():
    kept, report = limit_off_target_pool(candidates(), {'attribute': 'hits', 'max_probes_per_gene': 2}, SORTS)
    # A: p1, p2 take the two ACTG1 slots; p3 and p4 are dropped. B is limited on its own.
    assert kept.probe_id.tolist() == ['p6', 'p1', 'p2', 'p5', 'q1', 'q2']  # original order
    assert report == {'attributes': ['hits'], 'max_probes_per_gene': 2, 'removed': 3,
                      'genes_at_limit': {'ACTG1': 3}}


def test_limit_zero_and_one():
    kept, _ = limit_off_target_pool(candidates(), {'attribute': 'hits', 'max_probes_per_gene': 0}, SORTS)
    assert kept.probe_id.tolist() == ['p6']
    kept, report = limit_off_target_pool(candidates(), {'attribute': 'hits', 'max_probes_per_gene': 1}, SORTS)
    # p3 is dropped for ACTG1, so POTEF's slot goes to p5.
    assert kept.probe_id.tolist() == ['p6', 'p1', 'p5', 'q1']
    assert report['genes_at_limit'] == {'ACTG1': 5}


def test_without_sorts_earlier_rows_win():
    kept, _ = limit_off_target_pool(candidates(), {'attribute': 'hits', 'max_probes_per_gene': 1})
    assert kept.probe_id.tolist() == ['p6', 'p1', 'p5', 'q1']


def test_several_attributes_share_the_limit():
    # p6's readout hits POTEF, whose one slot p5 (better ranked) already holds.
    df = candidates().assign(readout_off_targets=['POTEF', '', '', '', '', '', '', '', ''])
    kept, _ = limit_off_target_pool(df, {'attribute': ['hits', 'readout'], 'max_probes_per_gene': 1}, SORTS)
    assert kept.probe_id.tolist() == ['p1', 'p5', 'q1']


def test_rank_order_follows_sorts():
    df = pd.DataFrame({'a': [1, 1, 0, None], 'b': [1, 2, 3, 4]})
    assert rank_order(df, {'is_ascending': ['a'], 'is_descending': ['b']}).tolist() == [2, 1, 0, 3]
    assert rank_order(df, {'is_ascending': ['missing']}).tolist() == [0, 1, 2, 3]


@pytest.mark.parametrize('config,message', [
    ({'attribute': 'hits'}, 'max_probes_per_gene must be an integer >= 0'),
    ({'attribute': 'hits', 'max_probes_per_gene': -1}, 'max_probes_per_gene must be'),
    ({'attribute': 'hits', 'max_probes_per_gene': 1.5}, 'max_probes_per_gene must be'),
    ({'max_probes_per_gene': 2}, 'attribute must name'),
    ({'attribute': 'hits', 'max_probes_per_gene': 2, 'k': 2}, 'unknown option'),
    ([], 'must be a mapping'),
])
def test_settings_are_checked(config, message):
    with pytest.raises(ValueError, match=message):
        pool_settings(config)


def test_missing_off_target_column():
    with pytest.raises(ValueError, match='needs other_off_targets'):
        limit_off_target_pool(candidates(), {'attribute': 'other', 'max_probes_per_gene': 2})


# -- composition with the other steps ------------------------------------------

def windows():
    """12 non-overlapping windows along target A; the even ones hit ACTG1."""
    return pd.DataFrame({
        'probe_id': [f'w{i}' for i in range(12)], 'target': ['A'] * 12,
        'start': [1 + 30 * i for i in range(12)], 'end': [20 + 30 * i for i in range(12)],
        'sub_region': [f'{1 + 30 * i}_{20 + 30 * i}' for i in range(12)],
        'exon_name': ['1_1_400_+'] * 12, 'transcript_names': ['T1'] * 12,
        'score': list(range(12)),
        'hits': [2 if i % 2 == 0 else 1 for i in range(12)],
        'hits_off_targets': ['ACTG1:1' if i % 2 == 0 else '' for i in range(12)],
    })


def test_pool_limit_runs_before_equal_space():
    assert STEP_ORDER.index('remove_overlap') < STEP_ORDER.index('off_target_pool') < STEP_ORDER.index('equal_space')
    config = {'post_process': {
        'equal_space': {'number_desired': 6},  # listed first: the protocol order does not matter
        'filters': {'hits': {'condition': 'hits <= 2'}},
        'off_target_pool': {'attribute': 'hits', 'max_probes_per_gene': 2},
        'remove_overlap': {'location_interval': 2},
        'sorts': {'is_descending': ['score']},
    }}
    result = post_process(windows(), config)
    # The pool keeps w10 and w8 (best-ranked ACTG1 hits) and the 6 clean windows;
    # equal_space then chooses 6 of these 8, so ACTG1 stays within the limit.
    assert len(result) == 6
    assert result.hits_off_targets.eq('ACTG1:1').sum() <= 2
    assert set(result.probe_id) <= {'w1', 'w3', 'w5', 'w7', 'w8', 'w9', 'w10', 'w11'}
    assert [s['stage'] for s in result.attrs['post_process_stages']] == [
        'Filters', 'Remove overlap', 'Off-target pool limit', 'Equal spacing', 'Sort']
    assert result.attrs['off_target_pool']['removed'] == 4
    # equal_space alone would have picked windows evenly, ACTG1 hits included.
    alone = post_process(windows(), {'post_process': {'equal_space': {'number_desired': 6}}})
    assert alone.hits_off_targets.eq('ACTG1:1').sum() > 2


def test_summary_reports_the_pool():
    result = post_process(windows(), {'post_process': {
        'off_target_pool': {'attribute': 'hits', 'max_probes_per_gene': 1}, 'sorts': {'is_descending': ['score']}}})
    report = {'attributes': {'hits': {'aligner': 'bowtie', 'max_mismatches': 2, 'strand': 'sense',
                                      'off_target_genes': 1}}}
    summary = build_run_summary({}, {'final': result, 'raw': windows(), 'off_targets': report})
    assert summary['off_targets']['pool'] == {'attributes': ['hits'], 'max_probes_per_gene': 1, 'removed': 5,
                                              'genes_at_limit': {'ACTG1': 5}}
    assert summary['stages'][0]['stage'] == 'Off-target pool limit'


# -- uprobe run / uprobe validate --------------------------------------------

@pytest.fixture
def rna():
    protocol = yaml.safe_load((DATA / 'RNA_format.yaml').read_text(encoding='utf-8'))
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    protocol['post_process']['off_target_pool'] = {'attribute': 'target_mappedGenes', 'max_probes_per_gene': 2}
    return protocol


def test_validate_accepts_and_checks_the_step(rna):
    assert validate_protocol(rna) == []
    pool = rna['post_process']['off_target_pool']
    pool.update(attribute='target_tm', max_probes_per_gene=-2, maximum=3)
    errors = {f.path: f for f in validate_protocol(rna) if f.severity == 'error'}
    assert 'not a mapped_genes attribute' in errors['post_process.off_target_pool.attribute'].message
    assert 'target_mappedGenes' in errors['post_process.off_target_pool.attribute'].suggestion
    assert 'max_probes_per_gene must be' in errors['post_process.off_target_pool.max_probes_per_gene'].message
    assert 'unknown option' in errors['post_process.off_target_pool.maximum'].message


def test_validate_warns_when_a_filter_leaves_nothing_to_limit(rna):
    rna['post_process']['filters']['target_mappedGenes'] = {'condition': 'target_mappedGenes <= 1'}
    warning = [f for f in validate_protocol(rna) if f.path == 'post_process.off_target_pool']
    assert warning and 'removes nothing' in warning[0].message
    rna['post_process']['filters']['target_mappedGenes'] = {'condition': 'target_mappedGenes <= 4'}
    assert validate_protocol(rna) == []


def test_run_rejects_a_pool_without_mapped_genes(rna):
    _validate_and_normalize_protocol(copy.deepcopy(rna))
    rna['post_process']['off_target_pool']['attribute'] = 'target_tm'
    with pytest.raises(ValueError, match="'target_tm' is not a mapped_genes attribute"):
        _validate_and_normalize_protocol(rna)
