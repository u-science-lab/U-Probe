"""probe_index: a deterministic per-probe number for expressions (e.g. MERFISH readout rotation)."""
import collections
import copy
import random

import pytest

from uprobe.core.api import UProbeAPI
from uprobe.core.gen.expression import (
    CONTEXT_NAMES, ExpressionError, check_expression, evaluate_expression, probe_index_nodes,
)
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import validate_protocol
from uprobe.templates import load_template

READOUTS = ['ATCCTCCTTCAATACATCCC', 'ACACTACCACCATTTCCTAT', 'ACTCCACTACTACTCACTCT', 'ACCCTCTAACTTCCATCACA']


def test_probe_index_and_integer_remainder():
    assert 'probe_index' in CONTEXT_NAMES
    namespace = {'encoding': {'G': {'readouts': READOUTS}}, 'target': 'G', 'probe_index': 6}
    expr = "encoding[target]['readouts'][(probe_index + 1) % 4]"
    assert check_expression(expr, CONTEXT_NAMES) == []
    assert evaluate_expression(expr, namespace) == READOUTS[3]
    assert evaluate_expression('probe_index % 4 + -1 % 4', namespace) == 5


@pytest.mark.parametrize('expr, message', [
    ("'%s' % target", "'%' takes two integers"),          # no printf-style string formatting
    ("target % 2", "'%' takes two integers"),
    ('probe_index % 0', 'remainder by zero'),
])
def test_remainder_takes_integers_only(expr, message):
    assert check_expression(expr) == []  # the grammar allows %; the values are checked
    with pytest.raises(ExpressionError, match=message):
        evaluate_expression(expr, {'target': 'G', 'probe_index': 3})


@pytest.mark.parametrize('expr', ['target_region * 2', 'probe_index // 2', 'probe_index ** 2'])
def test_other_operators_stay_rejected(expr):
    assert "is not allowed (only '+' and '%')" in check_expression(expr)[0]


def test_nodes_that_depend_on_probe_index():
    probes = load_template('merfish')['probes']
    assert probe_index_nodes(probes) == {'encoding_probe', 'encoding_probe.readout_1', 'encoding_probe.readout_2',
                                         'encoding_probe.readout_3'}
    probes['copy'] = {'template': '{a}', 'parts': {'a': {'expr': 'encoding_probe.readout_2'}}}
    assert {'copy', 'copy.a'} <= probe_index_nodes(probes)
    assert probe_index_nodes({'p': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}}) == set()


def merfish_api(tmp_path, sequence):
    """The MERFISH template on one direct sequence, with the attributes that need no genome."""
    protocol = copy.deepcopy(load_template('merfish'))
    protocol['targets'] = [{'TEST': sequence}]
    protocol['encoding'] = {'TEST': {'readouts': READOUTS}}
    del protocol['attributes']['target_mapped_genes']
    del protocol['post_process']['filters']['target_mapped_genes']
    protocol['summary'] = {'report_name': '', 'attributes': ['target_gc']}
    api = UProbeAPI(protocol, {}, tmp_path, require_genome=False)
    api.build_genome_index = lambda **kwargs: None
    return api


def test_merfish_rotates_three_of_four_readouts_over_the_final_probes(tmp_path):
    rng = random.Random(4)
    sequence = ''.join(rng.choice('ACGT') for _ in range(5000))
    api = merfish_api(tmp_path, sequence)
    final = api.run_workflow().sort_values('probe_index').reset_index(drop=True)
    count = len(final)
    assert count > 40
    # Numbered 0..n-1 along the RNA, after post-processing chose the probes.
    assert final['probe_index'].tolist() == list(range(count))
    assert final['start'].is_monotonic_increasing
    sites = [reverse_complement(r) for r in READOUTS]
    for k, row in final.iterrows():
        assert [row[f'encoding_probe.readout_{i}'] for i in (1, 2, 3)] == [sites[(k + i) % 4] for i in (1, 2, 3)]
        assert row['encoding_probe'] == ('N' * 20 + 'A' + sites[(k + 1) % 4] + reverse_complement(row['target_region'])
                                         + 'A' + sites[(k + 2) % 4] + sites[(k + 3) % 4] + 'N' * 20)
    # Every on-bit appears, each on 3 of every 4 probes; the 4 combinations are balanced.
    left_out = collections.Counter(k % 4 for k in range(count))
    carried = collections.Counter(site for i in (1, 2, 3) for site in final[f'encoding_probe.readout_{i}'])
    assert set(carried) == set(sites)
    assert max(left_out.values()) - min(left_out.values()) <= 1
    assert all(carried[site] == count - left_out[j] for j, site in enumerate(sites))
    tsv = (tmp_path / api.output_files['tsv']).read_text().splitlines()
    assert tsv[0].split('\t').count('probe_index') == 1


def test_validate_checks_every_probe_number(tmp_path):
    protocol = copy.deepcopy(load_template('merfish'))
    assert validate_protocol(protocol) == []
    # Three readouts cannot serve probe_index % 4.
    protocol['encoding']['GAPDH']['readouts'] = READOUTS[:3]
    errors = [f for f in validate_protocol(protocol) if f.severity == 'error']
    assert errors and 'index 3 out of range' in errors[0].message
    # An attribute on a rotated part is recomputed at the end; filters saw provisional numbers.
    protocol['encoding']['GAPDH']['readouts'] = READOUTS
    protocol['attributes']['probe_gc'] = {'target': 'encoding_probe', 'type': 'gc_content'}
    warning = [f for f in validate_protocol(protocol) if f.path == 'attributes.probe_gc.target']
    assert warning and 'depends on probe_index' in warning[0].message
