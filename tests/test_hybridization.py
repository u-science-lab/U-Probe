import copy
import json
from io import BytesIO
from zipfile import ZipFile

import pandas as pd
import primer3
import pytest
import RNA

from uprobe.core.attributes import add_attributes, cal_temp, cal_target_fold_score
from uprobe.core.cli import _validate_and_normalize_protocol
from uprobe.core.hybridization import (
    FORMAMIDE_MOLAR_PER_PCT, describe_conditions, effective_conditions,
    validate_hybridization, validate_protocol_hybridization,
)
from uprobe.core.report.compact import build_report_data
from uprobe.core.report.excel import workbook_bytes

SEQS = ['ACGTTGCAGGCTAGCTAGGATCGATCGTAGC', 'GTATGACAATGAA', 'GGGGCCCCAAAATTTTGGGGCCCC']
HAIRPIN = 'GGGGAAAACCCCTTTTGGGGAAAACCCC'


def tm_table(protocol, seqs=SEQS):
    return add_attributes(pd.DataFrame({'s': seqs}), protocol, {})


def tm_protocol(block=None, **attr_fields):
    attr = {'target': 's', 'type': 'annealing_temperature', **attr_fields}
    protocol = {'attributes': {'tm': attr}}
    if block is not None:
        protocol['hybridization'] = block
    return protocol


def legacy_tm(seq):
    return primer3.calc_tm(seq, dv_conc=0, formamide_conc=0)


def test_without_block_tm_is_unchanged():
    expected = [legacy_tm(seq) for seq in SEQS]
    assert [cal_temp(seq) for seq in SEQS] == expected
    assert tm_table(tm_protocol())['tm'].tolist() == expected
    # An empty block or one without Tm fields keeps every legacy value.
    assert tm_table(tm_protocol({}))['tm'].tolist() == expected
    assert tm_table(tm_protocol({'temperature_c': 42}))['tm'].tolist() == expected
    explicit = {'na_mM': 50, 'mg_mM': 0, 'dntp_mM': 0.6, 'probe_nM': 50, 'formamide_pct': 0}
    assert tm_table(tm_protocol(explicit))['tm'].tolist() == expected


def test_without_block_fold_score_is_unchanged():
    protocol = {'attributes': {'fold': {'target': 's', 'type': 'fold_score'}}}
    expected = [round(-RNA.fold_compound(s.replace('T', 'U')).mfe()[1], 2) for s in SEQS]
    assert tm_table(protocol)['fold'].tolist() == expected
    protocol['hybridization'] = {'temperature_c': 37}
    assert tm_table(protocol)['fold'].tolist() == expected


def test_tm_falls_with_formamide_and_rises_with_salt():
    seq = SEQS[0]
    by_formamide = [tm_table(tm_protocol({'formamide_pct': pct}), [seq])['tm'].iloc[0]
                    for pct in (0, 10, 30, 50)]
    assert by_formamide == sorted(by_formamide, reverse=True)
    # Blake & Delcourt via primer3: 0.61-0.72 °C per % (v/v) formamide.
    assert 0.6 < (by_formamide[0] - by_formamide[1]) / 10 < 0.73
    assert by_formamide[1] == primer3.calc_tm(
        seq, dv_conc=0, formamide_conc=10 * FORMAMIDE_MOLAR_PER_PCT)
    by_salt = [tm_table(tm_protocol({'na_mM': na}), [seq])['tm'].iloc[0] for na in (50, 390, 1000)]
    assert by_salt == sorted(by_salt) and by_salt[0] < by_salt[1]
    assert tm_table(tm_protocol({'mg_mM': 5, 'dntp_mM': 0}), [seq])['tm'].iloc[0] > legacy_tm(seq)
    assert tm_table(tm_protocol({'probe_nM': 500}), [seq])['tm'].iloc[0] > legacy_tm(seq)


def test_attribute_override_wins_over_protocol_block():
    protocol = tm_protocol({'na_mM': 390, 'formamide_pct': 30})
    protocol['attributes']['tm_no_formamide'] = {
        'target': 's', 'type': 'annealing_temperature', 'hybridization': {'formamide_pct': 0}}
    result = tm_table(protocol)
    assert result['tm'].tolist() == [cal_temp(s, mv_conc=390, formamide_conc=30 * FORMAMIDE_MOLAR_PER_PCT) for s in SEQS]
    assert result['tm_no_formamide'].tolist() == [cal_temp(s, mv_conc=390, formamide_conc=0) for s in SEQS]
    assert (result['tm_no_formamide'] > result['tm']).all()


def test_fold_score_uses_hybridization_temperature():
    protocol = {'hybridization': {'temperature_c': 60},
                'attributes': {'fold_hyb': {'target': 's', 'type': 'fold_score'},
                               'fold_37': {'target': 's', 'type': 'fold_score',
                                           'hybridization': {'temperature_c': 37}}}}
    result = tm_table(protocol, [HAIRPIN])
    assert result['fold_37'].iloc[0] == round(cal_target_fold_score(HAIRPIN), 2)
    assert result['fold_hyb'].iloc[0] == round(cal_target_fold_score(HAIRPIN, 60), 2)
    assert result['fold_hyb'].iloc[0] < result['fold_37'].iloc[0]


@pytest.mark.parametrize('block', [
    [], 'formamide 30%', {'formamide': 30}, {'na_mm': 390},
    {'formamide_pct': '30'}, {'formamide_pct': None}, {'formamide_pct': True},
    {'formamide_pct': float('nan')}, {'formamide_pct': float('inf')},
    {'formamide_pct': -1}, {'formamide_pct': 101}, {'temperature_c': -5},
    {'temperature_c': 120}, {'na_mM': 0}, {'probe_nM': 0}, {'mg_mM': -1}, {'dntp_mM': -0.1},
])
def test_invalid_values_rejected(block):
    with pytest.raises(ValueError, match='hybridization'):
        validate_hybridization(block)
    with pytest.raises(ValueError, match=r'attributes\.tm\.hybridization'):
        validate_protocol_hybridization(tm_protocol({'na_mM': 390}, hybridization=block))
    with pytest.raises(ValueError):
        tm_table(tm_protocol(block))


def test_valid_block_accepts_ints_and_floats():
    assert validate_hybridization({'temperature_c': 37, 'na_mM': 390.0, 'formamide_pct': 100}) == {
        'temperature_c': 37.0, 'na_mM': 390.0, 'formamide_pct': 100.0}
    assert validate_hybridization(None) == {}


def test_override_only_for_tm_and_fold_attributes():
    protocol = {'attributes': {'gc': {'target': 's', 'type': 'gc_content', 'hybridization': {'na_mM': 390}}}}
    with pytest.raises(ValueError, match='only supported'):
        validate_protocol_hybridization(protocol)


def minimal_protocol():
    return {'name': 'hyb', 'genome': 'test', 'targets': [{'x': 'ACGT' * 20}],
            'extracts': {'target_region': {'source': 'exon', 'length': 30}},
            'encoding': {}, 'probes': {'probe_1': {'template': '{part1}',
                                                   'parts': {'part1': {'expr': 'rc(target_region)'}}}}}


def test_cli_normalization_validates_and_keeps_defaults():
    protocol = minimal_protocol()
    protocol['hybridization'] = {'temperature_c': 37, 'na_mM': 390, 'formamide_pct': 30}
    normalized = _validate_and_normalize_protocol(copy.deepcopy(protocol))
    assert normalized['hybridization'] == protocol['hybridization']
    assert normalized['attributes']['target_tm']['type'] == 'annealing_temperature'
    assert 'target_tm' in normalized['post_process']['filters']
    protocol['hybridization']['formamide_pct'] = 150
    with pytest.raises(ValueError, match='formamide_pct'):
        _validate_and_normalize_protocol(copy.deepcopy(protocol))
    without = _validate_and_normalize_protocol(minimal_protocol())
    assert 'hybridization' not in without


def test_conditions_are_described_and_recorded():
    protocol = tm_protocol({'temperature_c': 37, 'na_mM': 390, 'dntp_mM': 0, 'formamide_pct': 30})
    protocol['attributes']['tm_readout'] = {'target': 's', 'type': 'annealing_temperature',
                                            'hybridization': {'formamide_pct': 10}}
    protocol['attributes']['fold'] = {'target': 's', 'type': 'fold_score'}
    protocol['attributes']['gc'] = {'target': 's', 'type': 'gc_content'}
    conditions = effective_conditions(protocol)
    assert conditions == {
        'tm': {'na_mM': 390, 'mg_mM': 0, 'dntp_mM': 0, 'probe_nM': 50, 'formamide_pct': 30},
        'tm_readout': {'na_mM': 390, 'mg_mM': 0, 'dntp_mM': 0, 'probe_nM': 50, 'formamide_pct': 10},
        'fold': {'temperature_c': 37}}
    assert describe_conditions(protocol) == [
        'Tm (tm): 390 mM Na+, 0 mM Mg2+, 0 mM dNTP, 50 nM probe, 30% formamide',
        'Tm (tm_readout): 390 mM Na+, 0 mM Mg2+, 0 mM dNTP, 50 nM probe, 10% formamide',
        'fold (fold): 37 °C']
    df = tm_table(protocol)
    assert build_report_data(df, protocol)['hybridization'] == describe_conditions(protocol)
    with ZipFile(BytesIO(workbook_bytes(df, protocol))) as book:
        custom = book.read('docProps/custom.xml').decode()
    assert 'hybridization:tm_readout' in custom
    assert json.dumps(conditions['tm_readout']) in custom


def test_legacy_outputs_record_nothing_new():
    protocol = tm_protocol()
    assert describe_conditions(protocol) == [
        'Tm (tm): 50 mM Na+, 0 mM Mg2+, 0.6 mM dNTP, 50 nM probe, 0% formamide']
    df = tm_table(protocol)
    assert build_report_data(df, protocol)['hybridization'] == []
    with ZipFile(BytesIO(workbook_bytes(df, protocol))) as book:
        assert 'docProps/custom.xml' not in book.namelist()
