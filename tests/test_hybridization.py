import copy
import json
import math
from io import BytesIO
from zipfile import ZipFile

import pandas as pd
import primer3
import pytest
import RNA

from uprobe.core.attributes import add_attributes, cal_temp, cal_target_fold_score
from uprobe.core.cli import _validate_and_normalize_protocol
from uprobe.core.hybridization import (
    FORMAMIDE_MOLAR_PER_PCT, default_tm_window, describe_conditions, effective_conditions,
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


def tm_filters(protocol):
    """Auto-generated Tm filter conditions of a protocol without post_process."""
    filters = _validate_and_normalize_protocol(copy.deepcopy(protocol))['post_process']['filters']
    return {name: spec['condition'] for name, spec in filters.items() if name.endswith('_tm')}


def windows(filters):
    return {name: condition.replace(f'{name} ', '') for name, condition in filters.items()}


@pytest.mark.parametrize('source', ['exon', 'genome'])
def test_default_tm_filters_follow_the_hybridization_temperature(source, monkeypatch):
    from unittest import mock
    from uprobe.core import cli as cli_module
    protocol = minimal_protocol()
    if source == 'genome':
        protocol['extracts']['target_region']['source'] = 'genome'
        protocol['targets'] = ['1:1-200']
    # RNA: target_tm, probe_1_tm and probe_1_part1_tm; DNA: only the part.
    names = ['target_tm', 'probe_1_tm', 'probe_1_part1_tm'] if source == 'exon' else ['probe_1_part1_tm']
    # Without declared conditions, or with an empty block: the legacy 50-90 °C, unchanged.
    legacy = {name: f'{name} >= 50 & {name} <= 90' for name in names}
    assert tm_filters(protocol) == legacy
    assert tm_filters({**protocol, 'hybridization': {}}) == legacy
    # Oligopaint-like: 42 °C, 2x SSC, 50% formamide -> 42-82 °C.
    log = mock.Mock()
    monkeypatch.setattr(cli_module, 'log', log)
    oligopaint = {**protocol, 'hybridization': {'temperature_c': 42, 'na_mM': 390, 'formamide_pct': 50}}
    assert windows(tm_filters(oligopaint)) == {name: '>= 42 & <= 82' for name in names}
    assert any('Tm filter 42-82 °C (hybridization.temperature_c' in str(call) for call in log.info.call_args_list)
    assert not log.warning.called
    assert windows(tm_filters({**protocol, 'hybridization': {'temperature_c': 37.5}})) == {
        name: '>= 37.5 & <= 77.5' for name in names}
    # Salt and formamide without a temperature: the default 37 °C, with a warning.
    assert windows(tm_filters({**protocol, 'hybridization': {'na_mM': 390, 'formamide_pct': 50}})) == {
        name: '>= 37 & <= 77' for name in names}
    assert 'temperature_c' in str(log.warning.call_args)


def test_default_tm_window_with_attribute_overrides():
    block = {'temperature_c': 37, 'na_mM': 390, 'formamide_pct': 30}
    attributes = {
        'plain': {},
        # The block's own values: still the hybridization conditions.
        'same': {'na_mM': 390},
        # Design conditions of a published window (MERFISH template): not at the hybridization temperature.
        'design': {'na_mM': 300, 'probe_nM': 5, 'formamide_pct': 0},
        'own': {'temperature_c': 45, 'formamide_pct': 10},
    }
    protocol = tm_protocol(block)
    for name, override in attributes.items():
        attr = {'target': 's', 'type': 'annealing_temperature'}
        protocol['attributes'][name] = {**attr, 'hybridization': override} if override else attr
    assert {name: default_tm_window(protocol, name) for name in attributes} == {
        'plain': (37, 77, 'protocol'), 'same': (37, 77, 'protocol'),
        'design': (50, 90, None), 'own': (45, 85, 'attribute')}
    # Without a protocol block only an attribute's own temperature moves its window.
    assert default_tm_window(tm_protocol(None, hybridization={'formamide_pct': 50}), 'tm') == (50, 90, None)
    assert default_tm_window(tm_protocol(None, hybridization={'temperature_c': 30}), 'tm') == (30, 70, 'attribute')
    assert default_tm_window(tm_protocol({'formamide_pct': 50}), 'tm') == (37, 77, 'default')
    # User-written attributes with post_process left out get one window each.
    rna = {**minimal_protocol(), 'hybridization': block, 'attributes': {
        'target_tm': {'target': 'target_region', 'type': 'annealing_temperature'},
        'design_tm': {'target': 'target_region', 'type': 'annealing_temperature',
                      'hybridization': attributes['design']}}}
    assert tm_filters(rna) == {'target_tm': 'target_tm >= 37 & target_tm <= 77',
                               'design_tm': 'design_tm >= 50 & design_tm <= 90'}


def test_validate_warns_when_generated_tm_filters_assume_37_c():
    from uprobe.core.validate import validate_protocol
    protocol = {**minimal_protocol(), 'hybridization': {'na_mM': 390, 'formamide_pct': 50}}
    warnings = [f for f in validate_protocol(protocol) if f.path == 'hybridization']
    assert len(warnings) == 1 and '37 °C' in warnings[0].message and 'temperature_c' in warnings[0].suggestion
    protocol['hybridization']['temperature_c'] = 42
    assert not [f for f in validate_protocol(protocol) if f.path == 'hybridization']


# Moldia/multi_padlock_design lib/screenseq.py calculatetm (HybISS design code), run on the
# published mouse Pvalb targets of Gyllborg et al. 2020 (Supplementary Table 1).
MPD_TM = {'CCCATCTCAGTCCTTCTCGCGGCCCTCCTG': 70.50057357343485,
          'CTCTATGGCCCTCGGATGACGCCATTCTTC': 67.16298334058791}
MPD_CONDITIONS = {'na_mM': 75, 'mg_mM': 10, 'probe_nM': 100, 'formamide_pct': 20}


def test_multi_padlock_design_tm_reproduces_the_hybiss_code():
    from uprobe.core.attributes import cal_temp_multi_padlock_design
    for seq, tm in MPD_TM.items():
        assert cal_temp_multi_padlock_design(seq, **MPD_CONDITIONS) == pytest.approx(tm, abs=1e-9)
        assert cal_temp_multi_padlock_design(seq.lower().replace('t', 'u'), **MPD_CONDITIONS) == pytest.approx(tm)
    assert math.isnan(cal_temp_multi_padlock_design('ACGTNACGT', **MPD_CONDITIONS))
    # The attribute takes the declared conditions; dNTP is not part of the formula.
    protocol = tm_protocol({'temperature_c': 45, 'dntp_mM': 1}, tm_method='multi_padlock_design',
                           hybridization=MPD_CONDITIONS)
    result = tm_table(protocol, list(MPD_TM))
    assert result['tm'].tolist() == pytest.approx(list(MPD_TM.values()), abs=1e-9)
    # primer3 under the same conditions is 3.3-4.7 °C lower for these two targets.
    primer3_tm = [cal_temp(s, mv_conc=75, dv_conc=10, dntp_conc=1, dna_conc=100,
                           formamide_conc=20 * FORMAMIDE_MOLAR_PER_PCT) for s in MPD_TM]
    assert ((result['tm'] - primer3_tm) > 3).all()
    assert describe_conditions(protocol) == [
        'Tm [multi_padlock_design] (tm): 75 mM Na+, 10 mM Mg2+, 100 nM probe, 20% formamide']
    # Unset fields keep the legacy defaults; salt, probe and formamide move Tm as in the repository.
    legacy = tm_table(tm_protocol(tm_method='multi_padlock_design'), [SEQS[0]])['tm'].iloc[0]
    for field, value in (('na_mM', 390), ('mg_mM', 10), ('probe_nM', 500)):
        block = {field: value}
        assert tm_table(tm_protocol(block, tm_method='multi_padlock_design'), [SEQS[0]])['tm'].iloc[0] > legacy
    with_formamide = tm_table(tm_protocol({'formamide_pct': 10}, tm_method='multi_padlock_design'), [SEQS[0]])
    assert with_formamide['tm'].iloc[0] == pytest.approx(legacy - 7.2)
    assert tm_table(tm_protocol(tm_method='primer3'), SEQS)['tm'].tolist() == [legacy_tm(s) for s in SEQS]


@pytest.mark.parametrize('attr, message', [
    ({'type': 'annealing_temperature', 'tm_method': 'sugimoto'}, 'must be one of primer3, multi_padlock_design'),
    ({'type': 'gc_content', 'tm_method': 'multi_padlock_design'}, 'is only supported for annealing_temperature'),
])
def test_invalid_tm_method_rejected(attr, message):
    from uprobe.core.validate import validate_protocol
    protocol = {**minimal_protocol(), 'attributes': {'x': {'target': 'target_region', **attr}}}
    with pytest.raises(ValueError, match=r'attributes\.x\.tm_method ' + message):
        validate_protocol_hybridization(protocol)
    with pytest.raises(ValueError, match='tm_method'):
        _validate_and_normalize_protocol(copy.deepcopy(protocol))
    errors = [f for f in validate_protocol(protocol) if f.severity == 'error']
    assert [f.path for f in errors] == ['attributes.x.tm_method'] and message in errors[0].message


# OligoMiner ExampleFiles/3_probes.bed (github.com/brianbeliveau/OligoMiner @ fcb124f, MIT): the Tm column of
# the shipped example output, made with the defaults of blockParse.py (-s 390, -c 25, -C 25, -F 50, DNA_NN3).
OLIGOMINER_BED = {
    'CCCCATCCTTATTCCCTTTAAGAACGTTCTGATGGGC': 42.85, 'ATAAGTCAGAGCTTTACAGAGGTGTCACCTAACAAAACGC': 42.27,
    'AATGACTAAGAATCATTCCAAGTGTCACCATCAAGACCACC': 42.24, 'AGGATTATGTCTGTTGCTACTGGGATATCACTGACAAGTCT': 42.12,
    'TTGAAGCTCATATGTCACACAATCTCCAGAAGGCCT': 42.22, 'ACTACATAGTGCAATTTTAATTCCAACTCGCTCTTTCCCCT': 42.14,
    'AGACTTTTGCCATTTTCTTCTCATCTTGCTGCAATCATGT': 42.08, 'TCTTACCACATGGCTAATTCAAATTGGGGTTACAACAGTGA': 42.26,
    'AAAAGATGGACACCCTGAGGTCAATGATTTCCTCCC': 42.63, 'TTACTACAACATGTTATGTTAACCAACCTCCCTAGTCCCGC': 42.87,
    'GCAATGAGCTCACAGAAGTCAGGATGTGCACAGGCT': 45.79, 'ATCCTTCGCTGTTTCCATTCTTGGTTTAAGATTTCCGTTCT': 42.49,
}
OLIGOMINER_CONDITIONS = {'na_mM': 390, 'probe_nM': 25, 'formamide_pct': 50}


def test_oligominer_tm_reproduces_the_published_example_output():
    from uprobe.core.attributes import cal_temp_oligominer
    for seq, tm in OLIGOMINER_BED.items():
        assert cal_temp_oligominer(seq, **OLIGOMINER_CONDITIONS) == pytest.approx(tm, abs=0.0051)
    assert math.isnan(cal_temp_oligominer('ACGTNACGT', **OLIGOMINER_CONDITIONS))
    # The attribute passes the declared conditions; Mg2+ and dNTP are not part of the formula.
    protocol = tm_protocol({'dntp_mM': 1, 'mg_mM': 5}, tm_method='oligominer', hybridization=OLIGOMINER_CONDITIONS)
    result = tm_table(protocol, list(OLIGOMINER_BED))
    assert result['tm'].tolist() == pytest.approx(list(OLIGOMINER_BED.values()), abs=0.0051)
    assert describe_conditions(protocol) == ['Tm [oligominer] (tm): 390 mM Na+, 25 nM probe, 50% formamide']


def test_oligominer_tm_agrees_with_biopython_when_available():
    mt = pytest.importorskip('Bio.SeqUtils.MeltingTemp')
    from uprobe.core.attributes import cal_temp_oligominer
    for seq in OLIGOMINER_BED:
        reference = mt.chem_correction(round(mt.Tm_NN(seq, nn_table=mt.DNA_NN3, Na=390, dnac1=25, dnac2=25), 2), fmd=50)
        assert cal_temp_oligominer(seq, **OLIGOMINER_CONDITIONS) == pytest.approx(reference, abs=0.0051)


def test_merfish_design_tm_formula():
    from uprobe.core.attributes import cal_temp_merfish_design
    seq = 'CCCATCTCAGTCCTTCTCGCGGCCCTCCTG'
    # Worked by hand from TRDesigner.m GetRegionTm: SantaLucia & Hicks 2004 sums, 0.2 / -5.7 initiation,
    # no terminal A/T (C ... G), 0.368 (N-1) ln[Na] at 0.3 M, ln(C) of 5 nM (no /4).
    nn = {'CC': (-8.0, -19.9), 'CA': (-8.5, -22.7), 'AT': (-7.2, -20.4), 'TC': (-8.2, -22.2),
          'CT': (-7.8, -21.0), 'TG': (-8.5, -22.7), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9),
          'CG': (-10.6, -27.2), 'AG': (-7.8, -21.0), 'GT': (-8.4, -22.4), 'TT': (-7.6, -21.3),
          'AA': (-7.6, -21.3), 'GA': (-8.2, -22.2), 'AC': (-8.4, -22.4), 'TA': (-7.2, -21.3)}
    dh = 0.2 + sum(nn[seq[i:i + 2]][0] for i in range(29))
    ds = -5.7 + sum(nn[seq[i:i + 2]][1] for i in range(29)) + 0.368 * 29 * math.log(0.3)
    expected = 1000 * dh / (ds + 1.9872 * math.log(5e-9)) - 273.15
    assert cal_temp_merfish_design(seq, na_mM=300, probe_nM=5) == pytest.approx(expected, abs=1e-9)
    assert cal_temp_merfish_design(seq.lower().replace('t', 'u'), na_mM=300, probe_nM=5) == pytest.approx(expected)
    # A terminal A or T adds 2.2 kcal/mol and 6.9 cal/K/mol at that end.
    a_end = 'A' + seq[1:]
    dh2 = dh + 2.2 + (nn['AC'][0] - nn['CC'][0])
    ds2 = ds + 6.9 + (nn['AC'][1] - nn['CC'][1])
    assert cal_temp_merfish_design(a_end, na_mM=300, probe_nM=5) == pytest.approx(
        1000 * dh2 / (ds2 + 1.9872 * math.log(5e-9)) - 273.15, abs=1e-9)
    # The probe concentration enters as C where Primer3 uses C/4, which makes this formula warmer by about 1.4 C
    # for a 30-mer at 5 nM (the rest is the SantaLucia & Hicks 2004 table, initiation and rounding).
    protocol = tm_protocol(tm_method='merfish_design', hybridization={'na_mM': 300, 'probe_nM': 5})
    assert tm_table(protocol, [seq])['tm'].iloc[0] == pytest.approx(expected, abs=1e-9)
    primer3_tm = cal_temp(seq, mv_conc=300, dv_conc=0, dntp_conc=0, dna_conc=5)
    assert 1.0 < expected - primer3_tm < 1.8
    assert describe_conditions(protocol) == ['Tm [merfish_design] (tm): 300 mM Na+, 5 nM probe']


def test_design_code_tm_methods_validate():
    from uprobe.core.validate import validate_protocol
    for method in ('merfish_design', 'oligominer'):
        protocol = {**minimal_protocol(), 'attributes': {'x': {'target': 'target_region', 'type': 'annealing_temperature',
                                                              'tm_method': method}}}
        assert not [f for f in validate_protocol(protocol) if f.severity == 'error' and f.path.startswith('attributes')]
    protocol['attributes']['x']['hybridization'] = {'formamide_pct': 30, 'mg_mM': 2, 'na_mM': 300}
    protocol['attributes']['x']['tm_method'] = 'merfish_design'
    warnings = [f.message for f in validate_protocol(protocol) if f.path == 'attributes.x.tm_method']
    assert len(warnings) == 2 and 'does not use formamide_pct, mg_mM' in warnings[0]
    assert 'probe_nM not declared' in warnings[1] and 'probe_nM 50' in warnings[1]


def test_design_code_tm_methods_without_conditions_warn_about_the_legacy_defaults():
    # Without a hybridization block these formulas run at U-Probe's PCR defaults (50 mM Na+, 50 nM, no formamide):
    # 64.86 C for OligoMiner's first example probe, which OligoMiner's own conditions put at 42.85 C.
    from uprobe.core.attributes import cal_temp_oligominer
    from uprobe.core.hybridization import TM_DEFAULTS, legacy_formula_warnings, tm_formula_arguments
    from uprobe.core.validate import validate_protocol
    seq = next(iter(OLIGOMINER_BED))
    assert cal_temp_oligominer(seq, **tm_formula_arguments('oligominer', {})) == pytest.approx(64.86, abs=0.01)
    for method, published in (('oligominer', 'na_mM 390, probe_nM 25, formamide_pct 50'),
                               ('merfish_design', 'na_mM 300, probe_nM 5')):
        protocol = {**minimal_protocol(), 'attributes': {'x': {'target': 'target_region', 'type': 'annealing_temperature',
                                                              'tm_method': method}}}
        [warning] = [f for f in validate_protocol(protocol) if f.path == 'attributes.x.tm_method']
        assert warning.severity == 'warning' and 'legacy PCR defaults' in warning.message
        assert published in warning.message and f'na_mM {TM_DEFAULTS["na_mM"]:g}' in warning.message
        assert 'hybridization: {' in warning.suggestion
        assert len(legacy_formula_warnings(protocol)) == 1
        # Declared conditions (here the published ones, protocol-wide) silence it.
        protocol['hybridization'] = dict(item.split(' ') for item in published.split(', '))
        protocol['hybridization'] = {key: float(value) for key, value in protocol['hybridization'].items()}
        assert not [f for f in validate_protocol(protocol) if f.path == 'attributes.x.tm_method']
        assert legacy_formula_warnings(protocol) == []
    # multi_padlock_design keeps its legacy defaults without a warning (its protocols always declared them).
    protocol = {**minimal_protocol(), 'attributes': {'x': {'target': 'target_region', 'type': 'annealing_temperature',
                                                          'tm_method': 'multi_padlock_design'}}}
    assert not [f for f in validate_protocol(protocol) if f.path == 'attributes.x.tm_method']
    protocol['attributes']['x']['tm_method'] = 'bogus'
    errors = [f for f in validate_protocol(protocol) if f.severity == 'error' and f.path == 'attributes.x.tm_method']
    assert errors and 'merfish_design, oligominer' in errors[0].message
