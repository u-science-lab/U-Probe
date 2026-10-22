import math

import pandas as pd
import pytest

from uprobe.core.attributes import add_attributes, cal_rna_dna_dg37
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import validate_protocol

# Probe sequences and dG37 (kcal/mol) from the output of the Oligostan script itself:
# Probes_XPO1_PNAS_1-2-4/Probes_XPO1_ALL.txt (column dG37) in Oligostan_Tutorial_Data_160527.zip,
# bitbucket.org/muellerflorian/fish_quant, Oligostan/Documentation (Tsanov et al. 2016, NAR 44:e165).
# Every one of the 68 probes of the file (and of the two sibling files) matches to 1e-9; ten are kept here.
OLIGOSTAN_OUTPUT = [
    ('TAATCACACATTTCTTCTGGAATCTCATGTGG', -32.078494),
    ('TTCTGGTTCTCTAGCAGCTGGGACATTTC', -32.678494),
    ('CTCCACTAATGGAGCCTATTGCCCAA', -31.478494),
    ('ATTTTCCGATATACACCTTTTCTTTCTCTACA', -31.378494),
    ('TCACAAAGAGCTTTACTTGAGCATCTTGTAGG', -31.278494),
    ('TGATAATCAATGAGAACTGCATCCAACAGAGG', -31.478494),
    ('TTAATGAAAGTATCACAAGCCATATCCTGGAC', -31.178494),
    ('TCTCTATAGAGTTCAGCAGCCAAATGATTCC', -32.778494),
    ('ATCACCATATTATTTTGACAGAGACTTTCGCT', -31.278494),
    ('GAACAAAATTTTCAGCGACCATCTGTGGATCA', -31.378494),
]


def test_dg37_reproduces_the_oligostan_output():
    for probe, published in OLIGOSTAN_OUTPUT:
        assert cal_rna_dna_dg37(probe) == pytest.approx(published, abs=1e-6)


def test_dg37_by_hand():
    # AA -0.2, AC -1.5, CC -2.2 plus the salt term -(ln 0.115 * -0.175 - 0.2) = -0.178494
    assert cal_rna_dna_dg37('AACC') == pytest.approx(-0.2 - 1.5 - 2.2 - 0.178494, abs=1e-6)
    assert -(math.log(0.115) * -0.175 - 0.2) == pytest.approx(-0.1784940514)
    # The salt term is a constant, not scaled by length; there is no initiation term.
    assert cal_rna_dna_dg37('AA') == pytest.approx(-0.2 - 0.178494, abs=1e-6)


def test_target_orientation_and_alphabet():
    probe = OLIGOSTAN_OUTPUT[0][0]
    target = reverse_complement(probe)  # the RNA in sense orientation
    assert cal_rna_dna_dg37(target, 'target') == pytest.approx(cal_rna_dna_dg37(probe))
    assert cal_rna_dna_dg37(probe.lower().replace('t', 'u')) == pytest.approx(cal_rna_dna_dg37(probe))
    # As in the script, a dinucleotide with an unknown base counts 0.
    assert cal_rna_dna_dg37('AANCC') == pytest.approx(-0.2 - 2.2 - 0.178494, abs=1e-6)
    assert math.isnan(cal_rna_dna_dg37('A')) and math.isnan(cal_rna_dna_dg37(None))


def protocol(**attr):
    return {'attributes': {'dg': {'target': 's', 'type': 'rna_dna_dg37', **attr}}}


def test_attribute_computes_dg37_of_the_probe_or_of_the_target():
    probes = [p for p, _ in OLIGOSTAN_OUTPUT]
    result = add_attributes(pd.DataFrame({'s': probes}), protocol(), {})
    assert result['dg'].tolist() == pytest.approx([v for _, v in OLIGOSTAN_OUTPUT], abs=1e-6)
    targets = [reverse_complement(p) for p in probes]
    by_target = add_attributes(pd.DataFrame({'s': targets}), protocol(sequence='target'), {})
    assert by_target['dg'].tolist() == pytest.approx(result['dg'].tolist())
    with pytest.raises(ValueError, match='sequence must be one of probe, target'):
        add_attributes(pd.DataFrame({'s': probes}), protocol(sequence='rna'), {})


def minimal(**attr):
    return {'name': 'dg', 'genome': 'test', 'targets': [{'x': 'ACGT' * 20}],
            'extracts': {'target_region': {'source': 'exon', 'length': 30}}, 'encoding': {},
            'probes': {'probe_1': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
            'attributes': {'dg': {'target': 'probe_1', 'type': 'rna_dna_dg37', **attr}},
            'post_process': {'filters': {'dg': {'condition': 'dg >= -33 & dg <= -31'}}}}


def test_validate_accepts_the_attribute_and_rejects_bad_options():
    assert not [f for f in validate_protocol(minimal()) if f.severity == 'error']
    assert not [f for f in validate_protocol(minimal(sequence='target')) if f.severity == 'error']
    errors = [f for f in validate_protocol(minimal(sequence='rna')) if f.severity == 'error']
    assert [f.path for f in errors] == ['attributes.dg.sequence'] and 'probe or target' in errors[0].message
    other = minimal()
    other['attributes']['dg'] = {'target': 'probe_1', 'type': 'gc_content', 'sequence': 'target'}
    warnings = [f for f in validate_protocol(other) if f.path == 'attributes.dg.sequence']
    assert len(warnings) == 1 and warnings[0].severity == 'warning'


def test_filters_accept_negative_numbers():
    from uprobe.core.process import filter_table, parse_condition
    assert parse_condition('dg >= -33 & dg <= -31.5') == [('dg', '>=', -33), ('dg', '<=', -31.5)]
    df = pd.DataFrame({'dg': [-34.0, -32.5, -31.2, -30.0]})
    kept = filter_table(df, {'dg': {'condition': 'dg >= -33 & dg <= -31'}})
    assert kept['dg'].tolist() == [-32.5, -31.2]
