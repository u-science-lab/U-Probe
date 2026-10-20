"""Region notation chrom:start-end is 1-based and inclusive everywhere."""
import random

import pytest

from uprobe.core.gen.fun import generate_target_seqs
from uprobe.core.process.otp import parse_region as otp_interval
from uprobe.core.utils import parse_region
from uprobe.core.validate import validate_protocol

GENOME = ''.join(random.Random(11).choice('ACGT') for _ in range(300))


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / 'g.fa'
    path.write_text(f'>chr1\n{GENOME}\n')
    return str(path)


def test_region_targets_are_one_based_and_inclusive(fasta):
    # chr1:101-160 is the 60 bases a genome browser shows from position 101 to 160.
    windows = generate_target_seqs('genome', ['chr1:101-160'], fasta, None, 20, step=20)
    assert windows['target_region'].tolist() == [GENOME[100:120], GENOME[120:140], GENOME[140:160]]
    assert windows['sub_region'].tolist() == ['1-20', '21-40', '41-60']
    # sub_region s-e is the window's position along the region: bases start + s - 1 to start + e - 1.
    for sub_region, sequence in zip(windows['sub_region'], windows['target_region']):
        s, e = map(int, sub_region.split('-'))
        assert GENOME[101 + s - 2:101 + e - 1] == sequence
    whole = generate_target_seqs('genome', ['chr1:1-300'], fasta, None, 300)
    assert whole['target_region'].tolist() == [GENOME]


def test_avoid_otp_reads_regions_like_targets(fasta):
    # The same string names the same bases as a target and as an avoid_otp region.
    chrom, start, end = otp_interval('chr1:101-160')  # 0-based, half-open
    region = generate_target_seqs('genome', ['chr1:101-160'], fasta, None, 60)['target_region'][0]
    assert (chrom, GENOME[start:end]) == ('chr1', region)


@pytest.mark.parametrize('region, message', [
    ('chr1:0-100', 'positions are 1-based'), ('chr1:200-100', 'ends before it starts'), ('chr1-100', 'not chrom'),
])
def test_invalid_regions(region, message):
    with pytest.raises(ValueError, match=message):
        parse_region(region)


def test_validate_rejects_a_zero_start():
    protocol = {'genome': 'g', 'targets': ['chr1:0-100', 'chr1:5-5'], 'encoding': {},
                'extracts': {'target_region': {'source': 'genome', 'length': 1}},
                'probes': {'p': {'template': '{a}', 'parts': {'a': {'expr': 'target_region'}}}}}
    errors = [f for f in validate_protocol(protocol) if f.severity == 'error']
    assert [(f.path, 'first base is 1' in f.message) for f in errors] == [('targets[0]', True)]
