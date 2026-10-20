"""Region notation chrom:start-end is 1-based and inclusive everywhere."""
import logging
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


@pytest.fixture
def log(caplog):
    """Records of uprobe.core.utils, whose logger does not propagate to caplog's root handler."""
    logger = logging.getLogger('uprobe.core.utils')
    logger.addHandler(caplog.handler)
    yield caplog
    logger.removeHandler(caplog.handler)


def test_masked_drop_skips_soft_masked_and_n_windows(tmp_path, log):
    # A soft-masked FASTA (repeats in lowercase), as UCSC hg38.fa or Ensembl dna_sm; N is hard-masked or a gap.
    genome = GENOME[:40] + GENOME[40:50].lower() + GENOME[50:80] + 'N' + GENOME[81:120]
    path = tmp_path / 'sm.fa'
    path.write_text(f'>chr1\n{genome}\n')
    kept = generate_target_seqs('genome', ['chr1:1-120'], str(path), None, 20, step=10)
    everything = generate_target_seqs('genome', ['chr1:1-120'], str(path), None, 20, step=10, masked='keep')
    dropped = generate_target_seqs('genome', ['chr1:1-120'], str(path), None, 20, step=10, masked='drop')
    # keep (the default) upper-cases every window, as before.
    assert kept.equals(everything)
    assert everything['target_region'].tolist()[2] == GENOME[20:40] and everything['target_region'].str.isupper().all()
    # drop skips the windows with a lowercase base (41-50) or the N (81); the others keep their probe_id.
    assert dropped['sub_region'].tolist() == ['1-20', '11-30', '21-40', '51-70', '61-80', '91-110', '101-120']
    assert dropped['probe_id'].tolist() == [f'chr1:1-120_{n}' for n in (1, 2, 3, 6, 7, 10, 11)]
    assert dropped['target_region'].tolist() == everything.set_index('probe_id').loc[
        dropped['probe_id'], 'target_region'].tolist()
    assert 'skipped 4 of 11 windows' in log.text


def test_masked_drop_warns_on_an_unmasked_fasta(fasta, log):
    windows = generate_target_seqs('genome', ['chr1:1-100'], fasta, None, 20, step=20, masked='drop')
    assert len(windows) == 5
    assert 'no soft-masked (lowercase) base' in log.text
    with pytest.raises(ValueError, match='masked must be one of keep, drop'):
        generate_target_seqs('genome', ['chr1:1-100'], fasta, None, 20, masked='hard')
