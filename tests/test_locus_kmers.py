"""kmer_count locus_only (jellyfish count --if) and extracts.target_region.mask."""
import gzip
import random
import shutil
import subprocess

import pandas as pd
import pytest

from uprobe.core import attributes, tools
from uprobe.core.attributes import add_attributes
from uprobe.core.utils import RepeatMask, extract_fasta, reverse_complement
from uprobe.core.validate import validate_protocol

K = 12


def genome(tmp_path, length=600, seed=3):
    rng = random.Random(seed)
    seq = ''.join(rng.choice('ACGT') for _ in range(length))
    seq = seq[:300] + seq[40:80] + seq[340:]  # an internal repeat: the 40-mer at 41-80 recurs at 301-340
    path = tmp_path / 'ref.fa'
    path.write_text(f'>chr1\n{seq}\n')
    return path, seq


# -- kmer_count locus_only ----------------------------------------------------------------------

@pytest.fixture
def fake_jellyfish(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    built, queried = [], []

    def build(fasta, k, out, threads, size, canonical=True, only_kmers_in=None):
        built.append({'fasta': fasta, 'k': k, 'threads': threads, 'size': size, 'canonical': canonical,
                      'listed': open(only_kmers_in).read() if only_kmers_in else None})
        open(out, 'w').write('index')

    def query(outdir, name, records, index, k, threads):
        queried.append((str(index), dict(records)))
        return {key: 0 for key in records}
    monkeypatch.setattr(tools, 'build_jf_index', build)
    monkeypatch.setattr(attributes, 'cal_kmer_count', query)
    return built, queried


def test_locus_only_counts_just_the_candidate_kmers(fake_jellyfish, tmp_path):
    built, queried = fake_jellyfish
    seqs = ['ACGTACGGTTCA' * 3, 'TTGACCAGTAGGCATCAGGAC', 'ACGTACGGTTCA' * 3]
    protocol = {'attributes': {'k': {'target': 's', 'type': 'kmer_count', 'aligner': 'jellyfish', 'kmer_len': K,
                                     'locus_only': True, 'strand': 'sense'}}}
    # No jellyfish: true needed and no index is built next to the FASTA.
    add_attributes(pd.DataFrame({'s': seqs}), protocol, {'fasta': str(tmp_path / 'ref.fa')}, {}, 3)
    assert len(built) == 1 and built[0]['k'] == K and built[0]['threads'] == 3
    assert built[0]['canonical'] is False and built[0]['fasta'] == str(tmp_path / 'ref.fa')
    # Each distinct sequence is listed once, with its reverse complement (a stranded index holds both strands).
    listed = built[0]['listed'].split()
    assert listed[1::2] == [seqs[0], reverse_complement(seqs[0]), seqs[1], reverse_complement(seqs[1])]
    assert int(built[0]['size']) >= 2 * sum(len(s) - K + 1 for s in seqs[:2])
    assert not (tmp_path / 'jf_genome').exists()
    assert queried[0][0].endswith('locus_kmers.jf') and len(queried[0][1]) == 3


def test_locus_only_size_and_canonical_options(fake_jellyfish, tmp_path):
    built, _ = fake_jellyfish
    protocol = {'attributes': {'k': {'target': 's', 'type': 'kmer_count', 'aligner': 'jellyfish', 'kmer_len': K,
                                     'locus_only': True, 'size': '50M'}}}
    add_attributes(pd.DataFrame({'s': ['ACGTACGGTTCA' * 3]}), protocol, {'fasta': str(tmp_path / 'ref.fa')}, {})
    assert built[0]['size'] == '50M' and built[0]['canonical'] is True


def kmer_protocol(**options):
    return {'name': 'k', 'genome': 'test', 'targets': [{'x': 'ACGT' * 20}], 'encoding': {},
            'extracts': {'target_region': {'source': 'exon', 'length': 30}},
            'probes': {'probe_1': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
            'attributes': {'k': {'target': 'probe_1', 'type': 'kmer_count', 'aligner': 'jellyfish', **options}}}


def test_validate_locus_only():
    assert not [f for f in validate_protocol(kmer_protocol(locus_only=True)) if f.severity == 'error']
    errors = [f for f in validate_protocol(kmer_protocol(locus_only='yes')) if f.severity == 'error']
    assert [f.path for f in errors] == ['attributes.k.locus_only'] and 'true or false' in errors[0].message
    other = kmer_protocol()
    other['attributes']['k'] = {'target': 'probe_1', 'type': 'gc_content', 'locus_only': True}
    warnings = [f for f in validate_protocol(other) if f.path == 'attributes.k.locus_only']
    assert len(warnings) == 1 and warnings[0].severity == 'warning'


@pytest.mark.skipif(shutil.which('jellyfish') is None, reason='jellyfish is not installed')
def test_locus_only_counts_equal_the_genome_index(tmp_path, monkeypatch):
    # Real Jellyfish (>= 2.0 for --if): the candidates' counts match a whole-genome index's, for both strands.
    monkeypatch.chdir(tmp_path)
    ref, seq = genome(tmp_path)
    windows = [seq[i:i + 30] for i in range(0, 560, 23)] + [reverse_complement(seq[100:130])]
    frame = pd.DataFrame({'s': windows})
    results = {}
    for strand in ('both', 'sense', 'antisense'):
        for locus_only in (False, True):
            protocol = {'attributes': {'k': {'target': 's', 'type': 'kmer_count', 'aligner': 'jellyfish',
                                             'kmer_len': K, 'size': '1M', 'locus_only': locus_only,
                                             'strand': strand}}}
            results[strand, locus_only] = add_attributes(frame.copy(), protocol,
                                                         {'fasta': str(ref), 'jellyfish': True}, {})['k'].tolist()
        assert results[strand, True] == results[strand, False], strand
    assert max(results['both', True]) >= 2  # the repeat gives counts above 1


# -- extracts.target_region.mask ---------------------------------------------------------------

def mask_bed(tmp_path, name='repeats.bed', lines=('chr1\t100\t140',)):
    path = tmp_path / name
    path.write_text('\n'.join(lines) + '\n')
    return path


def test_bed_mask_drops_windows_like_a_soft_masked_genome(tmp_path):
    ref, seq = genome(tmp_path)
    bed = RepeatMask(mask_bed(tmp_path))  # bases 101-140 (1-based)
    plain = extract_fasta(ref, 'chr1:1-400', 40, step=20)
    dropped = extract_fasta(ref, 'chr1:1-400', 40, step=20, masked='drop', mask=bed)
    overlapping = {w[0] for w in plain if int(w[2].split('-')[0]) <= 140 and int(w[2].split('-')[1]) >= 101}
    assert overlapping and {w[0] for w in plain} - {w[0] for w in dropped} == overlapping
    # The same windows as a genome copy soft-masked at those bases.
    soft = tmp_path / 'soft.fa'
    soft.write_text(f'>chr1\n{seq[:100]}{seq[100:140].lower()}{seq[140:]}\n')
    assert dropped == extract_fasta(soft, 'chr1:1-400', 40, step=20, masked='drop')
    # Without masked: drop the mask does nothing; a region starting inside the repeat is shifted correctly.
    assert extract_fasta(ref, 'chr1:1-400', 40, step=20, mask=bed) == plain
    inside = extract_fasta(ref, 'chr1:121-400', 40, step=20, masked='drop', mask=bed)
    assert inside[0][2] == '21-60'  # bases 121-140 are masked, the first clean window starts at 141


def test_soft_masked_locus_fasta_and_gzip_bed(tmp_path):
    ref, seq = genome(tmp_path)
    locus = tmp_path / 'locus.fa'
    locus.write_text(f'>chr1:51-250\n{seq[50:100]}{seq[100:140].lower()}{seq[140:250]}\n')
    from_fasta = extract_fasta(ref, 'chr1:1-300', 40, step=20, masked='drop', mask=RepeatMask(locus))
    gz = tmp_path / 'repeats.bed.gz'
    with gzip.open(gz, 'wt') as handle:
        handle.write('# comment\nchr1\t100\t140\trpt\t0\t+\n')
    assert from_fasta == extract_fasta(ref, 'chr1:1-300', 40, step=20, masked='drop', mask=RepeatMask(gz))
    with pytest.raises(ValueError, match='expected a BED line'):
        RepeatMask(mask_bed(tmp_path, 'bad.bed', ('chr1 a b',)))
    (tmp_path / 'wrong.fa').write_text('>chr1:1-10\nACGT\n')
    with pytest.raises(ValueError, match='its name says 10'):
        RepeatMask(tmp_path / 'wrong.fa')


def region_protocol(ref, **region):
    return {'name': 'm', 'genome': 'test', 'targets': ['chr1:1-400'], 'encoding': {},
            'extracts': {'target_region': {'source': 'genome', 'length': 40, **region}},
            'probes': {'probe_1': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}}}


def test_mask_through_the_api_and_validate(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ref, seq = genome(tmp_path)
    bed = mask_bed(tmp_path)
    from uprobe.core.api import UProbeAPI
    api = object.__new__(UProbeAPI)
    api.genome = {'fasta': str(ref), 'gtf': None}
    api.protocol = region_protocol(ref, step=20, masked='drop', mask=str(bed))
    kept = api.generate_target_seqs()
    api.protocol = region_protocol(ref, step=20, masked='drop')
    assert len(kept) < len(api.generate_target_seqs())

    def findings(**region):
        protocol = {**region_protocol(ref, **region), 'genome': 'test'}
        return [f for f in validate_protocol(protocol) if f.path.startswith('extracts.target_region.mask')]
    assert not [f for f in findings(masked='drop', mask=str(bed)) if f.severity == 'error']
    assert 'only has an effect with masked: drop' in findings(mask=str(bed))[0].message
    missing = findings(masked='drop', mask=str(tmp_path / 'nope.bed'))
    assert missing[0].severity == 'error' and 'cannot read the mask' in missing[0].message
    assert findings(masked='drop', mask=5)[0].severity == 'error'
