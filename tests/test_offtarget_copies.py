"""exempt_copies: off-target genes that copy most of the target's sequence, and are not protein-coding, are not counted."""
import logging
import os
import random
import shutil
import time
from pathlib import Path

import pandas as pd
import pytest

from uprobe.core import copies
from uprobe.core.attributes import add_attributes
from uprobe.core.copies import copies_settings, copy_fractions, kmer_set, read_transcripts, target_sequences
from uprobe.core.offtarget import OffTargetCounter
from uprobe.core.report.formatting import ordered_columns
from uprobe.core.report.run_summary import build_run_summary
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import validate_protocol

needs_bowtie = pytest.mark.skipif(not (shutil.which('bowtie') and shutil.which('bowtie-build')),
                                  reason='needs Bowtie 1 (bowtie, bowtie-build)')
SETTINGS = {'aligner': 'bowtie', 'strand': 'sense', 'max_mismatches': 2,
            'exempt_copies': {'kmer': 30, 'min_shared': 0.5}}


def random_sequence(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def substitute_every(sequence, period, rng, first=None):
    """One substitution every `period` bases."""
    bases = list(sequence)
    for i in range(period // 2 if first is None else first, len(bases), period):
        bases[i] = rng.choice([b for b in 'ACGT' if b != bases[i]])
    return ''.join(bases)


def brute_share(target, transcripts, k):
    """Highest share of the target's distinct k-mers that one transcript contains."""
    def kmers(seq):
        return {seq[i:i + k] for i in range(len(seq) - k + 1) if set(seq[i:i + k]) <= set('ACGT')}
    own = kmers(target.upper().replace('U', 'T'))
    return max((len(own & kmers(t.upper().replace('U', 'T'))) / len(own) for t in transcripts), default=0.0)


def windows(target, sequence, length=20, step=1, first=1):
    """Candidate rows of a target: windows placed by target_start/target_end, as generate_target_seqs writes them."""
    rows = [(target, sequence[i:i + length], first + i, first + i + length - 1)
            for i in range(0, len(sequence) - length + 1, step)]
    return pd.DataFrame(rows, columns=['target', 'target_region', 'target_start', 'target_end'])


# -- settings --------------------------------------------------------------------

def test_settings():
    assert copies_settings(None) is None and copies_settings(False) is None
    assert copies_settings(True) == {'kmer': 30, 'min_shared': 0.5}
    assert copies_settings({'min_shared': 0.8}) == {'kmer': 30, 'min_shared': 0.8}
    assert copies_settings({'kmer': 25, 'min_shared': 1}) == {'kmer': 25, 'min_shared': 1.0}


@pytest.mark.parametrize('value,message', [
    ('yes', 'must be true, false or a mapping'),
    ({'k': 30}, 'unknown option'),
    ({'kmer': 33}, 'kmer must be an integer from 12 to 32'),
    ({'kmer': 30.0}, 'kmer must be an integer'),
    ({'min_shared': 0}, 'min_shared must be a number above 0 and at most 1'),
    ({'min_shared': 1.5}, 'min_shared must be'),
    ({'min_shared': True}, 'min_shared must be'),
])
def test_settings_rejected(value, message):
    with pytest.raises(ValueError, match=message):
        copies_settings(value, 'attributes.hits.exempt_copies')


# -- k-mers and sequences ----------------------------------------------------------

@pytest.mark.parametrize('k', [12, 17, 30, 32])
def test_kmer_codes_match_a_plain_encoding(k):
    rng = random.Random(k)
    sequence = ''.join(rng.choice('ACGTacgtNU') for _ in range(400))
    expected = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k].upper().replace('U', 'T')
        if set(kmer) <= set('ACGT'):
            expected.add(int(''.join(str('ACGT'.index(b)) for b in kmer), 4))
    assert sorted(expected) == kmer_set(sequence, k).tolist()
    assert kmer_set('ACGT', k).size == 0


def test_target_sequences_from_windows():
    rng = random.Random(1)
    a, b = random_sequence(rng, 120), random_sequence(rng, 60)
    df = pd.concat([windows('A', a, step=5), windows('B', b, length=20, step=30)], ignore_index=True)
    sequences = target_sequences(df)
    assert sequences['A'] == a.encode()[:len(sequences['A'])] and len(sequences['A']) >= 115
    # B's windows at 1-20 and 31-50 leave 21-30 unknown.
    assert sequences['B'] == (b[:20] + 'N' * 10 + b[30:50]).encode()
    # Rows without target_start are placed by start/end or sub_region ("11-30"), as remove_overlap does.
    genome = pd.DataFrame({'target': ['chr1:1-100'] * 2, 'target_region': [a[10:30], a[30:50]],
                           'sub_region': ['11-30', '31-50']})
    assert target_sequences(genome)['chr1:1-100'] == ('N' * 10 + a[10:50]).encode()
    assert target_sequences(pd.DataFrame({'target_region': ['ACGT']})) == {}


def test_read_transcripts(tmp_path):
    path = tmp_path / 'transcript.fa'
    path.write_text('>G1_T1\nACGT\nTTTT\n>G2_T2 description\nGGGG\n>G1_T3\nCCCC\n>G3_T4\nAAAA\n')
    assert list(read_transcripts(path, {'G1', 'G3'})) == [('G1', b'ACGTTTTT'), ('G1', b'CCCC'), ('G3', b'AAAA')]
    assert list(read_transcripts(path, set())) == []


# -- shared k-mers -------------------------------------------------------------------

def test_copy_fractions_match_brute_force(monkeypatch):
    monkeypatch.setattr(copies, '_CHUNK_BASES', 700)  # several chunks
    rng = random.Random(3)
    targets = {f'X{i}': random_sequence(rng, 300 + 50 * i) for i in range(3)}
    targets['X1'] = targets['X1'][:100] + 'N' + targets['X1'][101:]
    genes = {}
    for g in range(8):
        transcripts = []
        for t in range(rng.randint(1, 3)):
            source = targets[rng.choice(sorted(targets))]
            start, end = sorted(rng.sample(range(len(source)), 2))
            piece = substitute_every(source[start:end], rng.choice([7, 40, 90, 10 ** 6]), rng)
            if rng.random() < 0.2:
                piece = reverse_complement(piece)
            transcripts.append(random_sequence(rng, rng.randint(0, 50)) + piece + random_sequence(rng, 20))
        genes[f'G{g}'] = transcripts
    candidates = {target: set(genes) for target in targets}
    records = [(gene, seq.encode()) for gene, transcripts in genes.items() for seq in transcripts]
    for k in (12, 30):
        shares = copy_fractions(targets, candidates, iter(records), k)
        for (target, gene), share in shares.items():
            assert share == pytest.approx(brute_share(targets[target], genes[gene], k)), (k, target, gene)


def test_near_identical_copy_is_a_copy_and_a_paralog_is_not():
    rng = random.Random(5)
    target = random_sequence(rng, 1500)
    genes = {
        'pseudogene': [substitute_every(target, 100, rng)],  # 99% identical, retrocopy-like
        'paralog': [substitute_every(target, 12, rng)[:700] + target[700:760] + substitute_every(target, 12, rng)[760:]],
        'antisense': [reverse_complement(target)],  # same sequence, other orientation
        'halves': [target[:600] + random_sequence(rng, 900), random_sequence(rng, 900) + target[900:]],
    }
    shares = copy_fractions({'T': target}, {'T': set(genes)},
                            ((gene, seq) for gene, seqs in genes.items() for seq in seqs), 30)
    assert shares[('T', 'pseudogene')] == pytest.approx(brute_share(target, genes['pseudogene'], 30))
    assert 0.65 < shares[('T', 'pseudogene')] < 0.75  # 70 of every 100 30-mers avoid the substitutions
    assert 0 < shares[('T', 'paralog')] < 0.03  # one identical 60-nt stretch
    assert shares[('T', 'antisense')] == 0
    # Each transcript holds about 40%: a copy needs the k-mers in one molecule.
    assert 0.35 < shares[('T', 'halves')] < 0.42


# -- counting ------------------------------------------------------------------------

GENES = {'ENSG1': 'ACTB', 'ENSG2': 'ACTG1', 'ENSG3': 'POTEF', 'ENSG4': 'ACTBP2'}
BIOTYPES = {'ENSG1': 'gene_biotype "protein_coding"; ', 'ENSG2': 'gene_biotype "protein_coding"; ',
            'ENSG3': 'gene_biotype "protein_coding"; ', 'ENSG4': 'gene_biotype "processed_pseudogene"; '}


def write_gtf(path, biotypes=BIOTYPES):
    path.write_text(''.join(f'1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "{gene}"; transcript_id "T{i}"; '
                            f'gene_name "{name}"; {biotypes[gene]}\n' for i, (gene, name) in enumerate(GENES.items())))


@pytest.fixture
def genome(tmp_path):
    """A GTF naming four genes, and a transcript.fa where ACTBP2 copies ACTB and ACTG1 is a paralog."""
    rng = random.Random(11)
    target = random_sequence(rng, 600)
    transcripts = {'ENSG1_T1': target, 'ENSG2_T2': substitute_every(target, 10, rng),
                   'ENSG3_T3': random_sequence(rng, 500), 'ENSG4_T4': substitute_every(target, 150, rng),
                   'ENSG4_T5': random_sequence(rng, 300)}
    (tmp_path / 'transcript.fa').write_text(''.join(f'>{name}\n{seq}\n' for name, seq in transcripts.items()))
    write_gtf(tmp_path / 'genome.gtf')
    return {'fasta': str(tmp_path / 'genome.fa'), 'gtf': str(tmp_path / 'genome.gtf'), 'align_index': ['bowtie'],
            'target': target}


def candidates(target):
    df = windows('ACTB', target)
    hits = {seq: {'ENSG1': 0} for seq in df.target_region}
    first, second = df.target_region[0], df.target_region[300]
    hits[first] = {'ENSG1': 0, 'ENSG2': 1, 'ENSG3': 2, 'ENSG4': 0}
    hits[second] = {'ENSG1': 0, 'ENSG4': 1}
    return df, hits, first, second


def test_copies_are_listed_and_not_counted(genome, caplog):
    df, hits, first, second = candidates(genome['target'])
    report = {}
    OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    rows = df.set_index('target_region')
    assert rows.hits[first] == 3 and rows.hits_off_targets[first] == 'ACTG1:1, POTEF:2'
    assert rows.hits_copies[first] == 'ACTBP2:0' and rows.hits_copies[second] == 'ACTBP2:1'
    assert rows.hits[second] == 1 and rows.hits_off_targets[second] == ''
    assert 'hits_unexpressed' not in df  # the exemption needs no expression table
    stats = report['attributes']['hits']
    assert stats['exempt_copies'] == {'kmer': 30, 'min_shared': 0.5} and stats['copy_genes'] == 1
    assert list(stats['copies']) == ['ACTB'] and 0.75 < stats['copies']['ACTB']['ACTBP2'] < 0.85
    assert stats['off_target_genes'] == 3 and 'targets_without_sequence' not in stats
    assert 'copies_not_exempt' not in stats
    assert '1 off-target gene(s) contain at least 0.5 of the 30-mers of a target' in caplog.text
    assert f'ACTB: ACTBP2 ({stats["copies"]["ACTB"]["ACTBP2"]:.2f})' in caplog.text
    summary = build_run_summary({'name': 'x'}, {'final': df.iloc[[0]], 'off_targets': report})
    entry = summary['off_targets']['attributes']['hits']
    assert (entry['final_off_target_genes'], entry['final_copy_genes']) == (2, 1)
    protocol = {'extracts': {'target_region': {}}, 'attributes': {'hits': {'target': 'target_region'}}}
    assert ordered_columns(['hits_copies', 'target_region', 'hits', 'hits_off_targets'], protocol) == [
        'target_region', 'hits', 'hits_copies', 'hits_off_targets']


def test_protein_coding_copy_is_not_exempt(genome, caplog):
    # A protein-coding gene holding most of the target is a real off-target (mouse Ifit3b of Ifit3).
    write_gtf(Path(genome['gtf']), {**BIOTYPES, 'ENSG4': 'gene_biotype "protein_coding"; '})
    df, hits, first, second = candidates(genome['target'])
    report = {}
    OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    rows = df.set_index('target_region')
    assert rows.hits[first] == 4 and rows.hits_off_targets[first] == 'ACTBP2:0, ACTG1:1, POTEF:2'
    assert rows.hits[second] == 2 and set(df.hits_copies) == {''}
    stats = report['attributes']['hits']
    assert stats['copies'] == {} and stats['copy_genes'] == 0
    entry = stats['copies_not_exempt']['ACTB']['ACTBP2']
    assert entry['not_exempt'] == 'protein_coding' and 0.75 < entry['share'] < 0.85
    assert 'stay counted' in caplog.text and 'ACTBP2 (0.' in caplog.text
    summary = build_run_summary({'name': 'x'}, {'final': df, 'off_targets': report})
    assert summary['off_targets']['attributes']['hits']['copies_not_exempt']['ACTB']['ACTBP2'] == entry


def test_pseudogene_copy_by_gencode_gene_type_is_exempt(genome):
    write_gtf(Path(genome['gtf']), {**BIOTYPES, 'ENSG4': 'gene_type "unprocessed_pseudogene"; '})
    df, hits, first, _ = candidates(genome['target'])
    report = {}
    OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    rows = df.set_index('target_region')
    assert rows.hits[first] == 3 and rows.hits_copies[first] == 'ACTBP2:0'
    assert list(report['attributes']['hits']['copies']['ACTB']) == ['ACTBP2']


def test_copy_of_unknown_biotype_is_not_exempt_and_warned(genome, caplog):
    write_gtf(Path(genome['gtf']), {**BIOTYPES, 'ENSG4': ''})
    df, hits, first, _ = candidates(genome['target'])
    report = {}
    with caplog.at_level(logging.WARNING):
        OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    assert df.set_index('target_region').hits[first] == 4 and set(df.hits_copies) == {''}
    stats = report['attributes']['hits']
    assert stats['copies_not_exempt']['ACTB']['ACTBP2']['not_exempt'] == 'unknown_biotype'
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert 'does not exempt 1 copy gene(s) whose biotype the GTF does not give' in warnings[0]


def test_copies_come_before_expression(genome, tmp_path):
    table = tmp_path / 'tpm.tsv'
    table.write_text('gene_id\tTPM\nENSG1\t1000\nENSG2\t500\nENSG4\t0\n')
    df, hits, first, _ = candidates(genome['target'])
    protocol = {'gene_expression': {'table': str(table), 'min_value': 1}}
    report = {}
    OffTargetCounter(protocol, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    rows = df.set_index('target_region')
    # ACTBP2 is unexpressed too, but listed once, as a copy; POTEF is not in the table and counts.
    assert (rows.hits_copies[first], rows.hits_unexpressed[first]) == ('ACTBP2:0', '')
    assert rows.hits_off_targets[first] == 'ACTG1:1, POTEF:2'
    stats = report['attributes']['hits']
    assert (stats['copy_genes'], stats['below_min_value'], stats['missing_from_table']) == (1, 0, 1)


def test_threshold_and_no_copies(genome):
    df, hits, first, _ = candidates(genome['target'])
    report = {}
    strict = {**SETTINGS, 'exempt_copies': {'kmer': 30, 'min_shared': 0.9}}
    OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, strict)
    assert df.set_index('target_region').hits[first] == 4 and set(df.hits_copies) == {''}
    assert report['attributes']['hits']['copies'] == {} and report['attributes']['hits']['copy_genes'] == 0
    # Without the option nothing is read and there is no _copies column.
    df, hits, _, _ = candidates(genome['target'])
    OffTargetCounter({}, {'gtf': genome['gtf']}).annotate(df, 'hits', 'target_region', hits,
                                                          {**SETTINGS, 'exempt_copies': None})
    assert 'hits_copies' not in df


def test_targets_without_sequence_are_reported(genome, caplog):
    df, hits, first, _ = candidates(genome['target'])
    df = df.drop(columns=['target_start', 'target_end'])  # nothing places the windows
    report = {}
    OffTargetCounter({}, genome, report).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    assert report['attributes']['hits']['targets_without_sequence'] == ['ACTB']
    assert df.set_index('target_region').hits[first] == 4
    assert 'exempt_copies cannot check 1 target(s)' in caplog.text
    summary = build_run_summary({'name': 'x'}, {'final': df, 'off_targets': report})
    warning = next(w for w in summary['warnings'] if w['code'] == 'copies_not_checked')
    assert warning['attributes'] == {'hits': ['ACTB']}


def test_transcript_fasta_is_extracted_when_missing(genome, monkeypatch, tmp_path):
    from uprobe.core.gen import fun
    (tmp_path / 'transcript.fa').rename(tmp_path / 'saved.fa')
    calls = []
    monkeypatch.setattr(fun, 'extract_trans_seqs',
                        lambda gtf, fasta, out: calls.append(out) or shutil.copy(tmp_path / 'saved.fa', out))
    df, hits, first, _ = candidates(genome['target'])
    OffTargetCounter({}, genome).annotate(df, 'hits', 'target_region', hits, SETTINGS)
    assert calls == [tmp_path / 'transcript.fa'] and df.set_index('target_region').hits_copies[first] == 'ACTBP2:0'


# -- uprobe validate -----------------------------------------------------------------

def test_validate_checks_exempt_copies():
    import yaml
    from pathlib import Path
    protocol = yaml.safe_load((Path(__file__).parent / 'data' / 'RNA_format.yaml').read_text(encoding='utf-8'))
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    attr = protocol['attributes']['target_mappedGenes']
    attr['exempt_copies'] = {'kmer': 30, 'min_shared': 0.5}
    assert validate_protocol(protocol) == []
    attr['exempt_copies'] = {'kmer': 31, 'min_shared': 2}
    found = [f for f in validate_protocol(protocol) if f.path == 'attributes.target_mappedGenes.exempt_copies']
    assert len(found) == 1 and found[0].severity == 'error' and 'min_shared must be' in found[0].message
    assert 'kmer: 30' in found[0].suggestion
    attr['exempt_copies'] = True
    protocol['extracts']['target_region'].update(step=100)
    found = [f for f in validate_protocol(protocol) if f.path == 'attributes.target_mappedGenes.exempt_copies']
    assert len(found) == 1 and found[0].severity == 'warning' and 'leaves gaps' in found[0].message


# -- with Bowtie 1 -------------------------------------------------------------------

@needs_bowtie
def test_add_attributes_exempts_a_planted_copy(tmp_path, monkeypatch):
    from uprobe.core.tools import build_transcripts_index
    rng = random.Random(21)
    target = random_sequence(rng, 800)
    transcripts = {'GT_T1': target,
                   'GCOPY_T2': random_sequence(rng, 50) + substitute_every(target, 120, rng) + random_sequence(rng, 50),
                   'GPARA_T3': substitute_every(target, 9, rng)[:400] + target[400:440] + random_sequence(rng, 300)}
    transcripts.update({f'GR{i}_T{i}': random_sequence(rng, 600) for i in range(5)})
    (tmp_path / 'transcript.fa').write_text(''.join(f'>{n}\n{s}\n' for n, s in transcripts.items()))
    (tmp_path / 'genome.gtf').write_text(''.join(
        f'1\tt\texon\t1\t9\t.\t+\t.\tgene_id "{n.split("_")[0]}"; transcript_id "x"; gene_name "{n.split("_")[0][1:]}"; '
        f'gene_biotype "{"processed_pseudogene" if n.startswith("GCOPY") else "protein_coding"}";\n'
        for n in transcripts))
    build_transcripts_index(tmp_path / 'genome.gtf', tmp_path / 'genome.fa', tmp_path / 'bowtie_transcript',
                            threads=1, aligner='bowtie')
    monkeypatch.chdir(tmp_path)
    df = windows('T', target)
    protocol = {'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                                        'max_mismatches': 2, 'strand': 'sense', 'exempt_copies': True}}}
    genome = {'fasta': str(tmp_path / 'genome.fa'), 'gtf': str(tmp_path / 'genome.gtf'), 'align_index': ['bowtie']}
    report = {}
    result = add_attributes(df, protocol, genome, report)
    stats = report['attributes']['hits']
    assert list(stats['copies']) == ['T'] and list(stats['copies']['T']) == ['COPY']
    assert stats['copies']['T']['COPY'] == pytest.approx(brute_share(target, [transcripts['GCOPY_T2']], 30), abs=1e-3)
    assert result.hits_copies.str.contains('COPY').mean() > 0.9  # nearly every window hits the copy
    assert not result.hits_off_targets.str.contains('COPY').any()
    assert result.hits_off_targets.str.contains('PARA').any()  # the paralog still counts
    assert (result.hits == 1 + result.hits_off_targets.str.contains('PARA')).all()


# -- Ensembl 113 transcriptomes (read only; skipped when absent) -----------------------

# Folders with records.tsv (transcript_id, gene_id, gene_symbol, ..., ensembl_canonical) and
# sequences.fa (Ensembl cDNA + ncRNA), e.g. prepared by the probe-design benchmark.
REFERENCE = Path(os.environ.get('UPROBE_REFERENCE_DATA', '~/Projects/probe-design-bench/resources/data')).expanduser()


def reference(name):
    folder = REFERENCE / name
    if not (folder / 'sequences.fa').is_file() or not (folder / 'records.tsv').is_file():
        pytest.skip(f'reference transcriptome not available: {folder}')
    return folder / 'sequences.fa', pd.read_csv(folder / 'records.tsv', sep='\t')


def ensembl_records(path, gene_of):
    """(gene id, sequence) of every record of an Ensembl cDNA/ncRNA FASTA."""
    gene, lines = None, []
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                if gene is not None:
                    yield gene, b''.join(lines)
                gene, lines = gene_of[line[1:].split(None, 1)[0].decode()], []
            else:
                lines.append(line.strip())
    if gene is not None:
        yield gene, b''.join(lines)


def canonical_shares(name, targets, candidate_symbols=None):
    """Highest share of each target's canonical mRNA 30-mers in the genes named (default: every other gene)."""
    fasta, records = reference(name)
    gene_of = dict(zip(records.transcript_id, records.gene_id))
    symbol_of = dict(zip(records.gene_id, records.gene_symbol))
    chosen = records[records.ensembl_canonical == 1]
    canonical = dict(zip(chosen.transcript_id, chosen.gene_symbol))
    sequences, current = {}, None
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                current = canonical.get(line[1:].split(None, 1)[0])
            elif current in targets:
                sequences[current] = sequences.get(current, '') + line.strip()
    genes = {gene for gene, symbol in symbol_of.items() if candidate_symbols is None or symbol in candidate_symbols}
    candidates = {target: {gene for gene in genes if symbol_of[gene] != target} for target in targets}
    started = time.perf_counter()
    shares = copy_fractions(sequences, candidates, ensembl_records(fasta, gene_of), 30)
    seconds = time.perf_counter() - started
    best = {}
    for (target, gene), share in shares.items():
        best[target, symbol_of[gene]] = max(share, best.get((target, symbol_of[gene]), 0.0))
    return best, len(genes), seconds


def test_reference_human_paralogs_and_pseudogenes_against_the_whole_transcriptome():
    targets = ['ACTB', 'GAPDH', 'EEF1A1', 'RPL13A']
    shares, genes, seconds = canonical_shares('human_grch38_ens113_transcriptome', targets)
    top = {target: max((s, g) for (t, g), s in shares.items() if t == target) for target in targets}
    print(f'\n{genes} genes x {len(targets)} targets: {seconds:.1f} s; ACTB vs ACTG1 {shares["ACTB", "ACTG1"]:.4f}; '
          + '; '.join(f'{t}: {g} {s:.4f}' for t, (s, g) in top.items()))
    assert 0.03 < shares['ACTB', 'ACTG1'] < 0.07  # a paralog: about 5% of ACTB's 30-mers
    assert top['ACTB'][0] < 0.2
    # GAPDH's closest processed pseudogene holds under a quarter of its mRNA's 30-mers: human
    # GAPDH pseudogenes are old copies, and none is exempt. EEF1A1P5 and RPL13AP5 hold most of
    # the coding sequence (see the README), but under half of these canonical mRNAs with UTRs.
    assert top['GAPDH'][1] == 'GAPDHP1' and 0.2 < top['GAPDH'][0] < 0.3
    assert top['EEF1A1'][1] == 'EEF1A1P5' and 0.3 < top['EEF1A1'][0] < 0.4
    assert top['RPL13A'][1] == 'RPL13AP5' and 0.4 < top['RPL13A'][0] < 0.5
    assert seconds < 120  # every gene of the 410,909 transcripts checked against 4 targets


def test_reference_mouse_gapdh_retrocopy_is_exempt():
    shares, _, seconds = canonical_shares('mouse_grcm39_ens113_transcriptome', ['Gapdh', 'Actb'],
                                          {'Gm13394', 'Actg1', 'Gapdh', 'Actb'})
    print(f'\nGapdh vs Gm13394 {shares["Gapdh", "Gm13394"]:.4f}; Actb vs Actg1 {shares["Actb", "Actg1"]:.4f}; '
          f'{seconds:.1f} s')
    assert 0.7 < shares['Gapdh', 'Gm13394'] < 0.8  # a retrocopy: 76% of Gapdh's 30-mers
    assert shares['Actb', 'Actg1'] < 0.1
