"""Off-target genes that copy the target: the ``exempt_copies`` option of mapped_genes.

A processed pseudogene or a retrocopy can contain most of the target's
sequence, so every probe of the target also hits it. Counting such a
copy as an off-target leaves the target with few or no probes, although
probe sets for these targets work in practice. A ``mapped_genes`` attribute
can exempt copies::

    target_mapped_genes:
      target: target_region
      type: mapped_genes
      aligner: bowtie
      exempt_copies:     # or `exempt_copies: true` for these defaults
        kmer: 30
        min_shared: 0.5

An off-target gene is a copy of a target when one of its transcripts contains,
in the same orientation, at least ``min_shared`` of the distinct k-mers
(``kmer`` nt) of the target's sequence. A copy is exempt only when its GTF
biotype (``gene_biotype``, or ``gene_type`` in GENCODE) is not
``protein_coding``: a processed or unprocessed pseudogene, a lncRNA retrocopy.
An exempt copy is not counted; the attribute lists it in ``<attribute>_copies``.
A protein-coding copy is a separately quantified transcript, and binding it is
a real off-target (mouse Ifit3b holds 62.5% of the 30-mers of Ifit3), so it
stays counted and the report lists it as not exempt; so does a copy whose
biotype the GTF does not give. Most paralogs share far fewer exact k-mers
(ACTG1 contains 5% of the 30-mers of ACTB's mRNA, 7% of those of its coding
sequence) and are no copies at all.

The target's sequence is what its candidate windows cover, placed by
``target_start``/``target_end``: for ``source: mRNA``, ``exon``, ``CDS`` or
``UTR`` the target's merged exon model (the extracted segments with the introns
removed, see gen.fun.segment_offsets; with one selected transcript its spliced
sequence, junctions included, see gen.fun.SPLICED_SOURCES), so the share is
that of the sequence the probes are drawn from (the whole mature mRNA for
``mRNA``, the coding sequence of a protein-coding gene for ``exon`` and
``CDS``); for a direct sequence the sequence itself. Bases that no
window covers are unknown, and k-mers across them are skipped. The off-target
genes' transcripts are read from ``transcript.fa``, the FASTA the
transcriptome index is built from; only genes the alignment hit are read.
"""
import numbers
import time

import numpy as np
import pandas as pd

DEFAULT_KMER = 30
DEFAULT_MIN_SHARED = 0.5
MIN_KMER, MAX_KMER = 12, 32  # a k-mer is packed into 64 bits, 2 per base
OPTIONS = ('kmer', 'min_shared')
NEVER_EXEMPT = 'protein_coding'  # the biotype of copies that stay counted
UNKNOWN_BIOTYPE = 'unknown_biotype'

# Bases (DNA or RNA, any case) -> 0-3; everything else -> 4, which ends a k-mer.
_CODE = np.full(256, 4, dtype=np.uint8)
for _code, _bases in enumerate(('Aa', 'Cc', 'Gg', 'TtUu')):
    for _base in _bases:
        _CODE[ord(_base)] = _code
_HASH_BITS = 24  # membership bitmap of 16 MB that screens transcript k-mers
_HASH_MULTIPLIER = np.uint64(0x9E3779B97F4A7C15)
_CHUNK_BASES = 1 << 21  # transcript bases encoded at a time


def copies_settings(value, where='exempt_copies'):
    """{kmer, min_shared} of an exempt_copies option, or None when it is off.

    ``true`` means the defaults; a mapping may set either option. Raises
    ValueError starting with ``where``.
    """
    if value is None or value is False:
        return None
    if value is True:
        value = {}
    if not isinstance(value, dict):
        raise ValueError(f'{where} must be true, false or a mapping with {" and/or ".join(OPTIONS)}, got {value!r}')
    unknown = [str(key) for key in value if key not in OPTIONS]
    if unknown:
        raise ValueError(f'{where}: unknown option(s) {", ".join(unknown)}; expected: {", ".join(OPTIONS)}')
    kmer = value.get('kmer', DEFAULT_KMER)
    if type(kmer) is not int or not MIN_KMER <= kmer <= MAX_KMER:
        raise ValueError(f'{where}.kmer must be an integer from {MIN_KMER} to {MAX_KMER}, got {kmer!r}')
    min_shared = value.get('min_shared', DEFAULT_MIN_SHARED)
    if isinstance(min_shared, bool) or not isinstance(min_shared, numbers.Real) or not 0 < min_shared <= 1:
        raise ValueError(f'{where}.min_shared must be a number above 0 and at most 1, got {min_shared!r}')
    return {'kmer': kmer, 'min_shared': float(min_shared)}


# -- k-mers ------------------------------------------------------------------

def _kmer_codes(bases: np.ndarray, k: int):
    """(code of every k-mer start, whether the k-mer has only A/C/G/T) for base codes 0-4."""
    n = len(bases) - k + 1
    if n <= 0:
        return np.empty(0, np.uint64), np.empty(0, bool)
    invalid = bases > 3
    if invalid.any():
        ends = np.concatenate(([0], np.cumsum(invalid, dtype=np.int64)))
        valid = ends[k:] == ends[:-k]
        bases = np.where(invalid, 0, bases)
    else:
        valid = np.ones(n, bool)
    # Codes of 1, 2, 4, ... -mers by doubling, then k as a sum of powers of two.
    pieces = {1: bases.astype(np.uint64)}
    size = 1
    while size * 2 <= k:
        prev = pieces[size]
        pieces[size * 2] = (prev[:-size] << np.uint64(2 * size)) | prev[size:]
        size *= 2
    codes, offset = None, 0
    for size in sorted(pieces, reverse=True):
        if k & size:
            part = pieces[size][offset:offset + n]
            codes = part if codes is None else (codes << np.uint64(2 * size)) | part
            offset += size
    return codes, valid


def kmer_set(sequence, k: int) -> np.ndarray:
    """Sorted distinct k-mer codes of a sequence; k-mers with bases other than A/C/G/T(U) are skipped."""
    data = sequence.encode() if isinstance(sequence, str) else bytes(sequence)
    codes, valid = _kmer_codes(_CODE[np.frombuffer(data, np.uint8)], k)
    return np.unique(codes[valid])


def _hash(codes: np.ndarray) -> np.ndarray:
    return (codes * _HASH_MULTIPLIER) >> np.uint64(64 - _HASH_BITS)


# -- sequences ---------------------------------------------------------------

def target_sequences(df: pd.DataFrame, column: str = 'target_region') -> dict:
    """Target -> its sequence as the candidate windows cover it ('N' where none does).

    Windows are placed by target_start/target_end (1-based, along the
    target's merged exon model), or by start/end and sub_region where those
    are missing. Targets whose windows cannot be placed are left out.
    """
    target_column = next((c for c in ('target', 'gene') if c in df.columns), None)
    if target_column is None or column not in df.columns:
        return {}
    from .process.regions import region_bounds
    starts = pd.Series(np.nan, index=df.index)
    if 'target_start' in df.columns:
        starts = pd.to_numeric(df['target_start'], errors='coerce')
    if starts.isna().any():
        try:
            starts = starts.fillna(region_bounds(df)['start'])
        except ValueError:
            pass
    sequences = {}
    for target, start, window in zip(df[target_column].astype(str), starts, df[column]):
        if pd.isna(start) or not isinstance(window, str) or not window:
            continue
        start = int(start) - 1
        bases = sequences.setdefault(target, bytearray())
        if len(bases) < start + len(window):
            bases.extend(b'N' * (start + len(window) - len(bases)))
        bases[start:start + len(window)] = window.encode()
    return {target: bytes(bases) for target, bases in sequences.items()}


def read_transcripts(path, genes):
    """(gene, sequence) of each transcript.fa record of the given genes.

    Records are named {gene_id}_{transcript_id}; the gene is the part before
    the first underscore, as attributes.transcript_gene reads it.
    """
    genes = set(genes)
    gene, keep, lines = None, False, []
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                if keep:
                    yield gene, b''.join(lines)
                name = line[1:].split(None, 1)[0].decode() if line[1:].strip() else ''
                gene = name.split('_')[0]
                keep, lines = gene in genes, []
            elif keep:
                lines.append(line.strip())
    if keep:
        yield gene, b''.join(lines)


# -- shared k-mers -------------------------------------------------------------

class _TargetKmers:
    """The distinct k-mers of several targets, sorted by code, for screening transcripts."""

    def __init__(self, sequences: dict, k: int):
        self.k = k
        self.names = list(sequences)
        sets = [kmer_set(sequences[name], k) for name in self.names]
        self.sizes = np.array([len(s) for s in sets], dtype=np.int64)
        codes = np.concatenate(sets) if sets else np.empty(0, np.uint64)
        owners = np.repeat(np.arange(len(sets), dtype=np.int64), self.sizes)
        order = np.argsort(codes, kind='stable')
        self.codes, self.owners = codes[order], owners[order]
        self.bitmap = np.zeros(1 << _HASH_BITS, dtype=bool)
        self.bitmap[_hash(self.codes)] = True

    def scan(self, sequences: list, first: int):
        """(record, entry) pairs: transcript first+i contains the target k-mer self.codes[entry]."""
        separator = b'N'
        data = separator.join(sequences)
        lengths = np.array([len(s) + 1 for s in sequences], dtype=np.int64)
        starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
        codes, valid = _kmer_codes(_CODE[np.frombuffer(data, np.uint8)], self.k)
        if not len(codes):
            return None
        positions = np.flatnonzero(valid & self.bitmap[_hash(codes)])
        if not len(positions):
            return None
        query = codes[positions]
        left = np.searchsorted(self.codes, query, 'left')
        right = np.searchsorted(self.codes, query, 'right')
        found = right > left
        positions, left, counts = positions[found], left[found], (right - left)[found]
        if not len(positions):
            return None
        # One pair per target holding the k-mer (usually one).
        entries = np.repeat(left, counts) + np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
        records = first + np.searchsorted(starts, np.repeat(positions, counts), 'right') - 1
        return records, entries


def copy_fractions(sequences: dict, candidates: dict, records, k: int = DEFAULT_KMER) -> dict:
    """Share of each target's distinct k-mers that one transcript of a candidate gene contains.

    ``sequences`` maps target -> sequence, ``candidates`` target -> genes to
    check, ``records`` yields (gene, transcript sequence), e.g. from
    read_transcripts. Returns {(target, gene): highest share over the gene's
    transcripts} for every candidate pair (0.0 when nothing is shared).
    """
    targets = {name: sequences[name] for name in sequences if candidates.get(name)}
    result = {(target, gene): 0.0 for target in targets for gene in candidates[target]}
    if not targets:
        return result
    index = _TargetKmers(targets, k)
    wanted = set().union(*(candidates[target] for target in targets))
    genes, batch, size, found = [], [], 0, []

    def flush():
        nonlocal batch, size
        if batch:
            pairs = index.scan(batch, len(genes) - len(batch))
            if pairs is not None:
                found.append(pairs)
        batch, size = [], 0

    for gene, sequence in records:
        if gene not in wanted:
            continue
        genes.append(gene)
        batch.append(sequence.encode() if isinstance(sequence, str) else sequence)
        size += len(sequence) + 1
        if size >= _CHUNK_BASES:
            flush()
    flush()
    if not found:
        return result
    record = np.concatenate([pairs[0] for pairs in found])
    entry = np.concatenate([pairs[1] for pairs in found])
    # Distinct k-mers per (transcript, target), then the best transcript of each gene.
    pair = np.unique(record * len(index.codes) + entry)
    record, owner = pair // len(index.codes), index.owners[pair % len(index.codes)]
    keys, shared = np.unique(record * len(index.names) + owner, return_counts=True)
    for key, count in zip(keys.tolist(), shared.tolist()):
        record, owner = divmod(key, len(index.names))
        target, gene = index.names[owner], genes[record]
        if gene in candidates[target]:
            share = count / index.sizes[owner]
            if share > result[(target, gene)]:
                result[(target, gene)] = share
    return result


def find_copies(df: pd.DataFrame, candidates: dict, transcript_fasta, settings: dict, biotypes: dict,
                column: str = 'target_region'):
    """({(target, gene): share} of the exempt copies among the candidates, report).

    ``candidates`` maps each target to the off-target genes its sequences
    hit; ``settings`` is copies_settings(); ``biotypes`` maps gene -> GTF
    biotype. Only copies whose biotype is known and not protein_coding are
    exempt. The report lists the other copies under ``not_exempt``
    ({(target, gene): (share, 'protein_coding' or 'unknown_biotype')}), the
    targets whose windows give no k-mer (they get no copies) and the seconds
    taken.
    """
    started = time.perf_counter()
    sequences = target_sequences(df, column)
    genes = set().union(*candidates.values()) if candidates else set()
    shares = copy_fractions(sequences, candidates, read_transcripts(transcript_fasta, genes), settings['kmer'])
    unknown = sorted(target for target, hit in candidates.items()
                     if hit and not len(kmer_set(sequences.get(target, b''), settings['kmer'])))
    copies, not_exempt = {}, {}
    for (target, gene), share in shares.items():
        if share < settings['min_shared']:
            continue
        biotype = biotypes.get(gene)
        if biotype and biotype != NEVER_EXEMPT:
            copies[(target, gene)] = share
        else:
            not_exempt[(target, gene)] = (share, biotype or UNKNOWN_BIOTYPE)
    return copies, {'not_exempt': not_exempt, 'targets_without_sequence': unknown,
                    'seconds': round(time.perf_counter() - started, 2)}
