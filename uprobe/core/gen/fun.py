import bisect
import itertools
import pandas as pd
import os
import tempfile
import re
import typing as t
from pyfaidx import Fasta

from uprobe.core.utils import *
from uprobe.core.sampling import normalize_sampling, iter_target_windows

log = get_logger(__name__)

GTF_FIELDS = [
    'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'score', 'info'
]


def read_gtf(
        gtf: str,
        filter_by_type: t.Optional[str]=None,
        extract_fields: t.List[str]=["gene_name"],
        get_length=True,
        basic_fields=GTF_FIELDS,
        ) -> pd.DataFrame:
    df = pd.read_csv(gtf, sep='\t', header=None, comment='#')
    df.columns = basic_fields
    if filter_by_type:
        df = df[df[basic_fields[2]] == filter_by_type]
    if get_length:
        df['length'] = df['end'] - df['start'] + 1
    def extract_attribute(key):
        # Match a complete GTF attribute key, not suffixes such as other_gene.
        pattern = rf'(?:^|;)\s*{re.escape(key)}\s+"([^"\r\n]*)"'
        return df[basic_fields[-1]].str.extract(pattern, expand=False)

    for f in extract_fields:
        df[f] = extract_attribute(f)
        if f == "gene_name":
            missing = df[f].isna() | df[f].str.strip().eq("")
            df.loc[missing, f] = extract_attribute("gene").loc[missing]
        elif f == "gene_id":
            # RefSeq annotations may omit gene_id but retain NCBI GeneID.
            for fallback in (
                df[basic_fields[-1]].str.extract(r'GeneID:(\d+)', expand=False),
                extract_attribute("gene_name"),
                extract_attribute("gene"),
            ):
                missing = df[f].isna() | df[f].str.strip().eq("")
                df.loc[missing, f] = fallback.loc[missing]
    chr_new = []
    for chr_ in df.chr:
        chr_new.append(str(chr_).replace('chr', ''))
    df.chr = chr_new
    return df

def process_gtf_inplace(filepath):
    """
    If gene_name does not exist, replace gene with gene_name.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
    updated_lines = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            updated_lines.append(line)
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            log.warning(f"not a valid GTF line: {line.strip()}")
            updated_lines.append(line)
            continue
        attributes = fields[8]
        if "gene_name" in attributes:
            updated_lines.append(line)
        else:
            updated_attributes = []
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                if attr.startswith("gene "):
                    updated_attributes.append(attr.replace("gene ", "gene_name ")) 
                else:
                    updated_attributes.append(attr)
            fields[8] = '; '.join(updated_attributes)
            updated_lines.append('\t'.join(fields) + '\n')
    with open(filepath, 'w') as file:
        file.writelines(updated_lines)

# What each RNA source of extracts.target_region extracts from the GTF, per gene:
#   exon  its CDS features when it has any, else its exon features: the coding
#         sequence of a protein-coding gene (no UTRs), the whole processed RNA of
#         a non-coding gene. Segments shorter than the window are left out.
#   CDS   the same as exon.
#   mRNA  all its exon features: the whole mature RNA (5' UTR, CDS and 3' UTR of
#         a protein-coding gene; the processed RNA of a non-coding gene). Exons
#         shorter than the window give no window but still count in positions.
#   UTR   its UTR features (GTF type UTR, or five_prime_utr and three_prime_utr
#         as Ensembl names them). Segments shorter than the window are left out.
# A target left without segments gets no windows (summary.json warns with
# no_candidate_windows); the other targets are extracted as usual.
# Each window lies within one segment, except for a target whose windows come
# from one transcript (see SPLICED_SOURCES).
RNA_SOURCES = ('exon', 'CDS', 'mRNA', 'UTR')
UTR_TYPES = ('UTR', 'five_prime_utr', 'three_prime_utr')

# Windows across exon-exon junctions. When transcripts (canonical, mane or a
# list) selects exactly one transcript of a target, mRNA tiles the spliced RNA
# of that transcript and exon and CDS its spliced coding sequence (its exons if
# it has no CDS): the segments are joined 5' to 3', every one however short, so
# windows also span the junctions between them and bind only the mature RNA.
#   target_start/end  the window on the spliced sequence (1-based, inclusive)
#   exon_name         the segment the window lies on or, across junctions, the
#                     segments it spans, 5' to 3', joined by ';'
#   start/end         the window on these segments joined, so on its segment
#                     for a window within one, as for the other windows
#   transcript_names  the transcript; n_trans 1
# transcripts: all (the merged exon model), several transcripts of a target,
# and source UTR keep one window per segment.
SPLICED_SOURCES = ('exon', 'CDS', 'mRNA')
SEGMENT_SEPARATOR = ';'


def tiles_spliced(source: str, entry: dict) -> bool:
    """Whether a target's windows tile one spliced transcript (see SPLICED_SOURCES).

    ``entry`` is the target's select_transcripts entry.
    """
    return (source in SPLICED_SOURCES and entry['selection'] != 'all'
            and entry['ids'] is not None and len(entry['ids']) == 1)


def gene_segments(df_gene: pd.DataFrame, source: str) -> pd.DataFrame:
    """The GTF rows of one gene that an RNA source extracts (before the length filter)."""
    if source == 'mRNA':
        return df_gene[df_gene.type == 'exon']
    if source == 'UTR':
        return df_gene[df_gene.type.isin(UTR_TYPES)]
    df_cds = df_gene[df_gene.type == 'CDS']
    return df_cds if df_cds.shape[0] else df_gene[df_gene.type == 'exon']


# extracts.target_region.transcripts: whose features an RNA source extracts.
#   all        every transcript of the gene (the union of their exons)
#   canonical  the transcript tagged Ensembl_canonical
#   mane       the transcript tagged MANE_Select; a gene without one (outside
#              human, most non-coding genes) uses its Ensembl_canonical transcript
#   [ids]      the listed transcripts; each target needs at least one of them
# Ensembl and GENCODE GTFs write the tags on every line of a transcript.
TRANSCRIPT_RULES = ('all', 'canonical', 'mane')
CANONICAL_TAG, MANE_TAG = 'Ensembl_canonical', 'MANE_Select'


class TranscriptSelectionError(ValueError):
    """The transcripts option does not fit the targets or the GTF: a protocol error, not a missing sequence."""

    def __init__(self, problems):
        self.problems = list(problems)
        super().__init__('extracts.target_region.transcripts: ' + '; '.join(self.problems))


def transcripts_setting(value):
    """'all', 'canonical', 'mane' or a list of transcript ids; ValueError for anything else."""
    if value is None:
        return 'all'
    if isinstance(value, str) and value in TRANSCRIPT_RULES:
        return value
    if isinstance(value, list) and value and all(isinstance(v, str) and v.strip() for v in value):
        return list(dict.fromkeys(v.strip() for v in value))
    raise ValueError(f'transcripts must be all, canonical, mane or a list of transcript ids, got {value!r}')


def _without_version(transcript_id: str) -> str:
    return re.sub(r'\.\d+$', '', transcript_id)


def _tagged(rows: pd.DataFrame, tag: str) -> t.List[str]:
    """Transcript ids of the rows with the GTF attribute tag "<tag>"."""
    found = rows['info'].str.contains(rf'(?:^|;)\s*tag\s+"{tag}"', regex=True, na=False)
    return sorted(set(rows.loc[found, 'transcript_id'].dropna()))


def _listed_transcripts(df_gtf: pd.DataFrame, listed, targets) -> t.Tuple[dict, t.List[str]]:
    """({target: listed ids of its transcripts as the GTF writes them}, problems)."""
    known = df_gtf[['transcript_id', 'gene_name']].dropna().drop_duplicates()
    genes = dict(zip(known.transcript_id, known.gene_name))
    by_base = {}
    for transcript_id in genes:
        by_base.setdefault(_without_version(transcript_id), []).append(transcript_id)
    chosen, problems = {}, []
    for name in listed:
        # Ensembl GTFs write ENST00000646664, GENCODE ENST00000646664.1: match either way.
        matches = [name] if name in genes else by_base.get(_without_version(name), [])
        if len(matches) != 1:
            problems.append(f'transcript {name} is not in the GTF' if not matches else
                            f'transcript {name} matches {", ".join(matches)} in the GTF; give the version')
        elif genes[matches[0]] not in targets:
            problems.append(f'transcript {name} belongs to {genes[matches[0]]}, which is not a target')
        else:
            chosen.setdefault(genes[matches[0]], []).append(matches[0])
    return chosen, problems


def select_transcripts(df_gtf: pd.DataFrame, targets, setting) -> t.Dict[str, dict]:
    """The transcripts each target's windows come from, per transcripts_setting().

    Returns {target: {'selection': all|canonical|mane|list, 'chosen_by':
    all|MANE_Select|Ensembl_canonical|listed, 'ids': [...] or None for all}}.
    Raises TranscriptSelectionError listing every problem.
    """
    setting = transcripts_setting(setting)
    targets = list(dict.fromkeys(targets))
    if setting == 'all':
        return {target: {'selection': 'all', 'chosen_by': 'all', 'ids': None} for target in targets}
    rows = df_gtf[df_gtf.gene_name.isin(targets) & df_gtf.transcript_id.notna()]
    rows_of = dict(tuple(rows.groupby('gene_name', sort=False)))
    problems, result = [], {}
    if isinstance(setting, list):
        chosen, problems = _listed_transcripts(df_gtf, setting, set(targets))
        for target in targets:
            if target in chosen:
                result[target] = {'selection': 'list', 'chosen_by': 'listed', 'ids': sorted(chosen[target])}
            elif target in rows_of:
                example = (_tagged(rows_of[target], MANE_TAG) or _tagged(rows_of[target], CANONICAL_TAG)
                           or sorted(set(rows_of[target].transcript_id)))[0]
                problems.append(f'no listed transcript belongs to target {target}; list one of its transcripts '
                                f'(e.g. {example}) or use canonical, mane or all')
    else:
        for target in targets:
            if target not in rows_of:
                continue  # not in the GTF: reported by target validation
            ids, chosen_by = [], MANE_TAG
            if setting == 'mane':
                ids = _tagged(rows_of[target], MANE_TAG)
            if not ids:
                ids, chosen_by = _tagged(rows_of[target], CANONICAL_TAG), CANONICAL_TAG
                if ids and setting == 'mane':
                    log.info(f'{target} has no {MANE_TAG} transcript; using its {CANONICAL_TAG} transcript '
                             f'{", ".join(ids)}')
            if not ids:
                wanted = f'{MANE_TAG} or {CANONICAL_TAG}' if setting == 'mane' else CANONICAL_TAG
                problems.append(f'no transcript of {target} is tagged {wanted} in the GTF; list its transcript '
                                'ids or use transcripts: all')
                continue
            result[target] = {'selection': setting, 'chosen_by': chosen_by, 'ids': ids}
    if problems:
        raise TranscriptSelectionError(problems)
    for target, entry in result.items():
        log.info(f'{target}: windows from transcript(s) {", ".join(entry["ids"])} ({entry["chosen_by"]})')
    return result


def restrict_to_transcripts(df_gtf: pd.DataFrame, selection: t.Dict[str, dict]) -> pd.DataFrame:
    """The GTF rows minus those of the targets' unselected transcripts (see select_transcripts)."""
    chosen = {target: set(entry['ids']) for target, entry in selection.items() if entry['ids'] is not None}
    if not chosen:
        return df_gtf
    pairs = set((target, transcript) for target, ids in chosen.items() for transcript in ids)
    selected = df_gtf.gene_name.isin(list(chosen))
    keep = [(gene, transcript) in pairs for gene, transcript in
            zip(df_gtf.gene_name[selected], df_gtf.transcript_id[selected])]
    mask = ~selected
    mask.loc[selected] = keep
    return df_gtf[mask]


def extract_exons_rca(df_gtf: pd.DataFrame, fa: Fasta,
                  genelist: pd.DataFrame, min_length: int=40, source: str = 'exon'
                  ) :
    """Extract the exons (CDS for source exon and CDS, see gene_segments) of each gene.

    A gene without such segments, or (source exon, CDS) with every segment
    shorter than ``min_length``, gets an empty list, and so no windows: one
    such target must not stop the extraction of the others.
    """
    gene2exons = {}
    for item in genelist.iterrows():
        gene = item[1]['geneID']
        df_gene = df_gtf[df_gtf.gene_name == gene]
        if df_gene.shape[0] <= 0:
            raise ValueError(f"Gene {gene} not exists in GTF file.")
        df_exons = gene_segments(df_gene, source).copy()
        segments = df_exons.shape[0]
        if source != 'mRNA':
            df_exons = df_exons[df_exons.length >= min_length]
        gene2exons[gene] = []
        if df_exons.shape[0] == 0:
            log.warning(_no_segment_message(gene, source, segments, min_length))
            continue
        df_exons['transcript_name'] = df_exons['info'].str.extract(r'transcript_id\s+"([^"]+)"')[0]
        exon_cnts = df_exons.groupby(
            by=['chr', 'start', 'end', 'length', 'strand'],
            as_index=False
        ).agg(
            count=('transcript_name', 'count'),  
            transcript_name=('transcript_name', lambda x: list(set(x))) 
        )
        for idx, row in exon_cnts.iterrows():
            chr_, start, end, strand, trans_name = str(row['chr']), row['start'], row['end'], row['strand'], row['transcript_name']
            exon_name = '_'.join([chr_, str(start), str(end), strand])
            n_trans = row['count']
            seq = _read_segment(fa, chr_, start, end, strand)
            exon = (exon_name, trans_name, seq, n_trans, (str(row['chr']), start, end, strand))
            gene2exons[gene].append(exon)
    return gene2exons


def _no_segment_message(gene: str, source: str, segments: int, min_length: int) -> str:
    """Why a gene target gets no windows from extract_exons_rca or extract_gene_features."""
    kind = {'UTR': 'UTR', 'mRNA': 'exon'}.get(source, 'CDS (or, without CDS, exon)')
    if not segments:
        return f"Gene {gene} has no {kind} features in the GTF; it gets no windows with source {source}."
    fixes = ("source: mRNA or a shorter window" if source == 'UTR' else
             "source: mRNA, one transcript per target (transcripts: mane or canonical, tiled across its exon "
             "junctions) or a shorter window")
    return (f"Gene {gene}: all its {segments} {kind} segment(s) are shorter than the {min_length}-nt window, "
            f"so it gets no windows with source {source}. Use {fixes}.")


def _read_segment(fa: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    """The sequence of a GTF interval (1-based, inclusive), 5' to 3' on its strand.

    The GTF chromosome is looked up in the FASTA as written, else without or
    with a 'chr' prefix.
    """
    chrom = str(chrom)
    if chrom not in fa.keys():
        for alias in (chrom.replace('chr', ''), 'chr' + chrom):
            if alias in fa.keys():
                chrom = alias
                break
    try:
        seq = fa[chrom][start - 1:end].seq.upper()
    except Exception:
        log.error(f"Cannot read {chrom}:{start}-{end} from FASTA (first keys: {list(fa.keys())[:5]})")
        raise
    return reverse_complement(seq) if strand == '-' else seq


SplicedTranscript = t.Tuple[str, t.List[t.Tuple[str, str]]]  # (transcript id, [(segment name, sequence)] 5' to 3')


def extract_spliced(df_gtf: pd.DataFrame, fa: Fasta, gene: str, source: str) -> SplicedTranscript:
    """The segments of a target's one selected transcript, 5' to 3' (see SPLICED_SOURCES).

    ``df_gtf`` holds only that transcript of the gene (restrict_to_transcripts).
    Segments are not filtered by length: all of them make up the spliced sequence.
    """
    df_gene = df_gtf[df_gtf.gene_name == gene]
    if df_gene.shape[0] <= 0:
        raise ValueError(f"Gene {gene} not exists in GTF file.")
    rows = gene_segments(df_gene, source)
    if rows.shape[0] == 0:
        raise ValueError(f"Gene {gene} can't found any exon records.")
    transcripts = sorted(set(rows.transcript_id.dropna()))
    intervals = sorted(set(zip(rows.chr.astype(str), rows.start, rows.end, rows.strand)), key=lambda iv: (iv[1], iv[2]))
    if len(transcripts) != 1 or len({(chrom, strand) for chrom, _, _, strand in intervals}) != 1:
        raise ValueError(f"Gene {gene}: expected the segments of one transcript on one strand, got "
                         f"{len(transcripts)} transcript(s) on {sorted({(c, s) for c, _, _, s in intervals})}")
    if intervals[0][3] == '-':
        intervals.reverse()
    segments = [('_'.join([chrom, str(start), str(end), strand]), _read_segment(fa, chrom, start, end, strand))
                for chrom, start, end, strand in intervals]
    return transcripts[0], segments

def get_exon_seq(genes, fa, gtf, min_length=40, source='exon'):
    genelist = pd.DataFrame(genes, columns=['geneID'])
    df_gtf = read_gtf(gtf, extract_fields=['gene_name', "transcript_id"], get_length=True)
    with Fasta(fa) as fasta:
        return extract_exons_rca(df_gtf, fasta, genelist, min_length, source)

def change_chrom_name(chrom):
    if chrom.startswith('chr'):
        return chrom[3:]
    else:
        return 'chr'+chrom
    
Interval = t.Tuple[str, int, int, str]  # (chrom, start, end, strand), 1-based inclusive
Exon = t.Tuple[str, str, str, int, Interval]  # (type, name, seq, n_trans, interval)
Utr = t.Tuple[str, str, str, int, Interval]  # (type, name, seq, n_trans, interval)


def segment_offsets(intervals: t.Iterable[Interval]) -> t.Dict[Interval, int]:
    """Offset of each segment of one gene on its merged exon model.

    The model is the union of the gene's extracted segments (exons, CDS or
    UTRs) with the introns removed, read 5' to 3'; for a gene with one
    transcript it is the spliced RNA. Windows are numbered 5' to 3' within
    their segment, so window s..e of a segment lies at offset+s..offset+e.
    Overlapping segments of different transcripts share model positions.
    """
    groups = {}
    for interval in intervals:
        groups.setdefault((interval[0], interval[3]), set()).add(interval)
    offsets, position = {}, 0
    for (_, strand), group in sorted(groups.items()):
        blocks = []  # [start, end, segments] of the merged intervals
        for interval in sorted(group, key=lambda iv: (iv[1], iv[2])):
            if blocks and interval[1] <= blocks[-1][1]:
                blocks[-1][1] = max(blocks[-1][1], interval[2])
                blocks[-1][2].append(interval)
            else:
                blocks.append([interval[1], interval[2], [interval]])
        if strand == '-':
            blocks.reverse()
        for start, end, segments in blocks:
            for interval in segments:
                # Bases between the 5' end of the block and that of the segment.
                offsets[interval] = position + (end - interval[2] if strand == '-' else interval[1] - start)
            position += end - start + 1
    return offsets

def extract_gene_features(df_gtf: pd.DataFrame, fa: Fasta,
                          genelist: pd.DataFrame, min_length: int=40
                          ) -> t.Mapping[str, t.List[t.Tuple[str, str]]]:
    """Extract all exons and UTRs of each gene (source UTR windows read only the UTRs).

    A gene without UTRs, or with every UTR shorter than ``min_length``, gets
    no windows; it does not stop the extraction of the other targets.
    """
    gene_features = {}
    for item in genelist.iterrows():
        gene = item[1]['geneID']
        df_gene = df_gtf[df_gtf.gene_name == gene]
        if df_gene.shape[0] <= 0:
            raise ValueError(f"Gene {gene} not exists in GTF file.")
        # extract CDS/exon (not windowed; source UTR reads the UTRs below)
        df_exons = df_gene[df_gene.type.isin(['CDS', 'exon'])].copy()
        df_exons = df_exons[df_exons.length >= min_length]
        gene_features[gene] = []
        df_exons['transcript_name'] = df_exons['info'].str.extract(r'transcript_id\s+"([^"]+)"')[0]
        df_exons = df_exons.groupby(
            by=['chr', 'start', 'end', 'length', 'strand'],
            as_index=False
        ).agg(
            count=('transcript_name', 'count'),  
            transcript_name=('transcript_name', lambda x: list(set(x))) 
        ) if df_exons.shape[0] else df_exons
        for idx, row in df_exons.iterrows():
            chr_, start, end, strand = str(row['chr']), row['start'], row['end'], row['strand']
            name = f"{gene}_{start}_{end}"
            n_trans = row['count']
            if chr_ not in fa.keys():
                chr_fixed = chr_.replace('chr', '')
                if chr_fixed in fa.keys():
                    chr_ = chr_fixed
                else:
                    chr_fixed = 'chr' + chr_
                    if chr_fixed in fa.keys():
                        chr_ = chr_fixed
            seq = fa[chr_][start - 1:end].seq.upper()
            if strand == '-':
                seq = reverse_complement(seq)
            gene_features[gene].append(("exon", name, seq, n_trans, (str(row['chr']), start, end, strand)))
        # extract utr (GENCODE: UTR; Ensembl: five_prime_utr, three_prime_utr)
        df_utrs = gene_segments(df_gene, 'UTR').copy()
        segments = df_utrs.shape[0]
        df_utrs = df_utrs[df_utrs.length >= min_length]
        if df_utrs.shape[0] == 0:
            log.warning(_no_segment_message(gene, 'UTR', segments, min_length))
        if df_utrs.shape[0] > 0:
            for idx, row in df_utrs.iterrows():
                chr_, start, end, strand = str(row['chr']), row['start'], row['end'], row['strand']
                name = f"{gene}_{start}_{end}"
                n_trans = 1
                if chr_ not in fa.keys():
                    chr_fixed = chr_.replace('chr', '')
                    if chr_fixed in fa.keys():
                        chr_ = chr_fixed
                    else:
                        chr_fixed = 'chr' + chr_
                        if chr_fixed in fa.keys():
                            chr_ = chr_fixed
                seq = fa[chr_][start - 1:end].seq.upper()
                if strand == '-':
                    seq = reverse_complement(seq) 
                if idx == 0:
                    utr_type = "5'UTR" if strand == '+' else "3'UTR"
                elif idx == df_utrs.shape[0] - 1:
                    utr_type = "3'UTR" if strand == '+' else "5'UTR"
                else:
                    utr_type = "unknown_UTR"
                gene_features[gene].append((utr_type, name, seq, n_trans, (str(row['chr']), start, end, strand)))
    return gene_features

def extract_trans_seqs(gtf_path, fa_path, output_fa_path):
    """
    Extract transcript's sequences.
    """
    log.info(f"extract transcript sequences from: {gtf_path}, {fa_path}")
    fa = Fasta(str(fa_path))
    try:
        # Match fisheye: identify transcripts by gene_id and transcript_id,
        # and write FASTA headers as >{gene_id}_{transcript_id}.
        exons_df = read_gtf(gtf_path, filter_by_type='exon', extract_fields=["gene_id", "transcript_id"])
        # NCBI primary chromosome names such as NC_000067.7 contain underscores.
        # Keep references actually present in FASTA, accepting chr aliases.
        fasta_keys = set(fa.keys())
        def resolve_chrom(chrom):
            chrom = str(chrom)
            alias = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            return next((name for name in (chrom, alias) if name in fasta_keys), None)
        resolved = exons_df["chr"].map(resolve_chrom)
        unmatched = resolved.isna()
        if unmatched.any():
            log.warning(f"Skipping {int(unmatched.sum())} exons on references absent from FASTA")
        exons_df = exons_df.loc[~unmatched].copy()
        exons_df["chr"] = resolved.loc[~unmatched]
        exons_df = exons_df[exons_df.start <= exons_df.end]
        for field in ("gene_id", "transcript_id"):
            invalid = exons_df[field].isna() | exons_df[field].astype(str).str.strip().eq("")
            if invalid.any():
                raise ValueError(f"Missing {field} for {int(invalid.sum())} matched exons in {gtf_path}")
        if exons_df.empty:
            raise ValueError(f"No usable transcript exons: check GTF/FASTA reference names in {gtf_path}, {fa_path}")
        exons_df = exons_df[['chr','start','end','strand','gene_id','transcript_id']]
        trans = {}  # (gene_id, trans_id) -> [chr, strand, exons],  exons: (start, end)
        for (_, row) in exons_df.iterrows():
            key_ = (row['gene_id'], row['transcript_id'])
            chrom, strand, left, right = str(row['chr']), row['strand'], row['start'], row['end']
            if key_ not in trans:
                trans[key_] = [chrom, strand, [[left, right]]]
            else:
                trans[key_][2].append([left, right])
        # Splice annotated exons exactly: even a short gap is not mature RNA.
        for chrom, strand, exons in trans.values():
            exons.sort()
        seq_dict = {}
        fasta_keys = set(fa.keys())
        for key_, [chrom, strand, exons] in list(trans.items()):
            # fisheye accepts both GTF/FASTA naming conventions (1 vs chr1).
            chrom = str(chrom)
            chrom_candidates = [chrom]
            if chrom.startswith('chr'):
                chrom_candidates.append(chrom[3:])
            else:
                chrom_candidates.append('chr' + chrom)
            fasta_chrom = next((candidate for candidate in chrom_candidates
                                if candidate in fasta_keys), None)
            if fasta_chrom is None:
                raise KeyError(
                    f"Chromosome {chrom!r} from GTF was not found in FASTA "
                    f"(also tried {chrom_candidates[1]!r})"
                )
            seq_lst = []
            for i in range(len(exons)):
                seq = fa[fasta_chrom][exons[i][0] - 1:exons[i][1]].seq
                if strand == '-':
                    seq = reverse_complement(seq)
                    seq_lst.append(seq)
                else:
                    seq_lst.append(seq)
            if strand == "-":
                seq_lst.reverse()
            seq = "".join(seq_lst)
            seq_dict[key_] = seq
        if not seq_dict or any(not seq for seq in seq_dict.values()):
            raise ValueError("Transcript extraction produced no sequences or an empty sequence")
        # Publish only a complete FASTA; failed extraction must not leave an empty file.
        output_dir = os.path.dirname(os.path.abspath(output_fa_path))
        fd, temp_path = tempfile.mkstemp(prefix=".transcript-", suffix=".fa", dir=output_dir)
        try:
            with os.fdopen(fd, 'w') as f:
                for (gene_id, tran_id), seq in seq_dict.items():
                    f.write(f">{gene_id}_{tran_id}\n{seq}\n")
            os.replace(temp_path, output_fa_path)
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        log.info(f"Saved {len(seq_dict)} transcripts to {output_fa_path}")
    finally:
        fa.close()

def generate_target_seqs(
                        source,
                        targets, 
                        fasta_path, 
                        gtf_path, 
                        min_length: int = 40, 
                        overlap: t.Optional[int] = None,
                        step: t.Optional[int] = None,
                        transcripts='all',
                         ):
    """Candidate windows of the targets; RNA_SOURCES describes what each RNA source extracts.

    For RNA sources, ``transcripts`` (see TRANSCRIPT_RULES) chooses the
    transcripts whose features are extracted; the windows' attrs['transcripts']
    record them per target (select_transcripts, with the ids of all
    transcripts that have such features for 'all'). A target with one
    selected transcript is tiled along its spliced sequence (SPLICED_SOURCES).
    """
    lower, _, step = normalize_sampling(min_length, step, overlap)
    if source in RNA_SOURCES:
        df_gtf = read_gtf(gtf_path, extract_fields=['gene_name', 'transcript_id'])
        selection = select_transcripts(df_gtf, targets, transcripts)
        df_gtf = restrict_to_transcripts(df_gtf, selection)
        single = [target for target, entry in selection.items() if tiles_spliced(source, entry)]
        genelist = pd.DataFrame([target for target in targets if target not in single], columns=['geneID'])
        with Fasta(fasta_path) as fa:
            if source == 'UTR':
                info = extract_gene_features(df_gtf, fa, genelist, lower)
            else:
                info = extract_exons_rca(df_gtf, fa, genelist, lower, source)
            spliced = {target: extract_spliced(df_gtf, fa, target, source) for target in single}
        data = _rna_windows(source, info, min_length, step, spliced, order=targets)
        for target, entry in selection.items():
            if entry['ids'] is None:
                rows = gene_segments(df_gtf[df_gtf.gene_name == target], source)
                entry['ids'] = sorted(set(rows.transcript_id.dropna()))
        data.attrs['transcripts'] = selection
        return data
    elif source == 'genome':
        data_list = []
        for target in targets:
            seq_list = extract_fasta(fasta_path, target, min_length, step=step)
            data_list.extend(seq_list)
        data = pd.DataFrame(data_list, columns=['probe_id', 'target', 'sub_region', 'target_region'])
        return data


def _spliced_windows(gene_name: str, transcript: SplicedTranscript, min_length, step) -> list:
    """Window rows (as _rna_windows writes them) along one spliced transcript; see SPLICED_SOURCES."""
    transcript_id, segments = transcript
    sequence = ''.join(seq for _, seq in segments)
    ends = list(itertools.accumulate(len(seq) for _, seq in segments))  # 0-based exclusive end of each segment
    rows = []
    for n, (i, window_end, tem) in enumerate(iter_target_windows(sequence, min_length, step), start=1):
        first, last = bisect.bisect_right(ends, i), bisect.bisect_right(ends, window_end - 1)
        start = i + 1 - (ends[first - 1] if first else 0)
        end = start + window_end - i - 1
        exon_name = SEGMENT_SEPARATOR.join(name for name, _ in segments[first:last + 1])
        rows.append([f"{gene_name}_{n}", gene_name, f"{start}_{end}", exon_name, [transcript_id], start, end,
                     i + 1, window_end, tem, 1])
    return rows


def _rna_windows(source, info, min_length, step, spliced=None, order=None) -> pd.DataFrame:
    """Windows of the segments extract_exons_rca (UTR: extract_gene_features) returned.

    ``spliced`` maps targets with one transcript to extract_spliced results,
    tiled along the spliced sequence; targets are written in ``order``.
    """
    # start/end count within the exon or UTR; target_start/target_end count
    # along the target (its merged exon model, see segment_offsets, or its
    # spliced transcript, see SPLICED_SOURCES).
    if source != 'UTR':
        exon_info = info
        rows = {}
        for gene_name, exon_list in exon_info.items():
            data_list = rows[gene_name] = []
            offsets = segment_offsets(exon[4] for exon in exon_list)
            n = 1
            for j, exon_data in enumerate(exon_list, start=1):
                exon_name, trans_name, seq, n_trans, interval = exon_data
                for i, window_end, tem in iter_target_windows(seq, min_length, step):
                    start = i + 1
                    end = window_end
                    gene_id = f"{gene_name}_{n}"
                    sub_region = f"{start}_{end}"
                    n += 1
                    data_list.append([gene_id, gene_name, sub_region, exon_name, trans_name, start, end,
                                      offsets[interval] + start, offsets[interval] + end, tem, n_trans])
        for gene_name, transcript in (spliced or {}).items():
            rows[gene_name] = _spliced_windows(gene_name, transcript, min_length, step)
        data_list = [row for gene_name in dict.fromkeys([*(order or []), *rows])
                     for row in rows.get(gene_name, [])]
        data = pd.DataFrame(data_list, columns=['probe_id', 'target', 'sub_region','exon_name', 'transcript_names','start', 
                                                'end', 'target_start', 'target_end', 'target_region', 'n_trans'])
        return data
    else:
        utr_info = info
        data_list = []
        for gene_name, utr_list in utr_info.items():
            utr_list = [utr for utr in utr_list if utr[0] != 'exon']
            offsets = segment_offsets(utr[4] for utr in utr_list)
            n = 1
            for j, utr_data in enumerate(utr_list, start=1):
                utr_type, utr_name, seq, n_trans, interval = utr_data
                trans_name = []
                # extract target region seqs
                for i, window_end, tem in iter_target_windows(seq, min_length, step):
                    start = i + 1
                    end = window_end
                    gene_id = f"{gene_name}_{n}"
                    sub_region = f"{start}_{end}"
                    n += 1
                    data_list.append([gene_id, gene_name, sub_region, utr_name, trans_name, start, end,
                                      offsets[interval] + start, offsets[interval] + end, tem, n_trans])
        data = pd.DataFrame(data_list, columns=['probe_id', 'target', 'sub_region', 'utr_name', 'transcript_names','start', 
                                                'end', 'target_start', 'target_end', 'target_region', 'n_trans'])
        return data

def validate_targets(targets, gtf_path, DTF_NAME_FIX=False):
    log.info(f"validating targets in gtf file: {targets}")
    df_gtf = read_gtf(gtf_path)
    # DTF_NAME_FIX is retained for caller compatibility; normalization now
    # happens in read_gtf without rewriting the input annotation.
    genome_genes = set(df_gtf["gene_name"].dropna().unique())
    valid_targets = []
    invalid_targets = []
    for target in targets:
        if target in genome_genes:
            valid_targets.append(target)
        else:
            invalid_targets.append(target)
    if invalid_targets:
        log.warning(f"invalid targets: {invalid_targets}")
        return False, valid_targets, invalid_targets
    else:
        return True, valid_targets, invalid_targets

