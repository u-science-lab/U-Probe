import math
from collections import defaultdict
import os
import typing as t
import primer3
import RNA
import subprocess as subp
from uprobe.core.utils import get_logger, reverse_complement, write_fastq
from uprobe.core.offtarget import STRAND_OPTIONS

log = get_logger(__name__)

Aln = t.Tuple[str, int, int]  # chr, start, end
Block = t.Tuple[str, str, t.List[Aln]]  # query_name, query_seq, alignments

def read_sam_align_blocks(
        sam_path: str,
        ) -> t.Iterable[Block]:
    """Yield all Bowtie2 alignments grouped by query name.

    No MAPQ filter is applied.  The
    mapped-gene metric is the number of genes among all alignments reported by
    Bowtie2 (up to the ``-k`` limit), not the number of uniquely mapped genes.
    """
    import pysam
    def yield_cond(old, rec, block, end=False):
        res = (old is not None)
        if res and not end:
            res &= rec.query_name != old.query_name
        return res
    with pysam.AlignmentFile(sam_path, mode='r') as sam:
        alns = []
        old = None
        rec = None
        for rec in sam.fetch():
            aln = rec.reference_name, rec.reference_start, rec.reference_end
            if yield_cond(old, rec, alns):
                yield old.query_name, old.query_sequence, alns
                alns = []
            if aln[0] is not None:
                alns.append(aln)
            old = rec
        if (rec is not None) and yield_cond(old, rec, alns, end=True):
            yield old.query_name, old.query_sequence, alns

def bowtie2_align_se_sen(
        fq_path: str,
        index: str,
        sam_path: str,
        threads: int = 10,
        log_file: t.Optional[str] = 'bowtie2.log',
        header: bool = True,
        strand: str = 'both',
        ) -> str:
    cmd = ["bowtie2", "-x", str(index), "-U", str(fq_path)]
    if not header:
        cmd.append("--no-hd")
    cmd += ["-t", "-k", "100", "--very-sensitive-local"] + STRAND_OPTIONS[strand]
    cmd += ["-p", str(threads)]
    cmd += ["-S", str(sam_path)]
    cmd_str = " ".join(cmd)
    if log_file:
        cmd_str += f" > {log_file} 2>&1"
    log.info(f"Call cmd: {cmd_str}")
    try:
        subp.check_call(cmd_str, shell=True)
    except subp.CalledProcessError as e:
        log.error(f"Error: {e}")
        if log_file:
            # print log file
            with open(log_file) as f:
                log.error(f.read())
        raise
    return sam_path

def parse_cigar(cigar: str) -> int:
    if cigar == "*":  
        return 0
    import re
    ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    ref_len = 0
    for length, op in ops:
        length = int(length)
        if op in "M=XDN":  
            ref_len += length
    return ref_len

def cal_mapped_sites(outdir: str, 
                target: str,
                recname2seq: t.Mapping[str, str],
                index_prefix: str, 
                threads: int = 10,
                ) -> t.Dict[str, t.List[t.Tuple[str, int, int]]]:
    """
    return:
        mapped_sites_dict: t.Dict[str, t.List[t.Tuple[str, int, int]]]
            - key: sequence name
            - value: list of (ref_name, start_pos, end_pos) tuples
    """
    fq_path = write_fastq(outdir, target, recname2seq)
    sam_path = f"{outdir}/{target}.mapped_sites.sam"
    if not os.path.exists(sam_path):
        bowtie2_align_se_sen(
            fq_path, index_prefix,
            sam_path, threads=threads,
            log_file=f"{outdir}/{target}.bowtie2.log")
    mapped_sites_dict: t.Dict[str, t.List[t.Tuple[str, int, int]]] = {}
    for seq_name in recname2seq.keys():
        mapped_sites_dict[seq_name] = []
    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):  
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:  
                continue
            flag = int(parts[1])
            if flag & 0x4:  
                continue
            query_name = parts[0]  
            ref_name = parts[2]  
            start_pos = int(parts[3])  
            cigar = parts[5] 
            ref_len = parse_cigar(cigar)
            end_pos = start_pos + ref_len - 1  
            if query_name in mapped_sites_dict:
                mapped_sites_dict[query_name].append((ref_name, start_pos, end_pos))
    os.remove(fq_path)
    return mapped_sites_dict

def cal_kmer_count(outdir: str,
                   target: str,
                   recname2seq: t.Mapping[str, str],
                   index_prefix: str,
                   kmer_len: int = 35,
                   threads: int = 10,
                ) -> t.Dict[str, int]:
    import numpy as np
    randomInt = np.random.randint(0, 1000000)
    temp_fasta = f"{outdir}/{target}_{kmer_len}_{randomInt}_temp.fa"
    with open(temp_fasta, 'w') as f:
        for rec_name, seq in recname2seq.items():
            f.write(f">{rec_name}\n{seq}\n")
    temp_output = f"{outdir}/{target}_{kmer_len}_{randomInt}_temp.txt"
    cmd = [
        'jellyfish', 'query', f'{index_prefix}',
        '-s', temp_fasta,
        '-o', temp_output
    ]
    try:
        subp.check_call(cmd, stderr=None, shell=False)
    except subp.CalledProcessError as e:
        log.error(f"Jellyfish query failed: {e}")
        raise
    with open(temp_output, 'r') as f:
        jf_lines = [line.strip() for line in f]
    seq_kmer_counts = {}
    line_idx = 0
    for rec_name, seq in recname2seq.items():
        seq_length = len(seq)
        num_kmers = seq_length - kmer_len + 1
        
        if num_kmers <= 0:
            seq_kmer_counts[rec_name] = 0
            continue
        max_count = 0
        for i in range(num_kmers):
            if line_idx < len(jf_lines):
                count_val = int(jf_lines[line_idx].split(' ')[1])
                max_count = max(max_count, count_val)
                line_idx += 1
        
        seq_kmer_counts[rec_name] = max_count
    try:
        os.remove(temp_fasta)
        os.remove(temp_output)
    except OSError:
        pass
    return seq_kmer_counts 

def transcript_gene(ref_name: str) -> str:
    """Gene of a ``transcript.fa`` record named ``{gene_id}_{transcript_id}``."""
    return ref_name.split("_")[0]

def gene_hits_of(transcript_hits: t.Mapping[str, t.Optional[int]]) -> t.Dict[str, t.Optional[int]]:
    """{gene: fewest mismatches} of one query's {transcript.fa record: fewest mismatches} (None for Bowtie2)."""
    genes = {}
    for ref_name, mismatches in transcript_hits.items():
        gene = transcript_gene(ref_name)
        if gene not in genes or (mismatches is not None and mismatches < genes[gene]):
            genes[gene] = mismatches
    return genes

def bowtie2_gene_hits(
        outdir: str,
        recname2seq: t.Mapping[str, str],
        name: str,
        index_prefix: str,
        threads: int = 10,
        strand: str = 'both',
        ) -> t.Dict[str, t.Dict[str, None]]:
    """Genes hit by each query among Bowtie2's local alignments (up to ``-k 100``).

    The Bowtie2 index must be built from ``transcript.fa`` whose record names
    follow the ``{gene_id}_{transcript_id}`` convention.  Multiple
    transcript hits from the same gene therefore contribute one mapped gene.
    Returns query -> {gene: None}: local alignments may be clipped, so they
    carry no mismatch count comparable to :func:`bowtie_gene_hits`.
    """
    return _bowtie2_hits(outdir, recname2seq, name, index_prefix, threads, strand, transcript_gene)

def bowtie2_transcript_hits(
        outdir: str,
        recname2seq: t.Mapping[str, str],
        name: str,
        index_prefix: str,
        threads: int = 10,
        strand: str = 'both',
        ) -> t.Dict[str, t.Dict[str, None]]:
    """Transcripts hit by each query, as :func:`bowtie2_gene_hits`: query -> {transcript.fa record: None}."""
    return _bowtie2_hits(outdir, recname2seq, name, index_prefix, threads, strand, lambda ref_name: ref_name)

def _bowtie2_hits(outdir, recname2seq, name, index_prefix, threads, strand, key):
    """query -> {key(transcript.fa record): None} of Bowtie2's local alignments."""
    fq_path = write_fastq(outdir, name, recname2seq)
    sam_path = f"{outdir}/{name}.sam"
    if not os.path.exists(sam_path):
        bowtie2_align_se_sen(
            fq_path, index_prefix,
            sam_path, threads=threads,
            log_file=f"{outdir}/{name}.bowtie2.log", strand=strand)
    # Include queries without an alignment explicitly.
    hits = {rec_name: {} for rec_name in recname2seq}
    for rec_name, seq, alns in read_sam_align_blocks(sam_path):
        hits[rec_name] = {key(ref_name): None for ref_name, start, end in alns}
    return hits

def count_n_bowtie2_aligned_genes(
        outdir: str,
        recname2seq: t.Mapping[str, str],
        name: str,
        index_prefix: str,
        threads: int = 10):
    """Count distinct genes hit by each query."""
    hits = bowtie2_gene_hits(outdir, recname2seq, name, index_prefix, threads)
    return {rec_name: len(genes) for rec_name, genes in hits.items()}

def bowtie_gene_hits(
        outdir: str,
        recname2seq: t.Mapping[str, str],
        name: str,
        index_prefix: str,
        max_mismatches: int = 2,
        threads: int = 10,
        strand: str = 'both',
        ) -> t.Dict[str, t.Dict[str, int]]:
    """Genes with a site within ``max_mismatches`` substitutions of each query.

    Bowtie 1 in ``-v`` mode ignores base qualities, and with ``-a`` it reports
    every end-to-end alignment with at most ``max_mismatches`` (0-3)
    mismatches, so no such site is missed; indels are not considered.
    Returns query -> {gene: fewest mismatches over the gene's sites}.
    """
    return _bowtie_hits(outdir, recname2seq, name, index_prefix, max_mismatches, threads, strand, transcript_gene)

def bowtie_transcript_hits(
        outdir: str,
        recname2seq: t.Mapping[str, str],
        name: str,
        index_prefix: str,
        max_mismatches: int = 2,
        threads: int = 10,
        strand: str = 'both',
        ) -> t.Dict[str, t.Dict[str, int]]:
    """Transcripts with a site within ``max_mismatches`` substitutions of each query, as :func:`bowtie_gene_hits`.

    Returns query -> {transcript.fa record: fewest mismatches over its sites}.
    """
    return _bowtie_hits(outdir, recname2seq, name, index_prefix, max_mismatches, threads, strand,
                        lambda ref_name: ref_name)

def _bowtie_alignments(outdir, recname2seq, name, index_prefix, max_mismatches, threads, strand):
    """Yield (query, strand, transcript.fa record, 0-based offset, mismatches) of Bowtie 1's alignments."""
    fa_path = f"{outdir}/{name}.fa"
    with open(fa_path, "w") as f:
        for rec_name, seq in recname2seq.items():
            f.write(f">{rec_name}\n{seq}\n")
    cmd = ["bowtie", "-f", "-v", str(max_mismatches), "-a", "-p", str(threads),
           "--suppress", "5,6,7"] + STRAND_OPTIONS[strand] + [str(index_prefix), fa_path]
    log_file = f"{outdir}/{name}.bowtie.log"
    log.info(f"Call cmd: {' '.join(cmd)} 2> {log_file}")
    with open(log_file, "w") as err:
        # Stream the output: with -a, repeat-like queries can have many alignments.
        proc = subp.Popen(cmd, stdout=subp.PIPE, stderr=err, text=True)
        for line in proc.stdout:
            # query, strand, reference, offset, mismatches (e.g. "3:G>T,10:A>C")
            fields = line.rstrip("\n").split("\t")
            mismatches = len(fields[4].split(",")) if len(fields) > 4 and fields[4] else 0
            yield fields[0], fields[1], fields[2], int(fields[3]), mismatches
        returncode = proc.wait()
    if returncode != 0:
        with open(log_file) as f:
            log.error(f.read())
        raise subp.CalledProcessError(returncode, cmd)

def _bowtie_hits(outdir, recname2seq, name, index_prefix, max_mismatches, threads, strand, key):
    """query -> {key(transcript.fa record): fewest mismatches} of Bowtie 1's alignments."""
    hits = {rec_name: {} for rec_name in recname2seq}
    for query, _, ref_name, _, mismatches in _bowtie_alignments(outdir, recname2seq, name, index_prefix,
                                                                max_mismatches, threads, strand):
        found = hits[query]
        hit = key(ref_name)
        if mismatches < found.get(hit, max_mismatches + 1):
            found[hit] = mismatches
    return hits

def bowtie_unit_hits(
        outdir: str,
        units: t.Iterable[t.Tuple[str, ...]],
        name: str,
        index_prefix: str,
        max_mismatches: int = 2,
        threads: int = 10,
        strand: str = 'both',
        max_gap: t.Optional[int] = None,
        ) -> t.Dict[t.Tuple[str, ...], t.Dict[str, int]]:
    """Transcripts on which every oligo of a unit has a site within ``max_mismatches`` substitutions.

    A unit is a tuple of binding sequences, one per oligo, all in the same
    orientation (e.g. layout parts of target_region); each is searched on its
    own with Bowtie 1 (``-v max_mismatches -a``). The longest oligo of a unit
    (the first of equal length) is its anchor; the other oligos keep only sites
    on the transcripts their anchors hit, so that short oligos with thousands
    of sites stay cheap. With ``max_gap``, a transcript counts only if each
    other oligo has a site on the anchor site's strand that does not overlap it
    and lies at most ``max_gap`` nt from it; without, any sites on the
    transcript count. Returns unit -> {transcript.fa record: mismatches of the
    unit's worst oligo there (the fewest over its placements)}.
    """
    units = list(dict.fromkeys(tuple(unit) for unit in units))
    anchors = {unit: max(range(len(unit)), key=lambda i: (len(unit[i]), -i)) for unit in units}
    # sequence -> {record: [(strand, offset, mismatches)]}
    sites = {}

    def search(sequences, label, keep=None):
        queries = {str(i): seq for i, seq in enumerate(sequences)}
        for seq in sequences:
            sites[seq] = {}
        for query, aln_strand, ref_name, offset, mismatches in _bowtie_alignments(
                outdir, queries, f"{name}_{label}", index_prefix, max_mismatches, threads, strand):
            seq = queries[query]
            if keep is None or ref_name in keep[seq]:
                sites[seq].setdefault(ref_name, []).append((aln_strand, offset, mismatches))

    search(list(dict.fromkeys(unit[anchors[unit]] for unit in units)), "anchors")
    keep = defaultdict(set)
    for unit in units:
        refs = sites[unit[anchors[unit]]].keys()
        for i, seq in enumerate(unit):
            if i != anchors[unit] and seq not in sites:
                keep[seq].update(refs)
    others = [seq for seq in keep if keep[seq]]
    for seq in keep:
        sites[seq] = {}
    if others:
        search(others, "oligos", keep)

    result = {}
    for unit in units:
        anchor = anchors[unit]
        anchor_sites = sites[unit[anchor]]
        found = {}
        for ref_name, placements in anchor_sites.items():
            other_sites = [sites[seq].get(ref_name) for i, seq in enumerate(unit) if i != anchor]
            if not all(other_sites):
                continue
            if max_gap is None:
                found[ref_name] = max(min(m for _, _, m in placed) for placed in [placements] + other_sites)
                continue
            length = len(unit[anchor])
            others_lengths = [len(seq) for i, seq in enumerate(unit) if i != anchor]
            best = None
            for aln_strand, start, mismatches in placements:
                worst = mismatches
                for placed, other_length in zip(other_sites, others_lengths):
                    near = [m for s, p, m in placed if s == aln_strand
                            and 0 <= max(p - (start + length), start - (p + other_length)) <= max_gap]
                    if not near:
                        break
                    worst = max(worst, min(near))
                else:
                    best = worst if best is None else min(best, worst)
            if best is not None:
                found[ref_name] = best
        result[unit] = found
    return result

def self_match(probe: str, min_match = 4):
    length = len(probe)
    probe_re = reverse_complement(probe)
    match_pairs = 0
    for i in range(0,length-min_match+1):
        tem = probe[i:min_match+i]
        for j in range(0,length-min_match+1):
            tem_re = probe_re[j:min_match+j]
            if tem == tem_re and i + j + min_match - 2 != length:
                match_pairs = match_pairs + 1
    return match_pairs

def cal_temp(seq: str, mv_conc=50, dv_conc=0, dntp_conc=0.6, dna_conc=50,
             formamide_conc=0): # Tm
    """Defaults are the legacy reaction conditions: 50 mM Na+, no Mg2+, 0.6 mM dNTP, 50 nM oligo, no formamide.

    Protocol ``hybridization`` conditions override them (mM, mM, mM, nM and
    mol/l formamide); see uprobe.core.hybridization.
    """
    return primer3.calc_tm(
        seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
        dna_conc=dna_conc, formamide_conc=formamide_conc, dmso_conc=0,
        tm_method="santalucia", salt_corrections_method="santalucia",
    )

# Sugimoto et al. 1996 (Nucleic Acids Res 24:4501-4505, doi:10.1093/nar/24.22.4501)
# DNA/DNA nearest-neighbor enthalpies (kcal/mol) and entropies (cal/K/mol), as
# tabulated in Moldia/multi_padlock_design lib/screenseq.py.
SUGIMOTO_1996_NN = {
    'AA': (-8.0, -21.9), 'AC': (-9.4, -25.5), 'AG': (-6.6, -16.4), 'AT': (-5.6, -15.2),
    'CA': (-8.2, -21.0), 'CC': (-10.9, -28.4), 'CG': (-11.8, -29.0), 'CT': (-6.6, -16.4),
    'GA': (-8.8, -23.5), 'GC': (-10.5, -26.4), 'GG': (-10.9, -28.4), 'GT': (-9.4, -25.5),
    'TA': (-6.6, -18.4), 'TC': (-8.8, -23.5), 'TG': (-8.2, -21.0), 'TT': (-8.0, -21.9),
}
SUGIMOTO_1996_INITIATION = (0.6, -9.0)

def cal_temp_multi_padlock_design(seq: str, na_mM=50, mg_mM=0, probe_nM=50, formamide_pct=0):
    """Tm (°C) by the formula of the HybISS padlock design code (Gyllborg et al. 2020).

    Moldia/multi_padlock_design lib/screenseq.py calculatetm: Sugimoto 1996
    nearest-neighbor sums plus one initiation term, Tm = 1000 dH / (dS + R ln C)
    - 273.15 with C the probe concentration (mol/l), + 16.6 log10([Na+] +
    3.795 sqrt([Mg2+])) (mol/l), - 0.72 °C per % formamide. The repository
    fixes 75 mM Na+, 10 mM Mg2+, 100 nM probe and 20% formamide; here they are
    arguments. NaN for a sequence with bases other than A, C, G, T/U.
    """
    seq = seq.upper().replace('U', 'T')
    if not seq or set(seq) - set('ACGT'):
        return float('nan')
    enthalpy, entropy = SUGIMOTO_1996_INITIATION
    for i in range(len(seq) - 1):
        h, s = SUGIMOTO_1996_NN[seq[i:i + 2]]
        enthalpy += h
        entropy += s
    tm = enthalpy * 1000 / (entropy + 1.9872 * math.log(probe_nM * 1e-9)) - 273.15
    tm += 16.6 * math.log10(na_mM / 1000 + 3.795 * math.sqrt(mg_mM / 1000))
    return tm - 0.72 * formamide_pct

# SantaLucia 1998 unified nearest-neighbor parameters (PNAS 95:1460-1465,
# doi:10.1073/pnas.95.4.1460): dH (kcal/mol), dS (cal/K/mol) per dinucleotide,
# the table Biopython calls DNA_NN3, which OligoMiner uses.  Complementary
# dinucleotides (AA/TT, AC/GT, ...) share a value.
SANTALUCIA_1998_NN = {
    'AA': (-7.9, -22.2), 'AC': (-8.4, -22.4), 'AG': (-7.8, -21.0), 'AT': (-7.2, -20.4),
    'CA': (-8.5, -22.7), 'CC': (-8.0, -19.9), 'CG': (-10.6, -27.2), 'CT': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9), 'GT': (-8.4, -22.4),
    'TA': (-7.2, -21.3), 'TC': (-8.2, -22.2), 'TG': (-8.5, -22.7), 'TT': (-7.9, -22.2),
}
SANTALUCIA_1998_TERMINAL = {'AT': (2.3, 4.1), 'GC': (0.1, -2.8)}  # per end base pair
# The nearest-neighbor table of the MERFISH design code (ZhuangLab/MERFISH_analysis,
# TRDesigner.m SantaLuciaNearestNeighbor): SantaLucia & Hicks 2004 (Annu Rev Biophys Biomol
# Struct 33:415-440), which differs from the 1998 values only for AA/TT (-7.6, -21.3).
SANTALUCIA_2004_NN = {**SANTALUCIA_1998_NN, 'AA': (-7.6, -21.3), 'TT': (-7.6, -21.3)}


def _nearest_neighbor_sums(seq, table):
    return (sum(table[seq[i:i + 2]][0] for i in range(len(seq) - 1)),
            sum(table[seq[i:i + 2]][1] for i in range(len(seq) - 1)))


def _clean_dna(seq):
    """Upper-case DNA, or None unless the sequence is at least 2 nt of A, C, G, T/U."""
    seq = seq.upper().replace('U', 'T') if isinstance(seq, str) else ''
    return seq if len(seq) >= 2 and not set(seq) - set('ACGT') else None


def cal_temp_oligominer(seq: str, na_mM=50, probe_nM=50, formamide_pct=0):
    """Tm (°C) by the formula of OligoMiner (Beliveau et al. 2018, PNAS 115:E2183).

    Biopython ``Tm_NN(seq, nn_table=DNA_NN3, Na=na_mM, dnac1=probe_nM, dnac2=probe_nM)``
    (blockParse.py / probeTm.py): SantaLucia 1998 nearest-neighbor sums with its
    terminal A/T and G/C terms, the entropy salt correction ``0.368 (N-1) ln[Na+]``
    (``saltcorr=5``), Tm = 1000 dH / (dS + R ln k) - 273.15 with R = 1.987 and the
    strand concentration k = dnac1 - dnac2 / 2 (nM), rounded to 2 decimals as OligoMiner
    does, then - 0.65 °C per % formamide (``chem_correction``).  OligoMiner's defaults
    are 390 mM Na+, 25 nM oligo and target (-c, -C) and 50% formamide.  Divalent
    cations and dNTPs are not part of the formula.  NaN for a sequence with bases other
    than A, C, G, T/U.
    """
    seq = _clean_dna(seq)
    if seq is None:
        return float('nan')
    enthalpy, entropy = _nearest_neighbor_sums(seq, SANTALUCIA_1998_NN)
    for end in (seq[0], seq[-1]):
        h, s = SANTALUCIA_1998_TERMINAL['AT' if end in 'AT' else 'GC']
        enthalpy += h
        entropy += s
    entropy += 0.368 * (len(seq) - 1) * math.log(na_mM / 1000)
    strand_nM = probe_nM - probe_nM / 2
    tm = enthalpy * 1000 / (entropy + 1.987 * math.log(strand_nM * 1e-9)) - 273.15
    return round(tm, 2) - 0.65 * formamide_pct


def cal_temp_merfish_design(seq: str, na_mM=50, probe_nM=50):
    """Tm (°C) by the formula of the MERFISH probe design code (Moffitt et al. 2016, PNAS 113:11046).

    ZhuangLab/MERFISH_analysis probe_construction/TRDesigner.m GetRegionTm:
    SantaLucia & Hicks 2004 nearest-neighbor sums, dH + 0.2 + 2.2 per terminal
    A/T and dS - 5.7 + 6.9 per terminal A/T, dS + 0.368 (N-1) ln[Na+], and
    Tm = 1000 dH / (dS + 1.9872 ln C) - 273.15 with C the probe concentration
    itself, where Primer3 uses C/4 (the code's own note says it should be
    probe minus half the target concentration).  The published window
    (66-76 °C for the 30-nt MERFISH targets) is computed at 300 mM Na+ and 5 nM
    probe, with no formamide, Mg2+ or dNTP term.  NaN for a sequence with bases
    other than A, C, G, T/U.
    """
    seq = _clean_dna(seq)
    if seq is None:
        return float('nan')
    enthalpy, entropy = _nearest_neighbor_sums(seq, SANTALUCIA_2004_NN)
    for end in (seq[0], seq[-1]):
        if end in 'AT':
            enthalpy += 2.2
            entropy += 6.9
    enthalpy += 0.2
    entropy += -5.7 + 0.368 * (len(seq) - 1) * math.log(na_mM / 1000)
    return enthalpy * 1000 / (entropy + 1.9872 * math.log(probe_nM * 1e-9)) - 273.15

# RNA/DNA nearest-neighbor free energies dG37 (kcal/mol) of the Oligostan script (Tsanov et al. 2016,
# NAR 44:e165, doi:10.1093/nar/gkw784; Oligostan.r, ConvertRNASeq2DeltaGat37, in
# bitbucket.org/muellerflorian/fish_quant).  Key: the dinucleotide of the DNA probe read 5'->3'.
# The script gives no reference for the numbers (the paper cites Wu et al. 2002, Eur J Biochem
# 269:2821, doi:10.1046/j.1432-1033.2002.02970.x, for the temperature dependence).
OLIGOSTAN_DG37 = {
    'AA': -0.2, 'AC': -1.5, 'AG': -0.9, 'AT': -1.0, 'CA': -1.0, 'CC': -2.2, 'CG': -1.2, 'CT': -1.4,
    'GA': -0.8, 'GC': -2.4, 'GG': -1.5, 'GT': -1.0, 'TA': -0.3, 'TC': -1.4, 'TG': -1.0, 'TT': -0.4,
}
OLIGOSTAN_SALT_M = 0.115
SEQUENCE_OF = ('probe', 'target')


def target_sense(column) -> bool:
    """True for an attribute target that holds the target RNA in sense orientation: target_region or a layout part
    of the window (target_parts.<name>), not a probe or its parts (rna_dna_dg37 needs sequence: target there)."""
    column = str(column or '').replace(':', '.')
    return column == 'target_region' or column.startswith('target_parts.')


def cal_rna_dna_dg37(seq: str, sequence_of='probe', salt_M=OLIGOSTAN_SALT_M):
    """RNA/DNA hybrid dG37 (kcal/mol, more negative is more stable) as the Oligostan script computes it.

    The sum of the nearest-neighbor values of the probe's dinucleotides (probe read 5'->3', so
    ``sequence_of='target'`` first reverse-complements an RNA target given in sense
    orientation, as the script does) plus the salt term ``-(ln[Na+] * -0.175 - 0.2)``, a constant
    of 0.1785 kcal/mol at its 0.115 M. There is no initiation or terminal term, and a
    dinucleotide with a base other than A, C, G, T/U contributes 0, as in the script.
    """
    if not isinstance(seq, str) or len(seq) < 2:
        return float('nan')
    probe = seq.upper().replace('U', 'T')
    if sequence_of == 'target':
        probe = reverse_complement(probe)
    total = sum(OLIGOSTAN_DG37.get(probe[i:i + 2], 0.0) for i in range(len(probe) - 1))
    return total - (math.log(salt_M) * -0.175 - 0.2)

def cal_fold(seq: str): # RNA fold
    return -RNA.fold_compound(seq).mfe()[1]

def cal_gc_content(seq: str):
    return (seq.count('G') + seq.count('C')) / len(seq)

def cal_target_fold_score(seq: str, temperature_c=None):
    if not seq or seq is None:
        return float('nan')
    clean_seq = seq.upper().replace('T', 'U')
    valid_chars = set('AUGC')
    clean_seq = ''.join(c if c in valid_chars else 'A' for c in clean_seq)
    if not clean_seq:
        return float('nan')
    try:
        if temperature_c is None:
            return -RNA.fold_compound(clean_seq).mfe()[1]
        # Fold at the hybridization temperature instead of ViennaRNA's 37 °C.
        md = RNA.md()
        md.temperature = temperature_c
        return -RNA.fold_compound(clean_seq, md).mfe()[1]
    except Exception as e:
        print(f"Warning: RNA folding failed for sequence '{seq[:20]}...': {e}")
        return float('nan')

def cal_target_blocks(seq: str, offset: float):
    whole_fold = RNA.fold_compound(seq).mfe()
    target_fold = whole_fold[0][offset:offset+len(seq)]
    target_blocks = len(target_fold) - target_fold.count('.')  # smaller is better
    return target_blocks

def cal_self_match(seq: str):
    return self_match(seq)
