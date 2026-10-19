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

    This intentionally mirrors fisheye: no MAPQ filter is applied.  The
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
    # Keep this extraction identical to fisheye/primer_design/seq_features.py.
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
    follow fisheye's ``{gene_id}_{transcript_id}`` convention.  Multiple
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
    """Count distinct genes hit by each query, matching fisheye semantics."""
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

def _bowtie_hits(outdir, recname2seq, name, index_prefix, max_mismatches, threads, strand, key):
    """query -> {key(transcript.fa record): fewest mismatches} of Bowtie 1's alignments."""
    fa_path = f"{outdir}/{name}.fa"
    with open(fa_path, "w") as f:
        for rec_name, seq in recname2seq.items():
            f.write(f">{rec_name}\n{seq}\n")
    cmd = ["bowtie", "-f", "-v", str(max_mismatches), "-a", "-p", str(threads),
           "--suppress", "5,6,7"] + STRAND_OPTIONS[strand] + [str(index_prefix), fa_path]
    log_file = f"{outdir}/{name}.bowtie.log"
    log.info(f"Call cmd: {' '.join(cmd)} 2> {log_file}")
    hits = {rec_name: {} for rec_name in recname2seq}
    with open(log_file, "w") as err:
        # Stream the output: with -a, repeat-like queries can have many alignments.
        proc = subp.Popen(cmd, stdout=subp.PIPE, stderr=err, text=True)
        for line in proc.stdout:
            # query, strand, reference, offset, mismatches (e.g. "3:G>T,10:A>C")
            fields = line.rstrip("\n").split("\t")
            mismatches = len(fields[4].split(",")) if len(fields) > 4 and fields[4] else 0
            found = hits[fields[0]]
            hit = key(fields[2])
            if mismatches < found.get(hit, max_mismatches + 1):
                found[hit] = mismatches
        returncode = proc.wait()
    if returncode != 0:
        with open(log_file) as f:
            log.error(f.read())
        raise subp.CalledProcessError(returncode, cmd)
    return hits

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
    """Defaults are the reaction conditions matching the fisheye reference output.

    Protocol ``hybridization`` conditions override them (mM, mM, mM, nM and
    mol/l formamide); see uprobe.core.hybridization.
    """
    return primer3.calc_tm(
        seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc,
        dna_conc=dna_conc, formamide_conc=formamide_conc, dmso_conc=0,
        tm_method="santalucia", salt_corrections_method="santalucia",
    )

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
