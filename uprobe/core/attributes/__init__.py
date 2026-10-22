import secrets
import os
import tempfile
from pathlib import Path
import pandas as pd
from ._attributes import *
from ..genome_paths import index_prefix as index_prefix_of, index_root
from ..hybridization import resolve_conditions, tm_arguments, tm_formula_arguments, tm_method
from ..offtarget import (BOWTIE2_LOCAL_MIN_LENGTH, STRAND_OPTIONS, OffTargetCounter, mapped_genes_settings,
                         short_for_bowtie2)
from ..utils import get_logger, reverse_complement

log = get_logger(__name__)


def locus_kmer_index(fasta, sequences, kmer_len, threads, size, canonical, directory):
    """A small Jellyfish index holding the genome counts of only the ``sequences``' k-mers.

    ``jellyfish count --if`` (Marcais & Kingsford 2011) reads the k-mers of a FASTA first and
    counts, over the genome, only those. Both strands of every sequence are listed, since a
    stranded index (strand: sense or antisense) holds each k-mer as written on the forward
    strand. The hash is sized for the k-mers listed unless ``size`` is given.
    """
    from uprobe.core.tools import build_jf_index
    unique = list(dict.fromkeys(seq for seq in sequences if isinstance(seq, str) and len(seq) >= kmer_len))
    locus_fasta = directory / "locus_kmers.fa"
    with open(locus_fasta, "w") as handle:
        for i, seq in enumerate(unique):
            handle.write(f">{i}\n{seq}\n>{i}r\n{reverse_complement(seq)}\n")
    kmers = sum(len(seq) - kmer_len + 1 for seq in unique) * 2
    index = directory / "locus_kmers.jf"
    log.info(f"Counting {kmers:,} candidate {kmer_len}-mers across {fasta} (jellyfish count --if)")
    build_jf_index(str(fasta), kmer_len, str(index), threads, size or str(max(1_000_000, 2 * kmers)),
                   canonical=canonical, only_kmers_in=str(locus_fasta))
    return index


TM_FORMULAS = {"multi_padlock_design": cal_temp_multi_padlock_design,
               "merfish_design": cal_temp_merfish_design,
               "oligominer": cal_temp_oligominer}

def add_attributes(
        df_probes: pd.DataFrame,
        protocol: dict,
        genome: dict,
        report: dict = None,
        threads: int = 10,
) -> pd.DataFrame:
    """Compute protocol['attributes']; ``report`` (a dict) receives off-target statistics.

    ``threads`` (``uprobe run --threads``) is used by every aligner and Jellyfish
    step of an attribute that sets no ``threads`` of its own.
    """
    # A private scratch directory per call: concurrent workflows sharing a
    # working directory must not share (and delete) each other's alignment files.
    with tempfile.TemporaryDirectory(prefix=".uprobe-tmp-", dir=".") as tmp_root:
        return _add_attributes(df_probes, protocol, genome, Path(tmp_root), report, threads)


def _add_attributes(
        df_probes: pd.DataFrame,
        protocol: dict,
        genome: dict,
        tmp_root: Path,
        report: dict = None,
        threads: int = 10,
) -> pd.DataFrame:
    attributes: dict = protocol['attributes']
    off_targets = None  # created by the first mapped_genes attribute
    def cached(series, function):
        values = {sequence: function(sequence) for sequence in series.unique()}
        return series.map(values)
    task_id = secrets.token_hex(6)  
    for attr_name, attr in attributes.items():
        if str(attr.get('type', '')).lower() == "difference":
            # |a - b| of two attributes computed before it, e.g. the Tm of two probe halves.
            first, second = attr['of']
            df_probes[attr_name] = (pd.to_numeric(df_probes[first], errors='coerce')
                                    - pd.to_numeric(df_probes[second], errors='coerce')).abs()
            continue
        target = attr.get('target')
        if target is None:
            print(f"Warning: 'target' key is missing in attribute '{attr_name}'")
            continue
        actual_target = target
        if target not in df_probes.columns:
            # Try replacing colon with dot
            dot_target = target.replace(':', '.')
            if dot_target in df_probes.columns:
                actual_target = dot_target
            else:
                print(f"Warning: '{target}' not found in csv columns.")
                continue
        attr_type: str = attr.get('type', '').lower() 
        attr_threads = attr.get('threads', threads)
        if attr_type == "mapped_sites":
            if attr.get('aligner') == "bowtie2":
                assert 'bowtie2' in genome.get('align_index', []), "Bowtie2 must be enabled in genome align_index" 
                fasta_path = Path(genome['fasta'])
                index_prefix = index_prefix_of(genome, 'bowtie2', 'genome')
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                if 'probe_id' in df_probes.columns and 'target' in df_probes.columns:
                    # DNA format (source: genome)
                    recname2seq = {f"{row['probe_id']}": row[actual_target] for _, row in df_probes.iterrows()}
                    mapped_sites_results = cal_mapped_sites(str(tmp_dir), attr_name, recname2seq, str(index_prefix), attr_threads)
                    mapped_sites = [mapped_sites_results.get(f"{row['probe_id']}", []) for _, row in df_probes.iterrows()]
                elif 'region' in df_probes.columns:
                    # RNA format (source: exon)
                    recname2seq = {f"{row['region']}": row[actual_target] for _, row in df_probes.iterrows()}
                    mapped_sites_results = cal_mapped_sites(str(tmp_dir), attr_name, recname2seq, str(index_prefix), attr_threads)
                    mapped_sites = [mapped_sites_results.get(f"{row['region']}", []) for _, row in df_probes.iterrows()]
                else:
                    raise ValueError(f"Unsupported DataFrame structure for mapped_sites attribute")
                df_probes[f"{attr_name}_num"] = [len(sites) for sites in mapped_sites]
                df_probes[attr_name] = [sites for sites in mapped_sites]
            else:
                raise NotImplementedError(
                    f"Aligner {attr['aligner']} is not implemented."
                )
        elif attr_type == "kmer_count":
            if attr.get('aligner') == "jellyfish":
                assert attr.get("locus_only", False) or genome.get('jellyfish', False) is True, \
                    "Jellyfish must be enabled in genome configuration" 
                kmer_len = attr.get("kmer_len", 35)
                # both: canonical k-mers (both strands); sense / antisense: the sequence's
                # k-mers / their reverse complements on the FASTA's forward strand.
                strand = attr.get("strand", "both")
                if strand not in STRAND_OPTIONS:
                    raise ValueError(f"attributes.{attr_name}.strand must be one of "
                                     f"{', '.join(STRAND_OPTIONS)}, got {strand!r}")
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                if attr.get("locus_only", False):
                    # Count only the candidates' k-mers across the genome, in one pass and a small temporary index.
                    if_path = locus_kmer_index(Path(genome['fasta']), df_probes[actual_target], kmer_len, attr_threads,
                                               attr.get("size"), strand == "both", tmp_dir)
                else:
                    # One index per k-mer length (and per strand setting), built when first needed.
                    from uprobe.core.tools import jellyfish_index
                    if_path = jellyfish_index(Path(genome['fasta']), kmer_len, attr_threads, attr.get("size", "1G"),
                                              canonical=strand == "both", directory=index_root(genome) / 'jf_genome')
                recname2seq = {f"{i}": row[actual_target] for i, (_, row) in enumerate(df_probes.iterrows())}
                if strand == "antisense":
                    recname2seq = {key: reverse_complement(seq) for key, seq in recname2seq.items()}
                kmer_counts = cal_kmer_count(
                    str(tmp_dir), 
                    f"{target}_{task_id}", 
                    recname2seq, 
                    str(if_path),
                    kmer_len, 
                    attr_threads
                )
                kmer_count_values = [kmer_counts[f"{i}"] for i in range(len(df_probes))]
                df_probes[attr_name] = kmer_count_values
            else:
                raise NotImplementedError(f"Aligner {attr['aligner']} is not implemented.")
        elif attr_type in {"mapped_genes"}:
            if attr.get('aligner') in {"bowtie2", "bowtie"}:
                settings = mapped_genes_settings(attr, f"attributes.{attr_name}")
                aligner = settings['aligner']
                assert aligner in genome.get('align_index', []), f"{aligner} must be enabled in genome align_index"
                fasta_path = Path(genome['fasta'])
                # fisheye calculates mapped genes against a transcriptome index,
                # not against the genomic chromosome index.
                index_prefix = index_prefix_of(genome, aligner, 'transcript')
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                # Row keys distinguish lengths at the same exon/start position.
                unique_sequences = list(dict.fromkeys(df_probes[actual_target]))
                recname2seq = {str(i): seq for i, seq in enumerate(unique_sequences)}
                if off_targets is None:
                    off_targets = OffTargetCounter(protocol, genome, report)
                # A table of transcripts (gene_expression.level: transcript) needs the transcripts hit.
                by_transcript = off_targets.transcript_level
                if aligner == "bowtie":
                    hits = (bowtie_transcript_hits if by_transcript else bowtie_gene_hits)(
                        str(tmp_dir), recname2seq, task_id, str(index_prefix),
                        settings['max_mismatches'], attr_threads, settings['strand']
                    )
                else:
                    short = sum(len(seq) < BOWTIE2_LOCAL_MIN_LENGTH for seq in unique_sequences)
                    if short:
                        log.warning(short_for_bowtie2(attr_name, short, len(unique_sequences)))
                    hits = (bowtie2_transcript_hits if by_transcript else bowtie2_gene_hits)(
                        str(tmp_dir), recname2seq, task_id, str(index_prefix),
                        attr_threads, settings['strand']
                    )
                sequence_hits = {seq: hits.get(str(i), {}) for i, seq in enumerate(unique_sequences)}
                transcript_hits = None
                if by_transcript:
                    transcript_hits = sequence_hits
                    sequence_hits = {seq: gene_hits_of(records) for seq, records in transcript_hits.items()}
                # Count, off-target list and, as configured, unexpressed genes and copies per row.
                off_targets.annotate(df_probes, attr_name, actual_target, sequence_hits, settings, transcript_hits)
            else:
                raise NotImplementedError(
                    f"Aligner {attr['aligner']} is not implemented."
                )
        elif attr_type == "annealing_temperature":
            # Only declared hybridization fields are passed; the rest keep the defaults.
            declared = resolve_conditions(protocol, attr_name)
            method = tm_method(attr, f"attributes.{attr_name}")
            if method != "primer3":
                conditions = tm_formula_arguments(method, declared)
                calculate = TM_FORMULAS[method]
            else:
                conditions = tm_arguments(declared)
                calculate = cal_temp
            vals = cached(df_probes[actual_target], lambda seq: calculate(seq, **conditions))
            df_probes[attr_name] = vals
        elif attr_type == "rna_dna_dg37":
            sequence_of = attr.get("sequence", "probe")
            if sequence_of not in SEQUENCE_OF:
                raise ValueError(f"attributes.{attr_name}.sequence must be one of {', '.join(SEQUENCE_OF)}, got {sequence_of!r}")
            df_probes[attr_name] = cached(df_probes[actual_target], lambda seq: cal_rna_dna_dg37(seq, sequence_of))
        elif attr_type == "gc_content":
            # Unrounded: filters compare it with the window bounds (19/30 = 0.633 > 0.63).
            vals = cached(df_probes[actual_target], cal_gc_content)
            df_probes[attr_name] = vals
        elif attr_type == "length":
            df_probes[attr_name] = df_probes[actual_target].str.len().astype("Int64")
        elif attr_type == "fold_score":
            temperature = resolve_conditions(protocol, attr_name).get('temperature_c')
            vals = cached(df_probes[actual_target], lambda seq: cal_target_fold_score(seq, temperature)).round(2)
            df_probes[attr_name] = vals
        elif attr_type == "self_match":
            vals = cached(df_probes[actual_target], cal_self_match).round(2)
            df_probes[attr_name] = vals
        elif attr_type == "blocks":
            if 'start' not in df_probes.columns:
                print(f"Warning: 'start' column not found for attribute '{attr_name}'")
                continue
            vals = df_probes.apply(lambda row: cal_target_blocks(row[actual_target], row['start']), axis=1)
            df_probes[attr_name] = vals
        else:
            raise NotImplementedError(
                f"Attribute type {attr_type} is not implemented."
            )
    return df_probes
