import secrets
import os
import tempfile
from pathlib import Path
import pandas as pd
from ._attributes import *
from ..hybridization import multi_padlock_design_arguments, resolve_conditions, tm_arguments, tm_method
from ..offtarget import BOWTIE2_LOCAL_MIN_LENGTH, OffTargetCounter, mapped_genes_settings, short_for_bowtie2
from ..utils import get_logger

log = get_logger(__name__)

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
                index_prefix = fasta_path.parent / 'bowtie2_genome' / fasta_path.stem
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
                assert genome.get('jellyfish', False) is True, "Jellyfish must be enabled in genome configuration" 
                kmer_len = attr.get("kmer_len", 35)
                # One index per k-mer length, built when first needed.
                from uprobe.core.tools import jellyfish_index
                if_path = jellyfish_index(Path(genome['fasta']), kmer_len, attr_threads, attr.get("size", "1G"))
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                recname2seq = {f"{i}": row[actual_target] for i, (_, row) in enumerate(df_probes.iterrows())}
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
                index_prefix = fasta_path.parent / f'{aligner}_transcript' / fasta_path.stem
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
            if tm_method(attr, f"attributes.{attr_name}") == "multi_padlock_design":
                conditions = multi_padlock_design_arguments(declared)
                calculate = cal_temp_multi_padlock_design
            else:
                conditions = tm_arguments(declared)
                calculate = cal_temp
            vals = cached(df_probes[actual_target], lambda seq: calculate(seq, **conditions))
            df_probes[attr_name] = vals
        elif attr_type == "gc_content":
            vals = cached(df_probes[actual_target], cal_gc_content).round(2)
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
