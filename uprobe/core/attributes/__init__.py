import secrets
import os
import tempfile
from pathlib import Path
import pandas as pd
from ._attributes import *

def add_attributes(
        df_probes: pd.DataFrame,
        protocol: dict,
        genome: dict,
) -> pd.DataFrame:
    # A private scratch directory per call: concurrent workflows sharing a
    # working directory must not share (and delete) each other's alignment files.
    with tempfile.TemporaryDirectory(prefix=".uprobe-tmp-", dir=".") as tmp_root:
        return _add_attributes(df_probes, protocol, genome, Path(tmp_root))


def _add_attributes(
        df_probes: pd.DataFrame,
        protocol: dict,
        genome: dict,
        tmp_root: Path,
) -> pd.DataFrame:
    attributes: dict = protocol['attributes']
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
                    mapped_sites_results = cal_mapped_sites(str(tmp_dir), attr_name, recname2seq, str(index_prefix), attr.get("threads", 10))
                    mapped_sites = [mapped_sites_results.get(f"{row['probe_id']}", []) for _, row in df_probes.iterrows()]
                elif 'region' in df_probes.columns:
                    # RNA format (source: exon)
                    recname2seq = {f"{row['region']}": row[actual_target] for _, row in df_probes.iterrows()}
                    mapped_sites_results = cal_mapped_sites(str(tmp_dir), attr_name, recname2seq, str(index_prefix), attr.get("threads", 10))
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
                fasta_path = Path(genome['fasta'])
                index_prefix = fasta_path.parent / 'jf_genome' / fasta_path.stem
                kmer_len = attr.get("kmer_len", 35)
                if_path = str(index_prefix) + '.jf'
                if not os.path.exists(if_path):
                    os.makedirs(str(index_prefix.parent), exist_ok=True)
                    from uprobe.core.tools.aligner import build_jf_index
                    build_jf_index(str(fasta_path), kmer_len, str(if_path), attr.get("threads", 10), attr.get("size", "1G"))
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                recname2seq = {f"{i}": row[actual_target] for i, (_, row) in enumerate(df_probes.iterrows())}
                kmer_counts = cal_kmer_count(
                    str(tmp_dir), 
                    f"{target}_{task_id}", 
                    recname2seq, 
                    str(if_path),
                    kmer_len, 
                    attr.get("threads", 10)
                )
                kmer_count_values = [kmer_counts[f"{i}"] for i in range(len(df_probes))]
                df_probes[attr_name] = kmer_count_values
            else:
                raise NotImplementedError(f"Aligner {attr['aligner']} is not implemented.")
        elif attr_type in {"mapped_genes"}:
            if attr.get('aligner') == "bowtie2":
                assert 'bowtie2' in genome.get('align_index', []), "bowtie2 must be enabled in genome align_index" 
                fasta_path = Path(genome['fasta'])
                # fisheye calculates mapped genes against a transcriptome index,
                # not against the genomic chromosome index.
                index_prefix = fasta_path.parent / 'bowtie2_transcript' / fasta_path.stem
                tmp_dir = tmp_root / attr_name
                tmp_dir.mkdir()
                # Row keys distinguish lengths at the same exon/start position.
                unique_sequences = list(dict.fromkeys(df_probes[actual_target]))
                recname2seq = {str(i): seq for i, seq in enumerate(unique_sequences)}
                n_mapped_genes = count_n_bowtie2_aligned_genes(
                    str(tmp_dir), recname2seq, task_id, str(index_prefix),
                    attr.get("threads", 10)
                )
                sequence_hits = {seq: n_mapped_genes.get(str(i), 0) for i, seq in enumerate(unique_sequences)}
                mapped_genes_values = [sequence_hits[seq] for seq in df_probes[actual_target]]
                df_probes[attr_name] = mapped_genes_values
            else:
                raise NotImplementedError(
                    f"Aligner {attr['aligner']} is not implemented."
                )
        elif attr_type == "annealing_temperature":
            vals = cached(df_probes[actual_target], cal_temp)
            df_probes[attr_name] = vals
        elif attr_type == "gc_content":
            vals = cached(df_probes[actual_target], cal_gc_content).round(2)
            df_probes[attr_name] = vals
        elif attr_type == "length":
            df_probes[attr_name] = df_probes[actual_target].str.len().astype("Int64")
        elif attr_type == "fold_score":
            vals = cached(df_probes[actual_target], cal_target_fold_score).round(2)
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
