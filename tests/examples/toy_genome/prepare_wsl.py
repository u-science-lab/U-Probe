"""Create and validate a small synthetic reference in the active WSL backend.

Run with the backend Python environment from its repository root.
"""
import json
import random
import shutil
from pathlib import Path

import yaml
from pyfaidx import Fasta
from uprobe.core.gen.fun import get_exon_seq, extract_trans_seqs, validate_targets
from uprobe.http.utils.paths import get_public_genomes_dir, get_genomes_yaml

NAME = "toy_uprobe"
out = get_public_genomes_dir() / NAME
if out.exists():
    raise SystemExit(f"Refusing to overwrite existing reference: {out}")
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(20261004)
sequences = {"1": "".join(rng.choices("ACGT", k=12000)),
             "2": "".join(rng.choices("ACGT", k=8000))}
fa_path = out / f"{NAME}.fa"
gtf_path = out / f"{NAME}.gtf"
with fa_path.open("w", newline="\n") as handle:
    for chrom, seq in sequences.items():
        handle.write(f">{chrom}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i:i + 80] + "\n")

# Standard GTF coordinates: 1-based, inclusive. Long exons support probe tiling.
genes = [
    ("1", "ToyGeneA", "+", [(1001, 1900), (2501, 3400), (4001, 4900)]),
    ("1", "ToyGeneB", "-", [(6501, 7400), (8001, 8900), (9501, 10400)]),
    ("2", "ToyGeneC", "+", [(1001, 2200), (3001, 4200)]),
]
rows = []
def row(chrom, feature, start, end, strand, attrs):
    rows.append("\t".join(map(str, [chrom, "synthetic", feature, start, end,
                                   ".", strand, ".", attrs])))

for chrom, gene, strand, exons in genes:
    attrs = f'gene_id "{gene}"; gene_name "{gene}";'
    transcript = attrs + f' transcript_id "{gene}.t1";'
    row(chrom, "gene", exons[0][0], exons[-1][1], strand, attrs)
    row(chrom, "transcript", exons[0][0], exons[-1][1], strand, transcript)
    for i, (start, end) in enumerate(exons, 1):
        number = i if strand == "+" else len(exons) + 1 - i
        row(chrom, "exon", start, end, strand,
            transcript + f' exon_number "{number}";')
    row(chrom, "UTR", exons[0][0], exons[0][0] + 149, strand, transcript)
    row(chrom, "UTR", exons[-1][1] - 149, exons[-1][1], strand, transcript)
gtf_path.write_text("# Synthetic test reference; not biological data.\n" +
                    "\n".join(rows) + "\n", encoding="utf-8")
(out / "targets.txt").write_text("ToyGeneA\nToyGeneB\nToyGeneC\n", encoding="utf-8")

targets = [gene[1] for gene in genes]
validate_targets(targets, gtf_path)
extracted = get_exon_seq(targets, str(fa_path), str(gtf_path))
assert {gene: len(exons) for gene, exons in extracted.items()} == {
    "ToyGeneA": 3, "ToyGeneB": 3, "ToyGeneC": 2}
extract_trans_seqs(str(gtf_path), str(fa_path), str(out / "transcripts.fa"))
with Fasta(str(out / "transcripts.fa")) as transcripts:
    lengths = {key: len(transcripts[key]) for key in transcripts.keys()}
# The current backend slices [start:end] without converting GTF start to zero
# based coordinates, so it drops one base per exon. Keep standard GTF data.
assert sorted(lengths.values()) == [2398, 2697, 2697], lengths

entry = {"description": "Synthetic 20 kb reference for UI and pipeline smoke tests",
         "species": "synthetic", "fasta": str(fa_path), "gtf": str(gtf_path),
         "out": str(out), "align_index": ["bowtie2"], "jellyfish": False}
config_path = get_genomes_yaml()
config = yaml.safe_load(config_path.read_text()) or {} if config_path.exists() else {}
if NAME in config:
    raise SystemExit(f"Existing configuration entry: {NAME}; no config changed")
if config_path.exists():
    backup = config_path.with_name(config_path.name + ".before_toy_uprobe.bak")
    if backup.exists():
        raise SystemExit(f"Backup already exists: {backup}; no config changed")
    shutil.copy2(config_path, backup)
config[NAME] = entry
config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
(out / "genomes.yaml").write_text(yaml.safe_dump({NAME: entry}, sort_keys=False), encoding="utf-8")
report = {"directory": str(out), "total_bases": 20000,
          "genes": targets, "exon_count": 8, "transcript_lengths": lengths,
          "standard_gtf_transcript_lengths": [2700, 2700, 2400],
          "backend_coordinate_note": "Current backend drops one base per exon (GTF start used directly as Python slice start).",
          "files": {p.name: p.stat().st_size for p in out.iterdir() if p.is_file()}}
(out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
