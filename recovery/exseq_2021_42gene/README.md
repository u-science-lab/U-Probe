# ExSeq 2021: Candidate Sequence Recovery

Candidate recovery benchmark for 334 padlocks from the 42-gene mouse visual-cortex panel in
[Alon et al., *Science* (2021)](https://doi.org/10.1126/science.aax2656),
using Supplementary Tables S9-S10 and mm10.

## Results

| Metric | Unfiltered candidates |
|---|---|
| Generated sequences | 2,246,118 |
| SRR | 99.7006% (333/334) |
| TRR | 100% (42/42) |
| MTR | 99.7024% |
| CPSR | N/A (single padlock) |
| Final selection | Not evaluated |

No reference records were excluded. See [summary](report/summary.csv),
[per-gene results](report/per_gene.csv) and [sequence matches](report/sequence_matches.csv).
Metric definitions are [shared across methods](../sequence_recovery_metrics.md).

## Design and limitations

[design.yaml](design.yaml) specifies 74-nt, 5'-phosphorylated padlocks with two
16-nt binding arms. Native U-Probe extraction and construction combine all
annotated 32-nt transcript windows with the public [codebook](codebook.csv).
Published binding arms and full probe sequences are used only for evaluation.

GENCODE M19 / Ensembl 94 covers all 42 genes and 271 transcripts, but does not
certify the S9 RefSeq accession versions. Binding regions map exactly for
333/334 padlocks; neither arm of FEZF2 `probeBC_01` maps to the same-gene
transcripts. The cause remains unresolved. See [annotation coverage](report/annotation_coverage.csv)
and [arm mappings](report/arm_mapping.csv).

Published quality filters and random final selection were not reproduced;
historical selection parameters are incomplete. Candidate recovery therefore
does not establish equivalent specificity or experimental performance.

## Reproduction

From the repository root, install U-Probe and [recovery dependencies](../requirements.txt).
Use the existing `mm10.fa`, `mm10.fa.fai` and `mm10.gtf` in `REFERENCE_DIR`.
Set `OUTPUT_DIR` outside the checkout; its `generated/` directory must be empty.

```bash
R=recovery/exseq_2021_42gene
PYTHONPATH=. python "$R/audit_reference.py" \
  --reference-dir "$REFERENCE_DIR" --outdir "$OUTPUT_DIR/audit"

PYTHONPATH=. python "$R/run_recovery.py" --config "$R/design.yaml" \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR/generated" --workers 8

python recovery/common/sequence_recovery.py \
  --reference-oligos "$R/reference_oligos.csv" \
  --method-version ExSeq_2021_mouse_visual_cortex_42gene \
  --generated "$OUTPUT_DIR/generated/candidates.csv.gz" \
  --output "$OUTPUT_DIR/report"
```

Add `--workbook "$SUPPLEMENT_TABLES_XLSX"` to the audit command to verify
[targets](targets.csv), [reference oligos](reference_oligos.csv) and codebook
against the original S9-S10 workbook.
