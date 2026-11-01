# ExSeq 2021: Sequence Recovery

Independent candidate recovery for the 42-gene mouse visual-cortex panel in
Alon et al., *Science* (2021), [doi:10.1126/science.aax2656](https://doi.org/10.1126/science.aax2656)
([paper and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC7900882/)).
[targets.csv](targets.csv) preserves Table S9 accessions and logical barcodes;
[reference_oligos.csv](reference_oligos.csv) preserves all 334 Table S10 padlocks,
5'-phosphate modifications and source rows.

## Results

| Measure | Candidate result |
|---|---|
| Generated padlocks | 2,246,118 |
| SRR | 99.7006% (333/334) |
| TRR | 100% (42/42) |
| MTR | 99.7024% |
| CPSR | N/A (single padlock) |
| Final selection | Not evaluated |

No reference records were excluded or duplicated. The sole unmatched full
sequence is FEZF2 probeBC_01, consistent with the independent arm-mapping audit.
See [summary](report/summary.csv), [per-gene](report/per_gene.csv) and
[sequence matches](report/sequence_matches.csv).
The completed run used Python 3.10.14 and the recorded recovery dependencies.

## Data compatibility

Existing mm10 FASTA and GENCODE M19 / Ensembl 94 GTF cover all 42 genes and
271 annotated transcripts. Each transcript is reconstructed separately.
GENCODE transcript IDs do not certify the exact S9 RefSeq accession versions:
matching uses gene names, retaining the published accessions for audit.
333/334 padlocks have an exact contiguous 32-nt match in same-gene transcripts.
FEZF2 probeBC_01 has neither 16-nt arm represented; these data cannot distinguish
accession/isoform differences from a published sequence discrepancy.
No reference is excluded. See [annotation coverage](report/annotation_coverage.csv)
and [arm mappings](report/arm_mapping.csv).

| Input | SHA-256 |
|---|---|
| mm10 FASTA | c7fe528da3fd081382ad65c1b1af3ea9c0b6251bc24f392b7b3ef4967eeed91c |
| mm10 GTF | 7ab6746c0417216e6fe299f91df8624e818ce419ae8f7d0f6118744e8cb36e6d |
| S9-S10 workbook | d34ea1cbd069c204435d2b26d1ea9ff04d582fda18a5d78c83fe2d93d91835a7 |

## Method and scope

[design.yaml](design.yaml) declares the 74-nt architecture:
`arm5 + AAA + SOLiD7 + TCTCGGGAACGCTGAAGACGGC + Illumina7 + AAA + arm3`,
with 16-nt arms and 5'-phosphate. The transcript-binding window is
`rc(arm3 + arm5)`. The runner uses native U-Probe extraction and construction,
scanning annotated transcripts in 32-nt windows at 1-nt steps.
Generation receives no reference binding arms or full oligos.
[codebook.csv](codebook.csv) contains only public gene/logical/physical-barcode
assignments, combined independently with generated target windows.

This is **unfiltered candidate recovery**, not reproduction of final selection.
The paper also applied homopolymer, GC, Tm, structure and junction-spanning BLAST
screening, advanced past accepted windows, and randomly selected at most eight
probes per gene and their barcode assignments. Historical selection state,
random seed and quantitative structure thresholds are unavailable.
These filters and final selection are not applied; candidate coverage does not
establish equivalent specificity or experimental performance.

The [common evaluator](../common/sequence_recovery.py) matches exact full DNA,
gene and probe role, recording modifications separately. SRR, TRR and MTR follow
the [shared definitions](../sequence_recovery_metrics.md). CPSR is N/A for
single-padlock ExSeq; final selection is not evaluated.

## Reproduction

From the repository root, install U-Probe and `recovery/requirements.txt`.
Keep reference data and generated candidates outside version control.
Set `REFERENCE_DIR` to the directory containing existing `mm10.fa`,
`mm10.fa.fai` and `mm10.gtf`; set `OUTPUT_DIR` to an external output root; its `generated` subdirectory must be empty.

```bash
PYTHONPATH=. python recovery/exseq_2021_42gene/audit_reference.py \
  --reference-dir "$REFERENCE_DIR" --outdir "$OUTPUT_DIR/audit"
```

To verify inputs against the original workbook, add
`--workbook "$SUPPLEMENT_TABLES_XLSX"` to the audit command.

```bash
PYTHONPATH=. python recovery/exseq_2021_42gene/run_recovery.py \
  --config recovery/exseq_2021_42gene/design.yaml \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR/generated" --workers 8

PYTHONPATH=. python recovery/common/sequence_recovery.py \
  --reference-oligos recovery/exseq_2021_42gene/reference_oligos.csv \
  --method-version ExSeq_2021_mouse_visual_cortex_42gene \
  --generated "$OUTPUT_DIR/generated/candidates.csv.gz" \
  --output "$OUTPUT_DIR/report"
```
