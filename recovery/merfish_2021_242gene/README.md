# MERFISH 2021: Candidate Sequence Recovery

Candidate recovery benchmark for 22,005 encoding probes from the 242-gene mouse
motor-cortex panel in [Zhang et al., Nature (2021)](https://doi.org/10.1038/s41586-021-03705-x).
Public inputs are Supplementary files 1-3 from
[Fang et al., eLife (2024)](https://doi.org/10.7554/eLife.90029), which reports using
the same encoding probes: codebook, complete encoding sequences and readouts,
respectively. The original 2021 supplementary workbooks were not directly
cross-checked. The 16 sequential two-color FISH genes are excluded.

## Results

| Metric | Unfiltered candidates |
|---|---|
| Generated sequences | 17,702,580 |
| SRR | 99.995456% (22,004/22,005) |
| TRR | 100.000000% (242/242) |
| MTR | 99.995508% |
| CPSR | N/A (single encoding oligo) |
| Final selection | Not evaluated |

No reference records were excluded. See [summary](report/summary.csv),
[per-gene results](report/per_gene.csv) and [sequence matches](report/sequence_matches.csv).
Metric definitions are [shared across methods](../sequence_recovery_metrics.md).

## Design and limitations

[design.yaml](design.yaml) specifies the published 112-nt encoding-template DNA
sequence, 5' to 3':
20-nt fixed flank, A, 20-nt reverse-complement readout, 30-nt target-binding region,
A, 20-nt reverse-complement readout and 20-nt fixed flank.
The fixed flanks are `CGCGGTGGTGGTACGACATG` and `AGCGCAACTAGGACGGACAA`.
Recovery evaluates these published DNA sequences; amplified or mature probe
products and T7-promoter incorporation are not evaluated.
Native U-Probe extraction and construction cross all annotated 30-nt transcript
windows with all 12 ordered distinct readout pairs permitted by each four-on-bit
barcode. Published binding sequences and complete probes are used only for
validation and evaluation, never candidate generation.

GENCODE M19 / Ensembl 94 covers all 242 genes and 1,212 transcripts. All 242
published unversioned transcript IDs are present; this does not establish
historical transcript-version identity. One THBS2 binding region is absent from
the same-gene annotated transcripts; its cause remains unresolved.
The Excel-converted March1 gene label is corrected by its transcript ID
`ENSMUST00000072482`.

Published quality filters and random final selection were not reproduced.
Candidate recovery does not establish equivalent specificity or experimental
performance. New candidate sequences are not counted as false positives.

## Reproduction

From the repository root, install U-Probe and [recovery dependencies](../requirements.txt).
Use `mm10.fa`, `mm10.fa.fai` and `mm10.gtf` in `REFERENCE_DIR`.
Set `OUTPUT_DIR` outside the checkout; its `generated/` directory must be empty.

```bash
R=recovery/merfish_2021_242gene
python "$R/prepare_reference.py" --raw-dir "$RAW_DATA_DIR" --output-dir "$R"

PYTHONPATH=. python "$R/audit_reference.py" \
  --reference-dir "$REFERENCE_DIR" --outdir "$OUTPUT_DIR/audit"

PYTHONPATH=. python "$R/run_recovery.py" --config "$R/design.yaml" \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR/generated" --workers 8

python recovery/common/sequence_recovery.py \
  --reference-oligos "$R/reference_oligos.csv" \
  --method-version MERFISH_2021_mouse_motor_cortex_242gene \
  --generated "$OUTPUT_DIR/generated/candidates.csv.gz" \
  --output "$OUTPUT_DIR/report"
```

`RAW_DATA_DIR` contains `elife-90029-supp1.xlsx`, `elife-90029-supp2.xlsx` and
`elife-90029-supp3.xlsx`; preparation verifies their SHA256 checksums.
The checked-in CSV and YAML inputs allow recovery without rerunning preparation.
