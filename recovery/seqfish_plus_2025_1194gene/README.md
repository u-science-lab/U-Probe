# seqFISH+ 2025: Probe-body Recovery

Unfiltered candidate recovery for the 1,194-gene mRNA panel in Takei et al.,
Nature (2025). Inputs are the cell-culture channel-1 probe, readout and barcode
files in the official [CaiGroup repository](https://github.com/CaiGroup/dna-seqfish-plus-multi-omics/tree/main/data/experimental-resource).

| Metric | Unfiltered candidates |
|---|---|
| Generated records | 28,851,260 |
| SRR | 99.089542% (37,548/37,893) |
| TRR | 99.832496% (1,192/1,194) |
| MTR | 99.050149% |
| CPSR | N/A (single probe body) |
| Final selection | Not evaluated |

No reference records were excluded or duplicated. Exact matching uses target,
oligo role and the complete 110-nt DNA probe body. See the four shared report
tables in [report/](report/) and [metric definitions](../sequence_recovery_metrics.md).

## Design and limitations

The published body is 5'-R1-R2-binding-R1-R3-R4-3', with a 35-nt binding
region and five 15-nt readout segments. Native U-Probe extraction and assembly
cross every valid annotated RNA window with the four cyclic arrangements of
each target's fixed codebook readouts. Controls are excluded. Published binding
and full-body sequences are used only as evaluation ground truth.

Reference sequences are mm10 with GENCODE M19 / Ensembl 94 annotation:
1,193 targets have annotated RNA candidates; DUX has none. Coverage includes
6,492 transcript IDs; transcripts shorter than 35 nt contribute no windows.
Of 345 unrecovered bodies, 29 belong to DUX and 316 have binding regions absent
from the same-target annotated RNA candidates. No readout-arrangement mismatch
was found. These observations do not establish the historical origin of the
missing binding regions. See [coverage](report/annotation_coverage.csv) and
[unmatched reasons](report/unmatched_reasons.csv).

This evaluates 110-nt probe bodies, not complete synthesis oligos. Published
quality filters and random final selection are not reproduced. New candidates
are not false positives; recovery does not establish experimental performance.

## Reproduction

From the repository root, install U-Probe and [recovery dependencies](../requirements.txt).
Place mm10.fa, mm10.fa.fai and mm10.gtf in REFERENCE_DIR. Set OUTPUT_DIR outside
the checkout; it must be empty. Checked-in design inputs are sufficient for generation.

```bash
R=recovery/seqfish_plus_2025_1194gene
PYTHONPATH=. python "$R/run_recovery.py" --config "$R/design.yaml" \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR" --workers 8
python recovery/common/sequence_recovery.py \
  --reference-oligos "$R/reference_oligos.csv" \
  --method-version seqfish_plus_2025_1194gene \
  --generated "$OUTPUT_DIR/candidates.csv.gz" --output "$OUTPUT_DIR/report"
python "$R/classify_unrecovered.py" \
  --matches "$OUTPUT_DIR/report/sequence_matches.csv" \
  --run-dir "$OUTPUT_DIR" --output "$OUTPUT_DIR/report/unmatched_reasons.csv"
```

To regenerate published input tables, place the three official channel-1 files
in RAW_DATA_DIR and run:

```bash
python "$R/prepare_reference.py" --raw-dir "$RAW_DATA_DIR" --output-dir "$R"
```

Both CSV and TXT filename suffixes are supported.
[Provenance](report/provenance.json) records input hashes.
