# MiP-seq 2024: Sequence Recovery

Independent candidate recovery of the 217-gene panel in Wu et al.,
*Nature Biomedical Engineering* (2024), [DOI](https://doi.org/10.1038/s41551-024-01205-7).
Only the `Fig 3` worksheet of the official
[supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-024-01205-7/MediaObjects/41551_2024_1205_MOESM5_ESM.xlsx)
is used. [Reference audit](reference_audit.json) records its SHA256 and source rows.

| Metric | Unfiltered candidates |
|---|---|
| Generated oligos | 20,696,144 |
| Padlock SRR | 98.847926% (429/434) |
| RCA initiator SRR | 98.847926% (429/434) |
| Overall SRR | 98.847926% (858/868) |
| TRR | 99.539171% (216/217 genes) |
| MTR | 98.847926% |
| CPSR | 98.617512% (428/434 reference pairs) |
| Final selection | Not evaluated |

[Shared evaluator reports](report/) retain full-sequence and same-generated-group
matching evidence; [metric definitions](../sequence_recovery_metrics.md) apply.
No reference sequence was excluded or duplicated. New candidate sequences are
not false positives.

## Design and limitations

The 58-nt padlock is `arm5-CC-B1-TGCGTCTATTT-GT-B2-TAGTGGAGCCT-arm3`,
with two 13-nt binding arms and two 3-nt code segments. The 25/26-nt RCA initiator
is `13-nt-binding-spacer-AGGCTCCACTA`. All published spacer forms
(`A`, `C`, `G`, `T`, `AC`, `CC`, `GC`, `TC`) are crossed with every site and each
gene's published B1/B2 combination. Internal short motifs are not treated as
unique gene identifiers.

Native U-Probe target extraction and assembly enumerate 40-nt windows at step 1
across all same-name annotated transcripts in mm10/GRCm38 with GENCODE M19
(Ensembl 94). A window contains two adjacent 13-nt padlock sites, a 1-nt gap,
and a 13-nt initiator site. Valid windows are deduplicated within each gene;
1,293,509 windows from 1,210 transcript IDs cover 216 targets. No final-count,
1,000-nt position, GC/Tm, specificity or ranking filter limits recovery.
Generation reads the YAML target list, fixed backbone and codebook; published
binding arms, full oligos and inferred coordinates are evaluation-only inputs.
Both output oligos retain their actual shared native assembly group ID.
[Provenance](report/provenance.json) records input and implementation hashes.

Four missing oligos belong to Ccn3, which has no same-name annotated transcript
in this GTF. Six missing oligos from Adgrl4-2, Kif11-2 and Pde1c-2 have binding
regions absent from legal same-gene generated windows; see
[unmatched reasons](report/unmatched_reasons.csv). Grpr-2's two oligos are
individually recovered at different windows but never in one generated group,
so this pair fails CPSR. These observations identify coverage/site differences
without establishing their historical origin. No spacer or backbone mismatch
was identified among the missing oligos. Final probe selection and experimental
performance are not evaluated.

All 868 raw oligos form 434 structurally valid reference pairs. Rows 20
(`Nts-1-PLP-2`) and 598 (`Junb-2-PLP-1`) are assigned the RCA role from their
25-nt structure, initiator suffix and adjacent matching padlock name. Syt1 rows
91/92 lack a pair index and retain their unique adjacent pairing. Five names
have outer whitespace ignored during parsing. Original names, rows and DNA are
preserved; no unresolved records remain.

## Reproduction

From the repository root, install U-Probe and [recovery dependencies](../requirements.txt).
Place `mm10.fa`, `mm10.fa.fai` and `mm10.gtf` in `REFERENCE_DIR`.
Set `OUTPUT_DIR` to an empty directory outside the checkout.

```bash
R=recovery/mipseq_2024_217gene
PYTHONPATH=. python "$R/run_recovery.py" --config "$R/design.yaml" \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR" --workers 8
python recovery/common/sequence_recovery.py \
  --reference-oligos "$R/reference_oligos.csv" \
  --reference-pairs "$R/reference_pairs.csv" \
  --method-version mipseq_2024_217gene \
  --generated "$OUTPUT_DIR/candidates.csv.gz" --output "$OUTPUT_DIR/report"
python "$R/classify_unrecovered.py" \
  --matches "$OUTPUT_DIR/report/sequence_matches.csv" \
  --run-dir "$OUTPUT_DIR" --output "$OUTPUT_DIR/report/unmatched_reasons.csv"
```

To reproduce reference extraction, download the linked workbook to `RAW_DATA_DIR`
and run (Python standard library):

```bash
python "$R/extract_reference.py" \
  "$RAW_DATA_DIR/41551_2024_1205_MOESM5_ESM.xlsx" --output-dir "$R"
```

The source workbook is not redistributed. Generation does not require it.
