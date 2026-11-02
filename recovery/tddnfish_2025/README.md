# TDDN-FISH 2025: primary sequence recovery

Source: Wang et al., *Nature Communications* (2025), [paper](https://doi.org/10.1038/s41467-025-64294-1), Supplementary Data 1 ([workbook](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-025-64294-1/MediaObjects/41467_2025_64294_MOESM3_ESM.xlsx)). Scope: the complete two-round HeLa panel, `Figure2!B82:B135`: nine human genes and 54 unique primary oligos. No records were excluded or deduplicated. `TGFB 1` is normalized to TGFB1; the six `NF-KB` bindings independently map to NFKB1, including ENST00000226574.9. Other panels are outside this benchmark.

Published primary DNA is **5′ readout(30)–binding(35) 3′**, totaling 65 nt. GRCh38 primary assembly and GENCODE v47 GTF supply independent spliced transcripts. Native `generate_target_seqs` enumerates every 35-nt window (step 1) of every annotated transcript; `construct_probes` executes the YAML reverse-complement and assembly DAG. Only non-ACGT and duplicate within-gene windows are removed. Each gene receives its two published fixed readouts. Binding DNA, full published primary DNA and mapped coordinates never enter generation. [Provenance](report/environment.json) records reference URLs, hashes and software versions.

| Measure | Candidate Recovery |
|---|---|
| Generated full primary oligos | 86,782 |
| SRR | 100% (54/54) |
| TRR | 100% (9/9) |
| MTR | 100% |
| CPSR | N/A |
| Final Recovery | Not evaluated |

The four [standard tables](report/) use the [common exact-match evaluator](../common/sequence_recovery.py) and [metric specification](../sequence_recovery_metrics.md). Generated IDs retain transcript and window origin. CPSR is not applied to single-primary units; complete physical nanostructure sets are not evaluated, and `pair_matches.csv` is header-only. No specificity, ranking or final-count selection was performed.

Structural Validation: [design.yaml](design.yaml) represents T0/T1/T2 as four strands each, with native DAG dependencies for primary→T0→T1→T2→fluorescent oligo. All 12 scaffold DNAs, ten external complementary interfaces, six round/color variants and four dependency perturbations pass. Native TSV/XLSX exports retain the explicit T0a disulfide after base 30 and 5′ Cy5 annotation; Cy3/FAM variant annotations are checked in memory. Fixed scaffolds are published design inputs, not independently recovered sequences.

Limitations: layer membership and intended 1:3:9:27 copy topology are metadata; native molecular-complex topology and folding are not evaluated. Pairwise core audits find longest contiguous complements of 16 nt (a–d) and 15 nt (b–c), rather than the stated 17-nt edge length; published DNA is preserved. Two hyb2 T0a entries lack an explicit SS marker, so no disulfide is inferred. [Structural evidence](report/structure.json) is separate from SRR. Recovery does not establish experimental performance.

## Reproduction

Install U-Probe, `recovery/requirements.txt`, and `openpyxl==3.1.5`. Download the workbook as `$WORKBOOK`. Download and decompress the GRCh38 FASTA and GENCODE v47 primary-assembly GTF from the URLs in `report/environment.json`; set `$FASTA`, `$GTF` and a fresh `$OUT` directory.

```bash
METHOD=recovery/tddnfish_2025
python "$METHOD/prepare_inputs.py" --workbook "$WORKBOOK" --output "$OUT/inputs"
python "$METHOD/generate.py" --codebook "$OUT/inputs/codebook.csv" \
  --fasta "$FASTA" --gtf "$GTF" --output "$OUT/candidates"
python recovery/common/sequence_recovery.py \
  --reference-oligos "$OUT/inputs/reference_oligos.csv" \
  --generated "$OUT/candidates/candidates.csv.gz" \
  --method-version tddnfish_2025 --output "$OUT/report"
python "$METHOD/assembly_validation.py" --workbook "$WORKBOOK" \
  --config "$METHOD/design.yaml" --output-dir "$OUT/structure"
```
