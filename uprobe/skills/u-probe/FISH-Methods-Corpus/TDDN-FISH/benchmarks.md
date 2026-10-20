# Benchmark Resources

## Benchmark Summary

TDDN-FISH has a public candidate dataset (Supplementary Data 1) containing primary probe sequences, TDDN component sequences, and fluorophore-labeled detection oligo sequences. However, the sequences were not fully examined in this research stage. The benchmark grade is C (Architecture-level only) because the probe set exists and provides architecture information, but full reconstruction of probe-design output requires examination of Supplementary Data 1 with reference context.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A. Consider target/reference context and relevant design metadata. Unassessed is not D.

## Candidate Datasets

### D01 — Nature Communications Supplementary Data 1

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | TDDN-FISH Supplementary Data 1 |
| Method version | V01 |
| Species | Homo sapiens (HeLa, A549, MDCK cells); Mus musculus (mouse brain tissue) |
| Cell / tissue / sample | HeLa cells, A549 cells, MDCK cells, mouse brain cryosections (10 μm), thick mouse brain sections (40 μm), influenza A virus-infected A549 cells |
| Genome build | Not reported |
| Transcriptome / annotation version | Not reported |
| Target genes | ACTB, POLR2A, HSP70, mTOR, B2M, MALAT1, NEAT1, miR-21, Satb2, Gad1, Slc17a6, Lmnb1, Slc17a7, PPIA, PA vRNA, PB1 vRNA, and others for 60-gene cell-type mapping (6 cell types × 10 marker genes) and 90-gene V1 mapping (9 cell types × 10 marker genes) |
| Target regions | 35 nt target-specific sequences per probe |
| Genomic/transcript coordinates | Not reported |
| Target-binding sequences | Available in Supplementary Data 1 (not examined in this research) |
| Full probe sequences | Available in Supplementary Data 1 (not examined in this research) |
| Barcode assignments | F^N combinatorial barcoding (3 channels × N rounds); experiment-specific assignments not documented in centralized codebook |
| Readout sequences | 30 nt readout sequences derived from MERFISH |
| Primer/linker sequences | TDDN component sequences (T0a-d, T1a-d, T2a-d) and fluorophore-labeled detection oligos |
| Design software | OligoArray 2.0 for target sequence design |
| Design parameters | Target region: 35 nt; Tm > 70°C; minimum off-target Tm 72°C; no internal secondary structures with Tm < 76°C; no contiguous runs > 6 |
| Source | S01, Supplementary Data 1 (Nature Communications article) |
| Availability | Available via article supplementary materials |
| Missing information | Exact sequences not examined; genome build not reported; transcriptome version not reported |
| Benchmark grade | C |
| Grade justification | Architecture and design parameters are reported, but exact sequences require examination of Supplementary Data 1. The probe set is useful for understanding TDDN-FISH probe structure and design principles, but full reconstruction for direct comparison requires examination of all sequences with reference context. |
| Comparison scope | Method comparison for understanding TDDN-FISH architecture and probe-design approach |
| Assessment date | 2026-10-15 |
| Evidence status | reported |
