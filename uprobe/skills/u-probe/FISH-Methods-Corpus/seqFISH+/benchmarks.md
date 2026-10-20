# Benchmark Resources

## Benchmark Summary

The highest justified V01 probe-design grade is **C — Architecture-level only** (D01). The public codebook is useful for encoding analysis, but it is not a sequence-level probe-design ground truth. Processed expression/location resources (D02) are **D — Insufficient information** for reconstructing oligo design. Neither grade evaluates the quality of the biological experiment or the utility of the data for other analyses. These scoped grades were retained after independent review and Fix reconciliation; they do not establish sequence-level ground truth. [S01, design/data-availability sections; S02, all sheets; S06, README/tree/archive; S07, dataset description; inferred assessment, 2026-10-15]

`benchmark_available: true` records a public candidate resource, not a reproducible design benchmark. `public_probe_sets_available: true` records the verified codebook subset of probe-related resources. Full synthesized primary/readout libraries remain unresolved. No Grade A or B candidate was established in this inspection.

## Candidate Datasets

### D01 — Original 10,000-gene codebook

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Supplementary Table 1: Codebook for 10,000 genes |
| Method version | V01, original RNA seqFISH+ (2019) |
| Species | Mus musculus; S01 experiment |
| Cell / tissue / sample | NIH/3T3; mouse cortex/SVZ and olfactory bulb; same panel reported in S01 |
| Genome build | Unresolved; S01 names UCSC masked genome but no assembly release |
| Transcriptome / annotation version | Unresolved; UCSC annotation and mouse transcriptome named without exact release |
| Target genes | 10,000 gene-name rows across three sheets; 3,334/3,333/3,333 for 647/561/488 nm, derived by counting |
| Target regions | Coding exons preferentially, UTR fallback per S01; per-probe regions not provided in table |
| Genomic/transcript coordinates | Not present in workbook; assembly/accession, strand and origin/interval conventions unresolved; gene names are not coordinates |
| Target-binding sequences | Not present in this workbook |
| Full probe sequences | Not present in this workbook |
| Barcode assignments | Present: four ordered integers 1–20 per gene; channel given by sheet |
| Readout sequences | Nucleotide strings and symbol-to-oligo map not present |
| Primer/linker sequences | Not present |
| Design software | Original-author sequence-design implementation unresolved; linked S06 software is an image-analysis pipeline |
| Design parameters | Reported filters in S01/probes.md; executable recipe, exact references and selection details unresolved |
| Source | S02: 41586_2019_1049_MOESM2_ESM.xlsx; S01, Primary-probe design and encoding strategy |
| Availability | Downloaded from official Caltech article deposit 2026-10-15; evidence copy retained |
| Missing information | Target sequences, final oligos, readout library, primer/linker strings, coordinates, exact reference/input versions and original implementation |
| Benchmark grade | C |
| Grade justification | Supports architecture and codeword comparison, including channel membership and code distance. Missing oligo-level outputs and references are major gaps, so cannot reproduce original sequence selection; B is not justified |
| Comparison scope | Coding scheme/codebook and reported rule comparisons only; not exact target-site or full-oligo overlap with U-Probe |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource contents; derived counts and distance; inferred benchmark grade |

Workbook locators: `Codebook-647nm!A2:F3335`, `Codebook-561nm!A2:F3334`, `Codebook-488nm!A2:F3334`; column B is blank and rounds are C–F. SHA-256: `45441bc0b0412a7899c686330562e786fb06b1d1c2de6744da52ea281bb8c2a6`. File identifiers and calculations are recorded in evidence/inspection.json and research_manifest.md. [S02; reported file layout, derived digest/counts; V01]

### D02 — Official processed RNA seqFISH+ measurements

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | CaiGroup/seqFISH-PLUS source data and NIH3T3 point-location deposit |
| Method version | V01 |
| Species | Mus musculus |
| Cell / tissue / sample | NIH/3T3 and mouse brain; Zenodo point locations specifically NIH/3T3 biological replicates 1 and 2 |
| Genome build | Not established by processed measurements |
| Transcriptome / annotation version | Not established by processed measurements |
| Target genes | Gene labels/count matrices; original panel mapping in D01 |
| Target regions | No oligo-level regions established by inspected archive members or deposit description |
| Genomic/transcript coordinates | Not established; spatial X/Y locations are not genomic/transcript target coordinates |
| Target-binding sequences | No sequence library identified in inspected repository archive |
| Full probe sequences | No sequence library identified in inspected repository archive |
| Barcode assignments | S06 states gene IDs correspond to codebook; use S02, not matrix column position assumptions |
| Readout sequences | Not identified in inspected archive/file metadata |
| Primer/linker sequences | Not identified in inspected archive/file metadata |
| Design software | S06 analysis code, not original oligo design software |
| Design parameters | Not reconstructed from decoded RNA counts or spatial point coordinates |
| Source | S06 README, recursive tree, sourcedata.zip; S07 Zenodo 2669683 Description/Files |
| Availability | S06 archive downloaded and inventories/CSV headers inspected; S07 landing page and file descriptions inspected, large payloads not downloaded; 2026-10-15 |
| Missing information | Oligo-level inputs, outputs, reference versions and sequence assignments; Zenodo payload content not independently audited |
| Benchmark grade | D |
| Grade justification | Post-decoding measurements cannot reconstruct the original probe-design outputs from the inspected content; useful biological/image-analysis data do not supply design truth |
| Comparison scope | Potential downstream expression/localization or decoding analyses, which need separate assessments; unsuitable as a strict U-Probe sequence-design benchmark |
| Assessment date | 2026-10-15 |
| Evidence status | reported repository/deposit scope; inferred design-specific grade, bounded by inventory inspection |

The downloaded S06 archive contains count CSVs, cell-centroid CSVs and `RNA_locations_cell_cortex.zip`; the nested archive contains per-cell RNA-location CSVs. No oligo sequence file is identified by these inventories. S07 describes DAPI, ROI and MATLAB point-location/gene-name files. Absence of original probe strings in these examined resources does not prove none exist elsewhere. [S06, archive inventories and CSV headers; S07, Description and Files; reported scope/inferred resource classification; V01]
