# Benchmark Resources

## Benchmark Summary

D01 is Grade C for architecture-level comparison; D02 is Grade D for biological probe-design reconstruction from the retrieved codebook. D03/D04 remain unassessed because bulk data were not inspected. No Grade A/B ground truth has been established. The metadata grade C summarizes the highest justified assessed candidate, not an average. Partial public probe-resource availability does not imply a complete biological probe set.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published sequence list is not automatically Grade A. These are Research-stage assessments, not independent validation.

## Candidate Datasets

### D01 — Author bridge and padlock reference library

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Author bridge and padlock reference library |
| Method version | V01; implementation/data-resource scope, not an independently established chemistry version |
| Species | Library sequences are synthetic; associated target-design example uses Mus musculus (taxonomy 10090). |
| Cell / tissue / sample | Mouse-brain background list provided; generic reusable library is not a tissue sample. |
| Genome build | Unresolved for biological probe design. mm39_tRNA.fa is a named blacklist input, insufficient to assign all designs to mm39. |
| Transcriptome / annotation version | RefSeq RNA database release unspecified; live Entrez retrieval is possible. |
| Target genes | Example_ProbeDesign.csv contains Gad2; this is an input template, not a completed probe output. |
| Target regions | Code selects two 30-nt binding regions with a 2-nt gap; no actual gene-specific output in this input file. |
| Genomic/transcript coordinates | No finalized coordinates. Code candidate start is zero-based index into selected input record; genomic mapping, accession version and feature offset need separate retention. |
| Target-binding sequences | Not supplied as finalized biological probe pairs in examined sequence/input files. |
| Full probe sequences | UsingPad.fa has 110 144-nt doubled references; not established as synthesized oligos. Biological A/B ordering list not established. |
| Barcode assignments | Padlock identifiers linked to bridge identifiers; complete experimental gene/RO/padlock mapping unresolved. |
| Readout sequences | Internal padlock segments and notebook example ROs; no verified fully annotated experimental readout ordering list. |
| Primer/linker sequences | Code literals ta; selectable padlock connectors ta/at/tt/aa; amplification primer and modifications unresolved. |
| Design software | S04 notebooks 01–05/openDesign; S03 companion implementation. |
| Design parameters | Code and notebook settings available; GC/P2 evalue defaults differ; final seeds and manual selections unresolved. |
| Source | S04 README, SEQUENCES/Bridge_sequences.csv, UsingPad.fa, Example_ProbeDesign.csv, notebooks 03/05 and openDesign/utils.py; S03 ProbeDesign tree. |
| Availability | Retrieved 2026-10-15 at S04 be490eb19503c742230665da1e665c6d6c3d97c8; local checksums preserved. |
| Missing information | Frozen biological reference records/database; final target probe outputs and manual selections; chemistry; experimental library annotation. |
| Benchmark grade | C |
| Grade justification | Useful for architecture, sequence assembly and library-consistency comparisons. Missing biological probe output/reference context prevents a strict target-design benchmark; gaps are too substantial for B or A. |
| Comparison scope | Modular architecture and reusable-library handling; not strict biological target-probe ground truth. |
| Assessment date | 2026-10-15 |
| Evidence status | inferred; sequence counts and duplicated lengths are derived by direct file inspection. |

### D02 — OpenFISH-MALDI protocol demo

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | OpenFISH-MALDI protocol demo |
| Method version | V01; implementation/data-resource scope, not an independently established chemistry version |
| Species | Not explicitly established by the retrieved metadata/codebook alone; mouse context appears in associated design example. |
| Cell / tissue / sample | Demo OpenFISH/MALDI data; actual tissue identity requires primary experimental metadata. |
| Genome build | Unresolved in inspected resources. |
| Transcriptome / annotation version | Unresolved in inspected resources. |
| Target genes | CodingBook.CSV lists Cck, Rgs5, P2ry12, Cbln1, Fn1, Dner, Pltp, Grm3, Synpr, Nfib (R1); Satb2, Pdgfra, Bcl11b, Shox2, Zfpm2, Otp, Gpr17, Mrc1, Nptxr, Agt (R2). |
| Target regions | Unresolved in inspected resources. |
| Genomic/transcript coordinates | Unresolved in inspected resources. |
| Target-binding sequences | Not in retrieved codebook. |
| Full probe sequences | Not in retrieved codebook. |
| Barcode assignments | 20 gene-to-two-color-to-round entries; exact public file retrieved. |
| Readout sequences | Codebook has channel labels rather than nucleotide strings. |
| Primer/linker sequences | Unresolved in inspected resources. |
| Design software | Associated author workflow S04; openDecode Version 2.0 notebook. |
| Design parameters | Demo decoding parameters in S04 notebook 06; not biological probe-design parameters. |
| Source | S05 10.5281/zenodo.20390362 metadata/file list and CodingBook.CSV; S04 notebook 06. |
| Availability | Metadata and 417-byte codebook retrieved 2026-10-15. Large CZI/imzML/IBD/Zarr/H5AD files listed but not downloaded. |
| Missing information | Probe sequences, reference release and design outputs; archive contents not examined. |
| Benchmark grade | D |
| Grade justification | Retrieved codebook alone cannot reconstruct biological probe design. Useful decoding demo status does not make it probe-design ground truth. |
| Comparison scope | Gene/round/channel codebook ingestion; possible image-decoding benchmark after separate bulk-data inspection. |
| Assessment date | 2026-10-15 |
| Evidence status | reported file contents; grade is inferred within retrieved scope. |

### D03 — Main study OpenFISH and MALDI datasets

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | Main study OpenFISH and MALDI datasets |
| Method version | V01; implementation/data-resource scope, not an independently established chemistry version |
| Species | Unresolved from retrieved metadata alone. |
| Cell / tissue / sample | Raw decoded OpenFISH, MALDI-MSI and processed annotations/registration data per deposit description. |
| Genome build | Unresolved in inspected resources. |
| Transcriptome / annotation version | Unresolved in inspected resources. |
| Target genes | Unresolved in inspected resources. |
| Target regions | Unresolved in inspected resources. |
| Genomic/transcript coordinates | Unresolved in inspected resources. |
| Target-binding sequences | Unresolved in inspected resources. |
| Full probe sequences | Unresolved in inspected resources. |
| Barcode assignments | Unresolved in inspected resources. |
| Readout sequences | Unresolved in inspected resources. |
| Primer/linker sequences | Unresolved in inspected resources. |
| Design software | S03 companion analysis repository (S01 Data and code availability). |
| Design parameters | Unresolved in inspected resources. |
| Source | S06 10.5281/zenodo.17659890 metadata and file listing; S01 Data and code availability. |
| Availability | Metadata/file list retrieved 2026-10-15: imzMLs.zip, spatialdatas.zip, processed.zip. Archives not downloaded or inspected. |
| Missing information | Probe-level contents and reference metadata inside archives unexamined. Deposit publication_date is 2027-08-19, later than inspection date; preserved as a source inconsistency rather than corrected by inference. |
| Benchmark grade | Unassessed. |
| Grade justification | Only metadata inspected. Do not infer absence of probe tables within uninspected archives or assign a default grade. |
| Comparison scope | Candidate data-analysis benchmark; probe-design suitability unresolved. |
| Assessment date | 2026-10-15 |
| Evidence status | reported metadata; suitability unresolved. |

### D04 — Adjacent-section OpenFISH versus Xenium V1 benchmark

| Field | Value |
|---|---|
| Dataset ID | D04 |
| Dataset name | Adjacent-section OpenFISH versus Xenium V1 benchmark |
| Method version | V01; implementation/data-resource scope, not an independently established chemistry version |
| Species | Mus musculus, reported as mouse brain in deposit description. |
| Cell / tissue / sample | Adjacent brain sections from matched tissue blocks; two replicates; reported RIN 8.7 and 8.4 (S07 description). |
| Genome build | Unresolved in inspected resources. |
| Transcriptome / annotation version | Unresolved in inspected resources. |
| Target genes | Unresolved in inspected resources. |
| Target regions | Unresolved in inspected resources. |
| Genomic/transcript coordinates | Unresolved in inspected resources. |
| Target-binding sequences | Unresolved in inspected resources. |
| Full probe sequences | Unresolved in inspected resources. |
| Barcode assignments | Unresolved in inspected resources. |
| Readout sequences | Unresolved in inspected resources. |
| Primer/linker sequences | Unresolved in inspected resources. |
| Design software | Unresolved in inspected resources. |
| Design parameters | Unresolved in inspected resources. |
| Source | S07 10.5281/zenodo.20917671 metadata/file list; S01 Data and code availability. |
| Availability | Metadata retrieved 2026-10-15; four Zarr ZIP archives listed (Slide1/Slide2 × OpenFISH/Xenium). Bulk contents not downloaded. |
| Missing information | Actual matched panel, coordinates, probe sequences/reference releases and full sample processing metadata require archive/supplement inspection. |
| Benchmark grade | Unassessed. |
| Grade justification | A platform-comparison experiment is not automatically probe-design ground truth. Metadata-only inspection does not establish resource completeness. |
| Comparison scope | Matched-panel measurement/platform comparison; exact target-design reproducibility unresolved. |
| Assessment date | 2026-10-15 |
| Evidence status | reported metadata; suitability unresolved. |
