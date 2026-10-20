# Benchmark Resources

## Benchmark Summary

D01–D02 are original supplementary-workbook candidates (S01/S02; V01). D03 is a ClusterMap image-resource candidate (S08). D04–D05 are STARmap PLUS deposits cited by S11 Data availability (V02). All remain **Unassessed**, and metadata benchmark.grade remains null. Direct workbook inspection and complete deposit inventories were not completed in this Fix round. Uninspected content is unknown, not unavailable; image/processed-data labels alone do not prove that no associated probe resources exist.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

No provisional grades or file-size-based content predictions are assigned. Later S03 design rules do not establish the exact reference/software settings for these datasets.

## Candidate Datasets

### D01 — Wang 2018 *Science* 160-gene panel (Supplementary Table S1)

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Wang 2018 STARmap 160-gene mouse V1 panel (Supplementary Table S1) |
| Method version | V01 (STARmap, Wang 2018) |
| Species | *Mus musculus* (mouse, C57BL/6) |
| Cell / tissue / sample | Mouse V1; S01 reports 8-µm imaged volumes for 160-gene work. Exact workbook-to-sample mapping not yet verified. |
| Genome build | Not yet verified for probe-design comparison; do not infer an assembly from publication year or species. |
| Transcriptome / annotation version | Not yet verified; scRNA-seq panel-selection provenance does not establish a design annotation release. |
| Target genes | 160 genes: 112 putative cell-type markers + 48 activity-regulated genes (ARGs) (per S01 main text) |
| Target regions | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Genomic/transcript coordinates | Not yet verified: assembly/accession, transcript version, strand, origin and interval convention remain unresolved. |
| Target-binding sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Full probe sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Barcode assignments | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Readout sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Primer/linker sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Design software | Exact implementation for this dataset unverified. Later S03 describes Picky 2.2 + Dedupe/BBMap + NUPACK + BLAST; S05 equivalence and original-panel use are unverified. |
| Design parameters | The V01 design rules are reported in S02 ("STARmap probe design", p. 3; Fig. S1A legend, p. 10; see probes.md): shortest isoform, coding region, Picky 2.2, 40–46 nt split into 20–25 nt halves with a 0–2 nt gap, designed Tm 60 °C, 4 pairs per gene. Whether every probe of this dataset followed them is unverified. Later S03 recommendations are separately scoped in probes.md; they are not original-panel ground truth. |
| Source | S01 supplementary-materials section; S02 workbook reference, Table S1 (V01); workbook contents not inspected. |
| Availability | Science article page; existence confirmed; content not directly retrieved for this draft |
| Missing information | Supplementary table content (probe sequences, barcode assignments, coordinates, per-gene target regions); the exact genome build and transcriptome annotation version used for design; the per-step design parameters actually used; the reference build for the BLAST specificity filter |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned: direct content/inventory and reference context have not been sufficiently inspected. File size is not interpreted as evidence of contents or reproducibility. |
| Comparison scope | Candidate for future assessment; strict probe-design comparison is not currently established. |
| Assessment date | 2026-10-15 (unassessed after Fix; no grade awarded) |
| Evidence status | reported |

### D02 — Wang 2018 *Science* 1020-gene panel (Supplementary Table S2)

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | Wang 2018 STARmap 1020-gene mouse V1 panel (Supplementary Table S2) |
| Method version | V01 (STARmap, Wang 2018) |
| Species | *Mus musculus* (mouse, C57BL/6) |
| Cell / tissue / sample | Primary visual cortex (V1), coronal sections, 8-µm thick (cryosections) |
| Genome build | Not yet verified for probe-design comparison; do not infer an assembly from publication year or species. |
| Transcriptome / annotation version | Not yet verified; scRNA-seq panel-selection provenance does not establish a design annotation release. |
| Target genes | 1,020 genes (per S01 main text) |
| Target regions | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Genomic/transcript coordinates | Not yet verified: assembly/accession, transcript version, strand, origin and interval convention remain unresolved. |
| Target-binding sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Full probe sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Barcode assignments | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Readout sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Primer/linker sequences | Not yet examined in the original workbook; do not infer completeness from the file title or size. |
| Design software | Exact implementation for this dataset unverified. Later S03 describes Picky 2.2 + Dedupe/BBMap + NUPACK + BLAST; S05 equivalence and original-panel use are unverified. |
| Design parameters | The V01 design rules are reported in S02 ("STARmap probe design", p. 3; Fig. S1A legend, p. 10; see probes.md): shortest isoform, coding region, Picky 2.2, 40–46 nt split into 20–25 nt halves with a 0–2 nt gap, designed Tm 60 °C, 4 pairs per gene. Whether every probe of this dataset followed them is unverified. Later S03 recommendations are separately scoped in probes.md; they are not original-panel ground truth. |
| Source | S01 supplementary-materials section; S02 workbook reference, Table S2 (V01); workbook contents not inspected. Note: S02 ("STARmap probe design", p. 3) says "All sequences were included in Table S1 (SNAIL probes) and Table S2 (SEDAL probes)", so Table S2 may hold the sequencing probes rather than a 1,020-gene panel; the identification of D02 with a 1,020-gene panel is unverified. |
| Availability | Science article page; existence confirmed; content not directly retrieved for this draft |
| Missing information | Same as D01 |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned: direct content/inventory and reference context have not been sufficiently inspected. File size is not interpreted as evidence of contents or reproducibility. |
| Comparison scope | Candidate for future assessment; strict probe-design comparison is not currently established. |
| Assessment date | 2026-10-15 (unassessed after Fix; no grade awarded) |
| Evidence status | reported |

### D03 — Wang 2018 STARmap raw image dataset (ClusterMap bundled)

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | Wang 2018 STARmap raw fluorescence image stacks (2D 1020-gene + 3D 28-gene mouse V1) |
| Method version | V01 (STARmap, Wang 2018) |
| Species | *Mus musculus* (mouse, C57BL/6) |
| Cell / tissue / sample | Primary visual cortex (V1), 2D 1020-gene and 3D 28-gene raw fluorescence image stacks |
| Genome build | Not yet verified for probe-design comparison; do not infer an assembly from publication year or species. |
| Transcriptome / annotation version | Not yet verified; scRNA-seq panel-selection provenance does not establish a design annotation release. |
| Target genes | 1,020 (2D) + 28 (3D) |
| Target regions | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Genomic/transcript coordinates | Not yet verified: assembly/accession, transcript version, strand, origin and interval convention remain unresolved. |
| Target-binding sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Full probe sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Barcode assignments | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Readout sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Primer/linker sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Design software | Exact implementation for this dataset unverified. Later S03 describes Picky 2.2 + Dedupe/BBMap + NUPACK + BLAST; S05 equivalence and original-panel use are unverified. |
| Design parameters | Exact parameters for this dataset unverified. Later S03 recommendations are separately scoped in probes.md; they are not original-panel ground truth. |
| Source | wanglab-broad/ClusterMap/datasets/RawImage (S08) |
| Availability | GitHub ClusterMap repository; existence confirmed; full file inventory not captured for this draft |
| Missing information | Probe sequences; barcode assignments; design parameters; the per-gene target regions; the per-gene coordinates; the design software parameters used |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned: direct content/inventory and reference context have not been sufficiently inspected. File size is not interpreted as evidence of contents or reproducibility. |
| Comparison scope | Candidate for future assessment; strict probe-design comparison is not currently established. |
| Assessment date | 2026-10-15 (unassessed after Fix; no grade awarded) |
| Evidence status | reported |

### D04 — Zeng 2023 STARmap PLUS processed data (SCP1375)

| Field | Value |
|---|---|
| Dataset ID | D04 |
| Dataset name | Zeng 2023 STARmap PLUS 2,766-gene mouse AD dataset (Broad Single-Cell Portal SCP1375) |
| Method version | V02 (STARmap PLUS, Zeng 2023) |
| Species | *Mus musculus* (mouse, TauPS2APP Alzheimer model) |
| Cell / tissue / sample | TauPS2APP/control brain cortex and hippocampus; S11 Methods: 20-µm coronal sections; exact deposit sample mapping unexamined |
| Genome build | Not yet verified for probe-design comparison; do not infer an assembly from publication year or species. |
| Transcriptome / annotation version | Not yet verified; scRNA-seq panel-selection provenance does not establish a design annotation release. |
| Target genes | 2,766 in main study; separate 64-gene validation (S11 Results and Methods); deposit contents unexamined |
| Target regions | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Genomic/transcript coordinates | Not yet verified: assembly/accession, transcript version, strand, origin and interval convention remain unresolved. |
| Target-binding sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Full probe sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Barcode assignments | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Readout sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Primer/linker sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Design software | Exact implementation for this dataset unverified. Later S03 describes Picky 2.2 + Dedupe/BBMap + NUPACK + BLAST; S05 equivalence and original-panel use are unverified. |
| Design parameters | Exact parameters for this dataset unverified. Later S03 recommendations are separately scoped in probes.md; they are not original-panel ground truth. |
| Source | Broad Single-Cell Portal SCP1375 (S09) |
| Availability | https://singlecell.broadinstitute.org/single_cell/study/SCP1375 ; existence confirmed; full file inventory not captured for this draft |
| Missing information | Probe sequences; barcode assignments; design parameters; the per-gene target regions; the per-gene coordinates; the design software parameters used |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned: direct content/inventory and reference context have not been sufficiently inspected. File size is not interpreted as evidence of contents or reproducibility. |
| Comparison scope | Candidate for future assessment; strict probe-design comparison is not currently established. |
| Assessment date | 2026-10-15 (unassessed after Fix; no grade awarded) |
| Evidence status | reported |

### D05 — Zeng 2023 STARmap PLUS processed data (Zenodo 10.5281/zenodo.7332091)

| Field | Value |
|---|---|
| Dataset ID | D05 |
| Dataset name | Zeng 2023 STARmap PLUS processed data (Zenodo mirror of SCP1375) |
| Method version | V02 (STARmap PLUS, Zeng 2023) |
| Species | *Mus musculus* (mouse, TauPS2APP Alzheimer model) |
| Cell / tissue / sample | Same as D04 |
| Genome build | Not yet verified for probe-design comparison; do not infer an assembly from publication year or species. |
| Transcriptome / annotation version | Not yet verified; scRNA-seq panel-selection provenance does not establish a design annotation release. |
| Target genes | 2,766 (same as D04) |
| Target regions | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Genomic/transcript coordinates | Not yet verified: assembly/accession, transcript version, strand, origin and interval convention remain unresolved. |
| Target-binding sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Full probe sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Barcode assignments | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Readout sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Primer/linker sequences | Complete deposit inventory not yet examined; presence or absence of associated sequence/codebook resources is unverified. |
| Design software | Exact implementation for this dataset unverified. Later S03 describes Picky 2.2 + Dedupe/BBMap + NUPACK + BLAST; S05 equivalence and original-panel use are unverified. |
| Design parameters | Exact parameters for this dataset unverified. Later S03 recommendations are separately scoped in probes.md; they are not original-panel ground truth. |
| Source | Zenodo 10.5281/zenodo.7332091 (S10) |
| Availability | https://zenodo.org/record/7332091 ; existence confirmed; full file inventory not captured for this draft |
| Missing information | Same as D04 |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned: direct content/inventory and reference context have not been sufficiently inspected. File size is not interpreted as evidence of contents or reproducibility. |
| Comparison scope | Candidate for future assessment; strict probe-design comparison is not currently established. |
| Assessment date | 2026-10-15 (unassessed after Fix; no grade awarded) |
| Evidence status | reported |
