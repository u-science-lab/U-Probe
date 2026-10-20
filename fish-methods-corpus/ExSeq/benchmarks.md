# Benchmark Resources

## Benchmark Summary

D01 is assessed **C** for full probe-selection/assembly comparison on 2026-10-15. Exact example output is useful, but historical off-target reference provenance remains material. This is the highest justified candidate grade, consistent with metadata.yaml. D02 remains unassessed pending original workbook inspection. D03-D07 are comparative tables, sequencing data or software/analysis resources, not established probe-design ground truth.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published sequence list is not automatically Grade A. Unassessed is not D. Original in situ chemistry and analysis references do not establish the historical design database.

## Candidate Datasets

### D01 — TargetedExSeqProbeGen worked example

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | TargetedExSeqProbeGen worked example |
| Method version | V04; pinned S04 revision |
| Species | Homo sapiens (design inputs) |
| Cell / tissue / sample | Not applicable (design-only example). |
| Genome build | Not established; inputs are transcript accessions, not a genome-coordinate design. |
| Transcriptome / annotation version | Input targets NM_001014432.1, NM_001799.4, NM_000422.3; complete historical off-target RefSeq release unresolved. |
| Target genes | AKT1, CDK7, KRT17; input Sheet1 A1:D3. |
| Target regions | 32-nt transcript segments; output has 16/13/12 probes respectively (41 total, derived). |
| Genomic/transcript coordinates | Homology FASTA headers record one-based forward-transcript start i; default interval [i,i+31] closed, inferred from Homology.m indexing. Exported DNA is reverse complement. Genomic mapping not reconstructed. |
| Target-binding sequences | Retrieved from example/example_output/HomologyRegionsExample/*.fa; already reverse-complemented DNA. |
| Full probe sequences | Retrieved ExampleOutputAllProbes.xlsx Sheet1 column D; 41 nonempty 74-nt sequences after excluding /5Phos/ notation (derived). |
| Barcode assignments | Input Sheet1 D1:D3: 3212010 / 2102033 / 0221001. Output names in column B specify embedded_len4. |
| Readout sequences | Physical barcodes in full sequences; Illumina AKT1 TGCGACA verified in D1 and by direct map; SOLiD conversion uses embedded_len4 specification. |
| Primer/linker sequences | Backbone and two AAA spacers in sequence cells and assembly line 482; RCA/SBS reagents separately in S13 Step 5 and S14 materials/Step 10. |
| Design software | S04 TargetedExSeqProbeGen at 2fe7748277df9d9f92ea4fe6a1d05e98796c527a. |
| Design parameters | V04 defaults and gene-specific Tm 55/53/52 °C; detailed rule/source table in probes.md. Search spacing differs from original brain rule. |
| Source | S04 input/output workbook, homology FASTAs, Homology.m and Assembly.m; S12 Steps 9-13. |
| Availability | Public files retrieved 2026-10-15; output workbook SHA-256 4e3de1bdcd60e9991a5a31e35a47386dca0c8bf0e0d5e4ade4a7a05396f43071. |
| Missing information | Historical complete off-target database release and environment needed to reproduce exact selection have not been reconstructed; no full MATLAB pipeline rerun performed. |
| Benchmark grade | C |
| Grade justification | Useful architecture, exact example outputs and code are available, but missing historical specificity-reference provenance is material to exact full selection reproduction. Public sequences alone do not support A/B. |
| Comparison scope | Full target-region selection plus probe assembly comparison; a narrower assembly-only benchmark has not been separately graded. |
| Assessment date | 2026-10-15 |
| Evidence status | reported file values; derived counts/lengths; inferred coordinate/code interpretation; assessment is an evidence-based corpus judgement. |

### D02 — Original targeted probe resources, Tables S9/S10

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | Original targeted probe resources, Tables S9/S10 |
| Method version | V02 brain and V03 cancer; retain experiment subdivisions. |
| Species | Mouse (brain) and human (cancer). |
| Cell / tissue / sample | 42-gene visual cortex; hippocampus initially 35 selected genes with Rgs5 excluded in silico (34 retained); 297-gene metastatic breast cancer panel. |
| Genome build | Not yet established for targeted design; do not reuse V01 mm10 analysis reference. |
| Transcriptome / annotation version | RefSeq GenBank inputs reported; individual workbook accessions/releases not yet examined. |
| Target genes | S9 genes/accessions/barcodes, explicitly directed from original brain and cancer methods; individual cells not yet examined. |
| Target regions | Original 32-mer recognition regions; original cortex maximum 8 probes/transcript, hippocampus maximum 16; cancer 9-16. |
| Genomic/transcript coordinates | Not yet examined in workbook; no genome-coordinate convention assumed. |
| Target-binding sequences | Probe resource reported in S10; separate column structure not yet examined. |
| Full probe sequences | S10 caption and methods report probe sequences; individual synthesis notation/cells not yet examined. |
| Barcode assignments | S9 reported; V02 four-symbol checksum; V03 seven-symbol distance-3; cell-level mapping pending. |
| Readout sequences | Physical barcode construction described in original methods; workbook parsing pending. |
| Primer/linker sequences | S02 pp. 19-20 report backbone and AAA linkers; full workbook notation pending. |
| Design software | Original implementation; later S04 public code is not assumed identical. |
| Design parameters | V02 cortex Tm initially 55 °C, lowered to minimum 52 °C if fewer than six probes/gene; hippocampus 50 °C; GC 40-65% per arm; next accepted search starts 5 nt past prior end. Preserve original ambiguous Tm inequality (probes.md). V03 inherits original design with stated encoding exceptions. |
| Source | S02 supplementary pp. 18-20, 23, 29-31 and p. 90 captions; NIHMS1670746-supplement-Supp_Tables.xlsx. |
| Availability | Public link/captions verified 2026-10-15; attempted workbook retrieval returned HTML challenge; individual cells not yet examined. |
| Missing information | Cell-level sequence/accession/coordinate inspection; historical database release; resolution of original Tm wording. |
| Benchmark grade | Unassessed. |
| Grade justification | Positive evidence for public probe resources, but workbook-level completeness is not yet assessed. No A/B claim. |
| Comparison scope | Potential experiment-specific target selection and full-probe comparison after workbook inspection. |
| Assessment date | 2026-10-15 (routing review only; no grade assigned). |
| Evidence status | reported original methods/table captions; workbook contents not yet examined. |

### D03 — Original comparative technology/yield tables S7/S8

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | Original comparative technology/yield tables S7/S8 |
| Method version | Original publication context (V02/V03). |
| Species | Varies across technologies; not a single probe collection. |
| Cell / tissue / sample | Comparative studies, not a unified design experiment. |
| Genome build | Not applicable (comparison tables, not probe designs). |
| Transcriptome / annotation version | Not applicable (comparison tables, not probe designs). |
| Target genes | Not applicable (comparison tables, not probe designs). |
| Target regions | Not applicable (comparison tables, not probe designs). |
| Genomic/transcript coordinates | Not applicable (comparison tables, not probe designs). |
| Target-binding sequences | Not applicable (comparison tables, not probe designs). |
| Full probe sequences | Not applicable (comparison tables, not probe designs). |
| Barcode assignments | Not applicable (comparison tables, not probe designs). |
| Readout sequences | Not applicable (comparison tables, not probe designs). |
| Primer/linker sequences | Not applicable (comparison tables, not probe designs). |
| Design software | Not applicable (comparison tables, not probe designs). |
| Design parameters | Not applicable (comparison tables, not probe designs). |
| Source | S02 supplementary pp. 86-89, Table S7/S8 captions and notes. |
| Availability | Captions/notes inspected in author-hosted supplement, 2026-10-15. |
| Missing information | No probe-design ground truth is established by these comparison tables; actual original probe resources route to D02. |
| Benchmark grade | Unassessed. |
| Grade justification | Reclassified resource; no probe-design grade applies to a technology/yield comparison table. |
| Comparison scope | Technology/yield context only. |
| Assessment date | 2026-10-15 (resource-type correction). |
| Evidence status | reported table captions; resource classification inferred. |

### D04 — Raw Illumina sequencing archive PRJNA663046

| Field | Value |
|---|---|
| Dataset ID | D04 |
| Dataset name | Raw Illumina sequencing archive PRJNA663046 |
| Method version | V01 ex situ sequencing; run-level coverage not yet examined. |
| Species | Not yet examined at run level. |
| Cell / tissue / sample | Untargeted ex situ archive cited by original paper; run/sample mapping pending. |
| Genome build | V01 downstream annotation uses GRCm38/mm10 (S02 p. 11); do not assume every deposited run has this context. |
| Transcriptome / annotation version | V01 downstream mouse RefSeq downloaded from UCSC 2018-08-14; archive run-level metadata pending. |
| Target genes | Not applicable (untargeted sequencing). |
| Target regions | Not applicable to this resource as probe-design ground truth. |
| Genomic/transcript coordinates | Not applicable to this resource as probe-design ground truth. |
| Target-binding sequences | Not applicable to this resource as probe-design ground truth. |
| Full probe sequences | Not applicable to this resource as probe-design ground truth. |
| Barcode assignments | Not applicable to this resource as probe-design ground truth. |
| Readout sequences | Not applicable to this resource as probe-design ground truth. |
| Primer/linker sequences | Not applicable to this resource as probe-design ground truth. |
| Design software | Not applicable (raw read archive). |
| Design parameters | Not applicable to this resource as probe-design ground truth. |
| Source | S01 Data and materials availability; S06 BioProject PRJNA663046. |
| Availability | Raw Illumina reads reported deposited; specific SRA runs/files not yet examined. |
| Missing information | Run inventory and sample mapping; this accession is not established as a microscope-image archive. |
| Benchmark grade | Unassessed. |
| Grade justification | Raw sequencing resource, not a published probe-design output set. |
| Comparison scope | Ex situ sequence-analysis validation after run-level inspection; not direct image-pipeline comparison. |
| Assessment date | 2026-10-15 (resource-type review). |
| Evidence status | reported deposition in S01; run-level content not yet examined. |

### D05 — ExSeqProcessing v1.0.0 software release

| Field | Value |
|---|---|
| Dataset ID | D05 |
| Dataset name | ExSeqProcessing v1.0.0 software release |
| Method version | Analysis software associated with original ExSeq. |
| Species | Not applicable (software). |
| Cell / tissue / sample | Not applicable (software). |
| Genome build | Not applicable (software). |
| Transcriptome / annotation version | Not applicable (software). |
| Target genes | Not applicable (software). |
| Target regions | Not applicable to this resource as probe-design ground truth. |
| Genomic/transcript coordinates | Not applicable to this resource as probe-design ground truth. |
| Target-binding sequences | Not applicable to this resource as probe-design ground truth. |
| Full probe sequences | Not applicable to this resource as probe-design ground truth. |
| Barcode assignments | Not applicable to this resource as probe-design ground truth. |
| Readout sequences | Not applicable to this resource as probe-design ground truth. |
| Primer/linker sequences | Not applicable to this resource as probe-design ground truth. |
| Design software | ExSeqProcessing is analysis software, not probe design. |
| Design parameters | Not applicable to this resource as probe-design ground truth. |
| Source | S07 Zenodo 10.5281/zenodo.4075515, description and Files. |
| Availability | Release landing page/file inventory inspected 2026-10-15; dgoodwin208/ExSeqProcessing-v1.0.0.zip, 6.5 MB. Tutorial/sample data linked externally to wiki. |
| Missing information | Tutorial data files and their independent provenance not yet inventoried; do not assert the software DOI contains the tutorial images. |
| Benchmark grade | Unassessed. |
| Grade justification | Software release, not a probe dataset. |
| Comparison scope | Software reproducibility; tutorial data require a separately verified resource record. |
| Assessment date | 2026-10-15 (resource-type review). |
| Evidence status | reported release description/file inventory. |

### D06 — ExSeq 3D visualization tool

| Field | Value |
|---|---|
| Dataset ID | D06 |
| Dataset name | ExSeq 3D visualization tool |
| Method version | Original ExSeq visualization. |
| Species | Not applicable (tool). |
| Cell / tissue / sample | Not applicable (tool). |
| Genome build | Not applicable (tool). |
| Transcriptome / annotation version | Not applicable (tool). |
| Target genes | Not applicable (tool). |
| Target regions | Not applicable to this resource as probe-design ground truth. |
| Genomic/transcript coordinates | Not applicable to this resource as probe-design ground truth. |
| Target-binding sequences | Not applicable to this resource as probe-design ground truth. |
| Full probe sequences | Not applicable to this resource as probe-design ground truth. |
| Barcode assignments | Not applicable to this resource as probe-design ground truth. |
| Readout sequences | Not applicable to this resource as probe-design ground truth. |
| Primer/linker sequences | Not applicable to this resource as probe-design ground truth. |
| Design software | Visualization tool; not probe design. |
| Design parameters | Not applicable to this resource as probe-design ground truth. |
| Source | S01 Data and materials availability; S09 Zenodo 10.5281/zenodo.4075649. |
| Availability | Deposition reported by S01; file inventory not yet examined. |
| Missing information | Archive file contents not yet examined. |
| Benchmark grade | Unassessed. |
| Grade justification | Visualization software is not probe-design ground truth. |
| Comparison scope | Visualization only. |
| Assessment date | 2026-10-15 (resource-type review). |
| Evidence status | reported deposition in S01; archive contents not yet examined. |

### D07 — Expansion Sequencing HTAPP Data Resource

| Field | Value |
|---|---|
| Dataset ID | D07 |
| Dataset name | Expansion Sequencing HTAPP Data Resource |
| Method version | V03 cancer analysis. |
| Species | Human. |
| Cell / tissue / sample | Metastatic breast cancer biopsy, Fig. 6. |
| Genome build | Not yet examined in deposited analysis files. |
| Transcriptome / annotation version | Not yet examined in deposited analysis files. |
| Target genes | GeneNames-1222.txt is listed; individual entries not yet examined. |
| Target regions | Not applicable to this resource as probe-design ground truth. |
| Genomic/transcript coordinates | Not applicable to this resource as probe-design ground truth. |
| Target-binding sequences | Not applicable to this resource as probe-design ground truth. |
| Full probe sequences | Not applicable to this resource as probe-design ground truth. |
| Barcode assignments | Not applicable to this resource as probe-design ground truth. |
| Readout sequences | Not applicable to this resource as probe-design ground truth. |
| Primer/linker sequences | Not applicable to this resource as probe-design ground truth. |
| Design software | Not a design pipeline; extracted reads and analysis scripts. |
| Design parameters | Not applicable to this resource as probe-design ground truth. |
| Source | S08 Zenodo 10.5281/zenodo.4479018, description and Files. |
| Availability | Description/file inventory inspected 2026-10-15; extracted transcript MAT files, expression matrices, scripts/notebook, segmentation and downsampled image files. Raw image terabytes available on request per deposit description. |
| Missing information | Large data files not downloaded; raw imaging is not an openly inventoried full archive here. |
| Benchmark grade | Unassessed. |
| Grade justification | Spatial read/analysis resource; not automatically full probe-design ground truth. |
| Comparison scope | Cancer spatial analysis reproduction; original probe resources separately route to D02. |
| Assessment date | 2026-10-15 (resource-type review). |
| Evidence status | reported deposit description/file inventory; large-file contents not yet examined. |
