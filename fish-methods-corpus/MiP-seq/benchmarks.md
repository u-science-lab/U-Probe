# Benchmark Resources

## Benchmark Summary

No candidate datasets have been formally graded. Public probe sequence data (D01) is reported as available in supplementary materials, but direct access and assessment of its completeness for ground-truth comparison was not performed during this research phase. Image data (D02) is publicly available for spatial analysis benchmarking.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A. Consider target/reference context and relevant design metadata. Unassessed is not D.

## Candidate Datasets

### D01 — Supplementary Probe Sequences

| Field | Value |
|-------|-------|
| Dataset ID | D01 |
| Dataset name | S01 Supplementary Data 1: Probe sequences (Table S1 in the preprint S02) |
| Method version | V01 |
| Species | Mouse (Mus musculus), Rat (Rattus norvegicus), Human (Homo sapiens), Zebrafish (Danio rerio) |
| Cell / tissue / sample | Brain tissue (hypothalamus, hippocampus, PVN, cerebellum), cultured cells |
| Genome build | Not explicitly reported |
| Transcriptome / annotation version | Not explicitly reported |
| Target genes | Multiple genes; examples include Cck, Penk, Crh, Oxt, Avp, Gad1, Npy, Kif5a, Snhg11, Tau, and 100+ genes in HDB experiments |
| Target regions | Exon regions near 5' end of transcripts |
| Genomic/transcript coordinates | Not provided in main text; expected in supplementary table |
| Target-binding sequences | Padlock arms and primer target parts, in S01 Supplementary Data 1 (41551_2024_1205_MOESM5_ESM.xlsx), sheets 'Fig 1' to 'Fig6'; see probes.md (examined 2026-10-20) |
| Full probe sequences | Padlock probes and initiator primers, in S01 Supplementary Data 1 (41551_2024_1205_MOESM5_ESM.xlsx), sheets 'Fig 1' to 'Fig6'; see probes.md (examined 2026-10-20) |
| Barcode assignments | Not listed separately; B1 and B2 (3 nt each) are read from the padlock sequences (derived; see probes.md) |
| Readout sequences | Fixed; not dataset-specific |
| Primer/linker sequences | Initiator primers in the per-figure sheets; anchors, query and detection probes in sheet 'In Situ Sequencing probe'; see probes.md |
| Design software | Self-developed probe searching program (specific name not reported) |
| Design parameters | Arm length: 12-16 bp; 2 target sites per gene; exon-proximal selection |
| Source | S01, Supplementary Data |
| Availability | Available as supplementary material with S01 |
| Missing information | Genome build, annotation version, target coordinates, design criteria, codebook mapping |
| Benchmark grade | Unassessed |
| Grade justification | Probe sequences are reported as available in supplementary materials (S01), but direct access to the supplementary table was not performed during this research phase. Without verifying the actual probe sequences, codebook assignments, and reference genome/annotation context, formal grade assignment is not possible. |
| Comparison scope | Probe design algorithm comparison, target region selection validation, barcode assignment verification |
| Assessment date | Not yet assessed |
| Evidence status | reported (existence of data); unavailable (actual content verification) |

### D02 — Source Image Data

| Field | Value |
|-------|-------|
| Dataset ID | D02 |
| Dataset name | Source image data from MiP-seq study |
| Method version | V01 |
| Species | Mouse, Rat, Zebrafish, Human (cell lines) |
| Cell / tissue / sample | Brain tissue, cultured cells, zebrafish brain |
| Genome build | Not applicable for imaging data |
| Transcriptome / annotation version | Not applicable for imaging data |
| Target genes | Multiple; 10-217 genes depending on experiment |
| Target regions | Subcellular localization (cytoplasmic/nuclear) |
| Genomic/transcript coordinates | Decoded from barcodes; original coordinates require codebook |
| Target-binding sequences | Not applicable (imaging data) |
| Full probe sequences | Not applicable (imaging data) |
| Barcode assignments | Encoded in decoded gene identities |
| Readout sequences | Fixed sequencing chemistry |
| Primer/linker sequences | Not applicable (imaging data) |
| Design software | InnerEye (https://github.com/GangCaoLab/InnerEye) for image processing |
| Design parameters | Not applicable for raw imaging data |
| Source | Zenodo record 10837475 |
| Availability | Publicly available at https://zenodo.org/records/10837475 |
| Missing information | Codebook linking barcodes to genes, raw decoded counts, cell segmentation masks |
| Benchmark grade | Unassessed |
| Grade justification | This dataset provides source imaging data but is not directly suitable as a probe-design benchmark. It would be valuable for benchmarking the image processing and decoding pipeline, but not for probe-design comparison. |
| Comparison scope | Image processing accuracy, decoding accuracy, spatial analysis validation |
| Assessment date | Not yet assessed |
| Evidence status | reported (existence and availability of data) |

### D03 — RNA-seq Validation Data (GSE190219)

| Field | Value |
|-------|-------|
| Dataset ID | D03 |
| Dataset name | RNA-seq data for validation |
| Method version | V01 |
| Species | Mouse, Rat |
| Cell / tissue / sample | C6 cells, offspring brain tissue |
| Genome build | Not reported |
| Transcriptome / annotation version | Not reported |
| Target genes | Multiple; used for identifying allele-specific SNPs |
| Target regions | Transcriptomic |
| Genomic/transcript coordinates | Not provided |
| Target-binding sequences | Not applicable (RNA-seq data) |
| Full probe sequences | Not applicable |
| Barcode assignments | Used to design allele-specific probes |
| Readout sequences | Not applicable |
| Primer/linker sequences | Not applicable |
| Design software | Standard RNA-seq analysis |
| Design parameters | Standard RNA-seq |
| Source | GEO GSE190219 |
| Availability | Publicly available at GEO |
| Missing information | Specific processing pipeline parameters, raw counts matrix |
| Benchmark grade | Unassessed |
| Grade justification | This RNA-seq data was used for identifying SNPs and validating MiP-seq detection, but is not itself a probe-design benchmark. |
| Comparison scope | Allele frequency validation, expression level comparison |
| Assessment date | Not yet assessed |
| Evidence status | reported (existence of data) |
