# Benchmark Resources

## Benchmark Summary

Highest justified probe-design candidate grade: **C**. Public sequence resources exist, but reference identity and historical design provenance are not sufficiently resolved for a strict U-Probe design-output benchmark. This grade concerns design reconstruction, not the quality of the imaging method or the accuracy of the published oligo strings. [inferred; S02/S05/S07; assessed 2026-10-15]

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A. Consider target/reference context and relevant design metadata. Unassessed is not D.

## Candidate Datasets

### D01 — Published 30/31-panel probe resources

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Published 30/31-panel probe resources |
| Method version | V01 |
| Species | Mus musculus; Homo sapiens; HBX viral target separately identified |
| Cell / tissue / sample | Mouse brain and embryo; human HCC |
| Genome build | Not resolved in examined sources |
| Transcriptome / annotation version | Not resolved; original transcript accession versions not supplied in inspected probe sheets |
| Target genes | 30 brain records; 30 embryo records; 91 HCC oligos representing a 31-panel set including HBX |
| Target regions | Preprint 40-nt sites; final HBX residual target-arm total 28 nt |
| Genomic/transcript coordinates | No source-confirmed assembly/accession, strand and interval convention recovered |
| Target-binding sequences | Exact residual terminal sequences can be separated by matching listed backbones; biological mapping unverified |
| Full probe sequences | 30+30 mouse 122-nt oligos; 90 HCC 122-nt + HBX 110-nt oligo |
| Barcode assignments | 31 codewords/backbones; sample padlock assignment recoverable by exact full-backbone match |
| Readout sequences | 17 20-nt imaging-probe identities; theoretical and empirical mixtures |
| Primer/linker sequences | Five modified RCA primers; 2-nt V01 prefix role unresolved |
| Design software | Current author-linked probe-designer public; original producing revision unknown |
| Design parameters | Preprint criteria available; original per-probe settings and reference snapshot unresolved |
| Source | S02, Imaging probe（31 panel）, RCA Primer (31 panel), three V01 padlock sheets; S03 dilution table |
| Availability | Public files downloaded 2026-10-15 |
| Missing information | Original reference/annotation and transcript versions; coordinates; producing code revision; one-versus-two mouse-target discrepancy |
| Benchmark grade | C |
| Grade justification | Many exact oligos and readout assignments are available, but missing biological reference and historical selection inputs are substantial; Grade B would understate these gaps. Suitable for architecture/sequence-handling comparisons, not strict selection-output ground truth. |
| Comparison scope | Parse and compare published oligo strings, backbone assignment and component handling; target-site selection requires independent remapping with explicit new reference provenance |
| Assessment date | 2026-10-15 |
| Evidence status | Resource contents reported; counts/segmentation derived; grade inferred from completeness audit |

### D02 — Published 64-panel mouse embryo probe resources

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | Published 64-panel mouse embryo probe resources |
| Method version | V02 |
| Species | Mus musculus |
| Cell / tissue / sample | Mouse embryo; worksheet panel; do not infer brain oligos from embryo list |
| Genome build | Not resolved in examined sources |
| Transcriptome / annotation version | Not resolved; original transcript accession versions not supplied in inspected probe sheets |
| Target genes | 64 listed padlocks; names in Padlock probe-mouse embryo(64) |
| Target regions | 20-nt residual flanks on each side of matched 80-nt backbone (derived) |
| Genomic/transcript coordinates | No source-confirmed assembly/accession, strand and interval convention recovered |
| Target-binding sequences | Exact residual terminal sequences can be separated by matching listed backbones; biological mapping unverified |
| Full probe sequences | 64 × 120-nt exact listed sequences |
| Barcode assignments | 64 backbones/codewords; high-G digit 5 in workbook versus high-level label 2 in article |
| Readout sequences | 38 sequence identities (two contain leading U+FEFF); separate high-green R/B/Y families; empirical mixtures |
| Primer/linker sequences | Twelve modified RCA primers; no separate universal linker established |
| Design software | Current author-linked probe-designer public; original producing revision unknown |
| Design parameters | Preprint criteria available; original per-probe settings and reference snapshot unresolved |
| Source | S02, Imaging probe (64 panel), RCA primer (64), Padlock probe-mouse embryo(64); S01 Extended Data Fig. 7 |
| Availability | Public files downloaded 2026-10-15 |
| Missing information | Original reference/annotation/transcript versions; coordinate conventions; exact producing code and final full-method execution; representation correspondence for high-G labels |
| Benchmark grade | C |
| Grade justification | The final workbook supports detailed architecture reconstruction, but lacks original biological reference context sufficient to reproduce target selection. A finite probe list alone is not Grade A/B design ground truth. |
| Comparison scope | Architecture and codebook/oligo parsing; preserve raw source labels and normalization history |
| Assessment date | 2026-10-15 |
| Evidence status | Resource contents reported; counts/segmentation derived; grade inferred from completeness audit |

### D03 — PRISM analysis and image resources

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | PRISM analysis and image resources |
| Method version | V01/V02, run-specific |
| Species | Mouse and human across paper-linked deposits |
| Cell / tissue / sample | Mouse brain, embryo, HCC |
| Genome build | Not resolved in examined sources |
| Transcriptome / annotation version | Not resolved; original transcript accession versions not supplied in inspected probe sheets |
| Target genes | Run-specific gene expression/spot data; not independently inventoried as probe-set membership |
| Target regions | Imaging spots; not target transcript coordinates |
| Genomic/transcript coordinates | No source-confirmed assembly/accession, strand and interval convention recovered |
| Target-binding sequences | Exact residual terminal sequences can be separated by matching listed backbones; biological mapping unverified |
| Full probe sequences | Not established from this image/analysis resource |
| Barcode assignments | Requires independent association to S02 and run metadata |
| Readout sequences | Not a replacement for S02 |
| Primer/linker sequences | Not a replacement for S02 |
| Design software | Current author-linked probe-designer public; original producing revision unknown |
| Design parameters | Preprint criteria available; original per-probe settings and reference snapshot unresolved |
| Source | S09, Zenodo 12755414 metadata and PRISM_data_publication_analysis.md; S01 Data availability |
| Availability | Public files downloaded 2026-10-15 |
| Missing information | Original design inputs/outputs and reference metadata; archive contents not downloaded/verified |
| Benchmark grade | D |
| Grade justification | For the explicitly assessed scope of probe-design reconstruction, metadata for image and downstream analysis archives cannot establish the design. This is not a negative grade for image-analysis benchmarking. |
| Comparison scope | Potential image/decoding benchmarking; unassessed for reproducibility of that separate scope |
| Assessment date | 2026-10-15 |
| Evidence status | Resource contents reported; counts/segmentation derived; grade inferred from completeness audit |
