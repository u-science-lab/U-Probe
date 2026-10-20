# Sources

## Primary Publications

### S01 — Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution

| Field | Value |
|---|---|
| Source ID | S01 |
| Source type | Original peer-reviewed publication; PMC author manuscript |
| Citation | Cheng Y. et al. Cell 188(24), 6953 onward (2025); online 2025-10-01, issue 2025-11-26. |
| Title | Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution |
| Authors | Yubao Cheng; Shengyuan Dang; Yuan Zhang; Yanbo Chen; Ruihuan Yu; Miao Liu; Shengyan Jin; Ailin Han; Samuel Katz; Siyuan Wang |
| Year | 2025 |
| DOI | 10.1016/j.cell.2025.09.006 |
| PMID | 41038164 |
| Official URL | https://doi.org/10.1016/j.cell.2025.09.006 |
| Repository / accession | PMC12662569; https://pmc.ncbi.nlm.nih.gov/articles/PMC12662569/ |
| Method version | V01/V02 |
| Relevant figures | Figure 1A–C; Figure 7; Figure S1 legend |
| Relevant tables | Key resources table; Table S1 link |
| Relevant supplementary files | PMC supplement-5.xlsx is Table S1; publisher equivalent is mmc1.xlsx. Figure supplements are separately numbered on PMC. |
| Information supported | Method identity, explicit 5′→3′ V01 architectures, library construction, design references, code schemes, sample-specific protocols, imaging, decoding and availability links. |
| Notes / limitations | PMC title omits the hyphen in whole genome; publisher/PubMed use whole-genome. STAR Methods has unresolved 94-bit/96-dimensional wording and splint residue-count wording. Manufacturer/template sequences must not be treated as final working probes. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Local evidence/article.html and article.txt; SHA-256 in evidence/checksums.json. |

### S06 — Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue

| Field | Value |
|---|---|
| Source ID | S06 |
| Source type | Original author preprint; discovery/version context |
| Citation | Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue |
| Title | Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue |
| Authors | Yubao Cheng; Shengyuan Dang; Yuan Zhang; Yanbo Chen; Ruihuan Yu; Miao Liu; Shengyan Jin; Ailin Han; Samuel Katz; Siyuan Wang |
| Year | 2025 |
| DOI | 10.1101/2025.03.06.641951 |
| PMID | 40161724 (preprint record, not final article PMID) |
| Official URL | https://www.biorxiv.org/content/10.1101/2025.03.06.641951v1 |
| Repository / accession | PMC11952344 |
| Method version | V01/V02 |
| Relevant figures | Not used as final-architecture evidence. |
| Relevant tables | Search excerpts identify separate human/mouse supplementary tables; not substituted for final S02. |
| Relevant supplementary files | Not applicable or not established; see limitations. |
| Information supported | Earlier publication stage/title and terminology; distinction from final article. |
| Notes / limitations | Only indexed abstract/results snippets and supplementary listing excerpts inspected in this research pass. Full preprint/supplement not independently downloaded. No final probe/protocol values are taken solely from S06. |
| Access date | 2026-10-15 |
| Release / commit / checksum | bioRxiv v1, posted March 2025; no local full-text snapshot. |

## Supplementary Materials

### S02 — Table S1. All oligos, related to all Figures and STAR Methods

| Field | Value |
|---|---|
| Source ID | S02 |
| Source type | Official supplementary spreadsheet |
| Citation | Supplement to S01; final publication Table S1. |
| Title | Table S1. All oligos, related to all Figures and STAR Methods |
| Authors | Yubao Cheng; Shengyuan Dang; Yuan Zhang; Yanbo Chen; Ruihuan Yu; Miao Liu; Shengyan Jin; Ailin Han; Samuel Katz; Siyuan Wang |
| Year | 2025 |
| DOI | Parent: 10.1016/j.cell.2025.09.006 |
| PMID | Parent: 41038164 |
| Official URL | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425010372-mmc1.xlsx |
| Repository / accession | PMC equivalent NIHMS2110109-supplement-5.xlsx; publisher mmc1.xlsx |
| Method version | V01/V02 |
| Relevant figures | All figures (publisher designation). |
| Relevant tables | All nine worksheets; row/column locators in probes.md and benchmarks.md. |
| Relevant supplementary files | publisher_mmc1.xlsx, 7,427,980 bytes. |
| Information supported | Exact ordered templates; endogenous and perturbation codebooks; adapters, dyes, blocking oligos; gRNA construction templates; library amplification primers. |
| Notes / limitations | Template—not uniformly working-probe—sequences. Numeric Ensembl_GeneID headers do not resolve stable accessions. Mouse paired probes may be duplicated. No verified reference-linked target coordinates. Upper/lower case and modification notation preserved. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Original workbook retained at evidence/publisher_mmc1.xlsx; SHA-256 in evidence/checksums.json; read-only inspection. |

## Official Protocols

S01 STAR Methods contains the examined official experimental protocols and probe-synthesis summary. No separate author-endorsed standalone RAEFISH protocol was established in the bounded search.

## Official Repositories and Software

### S03 — RAEFISH — Codes and data

| Field | Value |
|---|---|
| Source ID | S03 |
| Source type | Official author laboratory resource page and linked archive |
| Citation | RAEFISH — Codes and data |
| Title | RAEFISH — Codes and data |
| Authors | Siyuan Wang Lab |
| Year | Page date not specified. |
| DOI | Not applicable or not established; see limitations. |
| PMID | Not applicable or not established; see limitations. |
| Official URL | https://campuspress.yale.edu/wanglab/raefish/ |
| Repository / accession | http://192.99.37.111/downloads/RAEFISH/DataAndCodes.zip |
| Method version | V01/V02 |
| Relevant figures | Not applicable. |
| Relevant tables | ZIP central-directory inventory. |
| Relevant supplementary files | DataAndCodes.zip; total 696,162,447 bytes (HTTP Content-Range). |
| Information supported | Author endorsement of data/code resource; MATLAB analysis script; human transcript/target sequence MAT resources with numeric ENSG/ENST fields; scoped correspondence to final templates. |
| Notes / limitations | Page and 92-entry ZIP directory inspected; four Design240427/Design240516 GoodGeneList MAT files and RAEFISH/runM2_2_HPC.m retrieved by byte ranges and inspected. Full ZIP not downloaded. No probe-design implementation established. Initial sample retrieval approval timed out; retry succeeded. HTTP URL is author-linked. The S03 author archive contains four inspected GoodGeneList MAT files under Design240427/Design240516. They contain numeric ENSG/ENST/LongestENST fields, transcript Sequence strings and C/L binding sequences. Comparison by symbol plus both binding sequences matches 23,303 distinct published human sequence-pair records (derived: evidence/design_resource_inspection.json); this does not establish versioned stable transcript accessions, a reference build or a probe-selection algorithm. The human workbook has 23,312 non-control rows but 23,300 unique symbols; mouse has 21,955 rows but 21,954 unique symbols. Duplicate symbols must not be collapsed during extraction. The sampled mouse decoder RAEFISH/runM2_2_HPC.m sets RoundsOfHybs=94 and reshapes vectors with that size. This supports a 94-dimensional implementation for this archive example (inferred from code), while the publication's 96-dimensional wording remains an explicit source discrepancy. |
| Access date | 2026-10-15 |
| Release / commit / checksum | HTTP Last-Modified: Sat, 17 May 2025 19:35:16 GMT; no tag/commit established; local archive_inventory.json and yale.html hashed. |

## Official Datasets

### S04 — Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution

| Field | Value |
|---|---|
| Source ID | S04 |
| Source type | Author-deposited dataset |
| Citation | Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution |
| Title | Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution |
| Authors | Yubao Cheng; Shengyuan Dang; Yuan Zhang; Yanbo Chen; Ruihuan Yu; Miao Liu; Shengyan Jin; Ailin Han; Samuel Katz; Siyuan Wang |
| Year | 2025 |
| DOI | 10.17632/8kbv637pxh.2; paper-linked v1: 10.17632/8kbv637pxh.1 |
| PMID | Not applicable or not established; see limitations. |
| Official URL | https://data.mendeley.com/datasets/8kbv637pxh/2 |
| Repository / accession | 8kbv637pxh, versions 1 and 2; metadata endpoint https://data.mendeley.com/public-api/datasets/8kbv637pxh |
| Method version | V01/V02 |
| Relevant figures | Not applicable or not established; see limitations. |
| Relevant tables | API files list (23 entries in v2); codebooks and processed counts/cell lists. |
| Relevant supplementary files | Codebook_StevenRaeFISH23K_LabA549.mat; Codebook_RaeFISH.mat; other MAT files enumerated in evidence/mendeley.json. |
| Information supported | Versioned public data deposit, contributors, file names/sizes and provider SHA-256 values. |
| Notes / limitations | Version 1 landing page inspected through web; API returned version 2. Do not assume versions or Yale ZIP are identical. Binary MAT data not inspected; provider metadata checksum is not a locally verified file checksum. |
| Access date | 2026-10-15 |
| Release / commit / checksum | v1 landing page reports 2025-07-15; current API version 2. Raw API response saved and hashed; file checksums within API are provider supplied. |

### S05 — Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue

| Field | Value |
|---|---|
| Source ID | S05 |
| Source type | Author-deposited GEO series |
| Citation | Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue |
| Title | Sequencing-free whole genome spatial transcriptomics at molecular resolution in intact tissue |
| Authors | Yubao Cheng; Shengyuan Dang; Yuan Zhang; Yanbo Chen; Ruihuan Yu; Miao Liu; Shengyan Jin; Ailin Han; Samuel Katz; Siyuan Wang |
| Year | 2025 |
| DOI | Not applicable or not established; see limitations. |
| PMID | Not applicable or not established; see limitations. |
| Official URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302564 |
| Repository / accession | GSE302564 |
| Method version | V01/V02 |
| Relevant figures | Not applicable or not established; see limitations. |
| Relevant tables | Series summary and linked sample/resource listing. |
| Relevant supplementary files | Not applicable or not established; see limitations. |
| Information supported | Public human bulk RNA-seq comparator deposit; series publication status 2025-07-22. |
| Notes / limitations | This is not the RAEFISH probe library. Sequencing-based comparator data must not be described as the sequencing-free assay output. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Retrieved geo.html; SHA-256 in evidence/checksums.json. |

## Later Developer Publications

No later method-changing developer publication was established in the RAEFISH searches performed on 2026-10-15. This does not prove universal absence.

## Secondary Sources

Search results included news/review summaries; none is used as primary support for architecture, sequence, protocol or benchmark claims.
