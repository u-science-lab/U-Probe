# Sources

All access dates are 2026-10-15. Source records distinguish actual inspection from links and metadata; see research_manifest.md. Evidence states attach to claims in production files, not merely to a bibliography entry.

## Primary Publications

### S01 — Original Oligo-LiveFISH publication

| Field | Value |
|---|---|
| Source ID | S01 |
| Source type | Original peer-reviewed publication; full author manuscript on PMC plus PubMed bibliography |
| Citation | Zhu Y. et al. Cell 188(12), 3310–3328.e27 (2025); online April 15, issue June 12 |
| Title | High-resolution dynamic imaging of chromatin DNA communication using Oligo-LiveFISH |
| Authors | Yanyu Zhu; Ashwin Balaji; Mengting Han; Leonid Andronov; Anish R. Roy; Zheng Wei; Crystal Chen; Leanne Miles; Sa Cai; Zhengxi Gu; Ariana Tse; Betty Chentzu Yu; Takeshi Uenaka; Xueqiu Lin; Andrew J. Spakowitz; W.E. Moerner; Lei S. Qi |
| Year | 2025 |
| DOI | 10.1016/j.cell.2025.03.032 |
| PMID | 40239646 |
| Official URL | https://doi.org/10.1016/j.cell.2025.03.032 ; full text https://pmc.ncbi.nlm.nih.gov/articles/PMC12167157/ ; bibliography https://pubmed.ncbi.nlm.nih.gov/40239646/ |
| Repository / accession | PMC12167157; NIHMS2067673 |
| Method version | V01 |
| Relevant figures | Figures 1–2, 5; other experiment captions consulted for scope |
| Relevant tables | Key Resources Table; S1–S4 via S02/S03 |
| Relevant supplementary files | S02, S03; indexed supplementary-figure legends and movie descriptions |
| Information supported | Name/expansion; predecessor distinction; guide-pool architecture; hg38 and design filters; manufacturing, delivery, labeling and imaging conditions; resource availability |
| Notes / limitations | PMC is author manuscript (available June 12, 2026), not guaranteed identical to final typeset text. Paper/code differences are retained. Full microscopy dataset stated available on request; not a blanket open-data claim. Preflight DOI placeholder is superseded by verified DOI here. |
| Access date | 2026-10-15 |
| Release / commit / checksum | PMC NIHMS2067673 snapshot fetched; bibliographic final DOI above |

## Supplementary Materials

### S02 — Published supplemental information (Tables S1, S3, S4)

| Field | Value |
|---|---|
| Source ID | S02 |
| Source type | Official publisher supplementary PDF |
| Citation | Supplemental information accompanying Zhu et al. (2025), Cell |
| Title | High-resolution dynamic imaging of chromatin DNA communication using Oligo-LiveFISH — Supplemental information |
| Authors | Yanyu Zhu; Ashwin Balaji; Mengting Han; Leonid Andronov; Anish R. Roy; Zheng Wei; Crystal Chen; Leanne Miles; Sa Cai; Zhengxi Gu; Ariana Tse; Betty Chentzu Yu; Takeshi Uenaka; Xueqiu Lin; Andrew J. Spakowitz; W.E. Moerner; Lei S. Qi |
| Year | 2025 |
| DOI | Parent: 10.1016/j.cell.2025.03.032 |
| PMID | Parent: 40239646 |
| Official URL | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425003502-mmc1.pdf |
| Repository / accession | S0092867425003502-mmc1.pdf; local preserved copy resources/Supplemental_information.pdf |
| Method version | V01 |
| Relevant figures | Associated with Figures 1–7 and STAR Methods |
| Relevant tables | Table S1(1–7), PDF pp.2–4; Table S3 p.5; Table S4 p.6; PDF p.1 is cover |
| Relevant supplementary files | PMC Table S1 cross-reference points to NIHMS2067673-supplement-11.pdf; publisher PDF downloaded instead |
| Information supported | Exact scaffold/tracrRNA/primer/adapter strings; explicit RNA 5-prime/3-prime orientation; degenerate template definitions; target regions and probe numbers; microscope comparison |
| Notes / limitations | S1(1–4) DNA template orientation interpreted from transcription/primer roles; not every DNA string has explicit terminal labels. Degenerate template N placeholders do not define a synthesized oligo. PDF text extraction cross-checked visually on sequence pages. |
| Access date | 2026-10-15 |
| Release / commit / checksum | SHA-256 db2a1cda059cbcae8d3f9f0af9211a2f1373839d0936546e326a2d102a48e4df |

### S03 — Table S2 sequence workbook

| Field | Value |
|---|---|
| Source ID | S03 |
| Source type | Official publisher supplementary XLSX |
| Citation | Table S2 accompanying Zhu et al. (2025), Cell |
| Title | Sequences and coordinates of Oligo-LiveFISH, RNA FISH probes and qRT-PCR primers |
| Authors | Yanyu Zhu; Ashwin Balaji; Mengting Han; Leonid Andronov; Anish R. Roy; Zheng Wei; Crystal Chen; Leanne Miles; Sa Cai; Zhengxi Gu; Ariana Tse; Betty Chentzu Yu; Takeshi Uenaka; Xueqiu Lin; Andrew J. Spakowitz; W.E. Moerner; Lei S. Qi |
| Year | 2025 |
| DOI | Parent: 10.1016/j.cell.2025.03.032 |
| PMID | Parent: 40239646 |
| Official URL | https://ars.els-cdn.com/content/image/1-s2.0-S0092867425003502-mmc2.xlsx |
| Repository / accession | S0092867425003502-mmc2.xlsx; local preserved copy resources/Table_S2.xlsx |
| Method version | V01 |
| Relevant figures | Figures 1–7 |
| Relevant tables | Table S2; all 14 sheets inspected; genomic sheets A:E; RNA FISH and qRT-PCR sheets A:B |
| Relevant supplementary files | PMC indexed as NIHMS2067673-supplement-12.xlsx; actual publisher XLSX inspected |
| Information supported | Exact DNA-alphabet spacer strings, genomic intervals, primary RNA FISH sequences, qRT-PCR primers; derived nonempty-row counts and coordinate-length consistency checks |
| Notes / limitations | No genomic strand column or explicit coordinate-origin declaration; eight length discrepancies in short ACAP2 sheet (see benchmarks.md). Formatting inflates max_row/max_column on one sheet; do not treat worksheet dimensions as probe counts. No source data changed or formulas recalculated. |
| Access date | 2026-10-15 |
| Release / commit / checksum | SHA-256 2279267e301cf77caf7d0210b07dabbfd57b6b81812fd2954cdc39dc85c237e8 |

## Official Protocols

The inspected executable-method evidence is the original STAR Methods (S01) and construct definitions (S02). A separate developer-authored protocol publication was not established in the scoped search. Website instructions (S04) describe tool use and direct readers to the manuscript; they are not a replacement experimental protocol.

## Official Repositories and Software

### S04 — Developer website and design-tool instructions

| Field | Value |
|---|---|
| Source ID | S04 |
| Source type | Official author-provided website and tool documentation |
| Citation | Oligo-LiveFISH website; linked by S01 Data and code availability |
| Title | Oligo-LiveFISH; Instructions; Contact |
| Authors | Developers listed on contact page: Yanyu Zhu, Stanley Qi; website/tool contacts Betty Chentzu Yu, Xueqiu Lin |
| Year | Not dated on inspected pages |
| DOI | Not applicable or not separately provided. |
| PMID | Not applicable or not separately provided. |
| Official URL | https://oligo-livefish.org/ ; https://oligo-livefish.org/instructions.html ; https://oligo-livefish.org/contact.html |
| Repository / accession | Homepage GitHub link points to S05 |
| Method version | V01 |
| Relevant figures | Not applicable or not separately provided. |
| Relevant tables | Not applicable or not separately provided. |
| Relevant supplementary files | Not applicable or not separately provided. |
| Information supported | Official design tool endorsement; advertised 18–20-bp spacer design; cell-line/gene/region selection; direct endorsement of S05 |
| Notes / limitations | Landing/documentation pages reachable; no design job submitted. Source-code availability statement on contact page has academic/nonprofit/personal-use wording; S05 separately has MIT LICENSE. Do not merge those scopes or infer service uptime from a cached search result. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Mutable website; snapshot date only |

### S05 — Official guide-design source

| Field | Value |
|---|---|
| Source ID | S05 |
| Source type | Official software repository, endorsed by developer website S04 |
| Citation | Not applicable or not separately provided. |
| Title | Oligo-LiveFISH crRNA/sgRNA design |
| Authors | LinBiologyLab repository; Xueqiu Lin and colleagues associated through S04 |
| Year | Snapshot inspected 2026; release year not assigned |
| DOI | Not applicable or not separately provided. |
| PMID | Not applicable or not separately provided. |
| Official URL | https://github.com/LinBiologyLab/Oligo-LiveFISH_crRNA_sgRNA_design |
| Repository / accession | LinBiologyLab/Oligo-LiveFISH_crRNA_sgRNA_design |
| Method version | Current public code associated with V01; not asserted to be exact publication-generating revision |
| Relevant figures | Not applicable or not separately provided. |
| Relevant tables | Not applicable or not separately provided. |
| Relevant supplementary files | README.md; LICENSE; example_data/chr3_repeatregion_lf20K.bed; code/generate_possible_guides.py; 01–05 shell stages; OffTarget_Summary_3rd.py; Remove_Original_TargetSites.py; Remove_HighPresent_sgRNA.py; remove_overlap_3rd.py; run_all.sh |
| Information supported | Design input/output, dependencies, guide orientation transformation, fixed public-code rules and departures from paper wording; source availability |
| Notes / limitations | No end-to-end execution; reference FASTA/Cas-OFFinder/environment not reconstructed. One raw-file request failed; two_file_merge.py was not inspected. README output illustration is not a full experimental probe-set release. MIT LICENSE observed in fetched file. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Commit 210b34c2c54fd0f7e282a67e8eca144d96f44ad7 |

### S06 — Paper-linked lab analysis repository context

| Field | Value |
|---|---|
| Source ID | S06 |
| Source type | Repository on the lab account linked by S01; metadata/tree/README inspection only |
| Citation | Not applicable or not separately provided. |
| Title | Oligo-LiveFISH-ML |
| Authors | QilabGitHub |
| Year | Repository describes analysis for 2026 data; do not backdate to 2025 |
| DOI | Not applicable or not separately provided. |
| PMID | Not applicable or not separately provided. |
| Official URL | https://github.com/QilabGitHub/Oligo-LiveFISH-ML |
| Repository / accession | QilabGitHub/Oligo-LiveFISH-ML |
| Method version | 2026 analysis snapshot; not established as a separate method version |
| Relevant figures | Not applicable or not separately provided. |
| Relevant tables | Not applicable or not separately provided. |
| Relevant supplementary files | Not applicable or not separately provided. |
| Information supported | Distinction between image/trajectory/model analysis and probe design; README describes segmentation, tracking, feature extraction, and modeling |
| Notes / limitations | README/tree inspected; analysis source and data not exhaustively checked or executed. It does not establish publication design code or a complete probe dataset. The original 2025 paper links the account, not this exact repository. |
| Access date | 2026-10-15 |
| Release / commit / checksum | Commit 8d84763960c4e65a1919bbc5a93c84de9c554241 |

## Official Datasets

Public probe resources are S02/S03 and the limited S05 example (D01–D03 in benchmarks.md). S01 states microscopy data are available from the lead contact upon request. This does not establish a public complete imaging release; the later S06 example tree is a separate resource with unresolved experiment mapping.

## Later Developer Publications

No distinct later developer method publication/version was established in the scoped search. A software snapshot dated 2026 is not by itself a method version.

## Secondary Sources

Search results included the Conte and Nicodemi 2025 spotlight (DOI 10.1016/j.crmeth.2025.101079); it was used only for discovery. Technical corpus claims are traced to S01–S05, not that commentary.
