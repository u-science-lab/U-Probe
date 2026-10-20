# Sources

Source IDs are method-local. Access/inspection date: 2026-10-15. V01 denotes the original 2019 RNA implementation. A link does not imply that every payload was read; see research_manifest.md. Checksums use SHA-256 unless labelled otherwise.

## Primary Publications

### S01 — Original RNA seqFISH+ paper

| Field | Value |
|---|---|
| Source ID | S01 |
| Source type | Original peer-reviewed article, full paper and Methods |
| Citation | Eng, C.-H. L. et al. Nature 568, 235–239 (2019) |
| Title | Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+ |
| Authors | Chee-Huat Linus Eng; Michael Lawson; Qian Zhu; Ruben Dries; Noushin Koulena; Yodai Takei; Jina Yun; Christopher Cronin; Christoph Karp; Guo-Cheng Yuan; Long Cai |
| Year | 2019; online 25 March |
| DOI | 10.1038/s41586-019-1049-y |
| PMID | 30911168 |
| Official URL | [Publisher](https://www.nature.com/articles/s41586-019-1049-y); [institutional deposit](https://authors.library.caltech.edu/records/1ejwh-1ka64); [author-hosted published PDF inspected](https://qianzhulab.github.io/pdf/transcriptome.pdf) |
| Repository / accession | PMCID PMC6544023 |
| Method version | V01 |
| Relevant figures | Fig. 1, especially 1a architecture/encoding; Fig. 2d sensitivity; Fig. 3 tissue scope; Extended Data Figs. 1–3 stabilization/retention/channel scope |
| Relevant tables | Supplementary Table 1, separately recorded as S02 |
| Relevant supplementary files | S02–S04 |
| Information supported | Naming, primary/readout design, manufacturing, encoding, experimental parameters, imaging, decoding, source links |
| Notes / limitations | No exact primary/readout/primer library or UCSC release identified in examined design sections. Fig. 1 is schematic without 5′/3′ labels. PDF page 6 is first Methods page; page 7 continues tissue/imaging/decoding. Publisher and PMC access was intermittent, so full local author-hosted PDF was used |
| Access date | 2026-10-15 |
| Release / commit / checksum | Published-layout 24-page PDF; evidence/paper.pdf: 47d9fe62695a7ca8b7a4d1d24e1d0b1e3a6497ff4369842ff5bcf66e3f2fcb52 |

## Supplementary Materials

### S02 — Supplementary Table 1, codebook

| Field | Value |
|---|---|
| Source ID | S02 |
| Source type | Official supplementary XLSX |
| Citation | Supplementary Table 1 to S01 |
| Title | Codebook for 10,000 genes |
| Authors | S01 author group |
| Year | 2019 |
| DOI | Parent 10.1038/s41586-019-1049-y |
| PMID | Parent 30911168 |
| Official URL | [Institutional file](https://authors.library.caltech.edu/records/1ejwh-1ka64/files/41586_2019_1049_MOESM2_ESM.xlsx) |
| Repository / accession | CaltechAUTHORS 1ejwh-1ka64; filename 41586_2019_1049_MOESM2_ESM.xlsx |
| Method version | V01 |
| Relevant figures | Parent Fig. 1 |
| Relevant tables | All sheets: Codebook-647nm, Codebook-561nm, Codebook-488nm |
| Relevant supplementary files | MOESM2 XLSX, not MOESM1 |
| Information supported | Gene names, channel assignments and four ordered base-20 symbols; derived gene counts and minimum symbol distance |
| Notes / limitations | Rows 2–3335/3334/3334 contain genes; columns C–F are symbols; no nucleotide library, target coordinates or reference release. Four codeword symbols are not an oligo string |
| Access date | 2026-10-15; downloaded and every gene/codeword row processed |
| Release / commit / checksum | 45441bc0b0412a7899c686330562e786fb06b1d1c2de6744da52ea281bb8c2a6 |

### S03 — Supplementary Tables 2–4, biological analysis workbooks

| Field | Value |
|---|---|
| Source ID | S03 |
| Source type | Official supplementary workbook bundle |
| Citation | Supplementary Tables 2–4 to S01 |
| Title | Cell-cluster genes; cortical mRNA localization; ligand–receptor pairs and neighbouring-cell enrichments |
| Authors | S01 author group |
| Year | 2019 |
| DOI | Parent 10.1038/s41586-019-1049-y |
| PMID | Parent 30911168 |
| Official URL | [Table 2 / MOESM3](https://authors.library.caltech.edu/records/1ejwh-1ka64/files/41586_2019_1049_MOESM3_ESM.xlsx); [Table 3 / MOESM4](https://authors.library.caltech.edu/records/1ejwh-1ka64/files/41586_2019_1049_MOESM4_ESM.xlsx); [Table 4 / MOESM5](https://authors.library.caltech.edu/records/1ejwh-1ka64/files/41586_2019_1049_MOESM5_ESM.xlsx) |
| Repository / accession | CaltechAUTHORS 1ejwh-1ka64 |
| Method version | V01 |
| Relevant figures | Parent Figs. 3–4 |
| Relevant tables | MOESM3: Cortex_SVZ, Olfactory bulb; MOESM4: Sheet1, Sheet2; MOESM5: Olfactory bulb, Cortex, Enrichments |
| Relevant supplementary files | MOESM3/4/5 XLSX |
| Information supported | These files are biological result tables, not complete probe lists |
| Notes / limitations | All sheet dimensions, headers and content types examined; biological findings were not reanalysed |
| Access date | 2026-10-15 |
| Release / commit / checksum | MOESM3: 72411537f8a88aeb38d0046ed79c6e9aa3506654981f018c292cfc5470c212ed; MOESM4: 9b5f3719489a3eca5f3673110b6dbd7786b756f0e4a1ffbc073d546b22d0564b; MOESM5: 6403e13edc5c15d65bcb631a8e696be5bbea165a7685eb954388f15f7dfc8f9e |

### S04 — Reporting Summary

| Field | Value |
|---|---|
| Source ID | S04 |
| Source type | Official reporting-summary PDF |
| Citation | Reporting Summary to S01 |
| Title | Reporting Summary |
| Authors | S01 author group; corresponding author Long Cai |
| Year | 2019; form last updated 29 January 2019 |
| DOI | Parent 10.1038/s41586-019-1049-y |
| PMID | Parent 30911168 |
| Official URL | [MOESM1 PDF](https://authors.library.caltech.edu/records/1ejwh-1ka64/files/41586_2019_1049_MOESM1_ESM.pdf) |
| Repository / accession | CaltechAUTHORS 1ejwh-1ka64 |
| Method version | V01 |
| Relevant figures | Not applicable |
| Relevant tables | Software and code; Data; life-sciences reporting |
| Relevant supplementary files | 41586_2019_1049_MOESM1_ESM.pdf |
| Information supported | Data/code availability context; nature of analysis tools |
| Notes / limitations | Reporting summary, not supplementary probe-sequence methods. Does not independently supply oligo strings |
| Access date | 2026-10-15 |
| Release / commit / checksum | b2774820c15e7592f755c54d1e3add85f92d3fc4c0d866ba9a49f9e50cc11b18 |

## Official Protocols

### S05 — RNA seqFISH+ Supplementary Protocol

| Field | Value |
|---|---|
| Source ID | S05 |
| Source type | Author experimental protocol, cited by S01 reference 40 |
| Citation | Eng, C.-H. L. and Cai, L. Protocol Exchange (2019) |
| Title | RNA seqFISH+ Supplementary Protocol |
| Authors | Chee Huat (Linus) Eng; Long Cai; mirrored cover lists Cai then Eng |
| Year | 2019; mirrored cover posted 26 March |
| DOI | 10.1038/protex.2019.019 |
| PMID | Not established |
| Official URL | [DOI](https://doi.org/10.1038/protex.2019.019); [legacy protocol endpoint linked by Cai Lab](https://protocolexchange.researchsquare.com/article/nprot-7333/v1) |
| Repository / accession | Legacy nprot-7333/v1; current [protocols.io V.1 destination](https://www.protocols.io/view/rna-seqfish-supplementary-protocol-8epv52ny6v1b/v1) |
| Method version | V01 supplementary preparation guidance |
| Relevant figures | None used |
| Relevant tables | Procedure headings: Tissue section preparation; Primary probes hybridization/wash; Hydrogel preparation; Digestion/Clearing; Readouts hybridization |
| Relevant supplementary files | [Eight-page protocol copy inspected](https://scispace.com/pdf/rna-seqfish-supplementary-protocol-3iwcvg3wmm.pdf) |
| Information supported | Mirror-attributed anchoring identity, sectioning/wash/gel guidance, typical times and readout count; unique mirror claims remain inferred |
| Notes / limitations | Official legacy URL failed through browser fetch and direct request redirected to a generic protocols.io homepage. Fix independently opened the current protocols.io V.1 destination on 2026-10-15; it returned no readable procedure text. Only the identified protocol copy was read; publisher equivalence/checksum not established. Mirror-only claims are inferred, not independently authenticated official-body statements. Most numerical records use independently inspected S01. The protocol is author guidance, not a separate peer-reviewed method version |
| Access date | 2026-10-15 |
| Release / commit / checksum | Copy identifies original DOI and 2019 date; no independently verified checksum |

## Official Repositories and Software

### S06 — CaiGroup/seqFISH-PLUS

| Field | Value |
|---|---|
| Source ID | S06 |
| Source type | Original-author repository and processed source-data archive |
| Citation | CaiGroup, seqFISH-PLUS; endorsed by S01 Code availability |
| Title | seqFISH+ |
| Authors | CaiGroup; README credits Michael Lawson, Sheel Shah and Nico Pierson for code work |
| Year | Repository release year not assigned |
| DOI | Not established for repository |
| PMID | Associated paper 30911168 |
| Official URL | [Repository](https://github.com/CaiGroup/seqFISH-PLUS) |
| Repository / accession | CaiGroup/seqFISH-PLUS, master |
| Method version | V01 analysis resources; repository state inspected 2026-10-15 |
| Relevant figures | Not applicable |
| Relevant tables | README prerequisites, preprocessing/postprocessing and license |
| Relevant supplementary files | README.md, recursive tree, sourcedata.zip |
| Information supported | Analysis/decoding purpose, MATLAB dependencies, processed count/location scope, absence of a separately identified oligo-design file in inspected tree |
| Notes / limitations | Tree/README and archive inventory/CSV headers inspected, not execution of every script or a complete numeric-data audit. No public design pipeline confirmed; this is not proof of universal absence. README limits software to non-commercial/academic research under Caltech rules |
| Access date | 2026-10-15 |
| Release / commit / checksum | Git tree ee6c416006fe09262526eb364c6f66edd35adc88 (tree ID, not asserted commit); README SHA-256 dc23071d7f293ea5e51f8c1ad89e421b035e11fd7ee02e07734e18219f40c075; sourcedata.zip SHA-256 64244caf920e482d5920a13dbae8593a1d2e03316dca711d91dfeffba804aaff |

### S08 — Cai Lab official resource page

| Field | Value |
|---|---|
| Source ID | S08 |
| Source type | Original-developer resource portal |
| Citation | Cai Lab, seqFISH web page |
| Title | seqFISH |
| Authors | Cai Lab, Caltech |
| Year | Not established; site footer is not a publication date |
| DOI | Not applicable |
| PMID | Not applicable |
| Official URL | [Cai Lab seqFISH resources](https://spatial.caltech.edu/seqfish/) |
| Repository / accession | Not applicable |
| Method version | Broad method family; resource links include V01 |
| Relevant figures | No technical claims extracted from website figures |
| Relevant tables | seqFISH resources: Github, seqFISH+ data, seqFISH+ protocol |
| Relevant supplementary files | Links route to repository, Zenodo and legacy protocol |
| Information supported | Official resource discovery/endorsement |
| Notes / limitations | General/promotional claims are not used as V01 quantitative evidence |
| Access date | 2026-10-15 |
| Release / commit / checksum | Mutable webpage; not pinned |

### S10 — Later DNA seqFISH+ tissue repository

| Field | Value |
|---|---|
| Source ID | S10 |
| Source type | Developer repository for related DNA extension |
| Citation | CaiGroup, dna-seqfish-plus-tissue |
| Title | dna-seqfish-plus-tissue |
| Authors | CaiGroup; README credits Yodai Takei and analysis team |
| Year | README last-updated statement: 13 September 2021 |
| DOI | Not established from inspected repository page |
| PMID | Not assigned to this repository record |
| Official URL | [Repository](https://github.com/CaiGroup/dna-seqfish-plus-tissue) |
| Repository / accession | CaiGroup/dna-seqfish-plus-tissue |
| Method version | Later DNA extension; not V01 RNA probe evidence |
| Relevant figures | Not applicable |
| Relevant tables | README and root file list |
| Relevant supplementary files | probes/ directory discovered; individual files not downloaded |
| Information supported | Public DNA-extension probe resources exist separately; avoid accidental V01 substitution |
| Notes / limitations | Resource/lineage inspection only; no claim that this repository is the same release as S09's cell experiment; no full DNA research or benchmark assessment |
| Access date | 2026-10-15 |
| Release / commit / checksum | Mutable master page; no commit pinned |

## Official Datasets

S02 records the original codebook. S06 records the repository source archive.

### S07 — NIH3T3 point-location deposit

| Field | Value |
|---|---|
| Source ID | S07 |
| Source type | Author-deposited experimental data |
| Citation | Eng, C. H. (Linus) and Cai, L. Zenodo (2019) |
| Title | NIH3T3_point_locations for RNA seqFISH+ experiments |
| Authors | Chee Huat (Linus) Eng; Long Cai |
| Year | 2019; released 6 May |
| DOI | 10.5281/zenodo.2669683 |
| PMID | Associated paper 30911168 |
| Official URL | [Zenodo record](https://zenodo.org/records/2669683) |
| Repository / accession | Zenodo 2669683, v1 |
| Method version | V01 |
| Relevant figures | S01 NIH/3T3 results |
| Relevant tables | Description and Files |
| Relevant supplementary files | DAPI_experiment1/2.zip; ROIs_Experiment1/2_NIH3T3.zip; seqFISH+_NIH3T3_point_locations.zip |
| Information supported | Decoded mRNA locations, DAPI/bead images, ImageJ segmentation and MATLAB gene-name/replicate files as described by deposit |
| Notes / limitations | Landing metadata/file descriptions only; large archives not downloaded or audited; not evidence of oligo-level probe sequences |
| Access date | 2026-10-15 |
| Release / commit / checksum | v1; point-location ZIP landing-page MD5 a228f8e6f64b0506ce072a5e6d469038 (reported, not recomputed) |

## Later Developer Publications

### S09 — DNA seqFISH+ lineage boundary

| Field | Value |
|---|---|
| Source ID | S09 |
| Source type | Later developer peer-reviewed paper |
| Citation | Takei, Y. et al. Nature 590, 344–350 (2021) |
| Title | Integrated spatial genomics reveals global architecture of single nuclei |
| Authors | Yodai Takei and colleagues, including Long Cai |
| Year | 2021 |
| DOI | 10.1038/s41586-020-03126-2 |
| PMID | 33505024 |
| Official URL | [PubMed abstract and figure captions](https://pubmed.ncbi.nlm.nih.gov/33505024/); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7878433/) |
| Repository / accession | PMCID PMC7878433 |
| Method version | Related DNA seqFISH+ extension; no V01 replacement |
| Relevant figures | Fig. 1 and Extended Data Fig. 1 captions |
| Relevant tables | None extracted |
| Relevant supplementary files | No sequence files extracted in this scope |
| Information supported | DNA targets and distinct probe/encoding procedures justify separation from original RNA method |
| Notes / limitations | Abstract and caption-level lineage check only; not a complete DNA method corpus |
| Access date | 2026-10-15 |
| Release / commit / checksum | Published record; no local PDF checksum |

## Secondary Sources

No review or informal summary supplies V01 technical parameters. The following is developer-authored documentation of a third-party implementation, included to prevent software misattribution rather than to establish original experimental rules.

### S11 — Third-party SeqFISH+ Probe Designer documentation

| Field | Value |
|---|---|
| Source ID | S11 |
| Source type | Primary software documentation for a third-party pipeline |
| Citation | Oligo Designer Toolsuite, SeqFISH+ Probe Designer documentation |
| Title | SeqFISH+ Probe Designer |
| Authors | Oligo Designer Toolsuite developers; documentation references Mekki et al. and Kuemmerle et al. |
| Year | Mutable documentation; page year not assigned |
| DOI | Referenced software/paper DOIs are not a DOI for this documentation page |
| PMID | Not established |
| Official URL | [Project documentation](https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/seqfishplus_probe_designer.html) (official for that software, not for original seqFISH+) |
| Repository / accession | Oligo Designer Toolsuite |
| Method version | Contemporary third-party implementation; not original V01 execution evidence |
| Relevant figures | Pipeline schematic, not used for V01 component orientation |
| Relevant tables | API example; Pipeline Description; references |
| Relevant supplementary files | No original experimental probe library extracted |
| Information supported | A third-party seqFISH+ design pipeline exists; its configurable defaults cannot be attributed automatically to Eng et al. |
| Notes / limitations | Not executed; package version/dependency state not pinned; biological design rules above come from S01 instead |
| Access date | 2026-10-15 |
| Release / commit / checksum | latest documentation branch; no release assigned |
