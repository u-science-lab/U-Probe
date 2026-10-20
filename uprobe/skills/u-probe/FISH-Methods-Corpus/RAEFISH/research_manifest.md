# Research Manifest

Audit only. Initial Research stage, 2026-10-15. No independent review or validation claim. Repository preflight report Appendix E states its blocking template/schema errors were closed. Current schema, terminology, routing index and all Research templates were read. Existing unrelated edits were preserved. No remote U-Probe service access was required.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Primary article | PMC12662569; DOI in sources.md | Full HTML/text; Summary, Results, Figure 1/7 legends, Figure S1 legend, STAR Methods, Key resources, supplement links | Naming, 5′→3′ order, lengths, reference releases, code scheme, protocols, readout, data links | Figure pixel geometry not independently used; 94/96 decoding wording inconsistent | Exact filtering algorithm; splint residue count | 2026-10-15 | PMC retrieved HTML, checksums.json |
| S02 | Official supplement | Publisher mmc1.xlsx | All nine worksheet headers/sample rows; human/mouse complete sequence columns for lengths/counts; control-name counts; Perturb sheets; amplification/readout oligos | Exact ordered oligos, codebook columns, synthesis primers, variant lengths, counts | No alignment to reference FASTA; no exact reconstruction of final synthesized working libraries | Numeric gene IDs, target accessions/coordinates, bit mapping | 2026-10-15 | Original bytes retained/hash |
| S03 | Author page/archive | Yale RAEFISH page; author-linked ZIP | Page, 92-entry central directory, four GoodGeneList MAT files and mouse runM2_2_HPC.m | Transcript sequences/numeric identifiers/C-L binding pairs; scoped 94-dimensional decoder | Full ZIP not downloaded; initial approval timeout resolved on retry | Versioned accession reconstruction; complete dependencies; design-selection code | 2026-10-15 | Page snapshot; Last-Modified 2025-05-17 |
| S04 | Author dataset | Mendeley v1 page and public API | v1 description; v2 JSON file metadata | Release distinction; 23 file entries; codebook names and processed outputs | Binary MAT files not opened; metadata hashes not local verification | v1/v2/ZIP equivalence; target mapping within MAT codebooks | 2026-10-15 | API version 2, local raw JSON |
| S05 | GEO deposit | GSE302564 | Series summary/resource listing | Bulk RNA-seq comparator accession | Raw reads not downloaded | No probe-specific reference inferred | 2026-10-15 | Local HTML |
| S06 | Author preprint | bioRxiv/PMC11952344/PubMed indexed excerpts | Abstract/results and supplementary-listing search excerpts | Earlier title/stage; final-source disambiguation | Full preprint and supplements not downloaded | Detailed stage-to-stage changes not claimed | 2026-10-15 | bioRxiv v1 indexed excerpts |

### Discovery and access record

Queries included `RAEFISH RNA FISH`, `RAEFISH probe design github`, `RAEFISH Table S1`, and `RAEFISH protocol GENCODE GC`. Followed the final article's supplementary and data/code links, laboratory resource page, Mendeley dataset API and GEO record. Searches for repository/software/protocol resources established the author-linked analysis archive but did not establish a probe-design program. The laboratory menu's ProbeDealer link was not taken as proof of method-specific software.

PMC browser access sometimes returned reCAPTCHA; direct permitted HTTP download retrieved the complete article. The first PMC spreadsheet URL returned an HTML challenge rather than a workbook. Publisher `mmc1.xlsx` succeeded and was identified from actual worksheet contents, not inferred from file numbering. Publisher mmc5 is a differential-expression table, not Table S1. Only the verified Table S1 is used for sequence claims. Mendeley `/versions/1` endpoint returned 404; the dataset endpoint returned v2, explicitly retained as a distinct snapshot. No confidential resources or credentials were read.

### Research checklist and self-check scope

- Official method spelling and original peer-reviewed publication checked; V01/V02 separated.
- Primary, supplementary, author protocol, repository/software and dataset discovery performed; access limits recorded.
- V01 architecture/orientation checked against explicit STAR Methods wording; template vs final encoding orientation separated.
- Components, original terms, sequence classes, exact examples, chemical modifications and unresolved V02 assembly documented.
- Target-binding, barcode/readout, design rules, reference releases and public probe sets investigated.
- D01–D04 assessed separately; no Grade A/B or strict ground-truth claim.
- Major claims carry source IDs/locators; missing values remain explicit/null.
- Required sections, YAML structure, routing status, links, counts and exact copied sequence examples checked using self_check.py.
- This checklist is a Research self-check, not an independent Review or Fix stage. validation_report.md and fix_report.md were not created.

## Unresolved Research Questions

1. Which stable gene/transcript accessions correspond to S02's numeric Ensembl_GeneID fields, and which isoforms/sites were selected?
2. What reference FASTA/GTF file versions/checksums, GC/Tm bounds, repeat filters, off-target database and selection algorithm generated the libraries?
3. Is there official design code separate from the author-linked analysis archive?
4. What is the complete verified mapping from each bit to round, adapter and channel, and what synthesis transformations define the final V02 encoding strand? The archive decoder verifies 94 dimensions for one mouse example only.
5. How should splint “16-nt” wording be reconciled with 15 printed standard bases plus inverted dT and the 12+3 binding/overhang description?
6. Does endogenous image-analysis code implement 94 or 96 dimensions, and what is the intended meaning of the source's “94 rounds” phrase?
7. How do Mendeley v1, v2 and the 2025-05-17 Yale archive differ? Mendeley codebook binary contents remain uninspected. Four Yale design-associated MAT files and one decoder were inspected; The S03 author archive contains four inspected GoodGeneList MAT files under Design240427/Design240516. They contain numeric ENSG/ENST/LongestENST fields, transcript Sequence strings and C/L binding sequences. Comparison by symbol plus both binding sequences matches 23,303 distinct published human sequence-pair records (derived: evidence/design_resource_inspection.json); this does not establish versioned stable transcript accessions, a reference build or a probe-selection algorithm. The human workbook has 23,312 non-control rows but 23,300 unique symbols; mouse has 21,955 rows but 21,954 unique symbols. Duplicate symbols must not be collapsed during extraction. The sampled mouse decoder RAEFISH/runM2_2_HPC.m sets RoundsOfHybs=94 and reshapes vectors with that size. This supports a 94-dimensional implementation for this archive example (inferred from code), while the publication's 96-dimensional wording remains an explicit source discrepancy.
