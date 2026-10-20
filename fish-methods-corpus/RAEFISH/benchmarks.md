# Benchmark Resources

## Benchmark Summary

Research-stage assessment, 2026-10-15. D01–D03 are **C**, useful for architecture/codebook comparison but not strict probe-design ground truth. D04 is **C** after inspection of design-associated MAT files; the processed expression resources alone are insufficient for probe reconstruction. These are scoped evidence assessments, not experimental validation.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published sequence list is not automatically Grade A. Significant reference, target-mapping and selection-rule gaps prevent assigning B here.

## Candidate Datasets

### D01 — Human endogenous RAEFISH templates/codebook

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Human endogenous RAEFISH templates/codebook |
| Method version | V01 |
| Species | Homo sapiens |
| Cell / tissue / sample | A549-Cas9 cultured cells |
| Genome build | Not established for probe design. Do not import the GRCh38.p14 bulk RNA-seq comparator reference by assumption. |
| Transcriptome / annotation version | GENCODE v46 (S01 Probe design) |
| Target genes | 23,312 targets reported; 23,352 data rows minus 40 NegativeControl names = 23,312 (derived, S02 B3:B23354). |
| Target regions | Adjacent target RNA sites; transcript accession and isoform mappings unresolved. |
| Genomic/transcript coordinates | No assembly-linked positions, strand mapping, coordinate origin or interval convention established from inspected tables. S03 provides transcript sequence strings, not independently verified reference coordinates. |
| Target-binding sequences | 30-nt padlock-binding and priming-binding regions embedded in templates. S03 MAT C/L fields match all distinct published human symbol/sequence-pair records; versioned reference mapping unresolved. |
| Full probe sequences | 80-nt padlock, 107-nt priming and 150-nt encoding ordered templates. The latter two are not identical to working probes. |
| Barcode assignments | 94-bit MHD4_Code, S02 C3:C23354. |
| Readout sequences | S02 94-bit RAEFISH readout oligos; two common dyes, 94 adapters and blockers. |
| Primer/linker sequences | S02 Library amplification oligos; no unsupported linker asserted. |
| Design software | No official probe-design implementation established. Author-linked analysis scripts are not an equivalent. |
| Design parameters | Two adjacent 30-nt recognition regions; 0/1-nt gap; one set/target; 40 high-expression targets excluded. |
| Source | S01 Probe design; S02 Human RAEFISH probe templates A:F and shared readout/amplification sheets. |
| Availability | Official workbook downloaded and read on 2026-10-15; original bytes preserved with SHA-256 in evidence/checksums.json. |
| Missing information | Stable transcript IDs, site coordinates, exact design-reference files, off-target/GC/Tm/isoform rules and executable selection algorithm; final synthesis orientation must be reconciled. |
| Benchmark grade | C |
| Grade justification | Exact ordered templates and code assignments support architecture/readout checks. Substantial target/reference/algorithm gaps prevent a reproducible selection benchmark; these are not minor omissions suitable for B. |
| Comparison scope | Probe-component layout, template lengths, coding assignments and oligo resource extraction after explicitly resolving synthesis orientation. Not exact target-site-selection ground truth. |
| Assessment date | 2026-10-15 |
| Evidence status | reported for source contents; derived for counts/lengths; inferred for suitability grade. |

### D02 — Mouse endogenous RAEFISH templates/codebook

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | Mouse endogenous RAEFISH templates/codebook |
| Method version | V01 |
| Species | Mus musculus |
| Cell / tissue / sample | C57BL/6 liver, placenta and lymph node sections |
| Genome build | Not established for probe design. Do not import the GRCh38.p14 bulk RNA-seq comparator reference by assumption. |
| Transcriptome / annotation version | GENCODE vM35 (S01 Probe design) |
| Target genes | 21,955 targets reported; 21,981 data rows minus 26 NegativeControl names = 21,955 (derived, S02 B3:B21983). |
| Target regions | Adjacent target RNA sites; transcript accession and isoform mappings unresolved. |
| Genomic/transcript coordinates | No assembly-linked positions, strand mapping, coordinate origin or interval convention established from inspected tables. S03 provides transcript sequence strings, not independently verified reference coordinates. |
| Target-binding sequences | Two padlock/priming pairs per row; some pairs are duplicated rather than independently selected. |
| Full probe sequences | Two 80-nt padlocks, two 107-nt priming templates, two 150-nt encoding templates per row. |
| Barcode assignments | 94-bit MHD4_Code, S02 C3:C21983. |
| Readout sequences | Shared S02 94-bit RAEFISH readout oligos. |
| Primer/linker sequences | Shared S02 Library amplification oligos. |
| Design software | No official probe-design implementation established. Author-linked analysis scripts are not an equivalent. |
| Design parameters | Adjacent 30-nt sites; no overlap; 0/1-nt gap; most two sets/target, duplicate fallback; 26 high-expression targets excluded. |
| Source | S01 Probe design; S02 Mouse RAEFISH probe templates A:I and shared readout/amplification sheets. |
| Availability | Official workbook downloaded and read on 2026-10-15; original bytes preserved with SHA-256 in evidence/checksums.json. |
| Missing information | Stable transcript IDs, site coordinates, exact design-reference files, off-target/GC/Tm/isoform rules and executable selection algorithm; final synthesis orientation must be reconciled. |
| Benchmark grade | C |
| Grade justification | Exact ordered templates and code assignments support architecture/readout checks. Substantial target/reference/algorithm gaps prevent a reproducible selection benchmark; these are not minor omissions suitable for B. |
| Comparison scope | Probe-component layout, template lengths, coding assignments and oligo resource extraction after explicitly resolving synthesis orientation. Not exact target-site-selection ground truth. |
| Assessment date | 2026-10-15 |
| Evidence status | reported for source contents; derived for counts/lengths; inferred for suitability grade. |

### D03 — Perturb-RAEFISH gRNA templates/codebook

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | Perturb-RAEFISH gRNA templates/codebook |
| Method version | V02 |
| Species | Engineered gRNA library expressed in Homo sapiens cells |
| Cell / tissue / sample | A549-Cas9i CRISPR screen |
| Genome build | Synthetic readout targets; genome/assembly for original CRISPR guide selection not established. |
| Transcriptome / annotation version | Not applicable to synthetic gRNA scaffold; genome reference for gRNA targeting design unresolved. |
| Target genes | 574 perturbation genes with three gRNAs each reported (S01 CRISPR screen library construction); template table has 1,758 rows including controls (derived), not 1,758 endogenous genes. |
| Target regions | gRNA spacer/upstream constant RNA and downstream scaffold; do not equate these with the knockout gene mRNA. |
| Genomic/transcript coordinates | No assembly-linked positions, strand mapping, coordinate origin or interval convention established from inspected tables. S03 provides transcript sequence strings, not independently verified reference coordinates. |
| Target-binding sequences | 40-nt recognition in padlocks/encoding templates; shared 45-nt scaffold recognition. |
| Full probe sequences | 100-nt padlock and 160-nt encoding ordered templates; 55-nt shared priming oligo (derived counts); final working encoding strand unresolved. |
| Barcode assignments | 14-bit codebook in S02 E3:E1760, not V01 HD4. |
| Readout sequences | 14 directly labeled oligos in Perturb-RAEFISH readout oligos A2:B15. |
| Primer/linker sequences | Shared primer in Perturb-RAEFISH probe templates D3; library manufacturing oligos in final sheet. |
| Design software | No official probe-design implementation established. Author-linked analysis scripts are not an equivalent. |
| Design parameters | 20-nt spacer plus 20-nt U6-derived constant region; 45-nt shared scaffold recognition; weight-4 HD2; three guides/gene share code. |
| Source | S01 Probe design/CRISPR screen library construction; S02 Perturb-RAEFISH probe templates, readout oligos and gRNA template oligo library. |
| Availability | Official workbook downloaded and read on 2026-10-15; original bytes preserved with SHA-256 in evidence/checksums.json. |
| Missing information | Stable transcript IDs, site coordinates, exact design-reference files, off-target/GC/Tm/isoform rules and executable selection algorithm; final synthesis orientation must be reconciled. |
| Benchmark grade | C |
| Grade justification | Exact ordered templates and code assignments support architecture/readout checks. Substantial target/reference/algorithm gaps prevent a reproducible selection benchmark; these are not minor omissions suitable for B. |
| Comparison scope | Probe-component layout, template lengths, coding assignments and oligo resource extraction after explicitly resolving synthesis orientation. Not exact target-site-selection ground truth. |
| Assessment date | 2026-10-15 |
| Evidence status | reported for source contents; derived for counts/lengths; inferred for suitability grade. |

### D04 — Author-linked processed imaging/expression resources and analysis archive

| Field | Value |
|---|---|
| Dataset ID | D04 |
| Dataset name | Author-linked processed imaging/expression resources and analysis archive |
| Method version | V01/V02; ZIP last modified 2025-05-17; Mendeley API snapshot version 2. Snapshot equivalence not assumed. |
| Species | Human and mouse as described by S01/S04 |
| Cell / tissue / sample | A549 and mouse tissues; individual file mapping requires directory context. |
| Genome build | Not established for probe-design reconstruction. |
| Transcriptome / annotation version | Not established in inspected resource metadata. |
| Target genes | Four S03 human MAT files expose ShortName/ENSG/ENST fields and C/L binding sequences; cross-match to S02 is documented. Other codebook binaries remain uninspected. |
| Target regions | Not established. |
| Genomic/transcript coordinates | Not established; no coordinate convention assigned. |
| Target-binding sequences | Four S03 MAT files contain exact C/L binding strings and transcript Sequence strings; derived cross-matches documented in evidence/design_resource_inspection.json. |
| Full probe sequences | Not established in these inspected resource metadata; use S02 for ordered oligos. |
| Barcode assignments | API lists Codebook_StevenRaeFISH23K_LabA549.mat and Codebook_RaeFISH.mat; binary contents not inspected. |
| Readout sequences | Not established here; see S02. |
| Primer/linker sequences | Not established here; see S02. |
| Design software | ZIP directory lists decoding/single-cell analysis scripts, not an established probe-design workflow. |
| Design parameters | S03 MAT structures include ProbeCount/Tx/C/L and transcript sequences; numerical fields do not establish an executable selection algorithm. |
| Source | S03 Yale page/ZIP directory, four GoodGeneList MAT files and sampled decoder; S04 Mendeley API dataset 8kbv637pxh v2; original paper links v1. |
| Availability | Yale page and ZIP range access verified; API lists 23 files in v2. Full datasets/ZIP not downloaded; four human design-associated MAT resources and one decoder retrieved and inspected by byte range. |
| Missing information | Mendeley codebook binary inspection, release correspondence, versioned transcript reference files and design-generation algorithm. |
| Benchmark grade | C |
| Grade justification | Four human design-associated MAT files provide transcript sequences and C/L target-binding strings, with correspondence to published template records. Versioned accession/reference/selection metadata remain incomplete, so these resources support architecture and sequence comparison, not strict selection ground truth. The S03 author archive contains four inspected GoodGeneList MAT files under Design240427/Design240516. They contain numeric ENSG/ENST/LongestENST fields, transcript Sequence strings and C/L binding sequences. Comparison by symbol plus both binding sequences matches 23,303 distinct published human sequence-pair records (derived: evidence/design_resource_inspection.json); this does not establish versioned stable transcript accessions, a reference build or a probe-selection algorithm. The human workbook has 23,312 non-control rows but 23,300 unique symbols; mouse has 21,955 rows but 21,954 unique symbols. Duplicate symbols must not be collapsed during extraction. The sampled mouse decoder RAEFISH/runM2_2_HPC.m sets RoundsOfHybs=94 and reshapes vectors with that size. This supports a 94-dimensional implementation for this archive example (inferred from code), while the publication's 96-dimensional wording remains an explicit source discrepancy. |
| Comparison scope | Resource discovery and future image/decoding benchmark investigation; not strict probe-design ground truth. |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource descriptions; inferred resource-role classification and grade. |

GSE302564 (S05) contains bulk RNA-seq comparator data, not an official probe-design truth set; it is not separately graded as a probe-set candidate.
