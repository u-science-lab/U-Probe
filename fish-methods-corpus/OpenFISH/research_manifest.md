# Research Manifest

Audit only; exclude from ordinary retrieval. Research draft completed on 2026-10-15; no independent review or validation performed. This is a bounded initial draft with substantial experimental-source access gaps, not a complete extraction of the paper and supplements.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Peer-reviewed paper, partial access | Publisher indexed page and PubMed 42810347 | Title, bibliographic record, abstract; indexed Introduction and Data and code availability excerpts | Name, publication identity/date, application and deposits | Full publisher/Cell article and supplements could not be retrieved; figures/methods not read in full | Complete wet protocol, probe tables, experimental chemistry | 2026-10-15 | DOI 10.1016/j.neuron.2026.09.007; PII S0896627326006902 |
| S02 | Preprint discovery | bioRxiv 10.1101/2025.08.19.671030; bibliographic discovery through Sciety/author ORCID | Bibliographic record only | Earlier title/source lineage | Full HTML/PDF failed; no technical claims extracted | Changes from preprint to final | 2026-10-15 | v1 URL, uninspected full text |
| S03 | Author repository | evidence/paper-code | README; repository tree; ProbeDesign source and sequence files compared to S04 | Companion provenance; source equivalence | Reproduce figure notebooks inventoried, not fully analyzed; no code executed | Exact publication experiment settings | 2026-10-15 | 01ea661ef43e1c6c7bc19a3a6e3c214a420a262c |
| S04 | Author computational workflow | evidence/protocol | README; notebook 01–06 cell sources; openDesign utils/generator/generator_multi/filter/cmdtools; openDecode decoding source excerpt; env files; Bridge_sequences.csv, UsingPad.fa, Example_ProbeDesign.csv | Probe organization, strand selection, rules, codebook and software; sequence counts and duplication | No wet protocol; no package installed; no design or decoding run; blacklist biological provenance not independently verified | Synthesized modifications, initiation, complete target outputs, reference freezing | 2026-10-15 | be490eb19503c742230665da1e665c6d6c3d97c8 |
| S05 | Author demo deposit | evidence/zenodo.json; evidence/demo_CodingBook.CSV | API metadata, full file inventory, exact codebook | 20 gene/round/color assignments, demo availability | Large images/MSI/Zarr/H5AD not downloaded | Full experimental and probe metadata | 2026-10-15 | API revision 5; codebook MD5 verified |
| S06 | Author study deposit | evidence/zenodo_17659890.json | API metadata and file list | Main raw/processed archive provenance | Archive contents not inspected; future publication_date conflict | Probe tables in archives; chronology | 2026-10-15 | Record 17659890; version index 2 |
| S07 | Author comparison deposit | evidence/zenodo_20917671.json | API metadata and file list | Adjacent mouse-brain sections, two replicates, RIN 8.7/8.4 | Four Zarr archives not downloaded | Panel/probe completeness and reference context | 2026-10-15 | Record 20917671 |

## Unresolved Research Questions

1. Obtain the final paper's complete methods, supplementary inventory and sequence tables. Establish concentrations, times, temperatures, buffers, fixation/clearing and round-exchange procedures from those sources.
2. Verify chemically complete synthesized strands: padlock phosphorylation, readout labels, RCA initiation, blocking/anchoring chemistry and ligation junctions. Do not infer them from bare FASTA strings.
3. Locate final biological A/B probe sequences, target accession versions, genomic/transcript coordinates and manually retained pair indices. Clarify per-gene probe counts and reference releases.
4. Compare the 2025 preprint with the 2026 publication before defining separate experimental versions.
5. Inspect relevant bulk archives before assessing D03/D04 as probe-design benchmarks. Present grades are scoped to examined resources.
6. Resolve S06's API publication_date of 2027-08-19 relative to retrieval on 2026-10-15; do not silently revise the source metadata.
7. Resolve software ion-unit docstring discrepancy and notebook/constructor parameter differences before reproducing a design. Record random seeds, database snapshot and manual selection.

## Discovery and Access Log

- Broad searches for OpenFISH, Neuron, GitHub, padlocks, hybridization and supplementary tables found the 2026 article, earlier preprint, author repositories and linked datasets. Unrelated nanopore-basecalling, robotics, NFT and assistant projects sharing the name were excluded.
- Publisher/DOI/Cell and bioRxiv full-text attempts returned tool errors or HTTP 403. The correct publisher PII is S0896627326006902. An early guessed PII ending 00679-3 was not evidence and was discarded after bibliographic verification.
- Standard publisher supplementary endpoint patterns for mmc1.pdf, mmc2.xlsx and mmc3.xlsx were tried as discovery attempts; all returned HTTP 403. These names are **unconfirmed**, not asserted supplementary inventory. No table number was invented from these requests.
- Author Git repositories were cloned read-only at the listed commits. Local source snapshots and SHA-256 manifests are retained. openDesign utils/filter/generator/generator_multi/cmdtools and Bridge_sequences.csv/UsingPad.fa are byte-identical across the two snapshots.
- Source notebooks were parsed into cell-indexed text for inspection. No notebook cells or downloaded project code were executed. Only local extraction/checking scripts ran.
- Zenodo API metadata and the 417-byte demo codebook were retrieved. Large data archives were not needed for the scoped sequence/codebook assessment and remain unexamined.
- Existing unrelated edits in other methods and methods_index.yaml were preserved. Only the OpenFISH routing lifecycle entry was updated.

## Research Checklist

| Required check | Research-stage result |
|---|---|
| Official name / original publication | Verified against S01 bibliographic and abstract records; no unsupported expansion |
| Primary and supplementary discovery | Conducted; final full methods/supplements inaccessible; limits above |
| Version separation | V01 original 2026; preprint and software snapshots retained without invented chemistry versions |
| Architecture / orientation | Code-level A/B and padlock precursor organization recorded; biological target orientation checked through reverse-complement operations; experimental synthesis orientation remains bounded |
| Components / original terminology | Individual components and source terms recorded with inferred functional mappings |
| Exact vs schematic | Schematics, exact library/notebook instances and doubled alignment references distinguished |
| Target design rules | Code and notebook settings traced; divergent defaults retained |
| Barcode/readout | Round-plus-color-pair identities distinguished from RO nucleotide segments |
| Software | Two author repositories inspected at pinned commits; not executed |
| Public resources | Partial sequence libraries, exact demo codebook and data inventories recorded separately |
| Species / references | Mouse example identified; no genome/annotation release guessed from tRNA filename |
| Benchmarks | D01 C; D02 D for probe reconstruction; D03/D04 unassessed; no strict ground-truth claim |
| Source traceability / gaps | S01–S07, locators and access limits recorded |
| Internal consistency | Structural checks and arithmetic/file checks below; independent Review still required |

## Research Self-Checks

Structural checks cover required Research file existence, H1/H2 anchors, required table columns/benchmark fields, YAML parsing with duplicate-key rejection, metadata template keys, legal lifecycle/evidence/class values, local routing/link targets and source IDs. Sequence inspection confirms 110 bridges of 30 nt, 110 144-nt padlock references with identical 72-nt halves, and codebook MD5 2387da50df31bcc6e346eb8f2724e76e. These checks support the draft's internal consistency only; they do not confer biological validation.
