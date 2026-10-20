# Research Manifest

Audit only; exclude from ordinary biological retrieval. Research Agent, 2026-10-15. This is an initial evidence-grounded draft, not an independent review or a validated corpus.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Original publication | DOI, PubMed/Caltech record; author-hosted PDF retained as evidence/paper.pdf | Main text, Fig. 1a visually, Fig. 2d, tissue results, Methods pp. 6–9, Extended Data Fig. 1 visually and Fig. 2/3 captions | Official name, sample scope, design filters, probe lengths/manufacture, encoding, experimental parameters and decoder | Main publisher/PMC access intermittent; author-hosted published-layout PDF used | Full oligo strings/order and exact references absent from inspected design sections | 2026-10-15 | PDF SHA-256 in sources.md and evidence/inspection.json |
| S02 | Official supplement | Caltech deposit; evidence/MOESM2 workbook | All three sheets, all gene/codeword rows; dimensions and columns; unique projections computed | 10,000 assignments; channel counts; symbols 1–20; minimum symbol distance 2 derived | Codebook only, no molecular sequences | Target/oligo/readout sequence map | 2026-10-15 | 45441bc0b0412a7899c686330562e786fb06b1d1c2de6744da52ea281bb8c2a6 |
| S03 | Official supplements | Caltech MOESM3/4/5 XLSX | Every sheet enumerated; headers, dimensions and row content types inspected | Cluster genes, localization and ligand/receptor results; not probe tables | Biological results not independently reproduced | No probe library in these table types | 2026-10-15 | All three SHA-256 values in sources.md/evidence/inspection.json |
| S04 | Official reporting summary | Caltech MOESM1 PDF; evidence/reporting_summary.txt | Two-page extracted summary, software/code, data and life-science reporting fields | Supporting data/code context; sample-reporting scope | Not a sequence supplement | No additional sequence resource identified | 2026-10-15 | b2774820c15e7592f755c54d1e3add85f92d3fc4c0d866ba9a49f9e50cc11b18 |
| S05 | Author protocol | Original DOI and Cai Lab link; SciSpace copy through web PDF extraction | Eight-page copy, Procedure and Troubleshooting; cover identity; sample, primary hybridization, wash, hydrogel, digestion, restabilization, readout/imaging sections | Anchoring reagents, broader sectioning/time guidance, 240-readout statement | Official endpoint inaccessible/redirected; no publisher-to-copy checksum equivalence; paper used as principal numeric record | Recover an official archived protocol payload and verify copy identity | 2026-10-15 | Copy cover: DOI 10.1038/protex.2019.019, 26 March 2019 |
| S06 | Official repository/data | Github master README, recursive API tree and sourcedata.zip | README complete; tree file inventory; all six top-level archive members; five CSV headers; nested 913 location-CSV filenames | MATLAB analysis purpose, dependencies/license; counts/centroids/locations, no identified oligo design library | Scripts not executed; numeric matrices and every nested CSV not audited | Any original design scripts or full oligo files outside public tree | 2026-10-15 | Tree ee6c416006fe09262526eb364c6f66edd35adc88; archive checksum in sources.md |
| S07 | Official author dataset | Zenodo 2669683 | Landing page Description, Files, version and reported MD5 | Two NIH/3T3 replicate point-location deposits and segmentation/DAPI data | Payloads not downloaded; landing metadata only | Payload-level confirmation if used in downstream data benchmark | 2026-10-15 | v1, 6 May 2019 |
| S08 | Official lab portal | spatial.caltech.edu/seqfish | seqFISH resources and linked targets | Official repository/data/protocol discovery | General promotional statements not used as numerical evidence | Legacy protocol destination currently unavailable | 2026-10-15 | Mutable page |
| S09 | Later developer paper | PubMed 33505024; PMC search excerpt | Abstract and Fig. 1/Extended Data Fig. 1 captions | DNA extension differs in biological target and probe/encoding scheme | Lineage inspection, not full DNA research | Separate extension corpus needed before using DNA sequences/design parameters | 2026-10-15 | Published 2021 DOI |
| S10 | Later developer repository | CaiGroup/dna-seqfish-plus-tissue | README and root probes/ listing | Distinct DNA probe resource, not original RNA probe library | Individual probes not inspected; not assumed same release as S09 | Independent DNA resource audit if added to corpus | 2026-10-15 | Master webpage; README says updated 13 September 2021 |
| S11 | Third-party software documentation | Oligo Designer Toolsuite latest SeqFISH+ pipeline page | API example, Pipeline Description, references | Separate design implementation, inputs/outputs and third-party identity | Not executed; defaults excluded from original V01 rules | Version pin/dependency tests if comparing this implementation | 2026-10-15 | Mutable latest documentation |

## Unresolved Research Questions

1. Where are the complete original target-binding and synthesized primary probe strings, and the 240 readout sequences with their codeword-slot mapping? S02 is not such a library; S03 contains downstream biology; S06's inspected inventory contains processed measurements.
2. What are the exact 5′→3′ primary oligo order, spacer/residual bases, manufacturing primer strings and array-pool construct? Fig. 1 lacks strand labels; approximately 93 nt is not sufficient to reconstruct them.
3. Which mouse assembly, UCSC annotation release, transcriptome FASTA and ENCODE experiments were used? BLAST options, candidate stepping, overlap/spacing and isoform-selection implementation remain unresolved.
4. Is original-author probe-design code public outside the linked repository? This search did not establish it; metadata remains null. Third-party tools are not substitutes for provenance.
5. Can an official archived S05 payload be retrieved and matched to the inspected copy? The legacy link was verified from Cai Lab, but its current landing response was a generic protocols.io page.
6. Clarify printed poly-D-lysine units and readout-coupling solvent nomenclature before operational reuse; source wording is preserved. Resolve paper optical-plane versus Methods z-acquisition details if reproducing imaging exactly.
7. Later DNA seqFISH+ resources were identified for lineage separation, not fully researched or graded. They must receive their own version-scoped research before sequence-level retrieval.

## Discovery and Access Record

Fix-stage access update (2026-10-15): the current [official protocols.io V.1 destination](https://www.protocols.io/view/rna-seqfish-supplementary-protocol-8epv52ny6v1b/v1) was opened directly and identified by title, but returned no readable procedure. The eight-page SciSpace mirror was independently reopened; its cover and Procedure were checked for the review findings. Neither route discovery nor matching identifiers establishes official-body equivalence. The historical legacy-endpoint failures below remain valid. S05 authentication remains unresolved (R002); mirror-only claims are now qualified as inferred in production. Fix also inspected S01 Fig. 1a visually and the retained published-layout Methods (probe manufacture, readout synthesis, cell/tissue experiments and imaging), with live PMC full text available for source identity/context. See fix_report.md for corrections and final checks; the Research-stage statements below document the original stage, not the current lifecycle.

Repository specifications and all Research templates were read before authoring. Existing ExSeq edits and its index status were preserved. The preflight report records zero remaining blocking ERRORs after its fixes; no global schema changes were made. seqFISH+ directory spelling matches S01. No server access was required.

Searches covered the original title, supplementary resources, protocol DOI, CaiGroup repositories, probe sequences and probe-design tools. The Caltech record exposes MOESM1 (reporting summary), MOESM2 (codebook), MOESM3–5 (biological tables) and Extended Data image files. The 24-page author-hosted PDF includes those Extended Data figures. A PMC SI Guide was discovered by its listing but not inspected; no probe sequence claim relies on it. GSE98674 is comparison RNA-seq cited by S01, not a new seqFISH+ probe-library accession.

Initial network requests through the default proxy failed. Direct network requests succeeded for the paper, supplements and GitHub resources; evidence/retrieval_log.json records their status/checksums. Failed guessed nprot-7019 URLs were discarded; the actual nprot-7333 link came from S08. The attempted protocol PDF endpoint returned 403. The direct legacy protocol page returned HTTP 200 but its HTML canonical/title identified a generic protocols.io homepage; it was not counted as protocol-content inspection. The SciSpace copy was then read through web PDF extraction. No credentials or signed redirect query strings are stored as corpus sources.

The evidence/ directory is a construction-support archive, not normal retrieval content. It holds downloaded originals, checksums, inspection output and small reproducible inspection scripts. File hashes establish the downloaded bytes, not independent scientific validation. Full source-data download was used only for inventory/header inspection; Zenodo large payloads were not fetched. No sequences were reconstructed from diagrams or third-party defaults.

## Research Checklist

| Research requirement | Result |
|---|---|
| Official name and original publication | Checked against S01 title/abstract and PubMed; seqFISH+ preserved |
| Primary, supplementary and official resource discovery | S01–S08 examined within recorded boundaries; all four table workbooks inspected |
| Versions separated | Original RNA V01 only in quantitative metadata; DNA extension and smFISH comparison explicitly distinguished |
| Probe architecture and orientation | Every relevant molecular species identified; explicit unresolved 5′→3′ primary order retained; readout 5′ modification verified |
| Component classification and terminology | Original terms retained alongside corpus functional labels; sequence classes recorded |
| Exact versus schematic strings | Schematics labelled; no fabricated nucleotide sequences; codewords distinguished from DNA |
| Target-binding design | Length, GC, exon fallback, homopolymer, off-target expression burden and cross-hybridization filters recorded from S01 |
| Barcode/readout logic | 20 pseudocolours ×4 rounds ×3 channels; codebook counts and symbol distance computed; QC repeat distinguished |
| Official software | Image-analysis code verified; original design-code availability unresolved; third-party implementation identified separately |
| Public probe resources | Codebook verified; sequence-level completeness unresolved; processed data separately classified |
| Species/reference context | Mouse samples and cell/tissue scope recorded; assembly/annotation unresolved rather than inferred from year |
| Benchmark candidates | D01 C, D02 D for sequence-design purpose; no strict ground truth claimed |
| Evidence traceability and uncertainty | Stable S01–S11 records; reported/derived/inferred/unavailable distinctions and exact locators retained |
| Metadata/Markdown consistency | Structural and cross-file self-check performed; final check results recorded below |
| Research-stage boundary | No validation_report.md or fix_report.md created; no independent validation claimed |

## Research Self-Check

Completed on 2026-10-15. `evidence/self_check.py` checked all eight required Research files, all template-required headings, YAML duplicate keys/field structure and relevant types, null field-state consistency, S01–S11 references, Markdown table column counts, local routing links, method naming, lifecycle consistency and the absence of Review/Fix reports. Results: `evidence/self_check_results.json`. `git diff --check` also completed without findings for tracked edits.

Manual cross-file reconciliation checked 28-nt binding regions versus 15-nt readouts and approximately 93-nt manufactured primaries; 24 primary probes/gene; four barcoding rounds versus 80 serial hybridizations plus QC repeat; channel-specific 3,334/3,333/3,333 counts; C/D candidate grades versus aggregate C; public codebook availability versus unresolved complete sequences; and unresolved design-code availability versus confirmed analysis code. Metadata and routing both use `research_complete`. No unsupported sequence was added to reconcile missing order or length information.

The repository schema does not enumerate allowed biological category/application strings despite referring to terminology-governed classification lists; descriptive strings consistent with S01 were used without changing the global schema. `versions/` was not created because no separate fully researched version file is needed; DNA lineage remains an explicit boundary for later work.

These are Research Agent self-checks and cannot substitute for Stage 2 independent review. No validation status is claimed.
