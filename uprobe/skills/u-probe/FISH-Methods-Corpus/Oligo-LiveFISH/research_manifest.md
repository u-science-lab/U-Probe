# Research Manifest

Audit only; exclude from ordinary retrieval. This is initial Research-stage construction, not independent review or validation. Inspection date: 2026-10-15 (Asia/Shanghai).

## Repository Preflight and Scope

Read AGENTS.md, corpus_schema.md, terminology.md, methods_index.yaml and all Research templates. Reviewed the preflight report’s post-fix disposition; no unresolved global blocking schema ERROR was identified in that disposition. Existing unrelated ExSeq, MERFISH, seqFISH+ and other working changes were preserved. Initial Oligo-LiveFISH directory had no corpus files. No server access was needed.

The historical preflight note requests Oligopaint and Oligo-LiveFISH version files but supplies no primary evidence of that lineage. Its placeholder DOI is incorrect relative to the verified paper DOI. S01 explicitly identifies LiveFISH as earlier work and names the present method Oligo pool-based LiveFISH. It does not establish Oligopaint as a method version or alias. Following the governing schema’s evidence/version policy, this draft records only V01 and does not create the proposed unsupported lineage files. methods_index.yaml is updated only within this method entry: unsupported interchangeable aliases removed, lifecycle set to research_complete. No global schema or historical preflight record is rewritten.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Original article | PMC full author manuscript; PubMed | Introduction, design Results, Figures 1/2/5 captions, relevant variant/cell-type text; STAR Methods design, DNA pools, RNA preparation, dye labeling, delivery, sequential FISH, smiFISH, microscopy and image-analysis sections; availability and Key Resources | Name, architecture, rules, hg38, manufacturing and delivery parameters, assay separation | No fresh experiment or code reproduction; manuscript/typeset equivalence not assumed | Exact table-generating revision; publication-wide raw-data release | 2026-10-15 | NIHMS2067673 / PMID 40239646 |
| S02 | Official supplementary PDF | Publisher mmc1.pdf | All seven pages extracted; Table S1 pp.2–4; Table S3 p.5; Table S4 p.6; sequence pages 2/3 rendered and visually checked | Exact RNA sequences and terminal orientation; primers/adapters; degenerate template scope; regional target context | DNA orientation inferred from synthesis context; not every library instantiated | Adapter mapping and exact full synthesized pool | 2026-10-15 | SHA-256 in sources.md |
| S03 | Official supplementary workbook | Publisher mmc2.xlsx | All 14 worksheets, headers and populated rows; genomic A:E; RNA FISH/qPCR A:B; formulas retained | Spacer records and coordinate fields; counts; string lengths versus End−Start; eight mismatches | No external hg38 remapping; no strand column; no coordinate-origin declaration | Resolve discrepancies, strand, origin, missing short P1–P4 sheet | 2026-10-15 | SHA-256 in sources.md |
| S04 | Official website | Homepage, instructions, contact | Text and homepage links | Tool scope, workflow instructions, official GitHub endorsement | No completed online design job; backend not inspected | Functional service and release provenance | 2026-10-15 | Mutable snapshot |
| S05 | Official design repository | GitHub API tree and raw files at pinned commit | README, LICENSE, input BED, guide generator, wrappers, off-target summarizer, original-target/frequent-guide filters, overlap selector, run_all.sh | Inputs, outputs, dependencies, orientation code, special-case filters, deployment assumptions | Not executed; two_file_merge.py fetch failed and was not inspected; workflow image not used as sequence evidence | Exact publication linkage and full execution reproducibility | 2026-10-15 | 210b34c2c54fd0f7e282a67e8eca144d96f44ad7 |
| S06 | Lab analysis repository | GitHub API metadata/tree and README | README and recursive file inventory only | Analysis-only scope, examples distinct from guide-design resources | Source/data contents not independently validated | Link to original publication experiments | 2026-10-15 | 8d84763960c4e65a1919bbc5a93c84de9c554241 |

## Discovery and Access Record

Searched exact method name; primary-paper/PMC discovery; method + GitHub; method + Table S1/sequences; method + mmc; and protocol discovery including protocols.io. Followed the primary paper’s supplementary links, Data and code availability, author GitHub account, and the official website’s GitHub link. The search also surfaced later reviews and a spotlight; they were not used as primary technical evidence. No separate developer protocol publication or later version was established by this search; this is not proof of absence.

PMC attachment URLs returned HTTP-200 HTML download-challenge pages rather than PDF/XLSX payloads. Those responses were not treated as inspected supplements. Europe PMC’s supplementaryFiles endpoint returned “not open access”; publisher-hosted mmc1.pdf and mmc2.xlsx then returned genuine files and were inspected. The final retained resources are unmodified publisher files; rejected HTML payloads remain only in ignored temporary research work. A raw GitHub request failed for one helper; the manifest does not claim it was read. No credentials were requested or exposed.

## Research Checklist

- Official name, expansion, original DOI/PMID and V01 scope checked against S01.
- Primary methods, supplementary materials, official design software and public probe tables actively investigated.
- RNA orientation explicitly checked in S02; DNA orientation interpretations and spreadsheet strand gaps labeled.
- Probe species/components, exact fixed regions, degenerate templates and exact variable spacer examples separated.
- Target design, PAM/GC/motif/off-target/spacing rules extracted; chromatin context retained; barcode inapplicability explained.
- Public spacer resources and genomic assembly investigated; no transcriptome release invented.
- Candidate grades independently reasoned during Research; no A/B strict-ground-truth claim.
- Source IDs, locators, versions and evidence states provided; file links and source IDs checked structurally.
- Cross-file consistency and template headings checked; metadata and routing lifecycle agree.
- No validation_report.md or fix_report.md created; no independent review claimed.

## Unresolved Research Questions

1. What coordinate-origin/interval convention and genomic strand apply to each Table S2 record? Independent reference mapping is needed before strict strand-aware comparison.
2. What accounts for the eight short-spacer ACAP2 coordinate/string-length discrepancies? Published source data were preserved.
3. Are additional short P1–P4 data available elsewhere? No separately named sheet was found in the inspected workbook.
4. Which adapters and full DNA/RNA assemblies were used for each experimental pool? Exact dye occupancy/3′ extension lengths remain unresolved.
5. Which software revision and reference FASTA generated the published lists, and how should public-code/paper filter differences be reconciled?
6. Is the online tool fully functional and versioned? Documentation reachability is established; design execution was not tested.
7. Can the authors provide a complete public microscopy/probe-design reproducibility bundle? Availability on request is not equivalent to a public benchmark.
8. Is any additional publication sufficient to establish a later method version? No such version was established here.
