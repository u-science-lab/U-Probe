# Research Manifest

Audit only; exclude from ordinary methodological retrieval. Initial research performed 2026-10-15. Independent Review and Fix have not been performed.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Peer-reviewed publication | Publisher article and Fig. 1 page | Abstract; Fig. 1 caption; Extended Data 1–10 captions; data/code availability; supplement links | Naming, configurations, mechanism, public resource links | Final Methods subscription-restricted; no full body inspection | Final-method comparison to preprint | 2026-10-15 | Version of record 2025-10-30 |
| S02 | Official XLSX | evidence/supp4.xlsx | All 9 sheet dimensions/sequence inventory; imaging C/D/H/I columns; merged labels; actual/theoretical pooling; padlock sequences; primer footnotes | Exact sequence resources, component segmentation, modifications, lengths, mappings, 64-panel note | No transcript reference or original alignments established; validation sheet inventoried, not every smFISH design analyzed | Original target mapping; prefix function; HBX design exception | 2026-10-15 | SHA-256 c791073a5336e058ec97aa8c74ea4b11bec021dc884ff7752ca4c9099fb1c466 |
| S03 | Official XLSX | evidence/supp3.xlsx | Sheet1 dimensions, sample group headings and entries exported to JSON | Selective-amplification dilution resource | No competitor architecture determined | Irreplicable-padlock construction | 2026-10-15 | SHA-256 ec3c8951c9d20235880348cdcbc79fbaaeab4c08f564b29b291819920ed574bb |
| S04 | Official PDF | evidence/supp1.pdf | All captions extracted; captions 1–7 and 38–41 read; Fig. 2 PDF page 3 rendered/visually read | Barcode pseudocode; validation scope; empirical calibration | Not all image panels visually examined; PDF is figures rather than complete Methods | Complete final protocol | 2026-10-15 | SHA-256 bb4cd777df479fd344656c595af341392451d8be880128bb711795dd8537951f |
| S05 | Author preprint | evidence/preprint2.html and preprint2.txt; biorxiv_api.json | Fig. 1 text; conceptual passage; Materials and Methods Mice through cell assignment | Design criteria; direct RNA recognition; thin/thick protocols; imaging/decoding | v2 only; v1/v3 requests 429; no asserted equivalence to final paper | Mouse target-number discrepancy and final-method changes | 2026-10-15 | v2 2024-11-22; API also lists v1/v3 |
| S06 | Author repository | evidence/software/PRISM-Code | README; two codebook CSVs; pyproject | Official designer link; post-stitching scope; codebook conventions | Downloaded source; not executed; upstream spatial_img_core unavailable by README | Historical software revision and original codebook assignment | 2026-10-15 | ab7b72d |
| S07 | Author-linked repository | evidence/software/probe-designer | README; CHANGELOG first sections; pyproject/requirements; two config YAMLs; assembly and SingleSequenceStrategy source | Input/output, reverse complement/arm order, assembly, current config differences | Not executed; no original experiment replay; no unsupported iLock transfer | Original producing version; README/package/changelog version mismatch | 2026-10-15 | d7326be |
| S08 | Author repository | evidence/software/PRISM-Analysis | README; archive directory inventory | Downstream analysis purpose; conflicting expansion | No notebook/data execution; large images not used | README naming discrepancy | 2026-10-15 | 3b10d61 |
| S09 | Author dataset | evidence/zenodo.json; zenodo_readme.md | Metadata file inventory and readme | Archive/run layout; separates image/analysis data from design truth | Multi-GB archives not downloaded; other raw deposit IDs found only through paper/README | Raw-deposit version relationships; design-run inputs | 2026-10-15 | 12755414 revision 4; updated 2024-07-25 |
| S10 | Unrelated primary source discovery | PMC6763432 search-result abstract | Abstract naming only | Excludes protein PRISM homonym | No protein-method evidence imported | None within selected RNA-method scope | 2026-10-15 | Not recorded |

## Unresolved Research Questions

- Recover final 2025 Methods and compare all conditions with the preprint v2; confirm its DAPI concentration units before experimental reuse.
- Determine original transcript accessions/versions, reference assembly/annotation and coordinate conventions; do not infer them from current software configs.
- Resolve two mouse targets per gene in v2 versus one listed final mouse padlock per gene; explain the 28-nt HBX residual target region.
- Establish the purpose of the V01 two-nucleotide backbone prefix and verify target-pairing orientation against original reference sequences.
- Resolve original selective-amplification competitor construction and exact design-run settings.
- Preserve channel order and V02 high-G digit 5 versus level 2; establish explicit conversion before combining codebooks.
- Determine how the current author-linked designer relates to publication-time code; resolve package/changelog version labels.
- Verify standalone official protocol availability when spatialprism.org is accessible. No author was contacted.
- Resolve raw Zenodo identifier/version relationships before large-scale image benchmarking.

## Research Checklist

- Official name and original publication: checked against primary publisher abstract; homonym excluded; PRISM directory spelling retained.
- Primary, supplementary, code, dataset, protocol and public probe searches: performed; actual access successes/failures above. Search terms included PRISM RNA spatial imaging, full paper title, padlock/ligation, preprint DOI and spatialprism protocol. Secondary search results were discovery only.
- Versions: V01/V02 configurations separated; preprint v2 and 2026 software never assumed equivalent to the final paper.
- Probe architecture and orientation: explicit strand/topology schematics supplied with inferred scope; raw modification notation retained; target alignments remain unresolved.
- Components and terminology: individually mapped; schematic vs exact sequences separated; library strings not treated as universal constants.
- Target-binding rules, barcode/readout logic, software and public probe sets: investigated and source-scoped.
- Species/reference and benchmarks: species recorded; unknown original references remain explicit; D01/D02 graded C, D03 D for probe-design reconstruction only.
- Traceability: source IDs and section/sheet/cell locators supplied; downloaded primary inputs have checksums. Large repository archives were inspected via inventories and selected text, not executed.
- Cross-file consistency and template structure: checked by evidence/self_check.py; this is a Research self-check, not independent validation.
- Metadata status and routing: research_complete; no validation_report.md or fix_report.md created.
- Preflight: existing repository preflight report records closure of its four blocking global errors; schemas/templates were re-read. YAML/routing/type checks performed without modifying other methods. Existing ExSeq, MERFISH, seqFISH+ and routing changes were preserved.

## Evidence Preservation

Official workbook/PDF inputs, inspected HTML and extracted text, source-download inventories, selected pinned software files, raw workbook JSON exports and derived sequence inventory are under evidence/. Whole repository archives and unused extracted media were removed after retaining inspected files and archive hashes. No sequences were invented. Missing values reflect inspected-source boundaries or explicit access limitations.
