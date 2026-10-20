# Benchmark Resources

## Benchmark Summary

The highest justified grade in this draft is **C**. Actual spacer lists and construct definitions are public, but none of the examined candidates is accepted as strict probe-design ground truth. Grade C is conservative for the intended full design-output comparison, not a claim that only diagrams exist. V01; S02–S05; assessed 2026-10-15.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A. Unassessed is not D.

## Candidate Datasets

### D01 — Published genomic spacer lists

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Published genomic spacer lists in Table S2 |
| Method version | V01 (2025) |
| Species | Homo sapiens; S01 design methods |
| Cell / tissue / sample | U2OS/HeLa and human-cell applications in S01; per-experiment cell scope must be retained, not assigned from workbook names alone |
| Genome build | hg38, reported in S01 design methods and S02 Table S1(7) |
| Transcriptome / annotation version | Not applicable to genomic DNA spacers; enhancer/gene annotation release not established |
| Target genes | ACAP2-associated regions, TUBGCP3, PLCXD1, MUC4, FOS and MYC promoter/enhancer regions |
| Target regions | S02 Table S1(7); e.g. ACAP2 IDT region Chr3:195446441–195470879 and FOS promoter Chr14:75238783–75278649 (source notation retained) |
| Genomic/transcript coordinates | Chromosome/start/end in S03 columns A–C; assembly hg38 from S01/S02. No strand column or explicit 0/1-based and interval convention in workbook. End−Start usually equals spacer length; this does not establish coordinate origin. Code has both BED/TXT conventions and is not proof of the workbook convention. |
| Target-binding sequences | Exact spacer strings, DNA alphabet, column D on 11 genomic sheets; published RNA orientation is separately established in S02. Individual spreadsheet strand remains unresolved. |
| Full probe sequences | Not supplied as fully assembled, library-assigned, chemically labeled live-cell probes for every row. S02 provides construct/scaffold definitions. |
| Barcode assignments | Not applicable to the direct imaging architecture. Per-experiment color/pool identity is distinct from a combinatorial codebook. |
| Readout sequences | No separate sequential-readout library; tracrRNA/scaffold in S02 |
| Primer/linker sequences | PCR primers/adapters in S02 Table S1; per-library adapter assignment incomplete |
| Design software | Official S05 pipeline exists; the pinned public commit is not established as the exact generator of this workbook |
| Design parameters | Paper specifies length, GC, motifs, PAM, off-target and spacing filters (S01); code differences preserved in probes.md |
| Source | S03 Table S2; S02 Table S1; S01 design methods; S05 pinned code for context |
| Availability | Actual publisher XLSX downloaded 2026-10-15; the local copy is not included in this repository (Elsevier supplement; not redistributable); SHA-256 in sources.md |
| Missing information | Genomic strands; explicit workbook coordinate convention; full synthesized library assemblies/adapter mapping; exact generating software/environment; eight coordinate/string-length discrepancies in short ACAP2 sheet |
| Benchmark grade | C |
| Grade justification | Full design-output comparison is not strict ground truth: unresolved reference-strand/coordinate conventions and source discrepancies are material, rather than minor B-level gaps. Usable for architecture study and scoped spacer-list comparisons after independent mapping/QC. |
| Comparison scope | Assess genomic CRISPR guide selection/pool design, not RNA smFISH antisense probe selection. Do not compare strand-specific coordinates or full synthesized oligos as exact truth without resolving gaps. |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource contents; derived counts/consistency checks; inferred benchmark judgment |

The source workbook has 14 sheets. Eleven genomic spacer sheets contain 96, 96, 96, 136, 96, 96, 96, 237, 343, 340 and 375 sequence rows in the order listed in probes.md. Counts are derived from nonempty DNA strings in column D. RNA FISH and qRT-PCR sheets are separate validation reagents and excluded from those guide counts. S03.

The eight mismatches in `Chr3_ACAP2-nonrep_P5-P8 13-17bp` are preserved below. The calculation uses source columns B/C/D directly, without assuming an assembly coordinate convention:

| Excel row | End − Start | Spacer string length | Evidence |
|---|---:|---:|---|
| 6 | 17 | 16 | derived; S03 B6:D6 |
| 14 | 16 | 15 | derived; S03 B14:D14 |
| 21 | 17 | 14 | derived; S03 B21:D21 |
| 27 | 17 | 15 | derived; S03 B27:D27 |
| 45 | 14 | 16 | derived; S03 B45:D45 |
| 52 | 14 | 17 | derived; S03 B52:D52 |
| 57 | 15 | 14 | derived; S03 B57:D57 |
| 65 | 17 | 15 | derived; S03 B65:D65 |

No correction or cause is inferred. No independently remapped coordinate table is produced.

### D02 — Construct and repetitive-locus definitions

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | Table S1 construct/scaffold and repetitive-locus resource |
| Method version | V01 (2025) |
| Species | Homo sapiens; S01 design methods |
| Cell / tissue / sample | Human-cell imaging study; cell type depends on experiment |
| Genome build | hg38, reported in S01 design methods and S02 Table S1(7) |
| Transcriptome / annotation version | Not applicable to genomic DNA spacers; enhancer/gene annotation release not established |
| Target genes | ChrX-rep, Chr13-rep, ACAP2-rep, Chr3 sites 3/4, telomere; non-repetitive pool summaries also listed |
| Target regions | hg38 locus ranges and neighboring fixed-cell DNA FISH ranges in Table S1(7) |
| Genomic/transcript coordinates | Region-level hg38 ranges; no complete strand-resolved individual binding-site set for repeat copies; coordinate convention unresolved |
| Target-binding sequences | Published repetitive spacers (some shorter than the non-repetitive pool range), plus Table S2 links |
| Full probe sequences | Exact RNA common regions/tracrRNA; degenerate template schematics; not a fully instantiated synthesized pool |
| Barcode assignments | Not applicable to the direct imaging architecture. Per-experiment color/pool identity is distinct from a combinatorial codebook. |
| Readout sequences | No separate sequential-readout library; tracrRNA/scaffold in S02 |
| Primer/linker sequences | Exact primers, promoter and four adapter strings; see probes.md |
| Design software | Official S05 pipeline exists; the pinned public commit is not established as the exact generator of this workbook |
| Design parameters | Paper specifies length, GC, motifs, PAM, off-target and spacing filters (S01); code differences preserved in probes.md |
| Source | S02 Table S1(1–7), PDF pp.2–4 (printed pp.1–3) |
| Availability | Downloaded 2026-10-15; the local copy is not included in this repository (Elsevier supplement; not redistributable); SHA-256 in sources.md |
| Missing information | Per-copy repeat binding-site map; synthesized pool instantiation; unambiguous library-to-adapter mapping; exact chemical-label occupancy |
| Benchmark grade | C |
| Grade justification | Supports explicit molecular architecture and selected exact sequences. Does not reconstruct the full experimental design output for strict comparison. |
| Comparison scope | Architecture, scaffold/primer representation and repeat-guide examples; not whole-library design ground truth |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource contents; derived counts/consistency checks; inferred benchmark judgment |

### D03 — Official repository example

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | Official design repository example region and output illustrations |
| Method version | Current public design-code snapshot associated with V01; not demonstrated as the publication-generating revision |
| Species | Homo sapiens; S01 design methods |
| Cell / tissue / sample | Example region; no independently verified experiment/sample mapping |
| Genome build | hg38, reported in S01 design methods and S02 Table S1(7) |
| Transcriptome / annotation version | Not applicable to genomic DNA spacers; enhancer/gene annotation release not established |
| Target genes | Example Chr3 repeat-flanking region |
| Target regions | README/example BED: chr3 195458317 195478317; source names it Chr3Repeat_LFlanking20K |
| Genomic/transcript coordinates | Input BED and code-generated BED use 0-based intervals; code also writes 1-based-start TXT. End-exclusive BED interpretation is inferred from code slicing/writing. Genome source is hg38.fa in wrapper; file checksum not provided. |
| Target-binding sequences | README shows a few spacer examples, not a complete final experimental output set |
| Full probe sequences | Not supplied as complete experimentally synthesized probes |
| Barcode assignments | Not applicable to the direct imaging architecture. Per-experiment color/pool identity is distinct from a combinatorial codebook. |
| Readout sequences | No separate sequential-readout library; tracrRNA/scaffold in S02 |
| Primer/linker sequences | PCR primers/adapters in S02 Table S1; per-library adapter assignment incomplete |
| Design software | S05 pinned commit 210b34c2c54fd0f7e282a67e8eca144d96f44ad7 |
| Design parameters | Wrapper default 18–20 bp; paper/code differences and environment setup requirements in probes.md |
| Source | S05 README; example_data/chr3_repeatregion_lf20K.bed; code/generate_possible_guides.py and wrappers |
| Availability | Public source inspected 2026-10-15 |
| Missing information | Complete reference data/environment and execution record; full expected final output; exact connection to published pools |
| Benchmark grade | C |
| Grade justification | Useful pipeline/architecture example, but no complete expected design output was established for a strict reproducibility benchmark. |
| Comparison scope | Software inspection and future reconstruction work; not current strict ground truth |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource contents; derived counts/consistency checks; inferred benchmark judgment |
