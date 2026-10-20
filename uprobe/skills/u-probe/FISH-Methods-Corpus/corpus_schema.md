# Corpus Specification

## Scope and authority
This English-language corpus records FISH, in situ sequencing, and spatial molecular imaging methods for evidence-grounded retrieval and probe-design comparison, including comparison with U-Probe. This specification and terminology.md govern all method records. Markdown is the detailed evidence layer; YAML provides compact routing and high-level metadata. No field requires an unsupported biological value. Templates are unfilled forms, not evidence.

## Directory structure
Preserve existing directory spelling and official capitalization. The initial routing labels are user-supplied identifiers, not independently verified official names. Each researched method must have:

```text
METHOD/
├── index.md
├── metadata.yaml
├── description.md
├── protocols.md
├── probes.md
├── benchmarks.md
├── sources.md
├── research_manifest.md
├── validation_report.md
├── fix_report.md
└── versions/
    └── optional version-specific files when justified
```

The seven production files are index.md, metadata.yaml, description.md, protocols.md, probes.md, benchmarks.md, and sources.md. The three audit/construction files are research_manifest.md, validation_report.md, and fix_report.md. Ordinary retrieval excludes audit files; load them only for explicit provenance, inspection-history, or validation-history questions. sources.md remains part of normal evidence retrieval.

At specification setup, method directories may be empty. Index paths are planned destinations, not assertions that files exist. Do not copy empty templates into method directories merely to imply research progress. Create the required files and versions/ directory when method construction begins.

## Evidence and claim attribution
Use exactly these evidence states:
- reported: Explicitly stated in a primary or official source.
- derived: Directly calculated from reported information without biological assumptions. Record inputs, calculation, units, and source locator.
- inferred: Interpreted from figures, code, diagrams, or indirect evidence. State the interpretation and uncertainty; do not present it as reported.
- unavailable: Not reported by the examined official sources. State the examined scope; this does not establish universal absence.

Every substantive claim, parameter, sequence, design rule, software record, and grade needs a method-local Source ID, precise locator (page, section, figure, table, file, line, or accession), method version, and evidence state. These may be shared by a clearly delimited table or paragraph. Secondary-only evidence must be labeled inferred and attributed to the secondary source until verified against primary/official evidence. Preserve conflicting claims separately with their sources and versions; do not silently select or average them.

Approximate values remain approximate; ranges remain ranges. Preserve reported units and qualifiers; any conversion is a separately labeled derivation. Never silently reconstruct missing sequences or fill protocol values from common practice. Do not generalize experiment-specific examples. Do not merge later versions into the original. Unsupported fields may remain empty/null or explicitly unavailable. “Not reported” differs from “Not applicable.” Uninspected information is unknown, not unavailable evidence.

## Method classification (probe_architecture and readout_strategy)

Each method in this corpus is classified along two orthogonal, controlled dimensions. Classification values are short controlled labels that mirror the method's probe design and decoding logic; detailed evidence remains in `description.md`, `probes.md`, and `sources.md`. Classification is **informative, not normative**: a controlled label means a documented fit, not a quality judgment. These fields are not used to rank, score, or exclude methods.

### Primary classification: probe_architecture

A single controlled label describing the molecular architecture of the recognition-and-encoding reagent that physically contacts the target (or the immediate template derived from it).

Allowed values (controlled vocabulary, lowercase, snake_case):

| Value | Meaning | Required structural evidence (probes.md) |
|---|---|---|
| `linear` | One or more contiguous oligos that hybridize directly to the target (or to a cDNA) and carry encoding/readout sequence features inline on the same molecule; no mandatory circularization step. | Canonical Architecture shows a single-strand 5'→3' representation; no mandatory intramolecular ligation of a circular template. |
| `padlock` | At least one oligo whose 5' and 3' ends are brought together on a target (or a cDNA / splint) and joined by a ligation reaction to form a circular template that is amplified in situ (typically by RCA). | Canonical Architecture shows a circularized padlock; Probe Components include `padlock_arm` and a ligation site; `amplification` references RCA. |
| `crispr` | A Cas-family protein (commonly nuclease-dead Cas9, dCas9) programmed by a guide RNA (crRNA/tracrRNA or sgRNA) that confers sequence specificity, with a fluorophore or readout moiety on the protein or the guide. | Probe Components include guide-RNA species (sgRNA, crRNA, tracrRNA); no target-binding oligo with sequence-complementary arms. |
| `nanostructure` | DNA or other nucleic-acid nanostructures (tetrahedra, dendrimers, branched assemblies) that act as signal amplifiers attached to a primary recognition element; recognition does not depend on the nanostructure itself. | Probe Components reference predefined nanostructure monomers and a primary recognition probe; Canonical Architecture shows the primary probe plus amplifier hierarchy. |

These values are the **initial controlled set**. A method may justify an extension by documenting the structural distinction in `corpus_schema.md` and `terminology.md`; never invent a new label by silently overriding the schema.

### Secondary classification: readout_strategy

A multi-label set describing how encoded identity is converted into a digital readout. Multiple labels are allowed when the method supports several readout modes in different versions or experiments.

Allowed values:

| Value | Meaning | Required structural evidence (probes.md / protocols.md) |
|---|---|---|
| `sequential_fish` | Repeated cycles of fluorescent-oligo hybridization, imaging, and signal removal (chemical cleavage, photobleaching, displacement). Each cycle reads one or more bits/channels; identity is decoded from the temporal pattern. | Readout System records multiple imaging rounds and a stripping/cleavage chemistry between rounds. |
| `in_situ_sequencing` | A sequencing reaction (sequencing-by-ligation, SOLiD, SBS, or comparable) performed in situ on RCA amplicons or other amplified templates, read out via fluorescent query probes. | Readout System records a sequencing primer / query-probe ligation or incorporation step and reads nucleotide identity base-by-base. |
| `fluorescent_barcode` | Identity is encoded by the intensity ratio, color combination, or co-localization pattern of fluorophores within a single imaging round (no temporal sequencing of bases). | Readout System records single-round multichannel acquisition with intensity- or color-coded codewords. |
| `live_imaging` | The assay is compatible with live-cell imaging (e.g., time-lapse tracking of loci or transcripts) without requiring cell fixation or in situ chemistry beyond reagent delivery. | protocols.md records live-cell-compatible workflow; description.md Limitations / Strengths note live-cell constraints. |

These values describe **readout mechanisms**, not imaging modalities in general (e.g., widefield vs. confocal microscopy). Confocal microscopy used solely to image a sequential FISH readout is still `sequential_fish`. Likewise, "live imaging" applies only when the assay is designed for live cells; fixed-cell compatibility of a method that is otherwise `in_situ_sequencing` does not add the `live_imaging` label.

When a method supports more than one readout mode (e.g., MERFISH with sequential FISH plus an optional bDNA amplification still read by sequential FISH), list every documented label.

### Recording classification

`probe_architecture` is a single string. `readout_strategy` is a list of strings. Both are governed by this section; do not redefine the labels in `terminology.md` or method files. Extensions and amendments to the controlled vocabulary live here and in `terminology.md` only.

In `metadata.yaml`, add:

```yaml
classification:
  probe_architecture: null      # null | linear | padlock | crispr | nanostructure
  readout_strategy: []          # any of: sequential_fish, in_situ_sequencing, fluorescent_barcode, live_imaging
  evidence_status: unknown      # reported | derived | inferred | unavailable
  source: null                 # method-local Source ID supporting the classification
  notes: null                 # optional brief rationale; cross-reference description.md
```

The pre-existing `classification.categories` and `classification.applications` fields remain for free-form author terminology and are not replaced by the controlled labels. Adding the controlled labels is an extension, not a substitute.

In `description.md`, add a `## Classification` H2 section that explains the assignment in plain language, citing the version and source locator. This is the canonical prose site for the rationale; YAML holds the structured record only.

In `probes.md`, structural evidence for the `probe_architecture` label is found in `Canonical Architecture` and `Probe Components`. Structural evidence for `readout_strategy` is found in `Readout System`. `Classification` in description.md cites those sections.

### Version handling

Different versions of the same method may carry different labels (e.g., padlock chemistry plus rolling-circle amplification in one version, linear probes with branched-DNA amplification in another). Record the labels at the method level using the dominant or original-version architecture; record version-level overrides in the corresponding `## Major Versions or Variants` section of `description.md` and in `metadata.yaml.versions`.

### Validation rules

- `probe_architecture` must be one of `linear`, `padlock`, `crispr`, `nanostructure`, or `null`. Other values are ERROR.
- Every element of `readout_strategy` must be one of `sequential_fish`, `in_situ_sequencing`, `fluorescent_barcode`, `live_imaging`, or a value added by an explicit schema extension in this file. Other values are ERROR.
- A populated value must have a non-null `source` and a non-empty `description.md` Classification section. Missing source or rationale is a WARNING.
- If the labels do not match the evidence in `probes.md` Canonical Architecture / Readout System, the finding is an ERROR.
- `null` and `[]` are acceptable only for `not_researched` methods; once `metadata.yaml.validation.status` advances past `not_researched`, classification fields should be populated (use `unavailable` for `evidence_status` when no official source is accessible).
- Do not invent extensions silently. If a controlled value does not fit, either reuse the closest existing label and explain the limitation in `notes`, or propose a schema amendment by editing this file before using the new label anywhere.

## Missing values
In Markdown use: “Not yet examined.” for pending inspection; “Not reported.” or “Not available in the examined official sources.” for bounded unsuccessful inspection; “Not applicable.” only with a justification; “Exact sequence not provided.”; “Architecture reported, sequence unresolved.”; or “Version-specific.” with a link to the relevant scope. Templates use Not yet examined., not claims about source absence.

In YAML:
- Scalar null means unknown, unexamined, unresolved, or unavailable; Markdown explains which and why.
- Empty lists mean no entries recorded, not proven absence.
- Boolean true requires positive evidence; false requires explicit evidence of absence within a documented version/source scope. Failure to find a tool or dataset is null, not false.
- To distinguish applicability without mixing string tokens into typed values, an optional shallow field_states mapping uses dotted field paths as keys and values unknown, unavailable, or not_applicable. A null without a state defaults to unknown. unavailable requires examined-source documentation; not_applicable requires justification in Markdown. Never apply these states to a populated value.
- Do not use strings such as "unknown", "N/A", or "unavailable" in typed scalar fields. Controlled status fields use their defined enumerations.
- Metadata is a summary. Facts that differ by version remain null at aggregate level unless explicitly scoped in Markdown; do not manufacture a single combined value.

## File contracts

### index.md
Fast routing only. Required headings: Official Name; Full Name; Aliases; Method Category; Primary Application; Original Publication; Brief Summary; Known Versions; Corpus Files. Use one H1 for the method label and H2 for these sections. Keep summaries concise and source-linked. No detailed protocol or probe-design content. Corpus Files links to all seven production files; audit files are identified separately if linked.

### metadata.yaml
Use the exact top-level groups and fields in metadata.template.yaml. All groups/keys are structurally present, but unknown values remain null or empty lists. No biological field must be populated to pass structural checks.
- method: official_name and full_name are nullable strings; aliases is a string list.
- classification: categories and applications are string lists governed by terminology.md; preserve author wording in detailed Markdown. The controlled probe_architecture (nullable string) and readout_strategy (string list) fields are defined under "Method classification" above; record classification only via the schema-defined keys.
- publication: original_year is a nullable integer; original_reference is a nullable Source ID; doi and pmid are nullable strings (PMID must not become a numeric identifier).
- versions: list of shallow records with version_id, version_name, year, reference, status. version_id is required for a recorded version; other values may be null. year is an integer, reference is a Source ID, and status is null or original, variant, superseded, or unclear. These describe documented lineage, not corpus validation or universal adoption.
- probe: architecture_type, target_binding, barcode_system, readout_system, amplification, and sequencing_or_imaging are nullable short strings; modifications is a string list.
- design: official_design_software_available and public_probe_sets_available are nullable booleans.
- benchmark: benchmark_available is a nullable boolean indicating a public candidate resource is evidenced, not that it is reproducible; grade is null or A, B, C, D. Grade summarizes the highest justified candidate grade within documented scope, never an average; link all individual grades in benchmarks.md. No assessed candidates means null, not D.
- references: primary_source_ids is a list of method-local Source IDs.
- validation: status uses the corpus lifecycle below, not the reviewer verdict.
- field_states: optional mapping as defined above.

Biological metadata must be traceable to corresponding Markdown sections and sources; avoid per-field deeply nested evidence objects. Extensions require documentation here and a template update before use.

### description.md
H1 Method Overview. Required H2 headings are exactly those in description.template.md: Purpose; Core Principle; Overall Workflow; Target Recognition Strategy; Encoding Strategy; Amplification Strategy; Imaging / Readout Strategy; Decoding Strategy; Major Versions or Variants; Strengths; Limitations; Unresolved or Version-Specific Points. Explain concepts, not executable protocol details. Keep all headings; justified Not applicable. is allowed. Strengths and limitations require evidence, not promotional interpretation.

### protocols.md
H1 Experimental Protocol. Use all H2 headings in protocols.template.md as stable retrieval anchors. Protocol Scope and Unreported Parameters are required content; procedural sections may be marked Not applicable. with justification, or Not yet examined./Not reported. as appropriate. Never infer a missing experimental step from its heading.

Parameter records preserve reagent, concentration, temperature, duration, buffer, pH, cycle number, imaging round, and sequencing cycle exactly as reported. Use the standard parameter table in the template, one record per reported step/condition; include qualifiers, units, source locator, version/sample scope, and evidence. A blank value is not a zero. Separate alternative protocols, samples, and versions; do not combine conditions into a synthetic protocol.

### probes.md
H1 Probe Architecture. This is the highest-priority comparison file. All H2 sections in probes.template.md are required retrieval anchors; irrelevant sections may be justified Not applicable.

Canonical Architecture uses explicit 5' → 3' order for each distinct molecular species when evidenced. Record strand identity, orientation, pairing, topology, and chemical groups separately; do not force a circle, duplex, or fluorophore into a fictitious linear nucleotide sequence. Show unresolved positions as [Unresolved component] and unknown order explicitly; a neutral template diagram is not method evidence.

Probe Components uses exactly the minimum columns in the template: Component, Original term, Function, Position, Orientation, Length, Sequence, Sequence class, Required, Design rule, Evidence status, Source / locator. Component is a standardized functional identifier; give distinct instances stable local suffixes when needed. Original term preserves exact author terminology. Source / locator is the single canonical column name for the source reference across every table in probes.template.md; the locator must include page, section, figure, table, file, line, or accession, and the column always carries the method version alongside the locator. Length preserves exact/range/approximate values and units; unknown is explicit. Required is yes, no, conditional, or unclear. Sequence class is fixed, variable, barcode_library, sample_specific, version_specific, or unresolved. Class describes variation, not evidential precision.

Supported component concepts: target_binding_region, target_binding_arm, padlock_arm, barcode, readout_sequence, primer, linker, spacer, adapter, ligation_sequence, amplification_sequence, promoter, fluorophore, chemical_modification, cleavage_site, anchor_sequence, other_functional_sequence. Use only applicable components; roles may overlap when justified, but do not duplicate physical sequence length.

Sequence representation must separately identify exact sequence, sequence example, consensus, motif, schematic sequence, degenerate sequence, reverse complement, variable region, or fixed region using terminology.md. Record molecule type, direction, source notation, version, and provenance. Exact published sequences require a textual/table/file source; a diagram alone is never an exact sequence. Any reverse complement or other transformation is derived with its operation documented. Do not overwrite original notation. Large public sequence collections should be linked with accession/file/version/checksum when available, not silently transcribed or reconstructed.

Target-Binding Design supports target type, target region, probe length, GC constraints, Tm constraints, specificity filtering, repeat masking, exon selection, isoform handling, sequence complexity, spacing, overlap, and number of probes per target. Report each rule's scope and whether it is a requirement or an example.

Barcode System supports barcode length, number of barcodes, codebook, Hamming distance, error correction, round assignment, channel assignment, and fixed versus experiment-specific sets. Distinguish molecular sequence from abstract codeword.
Readout System supports readout sequence, readout probe, fluorophore, imaging round, sequencing primer, and decoding role.
Primer / Linker / Adapter Components, Amplification-Related Components, and Chemical Modifications cross-reference component records rather than inventing missing parts.
Probe-Design Workflow records evidenced inputs, reference/annotation versions, ordered operations, parameters, and outputs.
Official Probe-Design Software uses the template table: software name, source/repository, official status, version, programming language, input, output, design logic, dependencies, availability, maintenance status; add evidence and source. Official status requires author/developer endorsement, not merely a matching repository name. Availability and maintenance are time-scoped observations; an old commit alone does not establish abandonment.
Public Probe Sets records dataset ID, release/accession, version, contents, availability/access date, and source; benchmark assessment belongs in benchmarks.md.
Version Differences and Unresolved Information retain scoped differences and gaps explicitly.

### benchmarks.md
H1 Benchmark Resources; H2 Benchmark Summary and Candidate Datasets. Each candidate gets a stable Dataset ID and the complete field/value record in benchmarks.template.md. Unknown fields remain explicit. Reference context must include genome build, transcriptome/annotation release, coordinate system (assembly/accession, chromosome/transcript, 0/1-based, closed/half-open, strand), target identity/version, and relevant design metadata when available. Never guess coordinate conventions.

Grades are exactly:
- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A. Assess target/reference context and relevant design metadata. State the intended comparison scope, missing information, and why gaps are minor for B. Unassessed candidates have null/unassessed grade, never a guessed D. Record assessment date, evidence, source locator, and version; do not assign grades during specification setup.

### sources.md
H1 Sources; required H2: Primary Publications; Supplementary Materials; Official Protocols; Official Repositories and Software; Official Datasets; Later Developer Publications; Secondary Sources. Assign each source one canonical record; cross-link it from other categories when needed. Record fields: Source ID, Source type, Citation, Title, Authors, Year, DOI, PMID, Official URL, Repository / accession, Method version, Relevant figures, Relevant tables, Relevant supplementary files, Information supported, Notes / limitations. Also record access date and release/commit/file checksum when available for mutable resources. Unavailable bibliography fields are allowed. A source record does not prove it was inspected; the manifest records actual inspection. Distinguish official status from hosting location.

### research_manifest.md
Audit only. H1 Research Manifest. Use records containing Source ID, Source type, Location, Files/tables/figures examined, Information extracted, Limitations, Unresolved questions; also inspection date and inspected revision when available. Required H2 Unresolved Research Questions. Record actual inspection, not planned or assumed inspection.

### validation_report.md
Independent review output, never a location to edit production claims. H1 Validation Report. Each issue has Issue ID, Severity, File, Section, Current claim, Problem, Verified evidence, Source, Required action. Severity is ERROR, WARNING, or INFO. Final status is PASS or FAIL only after review; an unperformed review has null. PASS requires no unresolved ERROR and completed scoped checks; warnings remain visible. Report scope/date/revision, including limitations of independence. Reviewers record findings; production changes are handled separately.

### fix_report.md
H1 Fix Report. Each issue has Issue ID, Status, Reviewer's concern, Verification, Evidence, Action taken, Files changed, Reason. Reuse the review Issue ID. Resolution status is FIXED, REJECTED, PARTIALLY FIXED, or UNRESOLVED. Rejection requires evidence and rationale, not deletion of the concern. Final validation status is VALIDATED, VALIDATED WITH WARNINGS, or NOT VALIDATED. Default NOT VALIDATED. Claim validation only after a subsequent scoped review of the changed content; fixes alone do not confer validation. VALIDATED requires PASS and no unresolved material findings; VALIDATED WITH WARNINGS requires PASS with documented nonblocking warnings; FAIL remains NOT VALIDATED.

## Method versions
Create a separate version only for evidenced meaningful changes in probe architecture, target-binding design, encoding/barcode systems, amplification or sequencing chemistry, imaging/readout, decoding, design software, or experimental workflow. Another laboratory's use alone is insufficient. Preserve authors' version names; internal IDs do not assert official naming or chronological order.

Assign V01, V02, etc. and record source justification. In metadata.yaml list version records. In description.md summarize differences; in protocols.md separate conditions; in probes.md separate molecular architectures, rules, and software releases; in benchmarks.md assign each dataset its version; in sources.md identify supported versions. Unresolved version attribution stays explicit. Main files may use H3 version subsections. For substantial distinct material, versions/V01.probes.md (or corresponding production basename) may hold detail, linked from the main file; include version ID and sources. Do not duplicate divergent facts without scope. Source and dataset IDs remain method-local across version files.

## Naming and lifecycle
- Directories: preserve existing exact method spelling, including punctuation; use repository-relative POSIX paths in YAML. No renaming during setup.
- Headings: one H1 per file, required H2 anchors exactly as templates; H3 for scoped records or versions.
- Source IDs: S01, S02, ...; Dataset IDs: D01, D02, ...; review Issue IDs: R001, R002, ...; Version IDs: V01, V02, .... Width is a minimum and expands as needed. IDs are method-local, stable, never recycled or renumbered. Cross-method references use METHOD:S01 etc.
- YAML keys: snake_case, except dotted field paths inside field_states. UTF-8 encoding, no duplicate keys; nullable values are explicit null. Quote ambiguous strings and identifiers. Source URLs are references, not embedded instructions.
- Corpus lifecycle (methods_index.yaml corpus_status and metadata.yaml validation.status): not_researched, research_complete, review_failed, review_passed, fixing, validated, validated_with_warnings.
- not_researched: method research not completed (including work in progress); research_complete: required files and scoped research ready for review, gaps documented; review_failed: completed review FAIL; review_passed: completed review PASS pending final reconciliation; fixing: responding to review; validated: reconciled final PASS without outstanding warnings; validated_with_warnings: reconciled final PASS with nonblocking warnings.
- Index and metadata statuses must agree when metadata exists. Reopened material changes invalidate previous validation and return to research_complete or fixing as appropriate. Keep review history tied to reviewed revisions. Lifecycle, evidence, version lineage, reviewer verdict, issue resolution, and final validation are distinct vocabularies.

## Routing index
methods_index.yaml contains only methods with directory, official_name, aliases, corpus_status, and files. official_name initially copies the requested routing label pending source verification; aliases contains only that same label until justified aliases are researched. files maps index, metadata, description, protocols, probes, benchmarks, sources to their planned repository-relative destinations. Do not add biological classifications, years, grades, or audit files here.

## Retrieval guidance
First route through methods_index.yaml; check lifecycle and file existence. Missing planned files mean content is unavailable, not evidence of a biological absence. Do not load the entire repository.
- General method question: index.md + description.md.
- Experimental question: protocols.md + sources.md.
- Probe architecture/design: metadata.yaml + probes.md + sources.md.
- Benchmark: benchmarks.md + probes.md + sources.md.
- Version: metadata.yaml + relevant version sections/files + sources.md.
Follow exact source locators as needed. Preserve uncertainty and version scope in answers. Audit files are excluded unless inspection/provenance/validation history is requested.

## Internal consistency checks
Check all template headings and field contracts, terminology mappings, exact directory/path spelling, YAML parsing/types/duplicate keys, lifecycle enums, evidence enums, validation enums, and identical grade definitions. Check that production and audit routing remain separate, placeholders do not imply evidence, and no method claims, sequences, protocols, or grades were introduced during setup. Routing targets may be absent only for not_researched entries. This structural review does not validate any method.

## Required protocol H2 headings
- Protocol Scope
- Sample Preparation
- Fixation
- Permeabilization
- Reverse Transcription
- Probe Hybridization
- Washing
- Ligation
- Amplification
- Polymerization
- Sequencing / Readout
- Imaging
- Decoding
- Critical Parameters
- Version-Specific Differences
- Unreported Parameters

## Required probe H2 headings
- Overview
- Canonical Architecture
- Probe Components
- Target-Binding Design
- Barcode System
- Readout System
- Primer / Linker / Adapter Components
- Amplification-Related Components
- Chemical Modifications
- Probe-Design Workflow
- Official Probe-Design Software
- Public Probe Sets
- Version Differences
- Unresolved Information


