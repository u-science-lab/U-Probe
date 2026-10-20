# FISH Methods Corpus — Agent Instructions

## 1. Project Objective

This repository is a standardized, evidence-grounded corpus for FISH, in situ sequencing, spatial molecular imaging, and probe-design methods.

The corpus is intended to support:

- method retrieval and comparison
- experimental protocol retrieval
- probe architecture interpretation
- probe-design rule extraction
- barcode and readout system interpretation
- official probe-design software discovery
- public probe-set discovery
- benchmark construction
- downstream agent/skill retrieval
- comparison with probe-design software such as U-Probe

The corpus must prioritize:

1. factual accuracy
2. source traceability
3. reproducibility
4. version awareness
5. concise technical writing
6. machine-readable organization
7. explicit uncertainty

Completeness is never more important than correctness.

Never invent information merely to fill a schema field.

---

## 1A. Method classification (controlled vocabulary)

Every method in the corpus must record the controlled two-level classification defined in `corpus_schema.md` (`Method classification` section) and detailed in `terminology.md` (`Method classification vocabulary` section).

- Primary `probe_architecture` (one of `linear`, `padlock`, `crispr`, `nanostructure`, or `null`).
- Secondary `readout_strategy` (a subset of `sequential_fish`, `in_situ_sequencing`, `fluorescent_barcode`, `live_imaging`, possibly empty).

The labels are informative, not normative. They support retrieval; they do not rank or score methods. Do not invent a new label without amending `corpus_schema.md` first. Do not assign a label that the structural evidence in `probes.md` (Canonical Architecture, Probe Components, Readout System) does not support. When a controlled value does not fit, use the closest existing label and explain the limitation in `description.md` (`## Classification`) and in `metadata.yaml.classification.notes`. Multiple secondary labels are allowed when the method supports several readout modes; record every documented label.

`classification.categories` and `classification.applications` in `metadata.yaml` remain free-form and may hold author wording; the new keys are an addition, not a substitute.

Future agents must populate these fields for every new or migrated method, cite the supporting source by method-local Source ID, and add a `## Classification` section to `description.md` that cites the version and source locator.

---

# 2. Language

All production corpus content must be written in English.

This includes:

- index.md
- metadata.yaml
- description.md
- protocols.md
- probes.md
- benchmarks.md
- sources.md

Audit files should also preferably be written in English:

- research_manifest.md
- validation_report.md
- fix_report.md

Use concise, neutral, technical English.

Avoid:

- promotional language
- unnecessary historical narrative
- vague claims
- speculation
- exaggerated descriptions
- unsupported biological explanations

---

# 3. Repository Specification

Before performing method-specific work, read:

- AGENTS.md
- corpus_schema.md
- terminology.md
- methods_index.yaml
- relevant files in templates/

These files define the corpus standard.

Do not silently override them.

If a method exposes a limitation in the schema, document the problem separately rather than silently changing global rules during method research.

---

# 4. Method Directory Structure

Each method should eventually contain:

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
```

The `versions/` directory is optional and should only be used when meaningful method versions exist.

## Production corpus files

These files are intended for normal downstream Skill or Agent retrieval:

- index.md
- metadata.yaml
- description.md
- protocols.md
- probes.md
- benchmarks.md
- sources.md

## Audit files

These files document corpus construction and quality control:

- research_manifest.md
- validation_report.md
- fix_report.md

Audit files should normally not be loaded for ordinary biological or methodological questions.

---

# 5. Official Method Naming

Always use the method name used by the original authors.

Preserve:

- capitalization
- punctuation
- hyphenation
- plus signs
- spacing
- acronym style

Do not normalize a method name simply because another spelling appears more conventional.

Examples of potentially sensitive naming include:

- STARmap
- ExSeq
- seqFISH+
- MERFISH
- Oligo-LiveFISH
- RAEFISH
- TDDN-FISH
- PRISM
- OpenFISH
- MiP-seq

Do not invent an acronym expansion unless the original authors explicitly define it.

Do not treat:

- paper titles
- software names
- later aliases
- related techniques

as the official method name unless supported by primary evidence.

If the current directory name differs from the verified official short method name, report the discrepancy before renaming it.

---

# 6. Evidence Hierarchy

Prefer sources in the following order:

1. Original peer-reviewed publication
2. Official supplementary materials
3. Official experimental protocol
4. Official repository or software documentation
5. Author-hosted datasets
6. Later peer-reviewed publications from the original method developers
7. Independent peer-reviewed applications
8. Review articles

Blogs, commercial websites, copied protocol pages, informal summaries, and AI-generated summaries must not be used as primary evidence when original sources are available.

Secondary sources may be used for discovery, but technical claims should be traced back to primary or official evidence whenever possible.

---

# 7. Source Discovery Rules

Do not stop after finding only the main paper.

For every method, actively investigate whether the authors provide:

- supplementary methods
- supplementary tables
- supplementary spreadsheets
- supplementary sequence files
- probe lists
- barcode lists
- codebooks
- readout sequence libraries
- GitHub repositories
- design scripts
- notebooks
- official protocols
- public datasets
- repository source data
- author-hosted resources

A major goal of the corpus is to identify evidence that is often absent from the main article.

---

# 8. Evidence Status

Technical information must distinguish, where appropriate:

## reported

Explicitly stated in a primary or official source.

## derived

Directly calculated from reported information without biological assumptions.

## inferred

Interpreted from:

- figures
- diagrams
- source code
- sequence structure
- indirect evidence

but not explicitly stated by the source.

## unavailable

Not provided by the examined official sources.

Never present inferred information as reported fact.

Never present derived information as directly reported.

---

# 9. Missing Information

Missing information must remain missing.

Do not invent or guess values.

Use explicit wording such as:

- Not reported.
- Not available in the examined official sources.
- Not applicable.
- Exact sequence not provided.
- Architecture reported, sequence unresolved.
- Reported only for a later version.
- Unable to determine from available official evidence.

Do not treat these as equivalent:

- unknown
- unavailable
- not reported
- not applicable
- false

In YAML, prefer `null` for unresolved information unless a controlled status field exists.

For example:

```yaml
official_design_software_available: null
```

is preferable to:

```yaml
official_design_software_available: false
```

when the existence of official software has not yet been established.

---

# 10. Exact Values and Ranges

Preserve the precision of the source.

If a source reports:

```text
approximately 20 nt
```

do not write:

```text
20 nt
```

If a source reports:

```text
18–22 nt
```

do not write:

```text
20 nt
```

Ranges must remain ranges.

Approximate values must remain approximate.

Do not infer exact experimental conditions from standard laboratory practice.

---

# 11. Method Version Policy

Different versions of a method must not be silently merged.

A separate version may be justified when there are meaningful changes in:

- probe architecture
- target-binding strategy
- encoding strategy
- barcode system
- amplification chemistry
- sequencing chemistry
- imaging/readout strategy
- decoding strategy
- probe-design algorithm
- probe-design software
- experimental workflow

A later publication merely using the method does not automatically constitute a new version.

When version-specific information exists, identify:

- version name
- publication
- year
- source
- changes relative to earlier versions
- affected probe architecture
- affected protocol
- affected encoding
- affected benchmark resources

Never attribute later-version features to the original method without labeling them.

---

# 12. index.md

`index.md` is a routing document.

It should remain short.

Include:

- official method name
- full name, if officially defined
- aliases
- original publication
- method category
- primary application
- concise summary
- known versions
- available corpus files

Do not place detailed:

- protocol steps
- probe sequences
- barcode tables
- benchmark analysis

in `index.md`.

---

# 13. metadata.yaml

`metadata.yaml` provides machine-readable routing and high-level structured information.

It should remain:

- compact
- shallow
- readable
- conservative
- easy to query

Detailed evidence belongs in Markdown.

Do not duplicate complete:

- probe sequences
- protocols
- source annotations
- benchmark reasoning

inside YAML.

Unknown fields may remain `null`.

Do not populate fields merely because the schema contains them.

All YAML keys should use `snake_case`.

---

# 14. description.md

`description.md` explains the method conceptually.

It should cover, when applicable:

- purpose
- core principle
- overall workflow
- target recognition
- encoding strategy
- amplification strategy
- imaging/readout strategy
- decoding strategy
- meaningful variants
- strengths
- limitations

Do not turn `description.md` into a detailed experimental protocol.

Strengths and limitations must remain conservative and evidence-based.

Avoid promotional language such as:

- revolutionary
- superior
- extremely powerful
- best-in-class

unless reporting a clearly attributed source statement, which should still generally be paraphrased neutrally.

---

# 15. protocols.md

`protocols.md` describes experimental execution.

Preserve reported values such as:

- reagents
- concentrations
- temperatures
- incubation times
- buffers
- pH
- reaction volumes
- cycle numbers
- imaging rounds
- sequencing cycles

Do not fill missing values using:

- standard protocols
- similar methods
- common laboratory practice
- assumptions

If protocol values differ between method versions, record the version explicitly.

---

# 16. probes.md — Highest Priority File

Probe architecture and probe-design information are the highest-priority content in this corpus.

Every method should be investigated carefully for probe structure.

Whenever evidence permits, represent architecture explicitly in 5'→3' orientation.

Example:

```text
5' — [Component A] — [Component B] — [Component C] — 3'
```

If part of the architecture cannot be resolved:

```text
5' — [Target-binding region] — [Unresolved component] — [Barcode] — 3'
```

is preferable to guessing.

---

# 17. Probe Component Fields

For every relevant probe component, record whenever available:

- component name
- original author terminology
- standardized functional category
- function
- position
- orientation
- length
- length range
- exact nucleotide sequence
- sequence class
- required/optional status
- design rule
- evidence status
- source

Potential components include:

- target-binding region
- target-binding arm
- padlock arm
- barcode
- codeword-related sequence
- readout sequence
- readout probe
- primer
- sequencing primer
- amplification primer
- linker
- spacer
- adapter
- ligation sequence
- amplification sequence
- promoter
- anchor sequence
- cleavage site
- fluorophore
- chemical modification
- other functional sequences

Do not assume every method contains all components.

---

# 18. Sequence Classes

When appropriate, classify probe components as:

## fixed

An invariant sequence explicitly specified by the method.

## variable

Generated or selected according to the biological target.

## barcode_library

Selected from a predefined barcode or readout sequence library.

## sample_specific

Depends on the specific experiment or sample.

## version_specific

Depends on the method version.

## unresolved

Available evidence is insufficient to determine the sequence class confidently.

Do not classify a sequence as fixed simply because only one example is shown.

---

# 19. Exact Sequence vs Schematic Sequence

This distinction is mandatory.

An exact published nucleotide sequence must be distinguishable from:

- schematic sequence
- sequence example
- motif
- consensus
- degenerate sequence
- conceptual architecture
- inferred sequence
- reverse complement derived from another sequence

A figure showing:

```text
target arm — barcode — primer
```

does not provide an exact nucleotide sequence.

Never convert a schematic architecture diagram into a nucleotide sequence.

---

# 20. Probe Orientation

Always verify orientation before reporting probe architecture.

Check:

- 5'→3'
- target strand
- complementary strand
- reverse complement
- ligation orientation
- primer orientation

Do not assume that sequences printed in a supplementary table correspond directly to synthesized 5'→3' oligo orientation.

If orientation cannot be confidently established, mark it unresolved.

---

# 21. Probe-Design Rules

Actively investigate whether authors report rules for:

- target region selection
- target type
- probe length
- length range
- GC content
- melting temperature
- sequence complexity
- repeat masking
- off-target filtering
- genome alignment
- transcript alignment
- exon selection
- isoform handling
- transcript selection
- spacing
- overlap
- probes per target
- barcode assignment
- codebook design
- Hamming distance
- error correction
- readout assignment

Only record method-specific rules when supported by evidence.

Do not turn general oligo-design practice into a method-specific rule.

---

# 22. Barcode Terminology

The term `barcode` is method-dependent.

Do not assume it always means a physical nucleotide segment.

Distinguish when applicable:

- physical barcode sequence
- abstract codeword
- codebook entry
- imaging-round identity
- sequencing-round identity
- fluorescence-channel combination
- readout sequence
- readout probe
- barcode assignment

Use the original authors' terminology and separately record the standardized functional interpretation when helpful.

---

# 23. Original Terminology

Preserve the terminology used by the original publication.

For example, if authors use:

```text
homology arm
```

do not silently replace it everywhere with:

```text
target-binding region
```

Instead record both:

```text
Original term: homology arm
Standard functional category: target_binding_region
```

Consult `terminology.md` whenever terminology is ambiguous.

---

# 24. Probe-Design Software

For every method, explicitly investigate whether the original developers provide:

- official probe-design software
- command-line tools
- scripts
- notebooks
- web tools
- barcode-generation tools
- probe-selection tools
- oligo-generation pipelines

When available, record:

- software name
- official status
- repository
- version
- programming language
- inputs
- outputs
- design logic
- dependencies
- availability
- maintenance status

Do not classify third-party software as official software.

Do not state that no software exists unless sufficient evidence supports that conclusion.

---

# 25. Public Probe Sets

Actively search for official public probe resources.

Possible resources include:

- complete synthesized oligo sequences
- target-binding sequences
- barcode tables
- readout libraries
- primer lists
- probe design tables
- gene lists
- target coordinates
- supplementary spreadsheets
- repository sequence files

Do not describe all of these simply as:

```text
probe dataset available
```

The degree of completeness matters.

---

# 26. Benchmark Evaluation

Public probe resources must be evaluated separately from method descriptions.

Use the following grades.

## Grade A — Fully reproducible

Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.

## Grade B — Mostly reproducible

Most relevant information is available, but minor information remains missing.

## Grade C — Architecture-level only

Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.

## Grade D — Insufficient information

The design cannot be reliably reconstructed or compared.

A published probe sequence list is NOT automatically Grade A.

Benchmark assessment should consider:

- species
- genome build
- transcriptome/annotation version
- target genes
- target regions
- genomic/transcript coordinates
- target-binding sequences
- full synthesized probes
- barcode assignments
- readout sequences
- primer/linker sequences
- design parameters
- design software
- source provenance

Always explain the grade.

Never label incomplete probe information as strict ground truth.

---

# 27. Species and Reference Information

When relevant, identify:

- species
- genome assembly
- transcriptome reference
- annotation version
- cell line
- tissue
- sample type

Do not infer a genome build merely from the publication year.

If probe sequences are public but the reference genome or annotation is unknown, explicitly record that limitation.

---

# 28. sources.md

Every important technical claim must be traceable to a source.

Use stable source IDs:

```text
S01
S02
S03
```

For each source record when available:

- Source ID
- source type
- citation
- title
- authors
- year
- DOI
- PMID
- official URL
- repository/accession
- method version
- relevant figures
- relevant tables
- relevant supplementary files
- information supported
- limitations

A citation is not sufficient if the source does not actually support the associated claim.

---

# 29. Conflicting Sources

If official sources disagree:

1. preserve both values
2. identify the source for each value
3. identify the relevant method version
4. explain the discrepancy if possible
5. do not silently select one

Later sources do not automatically supersede original sources unless they explicitly describe a revised protocol or method version.

---

# 30. Research Manifest

The Research Agent must create:

```text
research_manifest.md
```

This file documents what was actually inspected.

For important sources record:

- Source ID
- source type
- location
- files/tables/figures examined
- information extracted
- important limitations
- unresolved questions

Include:

```text
## Unresolved Research Questions
```

This file supports independent review.

It is not part of normal downstream biological retrieval.

---

# 31. Research → Review → Fix Workflow

Method-specific corpus construction follows three independent stages.

## Stage 1 — Research

Research Agent:

- discovers sources
- extracts evidence
- creates initial production corpus
- creates research_manifest.md

The Research Agent must not claim that the corpus is validated.

## Stage 2 — Review

Review Agent:

- independently verifies important claims
- assumes the Research Agent may be wrong
- does not modify corpus content
- creates validation_report.md

## Stage 3 — Fix

Fix Agent:

- reads the corpus and validation report
- independently verifies review findings
- fixes supported errors
- rejects incorrect reviewer findings when appropriate
- creates fix_report.md
- performs final validation

No Agent should blindly trust the previous Agent.

---

# 32. Review Severity

Review findings use:

## ERROR

A factual, structural, provenance, versioning, architecture, or benchmark problem that must be corrected before validation.

## WARNING

Incomplete, ambiguous, weakly supported, or potentially misleading information.

## INFO

A non-critical improvement.

Use stable issue IDs:

```text
R001
R002
R003
```

Do not classify purely stylistic preferences as ERROR.

---

# 33. Fix Resolution Status

Each review issue must ultimately receive one of:

- FIXED
- REJECTED
- PARTIALLY FIXED
- UNRESOLVED

Definitions:

## FIXED

The reviewer was correct and the corpus was corrected.

## REJECTED

The reviewer was incorrect and primary evidence supports the existing corpus.

## PARTIALLY FIXED

The review identified a real issue, but the proposed resolution was incomplete or partly incorrect.

## UNRESOLVED

Available evidence is insufficient to resolve the issue confidently.

---

# 34. Final Validation Status

After Fix, assign one of:

## VALIDATED

No unresolved ERROR-level factual or structural problems remain.

## VALIDATED WITH WARNINGS

No ERROR-level issues remain, but meaningful uncertainty or incomplete evidence remains.

## NOT VALIDATED

One or more ERROR-level problems remain unresolved.

Do not mark a method validated simply because all files exist.

---

# 35. Corpus Status Vocabulary

Use the canonical method-level corpus status vocabulary defined in `corpus_schema.md` and `methods_index.yaml`.

Recommended lifecycle:

```text
not_researched
research_complete
review_failed
review_passed
fixing
validated
validated_with_warnings
```

Do not invent synonyms without updating the global schema.

---

# 36. Cross-File Consistency

Whenever a factual value changes, check all production files for dependent references.

For example, changing a barcode length may require updates to:

- metadata.yaml
- probes.md
- benchmarks.md
- description.md

The same fact must not have contradictory values in different files.

Before finishing any task, inspect for cross-file contradictions.

---

# 37. Template Compliance

Use files under:

```text
templates/
```

as the structural basis for corpus files.

Do not remove required sections simply because information is unavailable.

Instead use an appropriate explicit missing-value statement where necessary.

Sections that genuinely do not apply may be marked:

```text
Not applicable.
```

Do not confuse this with:

```text
Not reported.
```

---

# 38. Retrieval Design

The corpus should support selective retrieval.

Recommended routing:

## General method question

Read:

- index.md
- description.md

## Experimental protocol question

Read:

- protocols.md
- sources.md

## Probe architecture/design question

Read:

- metadata.yaml
- probes.md
- sources.md

## Benchmark question

Read:

- benchmarks.md
- probes.md
- sources.md

## Version question

Read:

- metadata.yaml
- relevant version sections/files
- sources.md

Do not design corpus files in a way that requires loading the entire repository for every question.

---

# 39. Production vs Audit Separation

Do not move research notes or reviewer commentary into production corpus files unless they represent validated knowledge.

Production corpus files should contain finalized knowledge.

Audit files should contain:

- search history
- uncertainty discussions
- review findings
- rejected interpretations
- correction history

This separation is intentional.

---

# 40. Prohibited Behaviors

Agents must NOT:

- invent nucleotide sequences
- reconstruct exact sequences without explicit evidence
- invent protocol parameters
- infer exact values from approximate figures
- silently normalize method names
- silently merge method versions
- silently resolve source conflicts
- present inferred claims as reported
- classify example sequences as universal without evidence
- classify third-party tools as official
- assume absence of evidence means false
- classify incomplete probe sets as Grade A ground truth
- use reviews as primary evidence when original sources are available
- populate unknown YAML fields merely to make the file look complete
- add speculative statements for completeness

---

# 41. Preflight Repository QA

Before beginning large-scale method research, repository-level QA should verify:

- directory names
- official method names
- methods_index.yaml routing
- YAML validity
- template completeness
- schema/template consistency
- terminology consistency
- evidence-status vocabulary
- benchmark-grade vocabulary
- method-status vocabulary
- missing-value conventions
- version rules
- retrieval design

If repository-level blocking errors remain, method-specific research should not proceed until they are resolved.

---

# 42. Final Checklist for Research Tasks

Before a Research Agent finishes a method, verify:

- official method name has been checked
- original publication is identified
- primary and supplementary sources were searched
- method versions are not mixed
- probe architecture is explicitly represented where possible
- probe orientation has been checked
- probe components are individually classified
- original terminology is preserved
- exact sequences are separated from schematic representations
- target-binding design rules were investigated
- barcode/readout logic was investigated
- probe-design software was investigated
- public probe sets were investigated
- species/reference information was investigated
- benchmark candidates were assessed
- benchmark grade is justified
- major claims are traceable to sources
- unresolved information is explicitly recorded
- metadata agrees with Markdown files
- research_manifest.md is complete

---

# 43. Final Checklist for Review Tasks

Before the Review Agent finishes:

- important claims were independently checked
- source citations were not automatically trusted
- probe architecture received detailed review
- 5'→3' orientation was checked
- version mixing was checked
- exact vs schematic sequences were checked
- reported vs inferred information was checked
- design software was checked
- public datasets were checked
- benchmark grades were independently reassessed
- cross-file contradictions were checked
- every issue has a stable issue ID
- no corpus file was modified
- PASS/FAIL reflects ERROR-level findings

---

# 44. Final Checklist for Fix Tasks

Before the Fix Agent finishes:

- every review issue was processed
- reviewer claims were independently verified
- corrections are evidence-supported
- incorrect reviewer findings were rejected when appropriate
- dependent files were updated consistently
- no new unsupported claims were introduced
- unresolved issues remain explicit
- fix_report.md records every resolution
- final validation was performed
- final status is justified

---

# 45. Current Corpus Scope

The initial corpus targets:

- STARmap
- ExSeq
- seqFISH+
- MERFISH
- Oligo-LiveFISH
- RAEFISH
- TDDN-FISH
- PRISM
- OpenFISH
- MiP-seq

This list may be expanded later.

Adding a new method must follow the same schema, evidence, versioning, and validation rules.

---

# 46. Governing Principle

When uncertain, prefer:

```text
accurate but incomplete
```

over:

```text
complete but inferred
```

The corpus is intended to function as a trusted technical knowledge base.

Every value that appears precise should be defensible from evidence.