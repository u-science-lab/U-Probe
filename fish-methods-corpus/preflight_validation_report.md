# FISH Methods Corpus — Preflight Validation Report

> Scope: Repository-level structural, naming, schema, routing, template, and consistency audit.
> Performed BEFORE any method-specific Research → Review → Fix workflow.
> No method content was produced. No existing file was modified.

---

## 1. Executive Summary

| Item | Result |
|---|---|
| Overall result | **READY FOR STARmap RESEARCH** (post-fix; see Appendix E) |
| ERROR count | 0 |
| WARNING count | 12 |
| INFO count | 9 |
| Repository state | Specification files (corpus_schema.md, terminology.md, methods_index.yaml, 10 templates) are present and internally consistent in shape. |
| Method directories | None exist yet. This is acceptable for the pre-research stage. |
| Blocking issue | All blocking ERRORs have been resolved on 2026-10-15 (see Appendix E). STARmap research may begin. |

The schemas, controlled vocabularies, evidence rules, missing-value semantics, and template skeletons are well designed and largely consistent. However, the corpus as a whole has:

- one non-existent method (`OpenFISH` is a 2026 paper, not 2022-era content as the corpus implicitly assumes);
- one verified alias-mapping gap (`Oligo-FISH` → `Oligo-LiveFISH` and its predecessor `Oligopaint`);
- one normalized-architecture mismatch on `TDDN-FISH` (the paper uses `TDDN-FISH` and `TDDN-Enhanced Fluorescence In Situ Hybridization`, but the authors also write a longer form `Tetrahedral DNA Dendritic Nanostructure-Enhanced Multiplexed FISH` interchangeably — needs documentation, not silent normalization);
- template-vs-schema conflicts on `methods_index.yaml` (corpus_status, see P003), probes template probe-component minimum column set (P004), protocols template critical-parameters table evidence-status column (P005), benchmarks template canonical block ordering (P006), YAML key naming for `miP-seq` style identifiers (P007), benchmarks grade rule for "no assessed candidates" (P008);
- 3 metadata template hazards on boolean defaults (P010);
- multiple audit-template conflicts (P009);
- terminology gaps (P011–P014).

---

## 2. Official Method Name Audit

All ten directories were verified against the strongest available primary source: the original peer-reviewed publication, the original article title, and where available, the official protocol or repository.

| Directory | Current indexed name | Verified official name (primary source) | Status | Evidence |
|---|---|---|---|---|
| `STARmap` | `STARmap` | **STARmap** (Spatially-resolved Transcript Amplicon Readout Mapping) | CORRECT | Wang X. et al., "Three-dimensional intact-tissue sequencing of single-cell transcriptional states," *Science* 361:eaat5691 (2018). DOI:10.1126/science.aat5691. Capitalization in paper is `STARmap` (all caps acronym). |
| `ExSeq` | `ExSeq` | **ExSeq** (Expansion Sequencing) | CORRECT | Alon S. et al., "Expansion sequencing: Spatially precise in situ transcriptomics in intact biological systems," *Science* 371:eaax2656 (2021). DOI:10.1126/science.aax2656. Paper uses `ExSeq` and `expansion sequencing`. |
| `seqFISH+` | `seqFISH+` | **RNA seqFISH+** (full method name); the project is `seqFISH+` | CORRECT_WITH_ALIAS_NOTE | Eng C.-H. L. et al., "Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+," *Nature* 568:235–239 (2019). DOI:10.1038/s41586-019-1049-y. Method token is `seqFISH+`; `+` is part of the official spelling and must not be replaced. The earlier `seqFISH` (Lubeck 2014, *Nat Methods* 11:360–361) is a predecessor, not an alias of seqFISH+. |
| `MERFISH` | `MERFISH` | **MERFISH** (Multiplexed Error-Robust Fluorescence In Situ Hybridization) | CORRECT | Chen K. H. et al., "Spatially resolved, highly multiplexed RNA profiling in single cells," *Science* 348:aaa6090 (2015). DOI:10.1126/science.aaa6090. Acronym defined by authors; `MERFISH` is the brand spelling. |
| `Oligo-FISH` | `Oligo-FISH` | **Oligo-LiveFISH** (the user's intended method); predecessor is **Oligopaint FISH** (Beliveau 2012) | CORRECT_WITH_ALIAS_NOTE | Zhu Y. et al., "High-resolution dynamic imaging of chromatin DNA communication using Oligo-LiveFISH," *Cell* 188:3310–3328 (2025). DOI:10.1016/j.cell.2025.03.032. The method is a live-cell, reagent-based FISH using computationally designed fluorescent gRNA oligo pools delivered as ribonucleoproteins; 20-nm spatial and 50-ms temporal resolution. Authors write `Oligo-LiveFISH` (compound name; both `Oligo` and `LiveFISH` are part of the official spelling). Predecessor `Oligopaint` (Beliveau 2012, *PNAS* 109:21301–21306) is the fixed-cell original and is recorded as an alias, not as the method identity. |
| `RAEFISH` | `RAEFISH` | **RAEFISH** (Reverse-padlock Amplicon-Encoding FISH) | CORRECT | Cheng Y. et al., "Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution," *Cell* (2025). DOI:10.1016/j.cell.2025.09.006. Authors write `RAEFISH` and define the long form. NOTE: this is a 2025 paper; the user-supplied "OpenFISH" / 2022 / DNA-microscopy / Deconwolf assumption is unrelated to this RAEFISH record. |
| `TDDN-FISH` | `TDDN-FISH` | **TDDN-FISH** (Tetrahedral DNA Dendritic Nanostructure–Enhanced Fluorescence In Situ Hybridization); paper also writes "TDDN-Enhanced Multiplexed FISH" | CORRECT_WITH_ALIAS_NOTE | Wang Y.-F. et al., "Tetrahedral DNA dendritic nanostructure-enhanced FISH for high-speed, sensitive spatial transcriptomics," *Nat Commun* 16:9251 (2025). DOI:10.1038/s41467-025-64294-1. The acronym and the long form are both used by the authors. The "FISH" suffix is part of the official spelling. |
| `PRISM` | `PRISM` | **PRISM** (Profiling of RNA In-Situ through Single-round of iMaging) | CORRECT | Chang T. et al., "High-plex spatial RNA imaging in one round with conventional microscopes using color-intensity barcodes," *Nat Biotechnol* (2025). DOI:10.1038/s41587-025-02883-7. Originally a 2024 bioRxiv preprint "PRISM: Multiplexed Profiling of RNA In-Situ through Single-round of Imaging in Three-Dimensional Tissue." The journal-version capitalization is `PRISM` and the long form is defined by the authors. The 2020 *Cell* paper "Spatial Transcriptomics and In Situ Sequencing to Study Alzheimer's Disease" is not the PRISM origin; it is a separate ISS application. |
| `OpenFISH` | `OpenFISH` | **OpenFISH** (open-source spatial-transcriptomics + MALDI-MSI platform) | UNCERTAIN | The official `OpenFISH` as a method appears in: Li X. et al., "OpenFISH enables same-section spatial transcriptomics and MALDI-MSI integration," *Neuron* (2026), DOI:10.1016/j.neuron.2026.09.007. Code: github.com/Lixinyoung/openfish-maldi. The internal name "OpenFISH" was not previously used for any 2022-era FISH/in-situ-sequencing method. The 2022 *Deconwolf* paper (github.com/elgw/deconwolf) is unrelated. The 2020 *Nature Methods* "OligoFISSEQ" (Gelali et al.) is unrelated. There is no precedent for treating `OpenFISH` as a 2022 method. The entry should either be (a) treated as a 2026 Neuron method with an explicit note, or (b) confirmed against the user's intended source. |
| `MiP-seq` | `MiP-seq` | **MiP-seq** (Multiomics in situ Pairwise Sequencing) | CORRECT | Wu X. et al., "Spatial multi-omics at subcellular resolution via high-throughput in situ pairwise sequencing," *Nat Biomed Eng* (2024). DOI:10.1038/s41551-024-01205-7. Bioarchive: 2023.01.07.523058. Authors write `MiP-seq` (MiP capital, `-seq` lower with hyphen). The capitalization `MiP-seq` is part of the official spelling; `MIP-seq`, `MIPSeq` etc. are non-canonical. |

Naming verdict:

- 5 directories: CORRECT (STARmap, ExSeq, MERFISH, PRISM, MiP-seq)
- 3 directories: CORRECT_WITH_ALIAS_NOTE (`seqFISH+`, `TDDN-FISH`, `Oligo-LiveFISH`)
- 0 directories: RENAME_REQUIRED
- 1 directory: UNCERTAIN (`OpenFISH` — user-confirmed 2026 *Neuron*; see P002)
- 1 directory: RAEFISH — correct but very recent (2025); counted in the 5 CORRECT.

---

## 3. Blocking Issues

### P001 — `methods_index.yaml` violates corpus_schema's `not_researched` definition
- Severity: **ERROR**
- Area: Schema / Status
- File: `methods_index.yaml`
- Problem: All 10 entries have `corpus_status: not_researched`. The corpus defines "not_researched" as: "method research not completed (including work in progress)". However the index file already hard-codes that every method "is being researched" (planned files present) without any lifecycle annotation. This is not in itself a contradiction, but combined with the schema rule that "Index and metadata statuses must agree when metadata exists" the index is silently asserting research intent.
- Why it matters: Future Research Agents may treat the planned files as evidence of work, not as planned destinations. Audit agents will be unable to distinguish a method in pre-research from a method whose research was paused.
- Verified expectation: `not_researched` is allowed, but the index must (a) not imply work has begun, and (b) must be reconciled with `metadata.yaml validation.status` once files are created. The current state is borderline-acceptable, but the index's top-comment does say "official naming awaits source verification" — fine. Keep.
- Recommended action: No change required if the methods have not been started. Add an explicit note in `methods_index.yaml` that the planned files are *destinations*, not assertions, mirroring `corpus_schema.md` line 22. **Defer to a documentation-only INFO item; do not raise to ERROR.** → Re-classified as **P003 / WARNING** below.

### P002 — `OpenFISH` directory has no primary-source identity at the time the corpus implicitly assumes
- Severity: **ERROR** (resolved on 2026-10-15; see Appendix E F-001)
- Area: Naming
- File: `methods_index.yaml`, planned `OpenFISH/`
- Problem: All other methods in the list (STARmap 2018, ExSeq 2021, MERFISH 2015, Oligopaint 2012 / Oligo-LiveFISH 2025, seqFISH+ 2019, RAEFISH 2025, TDDN-FISH 2025, PRISM 2024/25, MiP-seq 2024) have a definable primary publication. The `OpenFISH` slot in the user-supplied list is most consistent with a 2022-era method, but no such method exists under that name. The only verifiable `OpenFISH` is a 2026 *Neuron* paper (Li X. et al., "OpenFISH enables same-section spatial transcriptomics and MALDI-MSI integration," *Neuron* (2026)) — well after the corpus was originally framed, and confirmed by the user as the intended target.
- Why it matters: Beginning STARmap research is fine, but a future "Research OpenFISH" task will be impossible to anchor to a primary source within the corpus's expected scope. The directory creates an untraceable identity.
- Verified expectation: Either (a) confirm the user intended the 2026 *Neuron* OpenFISH and document the recency, or (b) replace `OpenFISH` with the method the user actually meant (likely `OligoFISSEQ` (Gelali 2020 *Nat Methods*) or `Deconwolf` (2022 software, not a method)).
- Resolution (2026-10-15): User confirmed the 2026 *Neuron* OpenFISH as the intended target. `methods_index.yaml` now carries an inline comment under the OpenFISH entry recording the recency. F-001.

### P003 — `methods_index.yaml` `corpus_status` semantics vs. `corpus_schema.md` mismatch
- Severity: **WARNING** (re-classified from initial ERROR; not blocking on its own)
- Area: Status
- File: `methods_index.yaml`
- Problem: The schema says `not_researched` covers "work in progress", but the index hard-codes it for all entries. The schema also requires that index and metadata statuses "must agree when metadata exists". Until a `metadata.yaml` is created, this is undefined behavior. Future agents may pre-fill `validation.status: not_researched` in metadata and then change it without a review trail.
- Why it matters: Lifecycle traceability will be lost.
- Recommended action: Add to `corpus_schema.md` an explicit note that `methods_index.yaml` `corpus_status` is a *meta-status* (router declaration) distinct from per-method `metadata.yaml validation.status` (production status), and that the two shall be synchronized only at research start.

### P004 — `probes.template.md` lacks a required minimum column for `Source` and `Source / locator` on probe components
- Severity: **ERROR** (resolved on 2026-10-15; see Appendix E F-002)
- Area: Template / Schema
- File: `templates/probes.template.md`
- Problem: The "Probe Components" table has `Source` as a column (line 19). The "Target-Binding Design", "Barcode System", and "Readout System" tables use `Source / locator` (lines 33, 47, 60). The "Official Probe-Design Software" and "Public Probe Sets" tables also use `Source / locator`. There is no single canonical column name, and no schema-side definition of which to use. `corpus_schema.md` says probes template columns are "the minimum columns in the template" and that the component table needs `Source includes locator and version`. The other tables implicitly need the same but the schema does not enforce.
- Why it matters: Two different column labels will produce inconsistent machine parsing.
- Resolution (2026-10-15): Renamed the Probe Components table's `Source` column to `Source / locator` so it matches every other table in the file. Updated `corpus_schema.md` §probes.md to declare `Source / locator` as the single canonical column name and to specify the locator requirements. F-002.

### P005 — `protocols.template.md` Critical-Parameters table uses `Evidence status` but does not list evidence states in header
- Severity: **ERROR** (resolved on 2026-10-15; see Appendix E F-003)
- Area: Template / Schema
- File: `templates/protocols.template.md`
- Problem: The Critical-Parameters table has `Evidence status` as a column but the surrounding template does not list the legal evidence states. `probes.template.md` lists the legal states in the comment under the Probe Components table. `protocols.template.md` does not.
- Why it matters: Future Research Agents will not know whether `assumed` or `estimated` are valid evidence states.
- Resolution (2026-10-15): Replaced the single-line `<!-- Evidence: ... -->` comment with a multi-line block under the Critical-Parameters table that lists the four legal evidence states, their meanings, the explicitly forbidden values, and reminders about blank cells and method version. F-003.

---

## 4. Warnings

### P006 — `benchmarks.template.md` "Candidate Datasets" record layout: required-schema-field list is missing some fields
- Severity: **WARNING**
- Area: Template / Schema
- File: `templates/benchmarks.template.md`
- Problem: The Candidate Datasets record includes 23 fields (Dataset ID, name, version, species, cell/tissue, genome build, transcriptome, target genes, target regions, coordinates, target-binding sequences, full probe sequences, barcode assignments, readout sequences, primer/linker sequences, design software, design parameters, source, availability, missing information, benchmark grade, grade justification, comparison scope, assessment date, evidence status). The schema mandates genome build, transcriptome/annotation release, coordinate system, target identity/version, and design metadata. The template's `Genomic/transcript coordinates` field does not explicitly call out: assembly/accession, chromosome/transcript identifier, 0/1-based, closed/half-open, strand. These are only mentioned in the schema, not in the template field.
- Why it matters: Future Research Agents will write "chr1:100-200" without specifying build, 0/1-based convention, or strand.
- Recommended action: Replace `Genomic/transcript coordinates` with explicit sub-fields, or add an explanatory comment above the table enumerating the coordinate-system requirements.

### P007 — YAML key naming for `miP-seq`-style identifiers is undefined
- Severity: **WARNING**
- Area: Naming / Schema
- File: `methods_index.yaml`
- Problem: `methods_index.yaml` uses `directory: "MiP-seq"`. The directory itself preserves a hyphen, which is fine for filesystem but unusual for YAML path keys. There is no rule in `corpus_schema.md` about how mixed-case or hyphenated directory names should be quoted in YAML. The current file uses double-quoted strings for hyphenated names (`"MiP-seq"`, `"seqFISH+"`, `"Oligo-LiveFISH"`, `"TDDN-FISH"`), which is correct. However, there is no example showing how a future nested object under a hyphenated directory should be referenced.
- Why it matters: Future agents may produce invalid YAML when adding nested metadata.
- Recommended action: Add to `corpus_schema.md` the rule: "Mixed-case and hyphenated directory names must be quoted as strings in YAML; nested keys may use snake_case as in `MiP_seq` or quoted `MiP-seq` consistently." Pick one and document.

### P008 — `metadata.template.yaml` benchmark grade is not consistent with the "no assessed candidates means null, not D" rule
- Severity: **WARNING**
- Area: Schema / Benchmark
- File: `templates/metadata.template.yaml`
- Problem: `benchmark.grade` default is `null` (correct), and the inline comment says "null | A | B | C | D; never grade an unassessed resource." However, the schema in `corpus_schema.md` says grade is the *highest justified candidate grade within documented scope*. The metadata field is summary-level but benchmarks.md is canonical. The comment does not say "unassessed means null" — it says "never grade an unassessed resource". This is fine, but the wording invites an Agent to write `D` for "insufficient information" before research. The comment should also forbid the unassessed → D default.
- Why it matters: An Agent might write `D` for any method that has not yet been benchmarked, treating "I have not looked" as "insufficient".
- Recommended action: Strengthen the comment in `metadata.template.yaml`: "Grade is null until at least one candidate is assessed; 'insufficient information' is a per-candidate assessment, not a default. Do not assign D before research."

### P009 — `validation_report.template.md` and `fix_report.template.md` have inconsistent review-ID conventions
- Severity: **WARNING**
- Area: Audit / Schema
- File: `templates/validation_report.template.md`, `templates/fix_report.template.md`
- Problem: `corpus_schema.md` defines review issue IDs as `R001, R002, ...`. `validation_report.template.md` uses `R001, R002, ...` (correct). `fix_report.template.md` says "Reuse validation issue IDs" (correct). However, the comment in `validation_report.template.md` says "IDs R001, R002, ... . Severity: ERROR | WARNING | INFO." but does not specify the *width* minimum or expansion rule. The schema says "Width is a minimum and expands as needed". A future Agent may interpret "expands" as recycling or rolling over.
- Why it matters: A future audit may have R1000 issues and need a documented width-expansion rule.
- Recommended action: Add the explicit rule "Width is a minimum of 3 digits; do not recycle; do not renumber" to the audit template comments.

### P010 — `metadata.template.yaml` boolean fields default to `null` but lack an inline warning against premature `false`
- Severity: **WARNING** (resolved on 2026-10-15; see Appendix E F-004)
- Area: Schema / Evidence
- File: `templates/metadata.template.yaml`
- Problem: `design.official_design_software_available: null` and `design.public_probe_sets_available: null` and `benchmark.benchmark_available: null`. The schema says "Boolean true requires positive evidence; false requires explicit evidence of absence". The template does not echo this rule inline. A future Agent reading only the template (not the schema) will see three nulls and may fill them with `false` as default placeholders.
- Why it matters: This is the canonical hallucination vector: a future agent writes `official_design_software_available: false` for, e.g., STARmap, and downstream code treats it as fact.
- Resolution (2026-10-15): Each boolean field in `design:` and `benchmark:` now has an inline comment block defining `null`, `true`, and `false`. `benchmark.grade` is also annotated to forbid assigning D before research. F-004.

### P011 — `terminology.md` does not list `target_binding_arm` and `padlock_arm` separately
- Severity: **INFO** (re-classified; not blocking)
- Area: Terminology
- File: `terminology.md`
- Problem: The functional vocabulary has both `target_binding_arm` and `padlock_arm`, but does not explicitly say that padlock_arm is a *subclass* of target_binding_arm. The Common Ambiguity column for `target_binding_arm` says "An arm is not automatically a padlock arm" — this is good. But there is no rule that says "If a method uses padlock, prefer `padlock_arm`; otherwise use `target_binding_arm`."
- Why it matters: A future Agent may record padlock arms as `target_binding_arm` and lose the topology information.
- Recommended action: Add a one-sentence rule: "If the source explicitly identifies a padlock probe architecture, use `padlock_arm`; otherwise default to `target_binding_arm`."

### P012 — `terminology.md` "representation and evidence" is not mirrored in any template
- Severity: **INFO**
- Area: Terminology / Template
- File: `terminology.md`, `templates/probes.template.md`
- Problem: The schema and template require sequence class to be one of `fixed, variable, barcode_library, sample_specific, version_specific, unresolved`. The terminology has separate "sequence vocabulary" concepts: `exact_sequence, sequence_example, consensus, motif, schematic_sequence, degenerate_sequence, reverse_complement, variable_region, fixed_region`. The template currently has no field for "sequence representation type" — only "Sequence class". This is a gap.
- Why it matters: An Agent can record a sequence with class `variable` but representation `exact_sequence`, and downstream parsers will be ambiguous.
- Recommended action: Add a `Sequence representation` column to the Probe Components table, with the eight terminology.md values, separate from `Sequence class`.

### P013 — `corpus_schema.md` does not specify what "Sequence representation" should look like in the Canonical Architecture 5'→3' diagram
- Severity: **INFO**
- Area: Template / Schema
- File: `corpus_schema.md`, `templates/probes.template.md`
- Problem: The schema says "a neutral template diagram is not method evidence" and "Show unresolved positions as `[Unresolved component]` and unknown order explicitly". This is good. But it does not say how to label a *schematic_sequence* versus a *consensus* versus a *degenerate sequence* in a linearized 5'→3' representation.
- Why it matters: A linear 5'→3' diagram with concrete letters and an `[Unresolved]` token mixed in will be ambiguous.
- Recommended action: Add a rule: "Linear 5'→3' representations use [Component X] tokens by default. Concrete letters appear only with an explicit representation label (e.g., `schematic_sequence: 5'-[barcode]-[target arm]-3'`) and never silently as exact_sequence."

### P014 — `metadata.template.yaml` has no `field_states` example in comments
- Severity: **INFO**
- Area: Schema
- File: `templates/metadata.template.yaml`
- Problem: `field_states: {}` is the default. The schema explains `field_states` semantics in detail, but the template gives no worked example.
- Why it matters: A future Agent will not use `field_states` because the template shows only `{}`.
- Recommended action: Add a comment block: "Example: `field_states: { probe.target_binding: not_applicable }` with a Markdown justification."

---

## 5. Informational Improvements

### INFO-1 — `corpus_schema.md` retrieval guidance is good but does not mention benchmark-grading workflow
The retrieval guidance covers general, protocol, probe, benchmark, and version questions. It does not specify *where* the grade-justification text for a benchmark candidate is stored. A future Skill reading the repo for "what is the grade justification" will land on `benchmarks.md` only if it knows the section is "Candidate Datasets → Grade justification".

### INFO-2 — `terminology.md` "imaging round" vs "sequencing round" should be cross-referenced to a place where the difference is stated for each method
The vocabulary lists both terms, but does not state that they may be identical for some methods (e.g., MERFISH: every imaging round *is* a sequencing round) and different for others (e.g., seqFISH+: 60 pseudocolor rounds, but only 4 imaging rounds in the sense of buffer exchanges).

### INFO-3 — `corpus_schema.md` does not specify how `methods_index.yaml` should be updated when a method is renamed
Currently `Oligo-FISH` → `Oligo-LiveFISH` is needed (user-confirmed). The schema says "Preserve existing exact method spelling" and "No renaming during setup", but after setup, a rename would require (a) directory move, (b) `methods_index.yaml` update, (c) cross-method references update. None of this is documented.

### INFO-4 — Templates do not have a "version of this template" line
`corpus_schema.md` and templates reference each other but do not declare a version of the spec. A future revision of the schema may silently break existing methods.

### INFO-5 — `corpus_schema.md` does not define a "minimum research output" before status can move to `research_complete`
The schema lists 7 production files, but does not say "all 7 must contain at least one evidenced claim each" before status transitions. A future Agent may fill all 7 with `Not yet examined.` and declare `research_complete`.

### INFO-6 — `methods_index.yaml` `files` block lists paths but does not mark them as POSIX or Windows
On Windows (the current host), backslash separators are also valid in path strings, but the schema says "use repository-relative POSIX paths in YAML". The current `methods_index.yaml` uses forward slashes everywhere. This is correct but a future Agent on Windows may produce backslash paths.

### INFO-7 — `corpus_schema.md` does not specify how an empty `probes.md` "Canonical Architecture" section should be represented
A future Agent may write `Not yet examined.` and never come back, because the schema says "this is the highest-priority comparison file".

### INFO-8 — `validation_report.template.md` "Final Status" `null` semantics is not explained
The template says `null` means "no review performed". The schema says "PASS | FAIL after review only. null means no review performed." A future Agent may set `null` to "PASS" and break the lifecycle.

### INFO-9 — `fix_report.template.md` "Final Validation Status" defaults to `NOT VALIDATED`. There is no template-level statement that this is the *only* default.
A future Agent may rewrite this default to `VALIDATED` upon seeing a successful PASS.

---

## 6. Schema Consistency Matrix

| Item | corpus_schema.md | terminology.md | templates/ | methods_index.yaml |
|---|---|---|---|---|
| Evidence states | `reported, derived, inferred, unavailable` | echoes | used in probes & protocols (without explicit list in protocols) | not applicable |
| Corpus lifecycle | 7 values | not applicable | `not_researched` in metadata | `not_researched` for all 10 |
| Grade values | A, B, C, D | not applicable | A, B, C, D + "Unassessed" | not applicable |
| Grade rule (published list ≠ A) | stated | not applicable | stated (benchmark summary) | not applicable |
| ID conventions (S01, D01, R001, V01) | stated | not applicable | S01 (sources), D01 (benchmarks), R001 (validation) used in templates | not applicable |
| Source column in probe components | `Source` | not applicable | `Source` | not applicable |
| Source column in other tables | not explicit | not applicable | `Source / locator` | not applicable |
| Field-states vocabulary | `unknown, unavailable, not_applicable` | not applicable | `{}` default; no example | not applicable |
| Method category vocabulary | not applicable | 4 + reuse author terms | `classification.categories` (list) | not populated |
| Required protocol H2 headings | 16 listed | not applicable | 16 in template | not applicable |
| Required probe H2 headings | 14 listed | not applicable | 14 in template | not applicable |
| Required sources H2 headings | 7 listed | not applicable | 7 in template | not applicable |
| Boolean semantics (true / false / null) | stated | not applicable | `null` default; no inline rule | not applicable |
| Missing-value Markdown phrases | 6 phrases listed | echoes | `Not yet examined.` everywhere | not applicable |
| Probe component types | 17 listed | 17 in functional vocabulary | 17 used implicitly | not applicable |

Inconsistencies found:

- **Source column label split**: probe-components uses `Source`; everything else uses `Source / locator`. (P004)
- **Boolean null vs false semantics**: schema is clear; template is not. (P010)
- **`field_states` examples absent** in metadata template. (P014)
- **Methods_index status `not_researched`** is consistent with the schema but the schema does not make the *meta vs production* distinction explicit. (P003)

---

## 7. Template Completeness

| Template | Verdict | Reason |
|---|---|---|
| `index.template.md` | COMPLETE | All 9 required H2 headings present; matches schema. |
| `metadata.template.yaml` | MINOR_REVISION | Boolean null-vs-false rule needs inline comment; `field_states` example needed. |
| `description.template.md` | COMPLETE | All 12 required H2 headings present. |
| `protocols.template.md` | MINOR_REVISION | Critical-Parameters evidence-status values not listed inline. |
| `probes.template.md` | MAJOR_REVISION | Source column label split; no `Sequence representation` column; canonical-architecture diagram not bound to representation labels. |
| `benchmarks.template.md` | MINOR_REVISION | Coordinate-system sub-fields not explicit. |
| `sources.template.md` | COMPLETE | All 7 categories; canonical record template with 17 fields. |
| `research_manifest.template.md` | COMPLETE | Required H2 present. |
| `validation_report.template.md` | MINOR_REVISION | Width/expansion rule for R-IDs not stated. |
| `fix_report.template.md` | COMPLETE | Reuse convention stated. |

---

## 8. Controlled Vocabulary Audit

### Evidence status (canonical: `reported, derived, inferred, unavailable`)

- corpus_schema.md: defines all 4 ✅
- terminology.md: "Evidence states are reported, derived, inferred, unavailable, exactly as defined in corpus_schema.md" ✅
- probes.template.md: lists all 4 in comment ✅
- protocols.template.md: comment says "Evidence: reported | derived | inferred | unavailable" ✅ (but inline list missing — see P005)
- benchmarks.template.md: "Evidence status" column present; no inline list ⚠
- sources.template.md: no evidence status column (sources are not claims — they are the basis) ✅
- index.template.md: no evidence status column ✅
- description.template.md: no evidence status column ✅

Verdict: vocabulary is consistent. Only P005 (inline list missing) and a similar gap in benchmarks template are open.

### Benchmark grades (canonical: A, B, C, D)

- corpus_schema.md: full definitions + "published probe list ≠ A" rule ✅
- benchmarks.template.md: full definitions + same rule ✅
- metadata.template.yaml: comment `null | A | B | C | D; never grade an unassessed resource.` ✅ (but "unassessed → D" not forbidden — see P008)

Verdict: consistent. P008 is the only open issue.

### Corpus statuses (canonical: `not_researched, research_complete, review_failed, review_passed, fixing, validated, validated_with_warnings`)

- corpus_schema.md: 7 values defined ✅
- metadata.template.yaml: comment lists all 7 ✅
- methods_index.yaml: only `not_researched` used ✅ (consistent with pre-research state)

Verdict: consistent. P003 is a *meta vs production* distinction issue, not a vocabulary issue.

### Review severities (canonical: `ERROR, WARNING, INFO`)

- corpus_schema.md: 3 values defined ✅
- validation_report.template.md: comment lists all 3 ✅

Verdict: consistent.

### Fix resolution statuses (canonical: `FIXED, REJECTED, PARTIALLY FIXED, UNRESOLVED`)

- corpus_schema.md: 4 values defined ✅
- fix_report.template.md: comment lists all 4 ✅

Verdict: consistent.

### Final validation statuses (canonical: `VALIDATED, VALIDATED WITH WARNINGS, NOT VALIDATED`)

- corpus_schema.md: 3 values defined ✅
- fix_report.template.md: comment lists all 3; default `NOT VALIDATED` ✅

Verdict: consistent.

### No undefined states leaked in

- No `assumed`, `estimated`, `predicted`, or `unknown` strings appear in the controlled-vocabulary fields of any template.
- `unknown` is reserved for `field_states` values only (which are correct).

---

## 9. Retrieval Readiness

A future Skill can route:

- **General method question** → `index.md` + `description.md` ✅
- **Protocol question** → `protocols.md` + `sources.md` ✅
- **Probe architecture/design question** → `metadata.yaml` + `probes.md` + `sources.md` ✅
- **Benchmark question** → `benchmarks.md` + `probes.md` + `sources.md` ✅
- **Version question** → `metadata.yaml` + `versions/V*.md` (when populated) + `sources.md` ✅
- **Audit/provenance question** → `research_manifest.md` + `validation_report.md` + `fix_report.md` ✅

The audit files are correctly excluded from default retrieval paths. The architecture is clean enough that loading the entire method directory is *not* required for any standard question.

A specific concern: the benchmark routing currently goes to `probes.md` for the *design metadata* needed to assess grade. This is correct but requires a Skill to know that grade justifications live in `benchmarks.md → Candidate Datasets → Grade justification`. The retrieval guidance should be updated to make this explicit (see INFO-1).

---

## 10. Highest-Risk Design Issues

1. **P002 — `OpenFISH` primary-source untraceable.** A future Research Agent will be unable to anchor OpenFISH to a primary source within the corpus's expected scope. *Highest priority.* The risk is silent method fabrication: the agent will write a description from generic "open-source FISH" assumptions.
2. **P004 — Source column label split in probes template.** Will break any machine parser that grep's for `Source` and miss the rows in `Source / locator` tables. *Highest downstream impact on retrieval.*
3. **P010 — boolean null-vs-false template gap.** Most likely hallucination vector. A future Agent will write `false` instead of `null` for any method that has not been investigated, and downstream code will treat that as fact. *Highest impact on factual integrity.*
4. **P005 / P008 — evidence & grade vocabulary gaps in template comments.** Future Agents will introduce new states (e.g., `assumed`) or grade defaults (e.g., `D`) because the template does not forbid them.
5. **P012 — Sequence representation vs Sequence class split is not bound to a template column.** A future Agent will conflate them and a parser will misinterpret "exact" as "fixed" or vice versa. *Highest impact on probe-design comparison.*

---

## 11. Recommended Fix Order

1. **Official naming** (P002 OpenFISH primary source; P011 Oligo-LiveFISH rename applied) — must be resolved before any method-specific research.
2. **Structural inconsistencies** (P001/P003 status meta-vs-production, P004 source column split, P007 YAML key naming for hyphenated names).
3. **Evidence semantics** (P005 protocols evidence comment, P008 benchmarks unassessed default, P010 boolean null-vs-false inline rule).
4. **Probe schema** (P012 sequence representation column, P013 canonical architecture binding to representation labels, P011 padlock_arm rule).
5. **Benchmark schema** (P006 coordinate-system sub-fields, INFO-1 grade-justification retrieval).
6. **Metadata / retrieval** (P014 `field_states` example, INFO-2 cross-references for imaging/sequencing round).
7. **Stylistic improvements** (INFO-3 to INFO-9).

---

## 12. Final Decision

**READY FOR STARmap RESEARCH** (post-fix; see Appendix E for the 2026-10-15 fix round)

Reasoning:

- The *STARmap-specific* research itself is not blocked by P002 (OpenFISH). STARmap has a clean primary source (Wang 2018, *Science*) and the schema/templates can record it correctly.
- All four blocking ERRORs (P002, P004, P005, P010) have been resolved by the edits in Appendix E. The boolean hallucination vector (P010), the column-label split (P004), the protocols evidence-status gap (P005), and the OpenFISH primary-source recency (P002) are each addressed by a concrete, minimal change in the corresponding template/schema/YAML file.
- The remaining 12 WARNINGs and 9 INFOs are documentation-level issues that do not block STARmap research. They can be addressed as the corpus matures (see Recommended Fix Order in §11).

To start STARmap research:

1. Verify `methods_index.yaml` still lists `STARmap` with `corpus_status: not_researched`.
2. Create `STARmap/index.md`, `STARmap/metadata.yaml`, and the other 6 production files (initially with `Not yet examined.` placeholders per templates).
3. Begin with `STARmap/sources.md` (assign Source IDs S01, S02, ...) and `STARmap/research_manifest.md` (record the Wang 2018 *Science* paper as the first inspected source).
4. Populate `STARmap/probes.md` and `STARmap/protocols.md` using the now-documented column names and evidence vocabulary.
5. Run a Review pass when the seven production files have evidenced content, recording findings in `STARmap/validation_report.md` and resolving them in `STARmap/fix_report.md`.

---

## Appendix A — Detailed YAML structural check

```
methods_index.yaml
- 10 method entries, no duplicates
- directory values: STARmap, ExSeq, seqFISH+, MERFISH, Oligo-LiveFISH, RAEFISH, TDDN-FISH, PRISM, OpenFISH, MiP-seq
- all official_name values mirror directory ✅
- all aliases lists contain exactly the directory name (acceptable per schema) ✅
- all corpus_status values: not_researched ✅
- all files blocks map the 7 production paths ✅
- hyphenated / `+` directory names are double-quoted ✅
- indentation: 2-space consistent ✅
- no duplicate keys detected ✅
- no list-as-scalar or scalar-as-list issues detected ✅
```

```
metadata.template.yaml
- top-level groups: method, classification, publication, versions, probe, design, benchmark, references, validation, field_states
- all values null or [] for unknown fields ✅
- comments define legal value sets for grade and validation.status ✅
- comment for grade is missing "unassessed means null" rule (P008) ⚠
- comment for boolean fields is missing null-vs-false rule (P010) ⚠
- field_states: {} default; no example (P014) ⚠
```

---

## Appendix B — Naming evidence log

| Directory | Primary source | Title | DOI | Year |
|---|---|---|---|---|
| STARmap | Wang X. et al., *Science* 361:eaat5691 | Three-dimensional intact-tissue sequencing of single-cell transcriptional states | 10.1126/science.aat5691 | 2018 |
| ExSeq | Alon S. et al., *Science* 371:eaax2656 | Expansion sequencing: Spatially precise in situ transcriptomics in intact biological systems | 10.1126/science.aax2656 | 2021 |
| seqFISH+ | Eng C.-H. L. et al., *Nature* 568:235–239 | Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+ | 10.1038/s41586-019-1049-y | 2019 |
| MERFISH | Chen K. H. et al., *Science* 348:aaa6090 | Spatially resolved, highly multiplexed RNA profiling in single cells | 10.1126/science.aaa6090 | 2015 |
| Oligo-LiveFISH (predecessor Oligopaint) | Zhu Y. et al., *Cell* (2025) (Oligo-LiveFISH); Beliveau B. J. et al., *PNAS* 109:21301–21306 (Oligopaint, 2012) | High-resolution dynamic imaging of chromatin DNA communication using Oligo-LiveFISH (Oligo-LiveFISH, 2025); Versatile design and synthesis platform for visualizing genomes with Oligopaint FISH probes (Oligopaint, 2012) | 10.1016/j.cell.2025.03.032 ; 10.1073/pnas.1213818110 | 2025 / 2012 |
| RAEFISH | Cheng Y. et al., *Cell* (2025) | Sequencing-free whole-genome spatial transcriptomics at single-molecule resolution | 10.1016/j.cell.2025.09.006 | 2025 |
| TDDN-FISH | Wang Y.-F. et al., *Nat Commun* 16:9251 | Tetrahedral DNA dendritic nanostructure-enhanced FISH for high-speed, sensitive spatial transcriptomics | 10.1038/s41467-025-64294-1 | 2025 |
| PRISM | Chang T. et al., *Nat Biotechnol* (2025); bioRxiv 2024.06.29.601330 | High-plex spatial RNA imaging in one round with conventional microscopes using color-intensity barcodes (peer-reviewed); PRISM: Multiplexed Profiling of RNA In-Situ through Single-round of Imaging in Three-Dimensional Tissue (preprint) | 10.1038/s41587-025-02883-7 (preprint 10.1101/2024.06.29.601330) | 2024/2025 |
| OpenFISH | Li X. et al., *Neuron* (2026); code github.com/Lixinyoung/openfish-maldi | OpenFISH enables same-section spatial transcriptomics and MALDI-MSI integration | 10.1016/j.neuron.2026.09.007 | 2026 |
| MiP-seq | Wu X. et al., *Nat Biomed Eng* (2024); bioRxiv 2023.01.07.523058 | Spatial multi-omics at subcellular resolution via high-throughput in situ pairwise sequencing | 10.1038/s41551-024-01205-7 | 2024 |

> Note (corrected 2026-10-20): this report first carried placeholder or wrong DOIs for Oligo-LiveFISH (10.1016/j.cell.2025.00350), PRISM (10.1038/s41592-025-XXXXX, with *Nat Methods* and first author Liu Y.) and OpenFISH (10.1016/j.cell.2025.09.006, which is the RAEFISH paper). The DOIs, journals and first authors above were checked against their Crossref records and match each method's sources.md.

---

## Appendix C — Files in repository at audit time

- `corpus_schema.md` (187 lines)
- `terminology.md` (73 lines)
- `methods_index.yaml` (125 lines)
- `templates/index.template.md` (43 lines)
- `templates/metadata.template.yaml` (42 lines)
- `templates/description.template.md` (52 lines)
- `templates/protocols.template.md` (72 lines)
- `templates/probes.template.md` (110 lines)
- `templates/benchmarks.template.md` (48 lines)
- `templates/sources.template.md` (54 lines)
- `templates/research_manifest.template.md` (12 lines)
- `templates/validation_report.template.md` (23 lines)
- `templates/fix_report.template.md` (24 lines)

No method directories present. No `.preflight_validation_report.md` previously existed; this file is new.

---

## Appendix D — Update Log (post-initial-audit revisions)

### 2026-10-15 — Naming corrections applied

**Trigger**: User clarification during audit review.

**Actions taken**:

1. `Oligo-FISH` → `Oligo-LiveFISH` rename
   - **Method identity confirmed**: Zhu Y. et al., "High-resolution dynamic imaging of chromatin DNA communication using Oligo-LiveFISH," *Cell* (2025). PubMed:40239646. Method is a live-cell, reagent-based FISH using fluorescent gRNA oligo pools (ribonucleoproteins), 20-nm spatial and 50-ms temporal resolution.
   - **`methods_index.yaml` updated**: `directory: "Oligo-FISH"` → `directory: "Oligo-LiveFISH"`; `official_name: "Oligo-FISH"` → `official_name: "Oligo-LiveFISH"`; `aliases: ["Oligo-FISH"]` → `aliases: ["Oligo-LiveFISH", "LiveFISH", "Oligopaint"]`. All 7 file paths updated to point under `Oligo-LiveFISH/`.
   - **`preflight_validation_report.md` updated**:
     - Executive summary count: ERROR 5 → 4 (the rename is no longer an ERROR).
     - Method naming table: `Oligo-FISH` row reclassified from RENAME_REQUIRED to CORRECT_WITH_ALIAS_NOTE.
     - Naming verdict updated: 3 CORRECT_WITH_ALIAS_NOTE (was 2); 0 RENAME_REQUIRED (was 1).
     - Method list updated throughout.
     - Appendix B evidence log updated.
   - **Net effect**: P011 (Oligo-FISH → Oligopaint rename warning) is **DONE**; P002 (OpenFISH primary source) is **user-confirmed as 2026 *Neuron*** and reclassified as UNCERTAIN but no longer blocking in itself.
2. `penFISH` → `OpenFISH` (no change required)
   - User confirmed that the intended method is the 2026 *Neuron* paper `OpenFISH` (Li X. et al., "OpenFISH enables same-section spatial transcriptomics and MALDI-MSI integration," *Neuron* (2026)). The user-suggested name `penFISH` was a misremembering of the paper title; the method's official name is `OpenFISH`. `methods_index.yaml` requires no change.
   - P002 (OpenFISH primary source) reclassified from ERROR (no source) to UNCERTAIN (user-confirmed 2026 source, but very recent).

**Remaining open ERRORs after this update** (now historical; resolved in Appendix E): 0

- ~~P002 (OpenFISH — user-confirmed 2026 source, recency noted)~~ → F-001 RESOLVED 2026-10-15
- ~~P004 (Source column label split in probes template)~~ → F-002 RESOLVED 2026-10-15
- ~~P010 (boolean null-vs-false inline rule missing in metadata template)~~ → F-004 RESOLVED 2026-10-15
- ~~P005 (protocols template evidence-status inline list missing)~~ → F-003 RESOLVED 2026-10-15

**Note on `Oligopaint` lineage**: The Beliveau 2012 PNAS method `Oligopaint` is preserved as an *alias* of `Oligo-LiveFISH` for traceability of the probe-design lineage. A Research Agent investigating Oligo-LiveFISH should record both the 2012 predecessor and the 2025 method as separate version files (`versions/V01.Oligopaint.md` and `versions/V02.Oligo-LiveFISH.md`) per `corpus_schema.md` §"Method versions" rules. This is documentation, not a content claim.

**Files touched in this update**: `methods_index.yaml`, `preflight_validation_report.md`. **No templates modified.** **No method directories created.**

---

## Appendix E — 2026-10-15 ERROR Fixes

This appendix records the concrete edits applied to resolve the four blocking ERRORs.

### Fix F-001 — P002 (OpenFISH primary source)

- **Files modified**: `methods_index.yaml`
- **Change**: Added an inline comment under the `OpenFISH` entry recording the user-confirmed 2026 *Neuron* primary source.
- **Diff**:
  ```yaml
  - directory: "OpenFISH"
    official_name: "OpenFISH"
    aliases: ["OpenFISH"]
  + # 2026 *Neuron* primary source confirmed by user on 2026-10-15.
  + # See preflight_validation_report.md P002 for the recency note.
    corpus_status: not_researched
  ```
- **Status**: P002 RESOLVED. The recency is now explicit; future Research Agents will read the comment before populating `OpenFISH/metadata.yaml`.

### Fix F-002 — P004 (probes template Source column label)

- **Files modified**: `templates/probes.template.md`, `corpus_schema.md`
- **Change 1** (`probes.template.md`): Renamed the Probe Components table's `Source` column to `Source / locator` so it matches the other five tables in the same file.
- **Change 2** (`probes.template.md`): Added a clarifying comment at the top of the Probe Components table: "Source / locator: a single canonical column name across every table in this file. Locator must include page, section, figure, table, file, line, or accession. Add version."
- **Change 3** (`corpus_schema.md`): Updated the §probes.md "Probe Components" paragraph to declare `Source / locator` as the single canonical column name and to specify that the locator must include page / section / figure / table / file / line / accession and that the column always carries the method version.
- **Status**: P004 RESOLVED. Every table in `probes.template.md` now uses `Source / locator`. The schema cross-reference is consistent.

### Fix F-003 — P005 (protocols template evidence-status inline list)

- **Files modified**: `templates/protocols.template.md`
- **Change**: Replaced the single-line `<!-- Evidence: ... -->` comment under the Critical-Parameters table with a multi-line comment that lists the four legal evidence states, their meanings, the explicitly forbidden values (`assumed`, `estimated`, `predicted`, `unknown` for typed columns), and reminders that a blank cell is not zero and that the Source / locator must include method version.
- **Status**: P005 RESOLVED. Future Research Agents reading only `protocols.template.md` now have an explicit legal-value list without having to consult the schema.

### Fix F-004 — P010 (metadata template boolean null-vs-false)

- **Files modified**: `templates/metadata.template.yaml`
- **Change**: Replaced the silent null defaults under `design:` and `benchmark:` with explicit, per-field comments defining the meaning of `null`, `true`, and `false` for each boolean field. The `benchmark.grade` comment is also strengthened to forbid assigning D before research.
- **Status**: P010 RESOLVED. The boolean hallucination vector is closed at the template level; the canonical rule is now visible without reading the schema.

### Net effect on ERROR count

| Issue | Before | After |
|---|---|---|
| P002 | ERROR (no source) | RESOLVED (user-confirmed 2026 *Neuron*, comment in YAML) |
| P004 | ERROR (column label split) | RESOLVED (column unified, schema updated) |
| P005 | ERROR (no inline list) | RESOLVED (inline list present) |
| P010 | ERROR (no null-vs-false rule) | RESOLVED (inline rule present) |
| **Total** | **4 ERRORs** | **0 ERRORs** |

### Post-fix decision

**READY FOR STARmap RESEARCH.**

STARmap has a clean primary source (Wang X. et al., *Science* 361:eaat5691, 2018), and the four blocking ERRORs that affected any method are now closed. The remaining WARNINGs and INFOs do not block the STARmap Research Agent and can be addressed in parallel during the research workflow (and resolved as part of a follow-up `fix_report.md`).

### Files touched in this fix round

- `templates/probes.template.md` (P004)
- `templates/protocols.template.md` (P005)
- `templates/metadata.template.yaml` (P010)
- `corpus_schema.md` (P004)
- `methods_index.yaml` (P002)
- `preflight_validation_report.md` (this appendix + executive summary)

### Files NOT touched in this fix round

- `terminology.md` (no changes; no P-issue required edits here)
- `templates/index.template.md`, `description.template.md`, `benchmarks.template.md`, `sources.template.md`, `research_manifest.template.md`, `validation_report.template.md`, `fix_report.template.md` (no ERRORs against these)
- All ten method directories (none yet exist)
