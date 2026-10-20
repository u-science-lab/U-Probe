# Fix Report

Audit only. Review findings were independently checked against primary/official evidence; reviewer recommendations were not treated as facts.

## Resolution Scope

- Method: STARmap; Fix date: 2026-10-15.
- Inputs: AGENTS.md, corpus_schema.md, terminology.md, methods_index.yaml, all templates, and all pre-existing STARmap files.
- Original independent review: validation_report.md, Final Status FAIL. Preserved unchanged as historical review output.
- Findings counted from the actual issue rows: 5 ERROR (R001–R005), 20 WARNING (R006–R025), 6 INFO (R026–R031). The review narrative miscounts/reclassifies these; this report uses the row severities.
- Primary verification: S01 PMC text/captions; S03 NIHMS-2117497 PDF text (design, procedure, Tables 2/4/5); S05 README/root listing; S11 PMC front matter/main text/Data availability. Precise locators are recorded below and in sources.md.
- Existing partial corrections were preserved where verified. Dependencies still contained stale values. No original supplementary-workbook contents or complete deposit inventories were inspected. No original reviewer report was overwritten.
- Subsequent review reference: the post-correction scoped review and checks recorded below; performed by the Fix Agent after editing, not a new independent Review Agent run.

## Issue Resolutions

| Issue ID | Status | Reviewer's concern | Verification | Evidence | Action taken | Files changed | Reason |
|---|---|---|---|---|---|---|---|
| R001 | FIXED | APS/TEMED vs VA-044 gelation | Valid source correction; some current protocol text already used APS, but dependencies remained stale. | S03 Steps 34–40 | Retained APS recipe, room temperature/1 h/nitrogen; removed conflicting dependent recipes; scoped to later protocol. | protocols.md; description.md; sources.md | No unsupported V01 attribution. |
| R002 | FIXED | Proteinase K concentration/time | Current protocol already contained 0.2 mg/mL/1 h; independently confirmed. | S03 Steps 42–43 | Retained correct values, recorded conditional omission for protein comapping and precise source. | protocols.md; sources.md | Correction includes version/sample applicability. |
| R003 | FIXED | MA-NHS 20 mM/overnight | Valid: S03 reports 25 mM/1 h. | S03 Steps 31–32 | Corrected all dependent production records to 25 mM in 100 mM bicarbonate, room temperature, 1 h. | protocols.md; probes.md; description.md; benchmarks.md; sources.md; metadata.yaml | Original chemistry scope retained; later recipe separated. |
| R004 | FIXED | PBSTR 3 × 30 min after clearing | Valid: post-clearing wash is PBST 3 × 5 min. | S03 Step 43; Reagent setup | Corrected reagent/time and conditional omission with proteinase K. | protocols.md; sources.md | PBST and PBSTR are distinct. |
| R005 | FIXED | RCA 8–14 h | Valid: 30 min at 4°C then 2 h at 30°C. | S03 Steps 24–26 | Corrected duration and BSA value; included >100-gene 2.5-fold enzyme adjustment. | protocols.md; description.md; sources.md | Later conditions do not reconstruct original experiments. |
| R006 | FIXED | SEDAL per-cycle concentrations missing | Valid; full mixture explicitly present. | S03 Step 50; Reagent setup | Recorded 10 µM reading, 5 µM total decoding mixture and 0.1875 U/µL ligase, minimum 3 h. | protocols.md; probes.md; description.md; sources.md | Mixture concentration is not concentration per decoding oligo. |
| R007 | FIXED | Permeabilization missing | Valid: methanol and fluorescence-preserving options reported. | S03 Step 10; Reagent setup | Added both options; retained quenching-buffer constituents rather than reducing buffer B to Triton/PBS. | protocols.md | Triton 1% is labeled a recipe-derived value. |
| R008 | FIXED | PFA conditions missing | Valid for S03; not evidence of original-panel fixation. | S03 Step 9A/B | Recorded perfusion, optional post-fixation, cell/tissue alternatives and 20-µm later cryosections. | protocols.md; benchmarks.md | Did not backfill S03 fixation into V01 samples. |
| R009 | FIXED | Hybridization conditions missing | Valid. | S03 Steps 15–18 | Added recipe, per-probe doses, denaturation and 40°C/12–36 h conditions. | protocols.md | Excluded RIBOmap splint from shared source table. |
| R010 | FIXED | Ligation conditions incomplete | Valid; review omitted >100-gene enzyme adjustment. | S03 Steps 21–22 | Added full mixture, room temperature/2 h, buffer handling and 2.5-fold adjustment. | protocols.md | Used complete source conditions. |
| R011 | FIXED | Post-hybridization washes missing | Concern valid; proposed 2× SSC high-salt composition is incorrect. | S03 Steps 19–20; Reagent setup | Added wash sequence; retained 2 mL 20× SSC + 8 mL PBSTR recipe, derived 4× SSC. | protocols.md | Rejected erroneous buffer concentration while fixing omission. |
| R012 | FIXED | Unsupported second embedding | Valid: no such step found in examined S03 procedure. | S03 Steps 31–48; S01 HTC main text | Removed second-embedding recipe and dependent workflow assertion. | protocols.md; description.md; sources.md | No uninspected supplementary attribution invented. |
| R013 | PARTIALLY FIXED | Cycle concentrations and absent color table | Concentrations concern valid; absence claim contradicted by Table 5. | S03 Overview; Tables 4–5; Steps 50, 61 | Added mixture and oligo-side dinucleotide/dye table; corrected disjoint-chunk interpretation. | protocols.md; probes.md; description.md | Rejected reviewer absence interpretation; experiment-specific codebooks remain unexamined, not claimed absent. |
| R014 | FIXED | Nissl vs DAPI/protein imaging | Source distinction needed: S11 uses PI, while later S03 uses DAPI during sequencing. | S01 160-gene Results; S03 Steps 52–53, 62–63 | Separated V01 Nissl, S11 PI and later S03 DAPI, plus post-sequencing secondary antibody. | protocols.md; probes.md; description.md; metadata.yaml | Stains removed from oligo-modification inventory. |
| R015 | FIXED | Aminoallyl-dUTP provenance | Values confirmed; source convention should not be silently interpreted as total dNTP. | S03 Step 24 | Cited Step 24; retained 250 µM dNTP mixture and 20 µM aminoallyl-dUTP notation. | protocols.md; probes.md; description.md; metadata.yaml; sources.md | No unsupported per-base/total conversion. |
| R016 | FIXED | Barcode capacity version scope | Valid scope concern; 1–5 nt is not a universal V02 ceiling. | S01 SEDAL Results; S03 design pp. 8–9 | Scoped 1,024 to five bases; included 4-nt and >5-nt/SeqA–SeqB configurations. | probes.md; description.md; metadata.yaml; benchmarks.md | The source explicitly discusses longer barcodes. |
| R017 | FIXED | Reading N-run orientation | Valid: not fixed N11 or universal 3′ N-run. | S01 Fig. 1E; S03 Table 4 and Step 61 | Recorded cycle-specific library, 5′ mixed-base growth in SeqA, 3′ inverted dT and source list/length caveats. | probes.md; description.md; metadata.yaml | No universal sequence reconstructed from a schematic. |
| R018 | REJECTED | Force primer barcode complement to internal padlock barcode | Reviewer proposed localization is not established; source permits barcode-containing paired terminal regions. | S03 Padlock and primer probe design, p. 8; Table 2 | Removed the existing overconfident localization and marked placement/overlap unresolved; did not apply proposed exclusive-central placement. | probes.md; metadata.yaml; sources.md | Do not invent complement placement or insist terminal arms lack barcode information. |
| R019 | PARTIALLY FIXED | Missing target-design rules | Most omissions valid; reviewer conflates heterodimer filtering with repeat masking and isoform extraction with exon selection. | S03 design pp. 8–9 | Added GC/Tm, spacing, isoforms, Dedupe, NUPACK, library-minus-strand BLAST and recommended pair counts; retained bounded gaps for repeat masking/exon rules. | probes.md; description.md; sources.md; benchmarks.md | Rejecting unsupported subclaims avoids new errors; original panel settings still unknown. |
| R020 | PARTIALLY FIXED | S05 represented as sole canonical pipeline | Concern valid; repository popularity does not establish authority or maintenance. | S03 design pp. 8–9; S05 README/root listing | Recorded multi-tool stack, third-party roles and S05 limited status; corrected README authors/date. | probes.md; description.md; metadata.yaml; sources.md; benchmarks.md | Did not copy stale stars/forks or equate developer-used tools with STARmap-authored software. |
| R021 | FIXED | Infer workbook content from file size | Valid. | S01 supplementary references; no direct workbook inspection | Removed size-based content predictions, speculative completeness and provisional grades. | probes.md; benchmarks.md; research_manifest.md | Uninspected content is unknown; all candidates remain unassessed. |
| R022 | FIXED | RNA-target interspace omitted | Valid for S03 design recommendation. | S03 design p. 9 | Added 1–2 nt gap; removed misleading assembled linear geometry and contiguous-footprint claim. | probes.md; description.md; metadata.yaml | Gap is RNA-site spacing, distinct from padlock/primer pairing. |
| R023 | FIXED | Reading 5′ phosphate attribution | Reported directly, not merely inferred from ligase. | S01 Fig. 1E; S03 Table 4 | Cited explicit sources and corrected decoding 3′ OH to reading 5′ phosphate ligation direction. | probes.md; description.md; metadata.yaml | Chemistry inference unnecessary where explicit notation exists. |
| R024 | FIXED | V01 acrylic-NHS absence scope | Valid scope refinement. | S01 HTC main text; S03 Steps 31–32 | Limited missing-parameter claim to examined main text; noted uninspected original supplements and later MA-NHS scope. | protocols.md; probes.md; sources.md | Not reported in inspected main text is not absent from all official evidence. |
| R025 | FIXED | Original article page range unresolved | Valid documentation cleanup; row severity is WARNING despite inconsistent review summary. | S01 PMC front matter | Retained Science 361(6400):eaat5691; removed unresolved print-page claim in current interpretation. | description.md; research_manifest.md | eaat5691 is an article locator. |
| R026 | FIXED | Explain V01/V02 routing rationale (INFO) | Version routing supported. | S03 Overview; S11 | Added explicit rationale based on protein detection and improved design. | description.md; probes.md | Meaningful workflow change supports a version record. |
| R027 | FIXED | Read-filter rates confirmation (INFO) | 20%/40% confirmed; no numerical correction required. | S01 Results, 1,020-gene experiment | Retained values with experimental scope. | protocols.md; description.md | Confirmation, not a new universal error rate. |
| R028 | FIXED | Full-name capitalization (INFO) | S01 abstract uses lowercase spelled-out form. | S01 abstract | Used original lowercase in metadata/index full-name field; retained capitalization as a display variant. | metadata.yaml; index.md | Preserve original naming. |
| R029 | FIXED | Optional authorship note (INFO) | S01 author contributions inspected; not necessary to probe-design retrieval. | S01 Author contributions | No optional authorship narrative added; corrected actual S05 credits separately under R020. | None for optional note | No substantive corpus finding requires a change. |
| R030 | FIXED | V02 deposit confirmation (INFO) | S11 Data availability confirms both accessions; Results/Methods contradict the reviewer’s 1,022-gene panel attribution. | S11 Data availability | Retained S09/S10 accession links; corrected parent-study panel to 2,766 genes, with separate 64-gene validation; deposit inventories left unexamined. | probes.md; benchmarks.md | Accession confirmation does not prove sequence completeness or absence. |
| R031 | FIXED | Long metadata probe strings (INFO) | Schema calls for compact short strings. | corpus_schema.md; metadata template | Shortened probe summaries and moved scoped detail to Markdown. | metadata.yaml | Metadata remains shallow and conservative. |

## Additional Dependent Corrections

The cross-file check found further errors connected to the reviewed claims: primer incorrectly called the RCA template; reading/decoding ligation ends reversed; readout described as three independent chunks; one probe pair confused with pairs per gene; source tables misnumbered; S03 bibliography conflated with S11; S11 pagination incorrect; circular RNA barcodes asserted without support; original imaged thickness treated as a cryosection recipe. Further S11 Methods/Results verification corrected a missed 1,022-gene claim to the reported 2,766-gene main panel (plus 64-gene validation), and separated S11 acrylic NHS/PI from later S03 MA-NHS/DAPI. S11 uses eight/four RNA rounds for these respective panels. S01 Fig. 6 linear seven-round readout is four genes per round for 28 genes, distinct from five-base decoding. These were corrected or qualified from S01/S03/S05/S11 without inventing original experimental settings. Historical research notes remain with a superseding inspection addendum.

## Post-Correction Scoped Review

The Fix Agent reread the changed files after edits and rechecked their assertions against the inspected sources. Outcome: **PASS within the documented scope**, with no remaining ERROR-level assertion from R001–R005. This is Fix-stage final validation, not a claim that another independent agent reviewed these edits. The original independent FAIL remains visible in validation_report.md; this subsequent scoped PASS reconciles its findings for the method lifecycle.

Nonblocking warnings remain: original supplementary content/reference context and dataset inventories are unexamined; strict benchmark candidates remain unassessed; internal barcode overlap/placement is unresolved; S03 Table 4/Step 61 library range and mixed-base chemistry wording need sequence-level confirmation; software maintenance and S05 equivalence remain unverified. These uncertainties are explicit, not reported facts or validated synthesis recipes.

## Final Fix Checklist

- Every review row processed: R001–R031, including all ERROR/WARNING rows; incorrect reviewer subclaims explicitly rejected with source rationale.
- Corrections independently checked against S01/S03/S05/S11; no fabricated sequences, reaction defaults, reference builds or exact codebooks introduced.
- Dependent values checked across all seven production files; original and later protocol scopes separated.
- Original terminology retained beside standardized component categories; exact-source table references distinguished from schematic/degenerate representations and calculations.
- Unknown/unexamined data distinguished from bounded unavailable values; benchmark grades remain null/unassessed rather than guessed A–D.
- Required template headings, tables, YAML structure/types/duplicate keys, source-ID references, production links, lifecycle agreement and issue coverage checked programmatically after correction.
- Subsequent scoped review completed after correction; source ambiguity and incomplete inspection documented above and in production files.
- Research/review history preserved; no claim of independent re-review or full benchmark reproducibility.

## Final Validation Status

VALIDATED WITH WARNINGS

Method lifecycle: validated_with_warnings in metadata.yaml and methods_index.yaml. This status applies to the conservative, scoped corpus statements; it does not certify unresolved sequences, original experimental reconstruction or strict benchmark ground truth.

### Recorded QA Result

Post-correction checks completed on 2026-10-15: eight check groups passed, zero structural/consistency errors; 31/31 issue rows covered exactly once (27 FIXED, 3 PARTIALLY FIXED, 1 REJECTED). The partial resolutions reject unsupported reviewer subclaims and retain explicit unknowns; they do not leave an ERROR-level assertion in production. Checked YAML duplicate keys/types and lifecycle agreement, required headings, table shapes and controlled vocabularies, source IDs, issue coverage, preservation of the original review, local index links, and targeted parameter/version regressions. `git diff --check` passed.

Preserved validation_report.md SHA-256: `d3281bcc0aa352d8efaa547741fbd265281fee0ef5620bd9c8751f8fbccd9ac7`.
