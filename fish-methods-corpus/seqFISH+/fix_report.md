# Fix Report

Audit only. Resolutions reconcile the independent review; missing evidence remains explicit.

## Resolution Scope

- Date: 2026-10-15 (Asia/Shanghai). Method: V01, original 2019 RNA seqFISH+.
- Review reference: validation_report.md, Stage 2 PASS, zero ERRORs, three WARNINGs and two INFO findings. The review is preserved unchanged.
- Revision context: HEAD 2b95f97c333ad7c1503fd78535cff03d5e84cf72; method files remain uncommitted. Existing ExSeq edits were preserved. Only the seqFISH+ lifecycle entry was changed in methods_index.yaml.
- Independently checked the governing instructions, schema, terminology, templates, all method production/audit text and supporting source records. Inspected S01 Fig. 1a visually and published-layout Methods (paper.txt, pp. 6–7); opened live PMC full text. Reopened the S05 eight-page mirror and current official V.1 destination. The latter returned no readable procedure; no authenticated official protocol payload or checksum comparison was recovered.
- Subsequent scoped validation: Fix-stage reconciliation and final checks below, after corrections. This is not a new independent Review Agent review and does not replace the preserved Stage 2 report.

## Issue Resolutions

| Issue ID | Status | Reviewer's concern | Verification | Evidence | Action taken | Files changed | Reason |
|---|---|---|---|---|---|---|---|
| R001 | FIXED | Component table loses reported versus inferred distinction | Fig. 1a has no 5′/3′ labels; internal/flanking placement is a schematic interpretation. Methods support lengths, roles and complement relationships | S01, Fig. 1a; Methods, Primary-probe design and Readout-probe design; V01 | Position cells explicitly say inferred; evidence cells distinguish reported roles/lengths from inferred position. Complete strand order remains unresolved | probes.md | Finding valid; no unsupported polarity was assigned |
| R002 | PARTIALLY FIXED | Mirror-only protocol claims lack official-body authentication | Mirror cover/Procedure agree with the attributed guidance; official V.1 route was opened but exposed no readable body. Identifiers cannot authenticate the procedure | S05 mirror pp. 3–4, 6–7; official protocols.io V.1 destination; S01 reference 40; V01 | Qualified all unique mirror claims and marked them inferred; separated derived 240 slot count from the mirror's unique-readout assertion. Retained S01-supported anchoring and protocol conditions independently. Authentication remains unresolved | sources.md; protocols.md; probes.md; description.md; research_manifest.md | Attribution corrected, but no authenticated payload obtained; warning remains open |
| R003 | FIXED | Important conditions retained only in prose | Re-read source Methods for Label-IT, PFA/Tris, AcX, label-X/re-embedding and dye coupling; checked remaining retained conditions | S01, Methods, Coverslip functionalization, Cell culture experiment, Tissue slices experiment, Primary-probe construction, Readout-probe synthesis, seqFISH+ imaging and smFISH; V01 | Added 35 scoped records for previously prose-only preparation, stabilization, manufacture, storage, acquisition and comparison conditions. Preserved approximate values, source wording and unavailable fields; did not combine alternative sample branches | protocols.md | Finding valid; structured retrieval now includes the omitted conditions |
| R004 | FIXED | Current official protocol destination should be recorded | Direct opening identified RNA seqFISH+ Supplementary Protocol V.1 but yielded zero readable procedure lines | S05, current protocols.io V.1 destination; access 2026-10-15 | Added current route alongside legacy accession; retained historical access failures and distinguished discovery from body inspection | sources.md; research_manifest.md | Finding valid; route recovery does not resolve R002 |
| R005 | FIXED | Unresolved length contribution assigned a functional component category | Approximate 93-nt total and 88-nt component sum do not establish another component or its role; terminology requires an evidenced function | S01, Fig. 1a and Primary-probe construction; terminology.md, other_functional_sequence; V01 | Removed residual component row; kept unresolved order and length-accounting notes in Canonical Architecture and Unresolved Information | probes.md | Finding valid; no spacer identity or exact residual length introduced |

## Final Validation Status

**VALIDATED WITH WARNINGS**

Subsequent scoped Fix verdict: **PASS**. No unresolved ERROR-level factual or structural problem was identified. R002 remains a nonblocking provenance warning: the procedure was read only in a mirror and its official-body equivalence is unresolved. Original oligo strings, synthesized strand order and exact reference/input releases remain explicitly unresolved; this verdict does not establish them.

The lifecycle is `validated_with_warnings` in metadata.yaml and methods_index.yaml. index.md now routes to the review and Fix reports. D01 remains C and D02 remains D for probe-design comparison; aggregate grade C and unresolved official design-code availability remain consistent across files.

### Full Fix checklist and final scoped checks

| Requirement | Final result |
|---|---|
| Every review issue processed | R001–R005 individually verified and resolved above; no findings silently discarded |
| Reviewer claims independently verified | Source figure/Methods and mirror re-examined; official route retried; reviewer conclusions were not assumed |
| Corrections evidence-supported | Layout explicitly inferred; added table records reproduce already retained source conditions; mirror-only claims demoted pending authentication |
| Incorrect findings rejected where appropriate | None required rejection after verification |
| Dependent files updated consistently | Architecture table/prose, S05 claims/source record, parameter prose/table, index and both lifecycle fields reconciled |
| No unsupported new claims | No nucleotide string, exact residual length, reference release, assumed temperature or manufactured reaction recipe introduced |
| Unresolved issues remain explicit | R002 authentication and original sequence/reference/design-code gaps retained; unusual source units/solvent wording preserved |
| Every resolution documented | Stable R001–R005 IDs and allowed resolution statuses used |
| Final validation performed | YAML duplicate-key/type/schema checks, all required headings and files, table widths, source IDs, local links, routing/status and targeted factual consistency checked; results in evidence/fix_check_results.json |
| Final status justified | Stage 2 PASS plus post-correction scoped PASS; no ERRORs, retained provenance warning prevents unqualified VALIDATED |

Validation limits: no wet-lab execution, original design pipeline reproduction, or new independent full scientific review. Existing downloaded evidence was used without claiming independently authenticated fresh downloads. Structural checks do not establish missing scientific facts. The historical Research self-check is preserved and was not rerun as if it were a Fix validator.
