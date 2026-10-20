# FISH methods corpus

A structured, source-annotated record of ten imaging-based spatial-omics methods (STARmap, ExSeq,
seqFISH+, MERFISH, Oligo-LiveFISH, RAEFISH, TDDN-FISH, PRISM, OpenFISH, MiP-seq): what each method
is, its protocols, its probe architecture and design rules, candidate benchmark probe sets, and the
sources behind every claim. It was researched and written by Qian Zhang (@jshn2022) and first
proposed in U-Probe pull request #24. Start at `methods_index.yaml`, `corpus_schema.md` and
`terminology.md`; in each method folder, start at `index.md`.

**Not part of the U-Probe Agent Skill.** Nothing here ships in the `uprobe` package or is copied by
`uprobe skill install`; the skill lives in `uprobe/skills/u-probe/`. Reviewed parts may later be
distilled into the skill's own technique cards.

`AGENTS.md`, `prompt.md`, `templates/`, the per-method `research_manifest.md`, `validation_report.md`
and `fix_report.md`, and `preflight_validation_report.md` record how the corpus was built (the
instructions given to the research, review and fix agents, and their reports). They are kept as
records, not as instructions for agents working in this repository.

## Review status

| Method | Author's Review/Fix pass | Spot-checked against the papers (#24 review) | Status |
| --- | --- | --- | --- |
| STARmap | yes | yes (incl. the 2025 STARmap PLUS protocol) | validated with warnings |
| ExSeq | yes | yes | validated with warnings |
| seqFISH+ | yes | yes | validated with warnings |
| MERFISH | no | yes | draft |
| RAEFISH | no | yes | draft |
| Oligo-LiveFISH | no | yes | draft |
| MiP-seq | no | the points raised in the review were corrected | draft |
| PRISM | no | no | draft |
| TDDN-FISH | no | no | draft |
| OpenFISH | no | no | draft |

"Validated with warnings" is the author's own Review/Fix verdict (`fix_report.md`); a spot check
covers selected claims, not every line. Treat drafts as unreviewed and check claims against the cited
source before using them.

## Not included

- `evidence/`: about 60 mentions (in `sources.md`, `research_manifest.md`, `probes.md`,
  `benchmarks.md` and the review reports) name files in the author's local research folder:
  downloaded papers, supplements, inventories and checksums. That folder is not in this repository;
  the mentions are kept as provenance, and the recorded checksums identify the files.
- `Oligo-LiveFISH/resources/`: the paper's supplementary PDF and Table S2 (Elsevier supplements to a
  Cell article that is not openly licensed) were removed because they cannot be redistributed. Get
  them from the publisher; their SHA-256 sums are in `Oligo-LiveFISH/sources.md`.

## Licence

The corpus falls under the repository's MIT licence (`../LICENSE`).
