# U-Probe Probe Sequence Recovery

## Overview

Evaluate whether U-Probe can independently regenerate experimentally published
probe sequences from reference target sequences under method-specific probe
architectures. This project records recovery of published molecules, not a
comparison of design tools. Methods share SRR, TRR and MTR, plus CPSR when a
functional unit contains multiple oligos.

[Metric definitions](sequence_recovery_metrics.md) and
[shared evaluation code](common/) apply across methods. Detailed configurations,
results and limitations belong in the method directories.

## Progress

| Method | Dataset | Reference | YAML | Candidate Recovery | Final Recovery |
|---|---|---|---|---|---|
| [STARmap PLUS](starmap_plus_2023_64gene/) | 2023 / 64 genes | Done | Done | Done: unfiltered | TODO |
| ExSeq | 2021 / 42 genes | Prepared externally | Initial specification | Partial external pilots; full panel TODO | TODO |
| MERFISH | TBD | TODO | TODO | TODO | TODO |
| seqFISH+ | TBD | Source collected | Engineering preset only | Published-sequence recovery not verified | TODO |
| MiP-seq | TBD | Source collected | Existing preset; recovery not audited | Not verified | TODO |
| PRISM | TBD | Source collected | TODO | Not verified | TODO |

The ExSeq reference archive includes the 42-gene design specification, targets,
oligos and barcode assignments. Separate external pilots cover different small
gene subsets; they do not establish complete 42-gene recovery. The seqFISH+
engineering pilot used synthetic readout placeholders, so it is not evidence of
published full-sequence recovery. Source-collected entries have no verified
recovery result in this record. Other methods' existing files are not migrated
or modified by the STARmap PLUS addition.

## TODO

### STARmap PLUS

- [x] Extract published probe sequences.
- [x] Prepare executable design YAML and portable candidate runner.
- [x] Generate unfiltered candidate probes.
- [x] Evaluate candidate sequence recovery and independently verify results.
- [ ] Evaluate recovery after predefined filtering.
- [ ] Evaluate recovery after final probe selection.

### ExSeq

- [x] Prepare published reference sequences in the external archive.
- [x] Prepare an initial probe architecture specification.
- [ ] Import and audit the intended 42-gene recovery package.
- [ ] Validate and finalize executable YAML.
- [ ] Run full-panel candidate sequence recovery.
- [ ] Validate and document full-panel recovery results.

### Other methods

- [ ] Prepare and normalize published reference probe sets.
- [ ] Implement method-specific YAML configurations.
- [ ] Run candidate sequence recovery.
- [ ] Document results and limitations.

### General

- [x] Standardize common recovery evaluation code.
- [x] Establish compact results and provenance formats for the first method.
- [x] Verify compact evidence, metric consistency and small-run equivalence.
- [ ] Generate a numerical cross-method recovery summary.
- [ ] Publish matching external reference-data acquisition artifacts.
- [ ] Verify complete fresh-machine reproduction with those public artifacts.
