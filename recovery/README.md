# U-Probe Probe Sequence Recovery

## Overview

Evaluate whether U-Probe can independently regenerate experimentally published
probe sequences from reference target sequences under method-specific probe
architectures. Methods share SRR, TRR, MTR and, where applicable, CPSR;
[metric definitions](sequence_recovery_metrics.md) specify matching and counting.

## Recovery Results & Progress

| Method | Dataset | Candidate SRR | CPSR | Final Recovery |
|---|---|---|---|---|
| [STARmap PLUS](starmap_plus_2023_64gene/) | 2023 / 64 genes | 100% | 100% | TODO |
| ExSeq | 2021 / 42 genes | - | - | TODO |
| MERFISH | TBD | - | - | TODO |
| seqFISH+ | TBD | - | - | TODO |
| MiP-seq | TBD | - | - | TODO |
| PRISM | TBD | - | - | TODO |

STARmap is **Unfiltered Candidate Recovery**, not final probe selection.
Dashes indicate no verified full-dataset result. ExSeq has external prepared
references and small-subset pilots; these do not establish 42-gene recovery.

## TODO

**STARmap PLUS**
- [ ] Evaluate predefined filtering and final four-pair selection.

**ExSeq**
- [ ] Import prepared references and finalize executable YAML.
- [ ] Evaluate and document full-panel recovery.

**Other Methods**
- [ ] Prepare reference sets, configurations and recovery records.

**General**
- [ ] Publish matching reference-data acquisition artifacts.
- [ ] Verify fresh-machine reproduction and extend cross-method summaries.

## Usage

Each method README provides two commands: generate candidates with the native
U-Probe API, then evaluate them with [common/sequence_recovery.py](common/sequence_recovery.py).
No method-specific evaluation wrapper is needed. Inputs use `gene`/`target_id`,
`oligo_type`, `full_sequence`, and `pair_id`/`probe_set_id` for functional groups.
Full DNA, target, role and same-generated-pair membership determine recovery;
chemical modifications are separate. Novel sequences are not false positives.

Keep the published references and all four compact report tables for auditing,
including full matched sequences and generated pair IDs. Do not commit large
generated candidates or reference databases. The small metric regression tests
are optional: `python -m unittest discover -s recovery/common`.
