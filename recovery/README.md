# U-Probe Probe Sequence Recovery

## Overview

Evaluate whether U-Probe can independently regenerate experimentally published
probe sequences from reference target sequences under method-specific probe
architectures. Methods share SRR, TRR, MTR and, where applicable, CPSR;
[metric definitions](sequence_recovery_metrics.md) specify matching and counting.

## Recovery Results & Progress

| Method | Dataset | Candidate SRR | CPSR | Final Recovery |
|---|---|---|---|---|
| [STARmap PLUS](starmap_plus_2023_64gene/) | 2023 / 64 genes | 100% | 100% | N/A |
| [ExSeq](exseq_2021_42gene/) | 2021 / 42 genes | 99.7006% | N/A | Not evaluated |
| [MERFISH](merfish_2021_242gene/) | 2021 / 242 genes | 99.995456% | N/A | Not evaluated |
| seqFISH+ | TBD | - | - | - |
| MiP-seq | TBD | - | - | - |
| PRISM | TBD | - | - | - |

STARmap, ExSeq and MERFISH report **unfiltered candidate recovery**, not final probe
selection. ExSeq recovers 333/334 padlocks and 42/42 targets using mm10;
TRR is 100% and MTR is 99.7024%. ExSeq final selection has not been evaluated.
MERFISH recovers 22,004/22,005 encoding templates and 242/242 targets; MTR is 99.995508%.
Dashes indicate no verified full-dataset result.

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
