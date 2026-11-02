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
| [seqFISH+](seqfish_plus_2025_1194gene/) | 2025 / 1,194 genes | 99.089542% | N/A | Not evaluated |
| [MiP-seq](mipseq_2024_217gene/) | 2024 / 217 genes | 99.308756% | 99.078341% | Not evaluated |
| [PRISM](prism_2025_64gene/) | 2025 / 64 genes | 100% | N/A | Not evaluated |
| [RAEFISH](raefish_2025_64gene/) | 2025 / 64 Human genes | 100% | 100% | Not evaluated |
| [TDDN-FISH](tddnfish_2025/) | 2025 / 9 HeLa genes | 100% | N/A | Not evaluated |

STARmap, ExSeq, MERFISH, seqFISH+ and MiP-seq report **unfiltered candidate recovery**, not final probe
selection. ExSeq recovers 333/334 padlocks and 42/42 targets using mm10;
TRR is 100% and MTR is 99.7024%. ExSeq final selection has not been evaluated.
MERFISH recovers 22,004/22,005 encoding templates and 242/242 targets; MTR is 99.995508%.
seqFISH+ recovers 37,548/37,893 110-nt probe bodies and 1,192/1,194 targets;
MTR is 99.050149%. Complete synthesis oligos are not evaluated.
MiP-seq recovers 862/868 full oligos and 217/217 targets; MTR is 99.308756%.
Same-generated-group recovery covers 430/434 padlock-RCA pairs.
PRISM recovers 64/64 complete 120-nt padlocks and 64/64 targets using mm10; MTR is 100%.
RAEFISH recovers 192/192 complete synthesis templates, 64/64 targets and 64/64
same-generated-set padlock/priming/encoding triplets using GRCh38/GENCODE v46; MTR is 100%.
TDDN-FISH recovers 54/54 complete 65-nt primary oligos and 9/9 targets using
GRCh38/GENCODE v47; TRR and MTR are 100%. Fixed tetrahedral components and
chemical modifications are checked separately; physical assembly is not evaluated.

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

## Planned Recoveries

| Priority | Method | Validation focus | Status |
|---|---|---|---|
| P1 | osmFISH | Linear RNA-FISH probes | Planned — reference audit required |
| P1 | split-FISH | Paired probes and bridge/readout assembly | Planned — reference audit required |
| P1 | DNA Oligopaint | Genomic-region probe pools | Planned — reference audit required |
| P1 | HCR v3.0 | Split-initiator probe pairs | Planned — reference audit required |
| P2 | SABER-FISH | Modular signal amplification | Planned — reference audit required |
| P2 | clampFISH 2.0 | Multi-component amplification probes | Planned — reference audit required |
| P2 | HT-smFISH | High-throughput probe libraries | Planned — reference audit required |
| P2 | FISHnCHIPs | Multiplexed RNA probe pools | Planned — reference audit required |
| P2 | π-FISH rainbow | Multi-level probe assembly | Planned — reference audit required |
| P2 | Oligo-LiveFISH | CRISPR-based live-cell DNA imaging | Planned — reference audit required |

Before execution, confirm complete original probe sequences, target mapping,
fixed functional segments, probe grouping and reference genome version.
