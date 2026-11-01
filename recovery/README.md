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

Follow each method README for candidate generation. The shared standard-library
evaluator supports plain CSV and gzip, arbitrary oligo roles, separate candidate
and final stages, and an optional wide pair table:

```bash
python recovery/common/sequence_recovery.py \
  --reference-oligos path/reference_oligos.csv --generated path/generated.csv.gz \
  --method-version example_v1 --output path/report
python -m unittest discover -s recovery/common -p 'test_*.py'
python -m unittest discover -s recovery/starmap_plus_2023_64gene -p test_evaluation.py
```

Required fields are `gene` (or `target_id`), `oligo_type`, and `full_sequence`.
Use `pair_id`/`probe_set_id` for functional groups and `probe_id`/`oligo_id` for
generated oligos. `stage` defaults to `candidate`; modifications are separate.
Pass `--reference-pairs` to cross-check reference groups. Outputs are the four
recovery CSV reports plus exclusion reasons. Generic pair reports use JSON
`member_matches`; method adapters may retain named role columns. New unpublished
oligos are not false positives. Generated data and reference databases are not
committed.
