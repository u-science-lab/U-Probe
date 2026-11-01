# Exact sequence recovery

`sequence_recovery.py` is a method-independent, standard-library evaluator. It reads a reference long table and a generated long table, compares complete 5-prime to 3-prime DNA sequences, and reports recovery without treating unpublished output as false positives.

```sh
python recovery/common/sequence_recovery.py \
  --reference-oligos path/reference_oligos.csv \
  --generated path/generated.csv.gz \
  --method-version example_v1 \
  --output path/report
python -m unittest discover -s recovery/common -v
```

Required fields are `gene` (or `target_id`), `oligo_type`, and `full_sequence`. The optional `method_version` field defaults to the CLI value, which is `unspecified` unless supplied. Reference and generated group IDs use `pair_id` or `probe_set_id`; generated oligo IDs use `probe_id` or `oligo_id`. Generated `stage` is `candidate` or `final`, defaulting to `candidate`. Candidate means the set before an output-number limit. A final-stage result is produced only when final-stage records are provided. Optional `dna_sequence` is checked against the normalized full sequence; `modification_5p` is recorded separately.

Whitespace is removed, DNA is uppercased, and explicit slash-delimited modification tokens are recorded separately. DNA-only slash tokens are rejected to avoid deleting nucleotides. The remaining complete sequence must contain only A/C/G/T. Matching never reverse-complements, aligns subsequences, permits mismatches, or substitutes other genes or roles. Reference deduplication uses `(method_version, casefold(target), casefold(role), normalized full DNA)`. Invalid references are excluded with reasons; invalid generated sequences fail rather than silently disappearing.

SRR is the fraction of unique reference sequences recovered. TRR is the fraction of reference targets with at least one recovered sequence. MTR is the arithmetic mean of per-target recovery fractions. CPSR requires every required reference role/sequence to coexist under one generated group ID at the same target and version. Generated group labels need not equal published labels; the report records the sequence-established mapping. Single-role groups are excluded from CPSR and have `cpsr_applicable=False`; blank CPSR means inapplicable. One group may have one member per role; conflicting different members assigned the same group and role are rejected.

An optional `--reference-pairs` wide table cross-checks reference group membership. It must contain the target and group IDs and, for each role in the long table, `<role>_full` or `<role>_dna`. Role names are unrestricted; no method-specific role names or default method version are embedded in the evaluator.

Plain CSV and gzip `.csv.gz` are supported. Generated rows are streamed; only exact reference hits and their group membership are retained. Output comprises `summary.csv`, `per_gene.csv`, `sequence_matches.csv`, `pair_matches.csv`, and `reference_exclusions.json`. The generic pair report uses a JSON `member_matches` field for arbitrary roles. Method-local adapters may translate this into established role columns. Rates use percent units. Stage-wide `generated_records` and `generated_unique_reference_hits` are repeated across role summary rows; the latter counts distinct reference keys recovered, not all unique new output sequences. No all-output uniqueness claim or precision/F1 score is made.

Tests cover explicit modifications and malformed DNA, complete-sequence direction, target/role identity, duplicates, uneven per-target denominators, stage separation, same-generated-group membership, reference-table consistency, gzip streaming, novel sequences, single-role CPSR inapplicability, and arbitrary functional role names.
