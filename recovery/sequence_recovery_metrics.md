# Probe Sequence Recovery metrics

Version: 1.0. These definitions preserve the existing shared sequence recovery
metric specification. This study evaluates recovery of published sequences;
it does not compare design tools or establish experimental probe performance.

## Core metrics

All rates are percentages. A reference record is a complete published oligo.

| Metric | Definition |
|---|---|
| SRR (Sequence Recovery Rate) | `100 * matched unique reference oligos / valid unique reference oligos` |
| TRR (Target Recovery Rate) | `100 * targets with at least one complete reference oligo recovered / reference targets` |
| MTR (Mean Target Recovery) | `100 * mean_g(matched reference oligos for g / reference oligos for g)` |
| CPSR (Complete Probe-Set Recovery Rate) | `100 * completely recovered reference functional sets / reference functional sets` |

For MTR, include only targets with at least one valid reference oligo, and give
each target equal weight. TRR requires at least one oligo, not a complete pair.
CPSR applies to functional sets with multiple required roles, such as padlock
plus primer. Single-oligo functional units use `N/A` for CPSR.

## Matching and counting

1. Compare complete 5'-to-3' nucleotide sequences with **100% exact equality**.
   Uppercase and remove whitespace. Do not reverse-complement, accept partial
   matches, or allow mismatches during evaluation.
2. Require the same target and oligo role. The implementation compares target
   and role labels case-insensitively; it does not infer gene aliases.
3. Fix published non-target functional sequences (barcodes, scaffold and other
   linkers) before evaluation. A different codebook changes full sequences.
4. Deduplicate references by `(method_version, target_id, oligo_type, full_DNA)`.
   Count each reference record at most once, regardless of generated duplicates.
   Identical DNA assigned to another target or role remains a separate record.
5. Record chemical modifications separately. Explicit slash-delimited tokens
   such as `/5Phos/` are removed before DNA comparison; arbitrary non-DNA text
   is not silently discarded. Modification export is a separate structural
   check and does not change SRR.
6. For CPSR, all required full-sequence members must occur in the **same
   generated functional set**, with the same target and corresponding roles.
   Generated IDs need not equal published IDs: report their exact-member-set
   mapping. Hits spread over unrelated generated pairs cannot recover a pair.
7. Main results use the complete candidate collection after its explicitly
   documented basic constraints and before final count limits. Report final
   selection as a separate stage. If no final selection was performed, use
   `N/A`, not zero and not the candidate rates.
8. Fix the method version, reference database, target list, codebook and design
   parameters; preserve code/configuration versions and their provenance.
9. Exclude incomplete or unassignable reference records before computing the
   denominators, and report exclusion counts and reasons. Invalid generated
   records are errors rather than silently excluded successes.

New generated oligos that were not published are **not false positives**.
Do not use an extra-output-as-error interpretation of precision or F1 here.
Recovery alone does not prove specificity, experimental validity or the
ability to select the original final probe set.

## Data and outputs

Canonical fields are `method_version,target_id,probe_set_id,oligo_type,full_sequence`.
The evaluator also accepts `gene,pair_id,probe_id` and optional
`stage,modification_5p,dna_sequence`. Stage defaults to `candidate`.

The compact outputs are `summary.csv`, `per_gene.csv`, `sequence_matches.csv`
and `pair_matches.csv`. The `all` role is the main oligo-level result; separate
role rows provide diagnostic rates. Empty CPSR fields denote non-applicability.
`generated_records` counts generated rows, including DNA repeated across pairs;
`generated_unique_reference_hits` counts distinct recovered reference keys, not
the number of distinct novel generated sequences.
