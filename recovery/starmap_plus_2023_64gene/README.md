# STARmap PLUS 2023 - Unfiltered Candidate Sequence Recovery

## 1. Reference publication and dataset

Zeng et al., *Nature Neuroscience* (2023),
[DOI: 10.1038/s41593-022-01251-x](https://doi.org/10.1038/s41593-022-01251-x).
[Author manuscript and supplementary material](https://pmc.ncbi.nlm.nih.gov/articles/PMC11332722/).
The 64-gene section of Supplementary Table 1 contains 512 oligos in 256
padlock/primer pairs, four pairs per gene. Source cells and checksums are in
[reference/source.json](reference/source.json).

## 2. Recovery objective

Independently regenerate published full oligo sequences from genomic reference
data under the panel's fixed architecture and gene-specific codebook:

`genes -> transcript extraction -> target-site enumeration -> probe assembly -> exact matching`.

Published target-binding sequences and full reference oligos are **evaluation
inputs only**. They are not read by the candidate generator. This is an
unfiltered candidate-space recovery study, not reconstruction of the original
ranking, specificity pipeline or final four-pair selection.

## 3. Reference probe sequences

- [reference_oligos.csv](reference/reference_oligos.csv): 512 complete records,
  with DNA, chemical modifications, role, target and published pair ID.
- [reference_pairs.csv](reference/reference_pairs.csv): 256 corresponding pairs.
- Reference deduplication leaves 512 oligos, 64 genes and 256 pairs; no records
  were excluded. The long and wide reference tables agree.

All 256 padlock arms and 256 primer arms map exactly to the assigned genes in
the recorded GRCm39 transcripts. Panel-wide geometry was inferred once as RNA
`padlock_site -> gap of 1 or 2 nt -> primer_site`. This supplies an architecture
constraint, not published target coordinates or binding sequences to generate.
Gap values are not claimed to be explicit published computational parameters.

## 4. U-Probe design configuration

The execution base is U-Probe commit
`a76baff9fa79c6ef45308b55877aa507e9c62a27`. No core algorithm was changed.
[config/design.yaml](config/design.yaml) uses the native extraction, layout and
probe-expression interfaces. Its portable runner iterates each annotated
transcript separately; invoking the default all-transcript extraction directly
is not an equivalent command.

| Setting | Value |
|---|---|
| Genome | GRCm39, assembly GCF_000001635.27 |
| Annotation | NCBI RefSeq GCF_000001635.27-RS_2024_02, as recorded in the GTF header |
| Transcript policy | All 308 annotated isoforms: 157 NM, 140 XM, 4 NR, 7 XR |
| Binding lengths | 19-25 nt independently for each arm |
| Gap | 1-2 nt |
| Whole target window | 39-52 nt, step 1; 98 geometries |
| Basic constraint | Uppercase A/C/G/T only |
| Candidate deduplication | Same gene, target-window DNA and layout |
| GC/Tm, off-target filtering and ranking | Not enabled; original panel-specific thresholds not confirmed |
| Final count limit | None; final selection not performed |

All 308 complete transcript sequences, including 113 negative-strand records,
agreed across independent raw-GTF splicing, prebuilt transcript FASTA and native
U-Probe splicing. The all-isoform policy includes the 148,183-nt predicted
Pcdh9 transcript XR_003950833.2; it was not removed to reduce computation.

Fixed DNA assembly (all oligos written 5'-to-3'):

```text
padlock = A + code + TA + rc(padlock_site) + AATTATTAC + code + CATACACTAAAGATA
primer  = rc(primer_site) + TA + rc(code) + TTATCTT
```

Padlock 5'-phosphate is exported separately. The initial architecture
specification was converted to native extraction/layout/probe fields. Reference
selection was resolved to all annotated GRCm39 isoforms, the complete 19-25 nt
arm Cartesian product replaced the initially observed total-length range, and
the gap was fixed to 1-2 nt following the panel-wide geometry audit. No recovery-
guided parameter retuning was performed. During repository curation, only the
protocol display name and unused execution-annotation key were renamed.

| Configuration artifact | SHA-256 |
|---|---|
| Historical executed YAML | `4bff8b65d58e2feb5d218f064b7279ea14d568f13c483b962cf4d72a64bcd608` |
| Current `config/design.yaml` | `db3d3555784580ac79d71d3082774ba492b926a28526e41b26d554f479aa7b81` |

The two YAML versions have identical genes, codebook, architecture and
generation parameters; their byte hashes differ because of the metadata renames.

## 5. Candidate generation

The actual complete run produced **44,271,276 candidate pair records** and
**88,542,552 oligo records**. The same oligo DNA may appear in multiple pairs;
these are record counts, not a claim about unique generated DNA molecules.
An independent transcript-window cardinality calculation agreed with every
gene's candidate count. All candidates were kept before any count restriction.

The historical run used native U-Probe extraction, layout expansion and
construction. Execution-only changes addressed a 64-GiB cgroup limit and
parallelized the long final gene, without changing design parameters. Small
complete Npy outputs were byte-identical after decompression across execution
versions. Design/export/merge took 3,356.25 seconds including recovery from
failed workers and execution optimization. The run used Python 3.10.14,
U-Probe 0.2.0, pandas 2.1.1, NumPy 1.26.0, PyYAML 6.0.1 and Biopython 1.83.
Additional dependency versions were not recorded exhaustively.

Large generated files, genome FASTA, annotation, indexes and caches are not
committed. The full candidate gzip has SHA-256
`24f6c47450b6e15bba359598e2efc33b7ff2d6c1d5f06e198136f556c23ac23f`.
The local historical candidate artifact is not a public download URL; public
large-artifact distribution is `not_recorded`.

## 6. Sequence recovery results

| Stage | SRR | TRR | MTR | CPSR |
|---|---|---|---|---|
| Unfiltered candidates | 512/512 = 100% | 64/64 = 100% | 100% | 256/256 = 100% |
| Final selection | N/A | N/A | N/A | N/A |

Padlock and primer recovery were each 256/256. No published reference oligo or
pair was unrecovered. Complete DNA equality requires the same gene and role;
CPSR additionally requires both members in the same generated pair. Reference
duplicates cannot inflate recovery; chemical modifications are separate.
All definitions are in [the shared metrics](../sequence_recovery_metrics.md).

The four files in [report/](report/) are byte-for-byte copies of the completed
run's compact outputs. Both the original evaluator and a separate implementation
fully scanned the generated collection. The independent implementation checked
88,542,552 rows and 44,271,276 pair records, all per-gene/role rates and each
reported sequence/pair match. No novel output was labelled a false positive.

## 7. Reproduction instructions

Run on Linux with Python 3.10 and the recorded U-Probe API. From the repository
root, install the repository and recorded recovery dependency versions:

```bash
python -m pip install -e .
python -m pip install -r recovery/requirements.txt
```

Supply genomic FASTA and GTF matching the versions above and these checksums:

| File | Bytes | SHA-256 |
|---|---:|---|
| `GRCm39.fa` | 2762330512 | `b1d94d96ef9ba5ec6ad957b7233a8ed18541d7da1a575f84cf97e31a9407e0dd` |
| `GRCm39.gtf` | 1115389575 | `a8183d575ccf62c08534a3df1e62339238c31f585dc8b62a839ef287703ab3ad` |

The exact historical database acquisition/preprocessing commands
were not recorded; a matching external reference artifact is required for an
exact candidate-count reproduction. Do not substitute a newer annotation and
assume its candidate count will be identical.

Generate candidates without any published binding sequences as input:

```bash
python recovery/starmap_plus_2023_64gene/scripts/run_recovery.py \
  --config recovery/starmap_plus_2023_64gene/config/design.yaml \
  --fasta /absolute/path/GRCm39.fa \
  --gtf /absolute/path/GRCm39.gtf \
  --output /absolute/path/generated/starmap_plus_2023_64gene \
  --workers 8
```

The output includes `candidates.csv.gz`, `candidates_by_gene/`, candidate-site
provenance and a run-details file. Mature-transcript source coordinates are
1-based inclusive. Cross-isoform window deduplication retains only the first
generating isoform, not every possible isoform mapping.

Evaluate the complete generated collection into a fresh output directory:

```bash
python recovery/starmap_plus_2023_64gene/scripts/evaluate_recovery.py \
  --generated /absolute/path/generated/starmap_plus_2023_64gene/candidates.csv.gz \
  --output /absolute/path/generated/starmap_plus_2023_64gene/report
```

Check committed references, compact evidence and metrics without generating or
scanning millions of candidates:

```bash
python recovery/starmap_plus_2023_64gene/scripts/verify_results.py
python -m unittest discover -s recovery/common -p 'test_*.py'
```

The compact verifier establishes report consistency, not a new full candidate
scan. Its optional `--generated` mode enables an independent full scan when the
large output is available. The historical full scan covered all 88,542,552
oligo rows and 44,271,276 pair records; repository curation did not repeat that
large calculation. The portable runner was separately validated on Npy:
51,107 pairs and 102,214 oligo rows, byte-identical after decompression to the
historical run. Its candidate CSV SHA-256 after decompression is
`9669fbeb4e217ad9c296597426d5ca930c85520ea69bc8c348d425b2516f3483`.
Source-site records were also identical; unchanged-input resume succeeded,
while modified-input resume was rejected.
The portable bounded-batch runner may take longer than the historical run,
which additionally parallelized the largest gene and compressed its merge.

## 8. Limitations and next steps

**100% candidate recovery is a permissive search-space coverage upper bound.**
It does not demonstrate that U-Probe selects the original four pairs, reproduce
GC/Tm or specificity filtering, or prove experimental performance of new probes.
The panel's original transcript accessions, computational thresholds and ranking
rules were not reconstructed. The recorded 2024 annotation is compatible with
the published sites, but is not claimed to be the original 2023 design database.

Next, retrieve and predefine original filtering/ranking conditions, freeze them
before comparison, and report filtered-candidate and final-selection recovery
separately. Publish an exact external reference-data package before claiming
fresh-machine full reproduction from publicly downloadable inputs alone.
