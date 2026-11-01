# STARmap PLUS 2023 - Unfiltered Candidate Recovery

## Overview

Zeng et al., *Nature Neuroscience* (2023),
[DOI 10.1038/s41593-022-01251-x](https://doi.org/10.1038/s41593-022-01251-x);
[article and supplementary material](https://pmc.ncbi.nlm.nih.gov/articles/PMC11332722/).

Independently regenerate published full padlock/primer sequences from GRCm39:
`genes -> transcript extraction -> target-site enumeration -> assembly -> exact matching`.
Published binding sequences and complete reference oligos are evaluation inputs
only; the generator does not read them.

## Recovery Results

| Item | Result |
|---|---|
| Reference genes | 64 |
| Published oligos | 512 |
| Published probe pairs | 256 |
| Generated candidate pairs | 44,271,276 |
| SRR | 100% (512/512) |
| TRR | 100% (64/64) |
| MTR | 100% |
| CPSR | 100% (256/256) |
| Final selection | Not evaluated |

**100% recovery refers to unfiltered candidate-space coverage, not successful
final probe selection.** Padlock and primer each recover 256/256 references.
The [four report files](report/) retain the original results unchanged.
The reported evaluation covered 88,542,552 oligo records and
44,271,276 pairs. Identical DNA can occur in different candidate pairs;
reference records are counted once. See [shared metrics](../sequence_recovery_metrics.md).

## Design Parameters & Reference Versions

Source: Supplementary Table 1 (worksheet name has a trailing space), genes
`P4:P67`, padlock names/sequences `Q5:Q260`/`R5:R260`, and primer names/sequences
`Y5:Y260`/`Z5:Z260`. Extracted [oligos](reference_oligos.csv) and
[pairs](reference_pairs.csv) agree: four pairs per gene, no exclusions or duplicates.

Original execution commit: `a76baff9fa79c6ef45308b55877aa507e9c62a27`.
No core algorithm changes. [design.yaml](design.yaml) and [run_recovery.py](run_recovery.py)
use native U-Probe extraction, layout and assembly, iterating isoforms separately.

| Parameter | Value |
|---|---|
| Genome / annotation | GRCm39, GCF_000001635.27 / NCBI RefSeq RS_2024_02 |
| Isoforms | All 308: 157 NM, 140 XM, 4 NR, 7 XR |
| Binding arms / gap | Each arm 19-25 nt / 1-2 nt |
| Target windows | 39-52 nt, step 1, 98 geometries |
| Constraints / deduplication | Uppercase ACGT / gene + window DNA + layout |
| GC/Tm, specificity, ranking, count limit | Not applied |

All isoforms include the 148,183-nt Pcdh9 XR_003950833.2 transcript. RNA geometry
is `padlock_site -> gap -> primer_site`, inferred once from the reference panel,
without supplying reference sites to generation. Gene-specific 3-nt codes and
fixed scaffold remain published values; [design.yaml](design.yaml) is the
assembly specification. Phosphorylation is separate metadata. The current YAML
preserves the executed design parameters; display-only metadata changes account
for the two configuration hashes below.

| Artifact | SHA-256 |
|---|---|
| Source workbook | `6388e90d4ab9e68b8ee72045310c7c69691ed3bef6f2c13dd121aa2813043b42` |
| reference_oligos.csv | `61da326e837e6107db14c7bc3a5c3140240fadcaab9b3d194467ca1d866b7ca0` |
| reference_pairs.csv | `0ca5588358067fc89196d29d2f76e692bd047e85a32c8ae3127dad6a29157eeb` |
| Executed design YAML | `4bff8b65d58e2feb5d218f064b7279ea14d568f13c483b962cf4d72a64bcd608` |
| design.yaml | `db3d3555784580ac79d71d3082774ba492b926a28526e41b26d554f479aa7b81` |
| GRCm39.fa (2,762,330,512 bytes) | `b1d94d96ef9ba5ec6ad957b7233a8ed18541d7da1a575f84cf97e31a9407e0dd` |
| GRCm39.gtf (1,115,389,575 bytes) | `a8183d575ccf62c08534a3df1e62339238c31f585dc8b62a839ef287703ab3ad` |
| Full candidate gzip (not committed) | `24f6c47450b6e15bba359598e2efc33b7ff2d6c1d5f06e198136f556c23ac23f` |

## Reproduction Commands

Use Linux, Python 3.10, the U-Probe API at the recorded execution base, and
matching FASTA/GTF files. The recorded environment used Python 3.10.14 and
U-Probe 0.2.0; additional pinned versions are in [requirements](../requirements.txt).
Install U-Probe and the recorded dependencies from the repository root:

```bash
python -m pip install -e .
python -m pip install -r recovery/requirements.txt
```

Replace `<REFERENCE_FASTA>`, `<REFERENCE_GTF>` and `<OUTPUT_DIR>` with your own
local inputs and output directory, then run these two commands:

```bash
python recovery/starmap_plus_2023_64gene/run_recovery.py \
  --config recovery/starmap_plus_2023_64gene/design.yaml \
  --fasta "<REFERENCE_FASTA>" --gtf "<REFERENCE_GTF>" \
  --output "<OUTPUT_DIR>" --workers 8
python recovery/common/sequence_recovery.py \
  --reference-oligos recovery/starmap_plus_2023_64gene/reference_oligos.csv \
  --reference-pairs recovery/starmap_plus_2023_64gene/reference_pairs.csv \
  --method-version STARmap_PLUS_2023_64gene \
  --generated "<OUTPUT_DIR>/candidates.csv.gz" \
  --output "<OUTPUT_DIR>/report"
```

Generation writes complete candidate pairs, source coordinates and run details.
Coordinates are mature-transcript, 1-based inclusive; window deduplication keeps
the first generating isoform. Evaluation streams all records and reports SRR,
TRR, MTR and CPSR. New generic pair reports use JSON `member_matches`; archived
pair reports retain the original `padlock_matched` and `primer_matched` columns.
The matching rule is identical; the archived results are not regenerated.

To audit a published pair, compare both full DNA/gene/role records with
`sequence_matches.csv`, then require a shared generated pair ID also listed in
`pair_matches.csv`. The committed evidence preserves all linked IDs. These compact reports support
auditing; reproducing them from candidates requires the full generated input.

**External inputs:** published probe references are included, with supplementary
table locations and checksums above. Matching genomic FASTA/GTF must be supplied
from the original data archive. The exact acquisition/preprocessing commands
and a public version-matched download URL are `not_recorded`. Consequently a
fresh-machine end-to-end reproduction cannot be claimed from public repository
files alone; supplying the exact hash-matched database artifacts is still required.
Do not silently substitute a newer annotation and expect identical candidate counts.

## Limitations

Original transcript accessions, GC/Tm and specificity thresholds, and ranking
rules were not reconstructed. The 2024 annotation is compatible with published
sites but is not claimed to be the original 2023 database. A public full-candidate download is also `not_recorded`.
New probes have not been experimentally validated.
Predefine filtering and final-selection rules before evaluating those stages.
