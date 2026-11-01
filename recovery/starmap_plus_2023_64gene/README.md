# STARmap PLUS 2023: Sequence Recovery

## Overview

Independent sequence recovery of the mouse (*Mus musculus*) STARmap PLUS
64-gene panel published by Zeng et al., *Nature Neuroscience* (2023)
([publication and supplementary data](https://doi.org/10.1038/s41593-022-01251-x)).

## Results

| Measure | Result |
|---|---|
| Published reference panel | 64 genes; 512 oligos; 256 padlock–primer pairs |
| Generated candidate pairs | 44,271,276 |
| SRR | 100% (512/512 oligos) |
| TRR | 100% (64/64 genes) |
| MTR | 100% |
| CPSR | 100% (256/256 pairs) |
| Final selection | Not evaluated |

**These results represent unfiltered candidate-space coverage, not final probe
selection.** [Summary and per-gene results](report/) include sequence-level and
pair-level matching evidence. Metric definitions are provided in the
[shared specification](../sequence_recovery_metrics.md).

## Data and Methods

- **Reference genome:** mouse GRCm39, assembly `GCF_000001635.27`;
  NCBI RefSeq annotation `RS_2024_02`.
- **Published probes:** Supplementary Table 1; included as
  [reference_oligos.csv](reference_oligos.csv) and
  [reference_pairs.csv](reference_pairs.csv), with four pairs per gene.
- **Design:** [design.yaml](design.yaml) specifies the published scaffold,
  gene-specific 3-nt barcodes and assembly order. All 308 annotated isoforms
  are processed separately. Binding arms span 19–25 nt each, separated by
  1–2 nt; target windows span 39–52 nt with a 1-nt step. Only ACGT windows
  are retained, deduplicated by gene, window sequence and layout.
- **Generation:** [run_recovery.py](run_recovery.py) uses U-Probe transcript
  extraction and probe assembly APIs. Published target-binding sequences and
  complete oligos are used only for evaluation, never as generation inputs.
- **Evaluation:** [sequence_recovery.py](../common/sequence_recovery.py)
  matches complete DNA sequences exactly by gene and oligo type, treats chemical
  modifications separately and counts each reference once. Pair recovery requires
  both corresponding oligos in the same generated pair. GC/Tm, specificity,
  ranking and final-count filters are not applied.

## Reproduction

From the repository root, install U-Probe and the recorded dependencies
(Python 3.10.14; U-Probe 0.6.2):

```bash
python -m pip install -e .
python -m pip install -r recovery/requirements.txt
```

Generate candidates from the reference genome and annotation, then evaluate:

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

## Reproducibility Notes

Results were generated using GRCm39 with NCBI RefSeq RS_2024_02 annotation. Exact reproduction requires matching reference FASTA/GTF files and a compatible U-Probe version. Reference data acquisition instructions are not yet provided. The reported recovery rates reflect unfiltered candidate-space coverage, not final probe selection.
