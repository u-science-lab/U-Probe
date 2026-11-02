# split-FISH 2020: 317-gene candidate recovery

Source: Goh et al., *Nature Methods* (2020), [paper](https://doi.org/10.1038/s41592-020-0858-0), Supplementary Tables 1–3 (`MOESM3`, `MOESM4`, `MOESM5`). The [bridge correction](https://doi.org/10.1038/s41592-020-01019-w) removes an erroneous extra G. All 26 corrected 18-nt cores agree with the [author design inputs](https://github.com/KHChenLab/split-fish/tree/master/probedesign) and published templates. The 26-bit weight-2 codebook contains 317 targets and eight blank controls without reference oligos.

Complete DNA templates, 5′→3′:

```text
P1: PCR-F(20)–GGTAC–rc(bridge[0:9])–TAT–rc(site[27:52])–GAATTC–PCR-R-tail(19)
P2: PCR-F(20)–GGTACC–rc(site[0:25])–TAT–rc(bridge[9:18])–AATTC–PCR-R-tail(19)
PCR-F: AACGAACGGAGGGTCATTGG
PCR-R-tail: CCCTAATGCGAGTAGAGCCT
```

Each template is 88 nt; target bindings are 25 nt with a 2-nt gap. Per gene and ascending on-bit, the reference table stores 36 P1 followed by 36 P2. Same-ordinal pairing was independently verified on the stated GENCODE M4 transcript for all 22,824 pairs; none were excluded. This audit is separate from candidate generation.

| Measure | Candidate Recovery |
|---|---|
| Generated complete templates | 12,367,668 |
| Ground truth | 317 genes / 45,648 oligos / 22,824 pairs |
| SRR | 100% (45,648/45,648) |
| TRR | 100% (317/317) |
| MTR | 100% |
| CPSR | 100% (22,824/22,824) |
| Final Recovery | Not evaluated |

Native `UProbeAPI.generate_target_seqs` enumerates all 52-nt windows (step 1) from every M4 transcript using mm10/GRCm38 genomic FASTA. Native layout extraction and `construct_probes` execute the [YAML/DAG](design.yaml), including reverse complements and complete assembly. Both published on-bits are applied to every site. Only non-ACGT and duplicate within-gene windows are removed. The generator never reads published binding DNA, full oligos, transcript selections or inferred coordinates. Pair IDs preserve native window IDs, namespaced by transcript and bit; both roles remain in the same native construction row.

The four [standard tables](report/) use the unchanged [common evaluator](../common/sequence_recovery.py) and [metric definitions](../sequence_recovery_metrics.md). [Provenance](report/provenance.json) records inputs, hashes, versions and generation counts. No oligos or confirmed pairs were unrecovered. Candidate recovery does not reproduce author specificity filters, random code assignment, alternating site selection, MATLAB ranking or final probe selection, and does not establish experimental performance. The existing genomic FASTA's original download URL is unrecorded; its SHA-256 is retained. Annotation is GENCODE M4, not the shared newer annotation.

## Reproduction

Install U-Probe and `recovery/requirements.txt`, plus `openpyxl`. Download the three workbooks and author `BridgeSequences.xlsx` using the provenance URLs. Set `$TABLE1`, `$TABLE2`, `$TABLE3`, `$AUTHOR_BRIDGES`, `$FASTA` (mm10 genomic FASTA), `$GTF` (decompressed GENCODE M4 annotation), and a fresh `$OUT`.

```bash
METHOD=recovery/splitfish_2020_317gene
python "$METHOD/prepare_inputs.py" --templates "$TABLE1" --codebook "$TABLE2" \
  --bridges "$TABLE3" --author-bridges "$AUTHOR_BRIDGES" \
  --fasta "$FASTA" --gtf "$GTF" --output "$OUT/inputs"
python "$METHOD/generate.py" --config "$OUT/inputs/design.yaml" \
  --fasta "$FASTA" --gtf "$GTF" --output "$OUT/candidates" --workers 8
python recovery/common/sequence_recovery.py \
  --reference-oligos "$OUT/inputs/reference_oligos.csv" \
  --reference-pairs "$OUT/inputs/reference_pairs.csv" \
  --generated "$OUT/candidates/candidates.csv.gz" \
  --method-version splitfish_2020_317gene --output "$OUT/report"
```
