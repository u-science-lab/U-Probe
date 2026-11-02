# PRISM 2025: sequence recovery

Unfiltered candidate recovery for the 64-gene mouse embryo panel of Chang et al., *Nature Biotechnology* (2025), [doi:10.1038/s41587-025-02883-7](https://doi.org/10.1038/s41587-025-02883-7). Final selection is not evaluated.

The official [Supplementary Table 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-025-02883-7/MediaObjects/41587_2025_2883_MOESM4_ESM.xlsx) supplies 64 unique 120-nt padlocks and their corresponding fixed 64-plex backbones. The structure is arm5–G–R–B–Y–arm3 (20 nt each); phosphate is a chemical modification. No 31-plex TC prefix is included.

Reference: mouse mm10/GRCm38, GENCODE M19 (Ensembl 94). Independent applicability checks found all 64 genes and exact antisense 40-nt transcript matches in both this reference and GRCm39/RefSeq. Coverage did not distinguish the assemblies; mm10 was selected with its existing GENCODE annotation. The authors' historical reference version could not be confirmed. These results describe candidate recovery using the specified reference.

The public fixed barcode data describe the required design scaffold; complete published padlocks are not bundled.

Generation uses U-Probe native mRNA extraction and assembly for every annotated isoform, 40-nt windows with step 1, and the fixed official barcode per gene. arm5 is the reverse complement of the first 20 target bases; arm3 is the reverse complement of the last 20. Non-ACGT windows are excluded and identical windows within a gene are deduplicated. No GC/Tm, specificity, ranking or final-count filters are applied. Published arms and coordinates are used only in the separate applicability audit and evaluation.

## Results

Python 3.10.14, pandas 2.1.1 and PyYAML 6.0.1 were used; reference and input checksums are recorded in [generation.json](report/generation.json).

389,656 native candidates recovered 64/64 complete padlocks: SRR 100%, TRR 100% (64/64), MTR 100%; CPSR N/A. No references were excluded or duplicated. There are no unrecovered sequences. [Public evaluation tables](report/) record exact matching by gene, oligo type and full DNA. Final selection: Not evaluated.

## Reproduction

Install U-Probe and recovery dependencies, plus openpyxl for workbook extraction. Download the official workbook as `$WORKBOOK`; set `REFERENCE_DIR` to the matching FASTA/GTF directory and `OUTPUT_DIR` to a new output directory.

```bash
python -m pip install -e .
python -m pip install -r recovery/requirements.txt openpyxl
python recovery/prism_2025_64gene/prepare_inputs.py "$WORKBOOK" "$OUTPUT_DIR/inputs"
python recovery/prism_2025_64gene/run_recovery.py \
  --config recovery/prism_2025_64gene/design.yaml \
  --fasta "$REFERENCE_DIR/mm10.fa" --gtf "$REFERENCE_DIR/mm10.gtf" \
  --output "$OUTPUT_DIR/candidates" --workers 8
python recovery/common/sequence_recovery.py \
  --reference-oligos "$OUTPUT_DIR/inputs/reference_oligos.csv" \
  --generated "$OUTPUT_DIR/candidates/candidates.csv.gz" \
  --method-version PRISM_2025_64gene --output "$OUTPUT_DIR/report"
```

The extractor normalizes Spp1 whitespace, TBX5→Tbx5, MEF2C→Mef2c and Hba-X→Hba-x. Complete published sequences are extracted locally from the workbook rather than redistributed. Full evaluator outputs are reproducible with the command above; committed sequence evidence omits the reference DNA column. Novel candidates are not counted as false positives. CPSR is not applicable to this single-oligo evaluation.
