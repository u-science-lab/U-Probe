# RAEFISH 2025: three-component sequence recovery

Unfiltered candidate recovery for a fixed 64-gene Human endogenous RAEFISH pilot from Cheng et al., *Cell* (2025), [doi:10.1016/j.cell.2025.09.006](https://doi.org/10.1016/j.cell.2025.09.006). Source: [Table S1 (`mmc1.xlsx`)](https://ars.els-cdn.com/content/image/1-s2.0-S0092867425010372-mmc1.xlsx), sheets `Human RAEFISH probe templates` and `94-bit RAEFISH readout oligos`. Perturb-RAEFISH and MERFISH are excluded.

The Human sheet contains 23,352 unique target records; 40 `NegativeControl`-suffix records are excluded, leaving 23,312 endogenous targets, all present in GENCODE v46. Before sequence matching, the pilot takes the first 32 sorted Ensembl IDs in each protein-coding/lncRNA stratum: 64 targets, 192 unique synthesis oligos and 64 three-member sets. The [input audit](report/input_audit.json) records source hashes and denominators.

## Architecture and generation

Reverse padlock synthesis templates are `P5(20)-RNA_binding(30)-P3(30)` (80 nt). Linear priming synthesis templates are `L5(20)-RNA_binding(30)-RCA_priming(10)-tail(47)` (107 nt); their functional prefix is 60 nt. The 150-nt encoding synthesis template is `E5(20)-R1(20)-R2(20)-rc(padlock_binding)(30)-R3(20)-R4(20)-E3(20)`. Readouts follow ascending active bits in the published 94-bit, weight-4 codebook. Functional encoding is the reverse complement of this synthesis template. Fixed regions are declared in [design.yaml](design.yaml).

Reference: official [GRCh38 FASTA](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz) and [GENCODE v46 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz), decompressed before generation. Stable Ensembl IDs select annotation records; version suffixes are stripped. U-Probe's native mRNA extraction enumerates spliced transcript windows with step 1. Predefined 60/61-nt layouts cover two adjacent nonoverlapping 30-nt sites, gap 0/1, and both padlock/priming orders. Only non-ACGT windows and duplicate complete windows within a gene are removed; the first sorted transcript origin is retained. No GC/Tm, specificity, ranking or final-count filters are applied.

Native `UProbeAPI.generate_target_seqs` and `construct_probes` execute the YAML DAG. Padlock and priming depend on one extracted window; encoding depends on padlock binding and the gene codebook. Native functional encoding is reverse-complemented by the DAG to produce the exported synthesis template. The three roles share an origin-specific ID containing target, transcript, interval and orientation. Published binding sequences and coordinates are absent from generator inputs.

## Results

The native DAG generated 1,596,454 candidate sets (4,789,362 oligos). Exact recovery was **192/192 synthesis templates: overall SRR 100%; padlock, priming and encoding SRR each 64/64 (100%)**. TRR was 64/64 (100%), MTR 100%, and CPSR **64/64 (100%)**. No reference records were excluded or duplicated, and no published pilot oligos or sets remained unrecovered. All generated sets passed independent structural checks; all four dependency perturbations passed. Native DAG: **PASS**; no U-Probe core changes were required.

The common evaluator compares complete synthesis templates by target and role, without reverse complementation or partial matching. CPSR requires all three exact members in one native set. [Report tables](report/) retain full DNA and generated set IDs. [Structural checks](report/structure.json) and [four DAG dependency tests](report/dag_checks.json) are separate from sequence recovery.

Final selection: **Not evaluated**. Full Human recovery: **Not evaluated**. These results do not measure experimental performance, specificity or the ability to select the original final probes. Novel generated sequences are not false positives.

## Reproduction

Download the workbook as `WORKBOOK`, and set `REFERENCE_DIR` and `OUTPUT_DIR`. Install U-Probe and recovery dependencies; [environment.json](report/environment.json) records the actual software versions. Create `OUTPUT_DIR` before running the dependency checks.

```bash
python -m pip install -e .
python -m pip install -r recovery/requirements.txt
METHOD_DIR=recovery/raefish_2025_64gene
mkdir -p "$OUTPUT_DIR"
python "$METHOD_DIR/test_dag.py" "$METHOD_DIR/design.yaml" "$OUTPUT_DIR/dag_checks.json"
python "$METHOD_DIR/prepare_inputs.py" --workbook "$WORKBOOK" \
  --gtf "$REFERENCE_DIR/gencode.v46.annotation.gtf" \
  --output-dir "$OUTPUT_DIR/inputs" --audit "$OUTPUT_DIR/input_audit.json"
python "$METHOD_DIR/generate.py" --config "$METHOD_DIR/design.yaml" \
  --targets "$OUTPUT_DIR/inputs/targets.csv" \
  --fasta "$REFERENCE_DIR/GRCh38.primary_assembly.genome.fa" \
  --gtf "$REFERENCE_DIR/gencode.v46.annotation.gtf" \
  --output "$OUTPUT_DIR/candidates" --workers 4
python recovery/common/sequence_recovery.py \
  --reference-oligos "$OUTPUT_DIR/inputs/reference_oligos.csv" \
  --generated "$OUTPUT_DIR/candidates/candidates.csv.gz" \
  --method-version raefish_2025_64gene --output "$OUTPUT_DIR/report"
python "$METHOD_DIR/audit_structure.py" \
  --generated "$OUTPUT_DIR/candidates/candidates.csv.gz" \
  --codebook "$OUTPUT_DIR/inputs/codebook.csv" \
  --readouts "$OUTPUT_DIR/inputs/readouts.csv" \
  --output "$OUTPUT_DIR/report/structure.json"
```
