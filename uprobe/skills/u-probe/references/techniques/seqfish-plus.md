# seqFISH+ primary probes

Read this when a user asks for seqFISH+ (Eng et al. 2019) or seqFISH-style barcoded primary probes:
28-nt RNA-binding regions carrying four 15-nt readout sites, 24 probes per gene, decoded by
pseudocolour rounds. Also covers how seqFISH, RNA SPOTs and later Cai-lab designs differ.

## 1. Method and when to choose it

- Primary probes (unlabelled) tile each mRNA; each carries the readout sites of the gene's barcode.
  Dye-labelled 15-nt readout probes are hybridized, imaged and stripped in serial rounds.
- seqFISH+ codebook: 60 pseudocolours = 3 channels x 20 pseudocolours; per channel 20^3 = 8,000 codes
  plus one parity round (4 barcoding rounds, 80 hybridizations). Genes per channel: 3,333-3,334
  (10,000 genes). Only 1/20 of a channel's genes is lit per image, which "dilute[s] mRNA molecules
  into 60 separate images" and allows ~35,000 transcripts per NIH/3T3 cell.
- Choose it for ~1,000-10,000 genes with sub-cellular resolution in cells or thin (5-20 um) tissue,
  with an automated fluidics setup (80+ hybridizations). For <~100 genes use smFISH/HCR rounds;
  for error-robust binary codes with 2-3 readouts per probe see `merfish.md`.

## 2. Probe architecture (5' -> 3')

Final primary probe after library amplification (~93 nt, Eng 2019):

```
5'-[ I 15 ]s[ II 15 ]s[ target-binding 28 = rc(mRNA) ]s[ III 15 ]s[ IV 15 ]CCC-3'
    rc(readout probe, round I..IV)                                         ^ T7 +1..+3 (GGG -> CCC)
s = 1-nt spacer (T in all 15,177 Lohoff 2022 probes, which reuse seqFISH+ readouts and primer;
    its Methods text says "2-nucleotide spacer"; Eng 2019 gives no spacer)
```

Oligo-pool oligo as ordered (before amplification; Lohoff 2022 Supplementary Table 1A, 135 nt with
30-nt target regions):

```
5'-[forward primer 21: GTCTTATCATGTGCCTTTCCT]-[probe body as above]-[rc(T7 promoter) 20: CCCTATAGTGAGTCGTATTA]-3'
   removed: RT primer carries a uracil, cut by USER enzyme
```

Readout probe: 15 nt, 5' amine coupled to NHS-ester dye (Alexa 647, Cy3B, Alexa 488).

- Order I-II-target-III-IV is drawn in Eng 2019 Fig. 1a; 5'/3' ends are not marked there. Lohoff 2022
  (135-nt pool oligos) confirms 5'-primer-R-T-R-T-target-T-R-T-R-rcT7-3' and rotates which round
  sits in which slot across a gene's probes (R1R2|R3R4, R2R3|R4R1, R3R4|R1R2, R4R1|R2R3, cycling in
  the row order of its table, which is not the order along the transcript).
- Length check: 4 x 15 + 28 = 88 nt; ~93 nt observed = spacers + the 3' CCC left by T7 transcription
  (inference; Eng 2019 does not explain the extra ~5 nt).

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Binding length | 28 nt | | Eng 2019 Methods, doi:10.1038/s41586-019-1049-y |
| Region | CDS exons; CDS + UTR exons if too few probes | masked UCSC genome | Eng 2019 |
| GC | 45-65 % | | Eng 2019 |
| Homopolymer | drop runs of >= 5 identical bases | | Eng 2019 |
| Tm window | none given | Eng 2019 filters GC only | Eng 2019 |
| Off-target (per probe) | BLAST 15-nt match to another gene = hit; drop probe if summed expected off-target copy number > 10,000 FPKM | ENCODE mouse RNA-seq | Eng 2019 |
| Probe-probe cross-hyb. | drop probes sharing >= 17 nt (from the larger probe set) | | Eng 2019 |
| Spacing | not stated; RNA SPOTs: no two probes within 2 nt | | Eng 2017, doi:10.1038/nmeth.4500 |
| Probes per gene | 24; genes with < 24 dropped | | Eng 2019; Cai protocol, doi:10.1038/protex.2019.019 |
| Readout sites | 4 x 15 nt (one per barcoding round) | rc of readout probes | Eng 2019 |
| Readout design | random 15-mers, GC 40-60 %, BLAST vs transcriptome, no 10-nt contiguous match between readouts; 240 readouts | | Eng 2019 |
| Primary hyb. | 40 % formamide, 2x SSC, 10 % dextran sulfate, ~1 nM per probe, 2 uM LNA polyT15/T30 | 37 C, 36-48 h cells, ~60 h tissue | Eng 2019; Cai protocol |
| Primary wash | 40 % formamide in 2x SSC | 37 C, 30 min | Eng 2019 |
| Readout hyb. | 50 nM (A647, Cy3B) / 100 nM (A488) in 10 % ethylene carbonate, 10 % dextran sulfate, 4x SSC, SUPERase In | RT, 17 min cells / 20 min tissue | Eng 2019 |
| Readout wash | 10 % formamide, 0.1 % Triton X-100, 2x SSC | RT | Eng 2019 |
| Strip | 55 % formamide, 0.1 % Triton X-100, 2x SSC | 1 min flow + 1 min | Eng 2019 |

Published seqFISH-family designs disagree; do not mix their numbers:

| | RNA SPOTs 2017 | seqFISH+ 2019 | Shah 2016 / 2018 (seqFISH, intron) | Lohoff 2022 (embryo seqFISH) |
| --- | --- | --- | --- | --- |
| Binding length | 25 nt | 28 nt | 35 nt | 30 nt |
| GC | 45-70 % | 45-65 % | ~55 % | 45-65 % |
| Probes per gene | 28-32 | 24 | 24 / 25 | >= 28 (barcoded) |
| Off-target limit | 10,000 FPKM, <= 6 probes per off-target gene | 10,000 FPKM | 5,000 copies / <= 7 hits | <= 7 hits |
| Readouts | 20 nt, 4 of 5 rounds per probe | 15 nt, all 4 rounds | HCR (2016); 15 nt + 3-nt spacer (2018) | 15 nt, all 4 rounds |
| Primary formamide | 30 % | 40 % | 50 % (2018) | 40 % |
| Primers in final probe | kept (149-mer) | removed (~93 nt) | | |

## 4. Target region

- CDS exons first (Eng 2019); U-Probe template uses the CDS of the MANE Select transcript.
- Isoforms: Eng 2019 says "designed as previously described" (RNA SPOTs: "consensus regions of all
  spliced isoforms"). U-Probe has no shared-region option; the reference transcript's CDS contains
  the constitutive CDS but may include isoform-specific exons. For strict isoform-agnostic counting,
  drop probes whose `transcript_names` do not cover all isoforms (run with `transcripts: all` to see).
- Genes short of 24 CDS probes: Eng 2019 added UTR exons. U-Probe's `source` is global, so run those
  genes in a second protocol with `source: mRNA`.
- Introns (nascent transcription, intron seqFISH, Shah 2018): 35-nt probes on constitutive introns
  near the 5' end; needs `source: genome` regions, not covered by this template.

## 5. Specificity

- A spot needs several of the 24 probes; a single off-target probe adds diffuse background, not a
  decoded barcode. But in barcoded methods a false spot also needs the right readout pattern, so the
  real risk is an off-target RNA bound by many probes of one gene, or a highly expressed off-target
  RNA that raises background in every round. Hence Eng 2019 weights off-targets by expression
  (summed FPKM > 10,000 -> drop), which also removes probes against housekeeping/ribosomal families.
- Readout orthogonality matters more than in smFISH: a readout that binds another gene's overhang or
  the transcriptome corrupts every barcode using it.
- Probe-probe complementarity (>= 17 nt) between genes would bridge barcodes; Eng 2019 removed it.
- Copies of the target: the template exempts non-protein-coding genes whose transcript holds half the target's 30-mers
  (`exempt_copies: {kmer: 30, min_shared: 0.5}`, as `smfish_stellaris`); they are listed in
  `target_mapped_genes_copies`. Mouse Gapdh, whose retrocopies (Gm13394 and others) are expressed in
  NIH/3T3, kept 2 of 24 probes without it and 24 with it.

## 6. U-Probe template: `seqfish_plus`

```bash
uprobe templates show seqfish_plus -o protocol.yaml
```

```yaml
encoding:
  GAPDH: {readouts: [NNNNNNNNNNNNNNN, NNNNNNNNNNNNNNN, NNNNNNNNNNNNNNN, NNNNNNNNNNNNNNN]}   # rounds I-IV
extracts:
  target_region: {source: CDS, transcripts: mane, length: 28, step: 1}
probes:
  primary_probe:   # 133-nt pool oligo, Lohoff 2022 layout
    template: "{primer_fwd}{overhang_1}T{overhang_2}T{binding}T{overhang_3}T{overhang_4}{t7_rc}"
    # primer_fwd GTCTTATCATGTGCCTTTCCT, t7_rc CCCTATAGTGAGTCGTATTA;
    # overhang_j = rc(encoding[target]['readouts'][(probe_index + j - 1) % 4])
hybridization: {temperature_c: 37, na_mM: 390, probe_nM: 1, formamide_pct: 40, dntp_mM: 0}
post_process:
  filters:
    target_gc: {condition: target_gc >= 0.45 & target_gc <= 0.65}
    homopolymers: {type: sequence_pattern, target: target_region, exclude_patterns: [AAAAA, TTTTT, CCCCC, GGGGG]}
    target_mapped_genes: {condition: target_mapped_genes <= 1}   # Bowtie 1, <= 2 mismatches
  remove_overlap: {location_interval: 0}
  equal_space: {number_desired: 24}
```

What to edit:
- `encoding`: per gene, `readouts`: the readout PROBE sequences of rounds I-IV, in round order, from
  your codebook (the template inserts their reverse complements). Placeholders `N` must be replaced.
- `primer_fwd`: the template uses the seqFISH+ forward primer of Lohoff 2022; use your library's.
  Keep the rc(T7) tail if you amplify by T7 transcription (else replace it with your reverse primer
  site). The template reproduces the first two published Lohoff 2022 Abcc4 probes exactly (with
  `length: 30`).
- Slot rotation: `probe_index` numbers a gene's final probes 0, 1, 2, ... along the RNA, and the
  rounds shift by one slot per probe, as in Lohoff 2022: 12|34, 23|41, 34|12, 41|23, then again; every
  round sits in every slot on a quarter of the probes. The final TSV has a `probe_index` column.
- Add `gene_expression` (your sample's table, e.g. TPM >= 1) so unexpressed off-target genes do not
  count; this is the closest U-Probe gets to the expression-weighted 10,000-FPKM rule.
- Optional: `remove_overlap: {location_interval: 2}` to follow the RNA SPOTs spacing.
- `genome` needs `bowtie` in `align_index` (20-30-nt windows; Bowtie2 misses many hits).

Template departures from Eng 2019 (from the template `meta.limitations` and the literature):
- Off-targets: Bowtie 1 within 2 substitutions, unweighted by expression, vs BLAST 15-nt matches
  weighted by FPKM. Probe-probe 17-nt filter: not modelled.
- Spacers, primer, tail and slot rotation follow Lohoff 2022, which reuses the seqFISH+ readouts and
  primer; the Eng 2019 library sequences are not public. The rotation follows the RNA, Lohoff's the
  row order of its table. Binding regions are 28 nt (Eng 2019), not Lohoff's 30.
- CDS only, no automatic UTR fallback; spacing 0 nt.
- `binding_tm` is reported only (no published Tm window).

## 7. Pitfalls

- Strand: probes = rc(mRNA); overhangs = rc(readout probe). Do not insert readout probe sequences
  directly.
- Highly expressed genes: Eng 2019 picked 10,000 genes "while avoiding highly abundant housekeeping
  genes, such as ribosomal proteins". Remove or image separately genes that dominate the counts;
  Shah 2016 imaged its 25 most abundant genes individually, Lohoff 2022 used non-barcoded smFISH for 36.
- Paralogs with high expression lose most probes under the off-target filter; check `summary.json`
  `final` per gene and drop genes with < 24 probes (Eng 2019) rather than accepting thin sets.
- Codebook: Eng 2019 round 4 = (R1 + R2 - R3) mod 20 (0 written as 20), minimum Hamming distance 2,
  decoded when >= 3 of 4 rounds agree. U-Probe has no codebook generator; build it externally.
- Detection efficiency is ~49 % of smFISH (Eng 2019); do not cut probes per gene below 24.

## 8. Checklist

- [ ] Every gene has 24 probes (or is dropped / redesigned with UTRs); `summary.json` status `ok`.
- [ ] Binding regions 28 nt, GC 45-65 %, no 5-nt homopolymer; `target_mapped_genes` = 1 or explained
      by `_unexpressed` / `_copies`.
- [ ] Overhangs are rc(readout probes); each gene has exactly its codebook's 4 readouts; each channel's
      readouts unique; readouts share no >= 10-nt contiguous match with each other.
- [ ] No two probes of different genes share >= 17 nt (k-mer check on the final table).
- [ ] Pool oligo = primer + 4 overhangs with T spacers + target + rc(T7) (133 nt with 28-nt targets);
      each gene's rounds rotate through the 4 slots (`probe_index`).
- [ ] Highly expressed genes excluded or routed to non-barcoded rounds.
- [ ] Run last, with the published values (Eng 2019 filters GC only; spacing is not stated, windows must not overlap):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --summary results/run1/summary.json \
    --transcripts <genome dir>/transcript.fa --gc 0.45 0.65 --min-gap 0 --min-probes 24
```

## 9. References

- Eng C-HL et al. 2019. Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+. Nature 568:235-239. doi:10.1038/s41586-019-1049-y
- Eng C-HL, Cai L. 2019. RNA seqFISH+ supplementary protocol. Protocol Exchange. doi:10.1038/protex.2019.019
- Eng C-HL et al. 2017. Profiling the transcriptome with RNA SPOTs. Nat Methods 14:1153-1155. doi:10.1038/nmeth.4500
- Lubeck E et al. 2014. Single-cell in situ RNA profiling by sequential hybridization. Nat Methods 11:360-361. doi:10.1038/nmeth.2892
- Shah S et al. 2016. In situ transcription profiling of single cells reveals spatial organization of cells in the mouse hippocampus. Neuron 92:342-357. doi:10.1016/j.neuron.2016.10.001
- Shah S et al. 2018. Dynamics and spatial genomics of the nascent transcriptome by intron seqFISH. Cell 174:363-376. doi:10.1016/j.cell.2018.05.035
- Lohoff T et al. 2022. Integration of spatial and single-cell transcriptomic data elucidates mouse organogenesis. Nat Biotechnol 40:74-85. doi:10.1038/s41587-021-01006-2
