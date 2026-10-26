# MERFISH encoding probes, readouts and codebooks

Read this when a user asks for MERFISH (or MERSCOPE-style) probes: encoding probes with 20/30-nt
readouts, Hamming codebooks (MHD4, MHD2), primers for oligo-pool amplification, probes per gene,
and which genes a panel can hold before optical crowding breaks decoding.

## 1. Method and when to choose it

- Unlabelled encoding probes tile each RNA; each carries 2-3 readout sequences. Every gene gets a
  binary code word of N bits; bit i is "on" for a gene when its encoding probes carry readout i.
  Round i hybridizes the dye-labelled readout probe i, images, then removes the dye.
- Codes: constant Hamming weight 4 (four "on" bits per gene, equal brightness and equal 1->0 risk).
  MHD4 = minimum Hamming distance 4: 1-bit errors corrected, 2-bit errors detected.
  Unused code words are kept as blank controls (false-positive rate).
- Choose it for 100-10,000 genes, cell typing in tissue, subcellular RNA maps. Do not use it for
  1-5 genes (see `smfish.md`, `hcr-v3.md`) or for genes expressed at hundreds of copies per cell
  in the combinatorial set (section 5).

## 2. Probe architecture (5' -> 3')

Moffitt 2016 / ZhuangLab `library_design_example.m` (the U-Probe template, 132 nt):

```
5'-[primer 1, 20][A][rc(R1), 20][rc(RNA), 30][A][rc(R2), 20][rc(R3), 20][rc(primer 2), 20]-3'
alternative placement:  [P1][A][rc(R1)][rc(R2)][rc(RNA) 30][A][rc(R3)][rc(P2)]
```

- R1-R3 are 3 of the gene's 4 readout PROBE sequences (dye-labelled 20-mers), chosen at random per
  probe, with 1 or 2 readouts 5' of the target region (random) (Moffitt 2016 SI: "three readout
  sequences that were randomly selected from the four readout sequences associated with the gene").
  The probe carries their reverse complements. Single A spacers avoid terminal G triplets / G
  quadruplets. All 4 readouts must appear across a gene's probes: with only 3, every molecule reads
  as a 1-bit error and loses its error correction.
- Primer 1 is the PCR forward and RT primer; primer 2 carries a T7 promoter (TAATACGACTCACTATAGGG) at
  its 5' end, and the probe holds its reverse complement at the 3' end (library: PCR -> in vitro
  transcription -> RT -> alkaline hydrolysis).

Other published layouts:

| Paper | Readouts per probe | Layout | Probe |
| --- | --- | --- | --- |
| Chen 2015 Science | 2 x 30 nt, one each side | [P1 20][R 30][target 30][R 30][rc P2 20] | 130 nt |
| Moffitt 2016 PNAS | 3 x 20 nt (3 of 4) | above | 132 nt |
| Wang 2018 Sci Rep | 2 x 20 nt | [P1][R 20][target 30][R 20][rc T7]; P1 removed by USER | 72 nt (final) |
| Zhang 2021 Nature (22-bit) | 2 x 20 nt | readout on each side of the target | |
| Xia 2019 PNAS (10,050 genes) | 3 x 20 nt | [20][A][R][A][R][A][target 30][A][R][20] | 134 nt |
| Xia 2019 Sci Rep (bDNA) | 3 x 20 nt (3 of 4) | Moffitt 2016 library unchanged (92 or 16 per RNA); readout sites bind amplifiers (section 6) | as Moffitt 2016 |

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Target region | 30 nt | all MERFISH papers | Chen 2015; Moffitt 2016 |
| GC | 43-63 % | Xia 2019: 33-73 %; Blosser 2021: 35-75 % | Moffitt 2016 SI, doi:10.1073/pnas.1612826113 |
| Tm | 66-76 C | SantaLucia NN, 300 mM monovalent salt, 5 nM probe, no formamide, by the design code's formula (`TRDesigner.m` `GetRegionTm`: strand term ln C_T with C_T = 5 nM, not C_T/4), which is U-Probe's `tm_method: merfish_design`; Primer3 at `probe_nM: 5` reads 1.36 C lower (+1.27 to +1.50, random 30-mers) | Moffitt 2016 SI; TRDesigner.m |
| Tm (other papers) | Chen: > 70 C, off-target < 72 C, no secondary structure > 76 C (OligoArray); Xia: 61-81 C | | Chen 2015; Xia 2019 |
| Homopolymer | none in Moffitt 2016; Chen: no runs >= 6 | | Chen 2015 |
| Specificity | isoform- and gene-specificity index 0.75-1 (17-mer tables, expression-weighted) | Blosser 2021: isoform 0-1 | Moffitt 2016 |
| rRNA / tRNA | no homology >= 15 nt (code; SI says "> 15") | Chen: > 14 nt to rRNA/tRNA, > 17 nt to genes > 10,000 FPKM | library_design_example.m L195 |
| Tiling | non-overlapping, from the 5' end of each isoform | Wang 2018: up to 20-nt overlap, 92 regions picked at random per RNA; Xia 2019: 20-nt overlap only for the 1,000 genes < 1,440 nt or short of targetable regions | Moffitt 2016; Wang 2018 Methods; Xia 2019 PNAS Results |
| Probes per gene | 92 | Chen ~192 (140 genes) / ~94 (1001 genes); MIE >= 50; Blosser 64; Xia <= 48 | Moffitt 2016; Zhang 2021 |
| Readouts | 20 nt; A/T/G only (G 50 %), GC 40-50 %, no G run > 3, <= 11 nt cross-homology and <= 11 nt homology to the transcriptome | Chen 2015: 30 nt, 4-letter | Moffitt 2016 SI |
| Primers | 20 nt, Tm 70-72 C, GC 50-65 %, no run >= 4, self-compl. <= 6 nt, cross-homology <= 8 nt | Chen/MIE: Tm 70-80 C, GC clamp | Moffitt 2016; PrimerDesigner.m |
| Encoding hyb. | 2x SSC, 30 % formamide, 10 % dextran sulfate, 0.1 % yeast tRNA, 1 % murine RNase inhibitor | 37 C, 36-48 h | Moffitt 2016 |
| Encoding probe conc. | ~50-200 uM total (2015-2018): 100 uM (140 genes) / 200 uM (1,001 genes) (Chen 2015); ~200 uM (Moffitt 2016, PNAS 113:11046); 100 uM (Moffitt 2016, PNAS 113:14456); ~50 uM (Wang 2018, corrected from ~300 uM). Later ~1 nM per probe (Xia 2019, Zhang 2021) | titrate per probe batch (Moffitt 2016); MIE: 10-200 uM by probe-set size, titrate over 5-200 uM. PNAS 113:14456 used 300 uM only as a 3x high-background test | Chen 2015; Moffitt 2016 (both); Wang 2018 Author Correction, doi:10.1038/s41598-018-24844-8; MIE; Xia 2019 SI; Zhang 2021 |
| Encoding wash | 2x SSC, 30 % formamide | 47 C, 2 x 30 min | Moffitt 2016 |
| Readout hyb. | 1 nM each, 2x SSC, 1 % ethylene carbonate, 10 % dextran sulfate, 2 mM VRC; wash 2x SSC + 10 % EC | room temperature | Moffitt 2016 |
| Readout hyb. (later) | 3 nM, 2x SSC, 10 % EC, 0.1 % Triton (Zhang 2021); 5 % EC (Wang 2018); 10 nM in 10 % formamide at 37 C (Chen 2015) | | as cited |
| Dye removal | 50 mM TCEP cleaves the disulfide-linked dye, 15 min | Chen 2015 used photobleaching | Moffitt 2016 |

Published readout probes: Moffitt 2016 Table S1 (16 x 20 nt, RS0015 ... RS0548); Xia 2019 Dataset S3
(69 x 20 nt; bits 1-16 = Moffitt 2016); Chen 2015 / MIE Table 1 (16 x 30 nt).

## 4. Codebook

| Code | Bits | Code words | Used for |
| --- | --- | --- | --- |
| MHD4, weight 4 | 16 | 140 | 130 genes + 10 blanks (Chen 2015, Moffitt 2016) |
| MHD4, weight 4 | 22 | 385 | 242 genes + 10 blanks (Zhang 2021); Blosser 2021 |
| MHD4 (HD4), weight 4 | 69 | 12,903 | 10,050 genes (Xia 2019, with expansion) |
| MHD2, weight 4 | 14 | 1,001 | 985 genes + 16 controls (Chen 2015; detects, cannot correct errors) |

- Bound: a weight-4, distance-4 code has at most C(N,3)/4 words (no two words may share 3 on-bits);
  16 and 22 bits reach it (140, 385).
- Leave 5-10 % of code words unassigned as blanks (Moffitt & Zhuang 2016, MIE).
- Readouts: N distinct readout probes, one per bit; channels split bits across 2-3 dyes per round.
- U-Probe has no codebook generator (`uprobe generate-barcodes` makes DNA sequences, not binary
  codes). Use a published codebook or build one externally, then map each gene's 4 on-bits to their
  readout probe sequences.

## 5. Abundance limits and optical crowding

- Each round lights all genes whose bit is on, so the spot density per image scales with the summed
  expression of the combinatorial panel. Chen 2015: ~2-3 molecules/um^3 resolvable per round, ~20
  molecules/um^3 with MERFISH, enough for "all but the top 10 % most expressed genes" of IMR90.
- Wang 2018: ~3-4 molecules/um^2 halves detection efficiency; a library 14x more abundant than earlier
  panels was detected at 16-21 % without expansion and ~100 % after 12-fold volume expansion.
- Practical thresholds: Moffitt 2016 genes 0.1-1,000 FPKM; Xia 2019 0.01-500 FPKM; Zhang 2021 routed
  genes with >= 200 counts per cell in any cluster (or < 48 probes) to sequential rounds; Moffitt 2018
  imaged 20 short and/or highly expressed genes sequentially. Vizgen guidance: each gene < 800 FPKM,
  panel total ~9,000 FPKM, abundant genes by sequential imaging (a 2025 benchmark quotes 900).
- Agent procedure: get a reference expression table for the tissue, flag genes above the per-gene
  limit, sum the rest, and move flagged genes to sequential (non-barcoded) rounds or drop them.

## 6. Target region and probes per gene

- Whole mRNA incl. UTRs on the most abundant isoform (Chen 2015; Ensembl cDNA in MIE). U-Probe template:
  MANE Select mRNA. Probes may hit other isoforms; isoform specificity is not modelled by U-Probe.
- 92 non-overlapping 30-mers need ~2.8 kb ("92 encoding probes for RNAs that are 3 kb or longer").
  With up to 20-nt overlap, RNAs down to ~1.2 kb carry 92 (Wang 2018), ~500 nt carry 48 (Xia 2019
  PNAS). Xia 2019 PNAS kept 48 non-overlapping probes for the 9,050 genes >= 1,440 nt and overlapped
  only the other 1,000; these read at a median 68 % of the FPKM-normalized copy number of the rest.
  In U-Probe, `remove_overlap: {location_interval: -20}` lets kept windows overlap by up to 20 nt (a
  start every >= 10 nt) for every target of the protocol: run short genes in a second protocol with
  -20 and `equal_space: {number_desired: 48}`, the rest with 0.
- Overlap is the exception: use the published non-overlap rule (`location_interval: 0`) and allow overlap
  only for a transcript too short to carry the requested probes without it, as Xia 2019 did for genes
  < 1,440 nt (Wang 2018 overlapped up to 20 nt for every RNA). Report the departure and keep the
  non-overlapping probes of such a gene where the count is reached.
- Fewer probes lower brightness and detection; genes with very few probes go to sequential rounds.
- bDNA amplification (Xia 2019 Sci Rep) changes no encoding probe: a primary amplifier ends (3') in the
  readout probe's sequence (bits 1-4 of its Table S1 are the template's readouts) after 5 x (unique
  20-mer + 2 nt); a secondary ends in rc(that 20-mer) after 5 repeats of its own 20-mer, which a new
  dye probe reads (5 x 5, 130-nt oligos; 9 x 9, 200 nt). Repeats use 3 letters (primary A/T/C,
  secondary A/T/G, C or G at 50 %), screened for <= 11-nt transcriptome homology and then in cells:
  16 of 20 pairs worked (G-rich secondaries bound cellular RNA, some pairs amplified poorly). Signal
  5.5x / 10.5x / 30.6x (4 x 4 / 5 x 5 / 9 x 9); with 16 probes per RNA, 61 % of the 92-probe counts
  vs 33 % unamplified, so RNAs >= 480 nt fit. U-Probe designs the encoding probes (`equal_space:
  {number_desired: 16}`), not the amplifiers; `uprobe generate-barcodes --alphabet ACT --length 20`
  only drafts orthogonal 3-letter 20-mers, without that composition or a transcriptome screen.

## 7. Specificity

- A false call needs the off-target RNA to show the gene's whole barcode, i.e. many probes of one gene
  on one off-target RNA, so a few promiscuous probes cost background, not false genes. Expression-
  weighted specificity (17-mers) and the rRNA/tRNA screen matter most: abundant RNAs bound by a few
  probes add background in every round that uses those readouts.
- Readout orthogonality is critical: a readout with homology to the transcriptome, to another readout
  or to a target region lights wrong bits for every gene. Keep cross-homology <= 11 nt.
- Copies of the target (processed pseudogenes, retrocopies) hold most of its sequence,
  so no probe set avoids them; the template exempts them unless they are protein_coding (`exempt_copies: {kmer: 30, min_shared: 0.5}`,
  as `smfish_stellaris`): a gene whose transcript holds half the target's 30-mers is listed in
  `target_mapped_genes_copies` instead of counted. Name the copies in the report; a copy that is
  expressed contributes signal.

## 8. U-Probe template: `merfish`

```bash
uprobe templates show merfish -o protocol.yaml
```

```yaml
encoding:   # readout PROBES of the gene's 4 on-bits (Moffitt 2016 Table S1, bits 1-4: RS0015, RS0083, RS0095, RS0109)
  GAPDH: {readouts: [ATCCTCCTTCAATACATCCC, ACACTACCACCATTTCCTAT, ACTCCACTACTACTCACTCT, ACCCTCTAACTTCCATCACA]}
extracts:
  target_region: {source: mRNA, transcripts: mane, length: 30, step: 1}
probes:
  encoding_probe:
    template: "{primer_1}A{readout_1}{binding}A{readout_2}{readout_3}{primer_2_rc}"
    # readout_i: rc(encoding[target]['readouts'][(probe_index + i) % 4]), i = 1, 2, 3
hybridization: {temperature_c: 37, na_mM: 390, formamide_pct: 30, dntp_mM: 0}
attributes:
  binding_tm: {target: encoding_probe.binding, type: annealing_temperature, tm_method: merfish_design,
               hybridization: {na_mM: 300, probe_nM: 5}}   # the design code's formula and conditions of 66-76 C
post_process:
  filters:
    target_gc: {condition: target_gc >= 0.43 & target_gc <= 0.63}
    binding_tm: {condition: binding_tm >= 66 & binding_tm <= 76}
    target_mapped_genes: {condition: target_mapped_genes <= 1}   # Bowtie 1, <= 2 mismatches;
                                                                 # exempt_copies: {kmer: 30, min_shared: 0.5}
  remove_overlap: {location_interval: 0}
  equal_space: {number_desired: 92}
```

Tm under the source formula: `tm_method: merfish_design` is the design code's `GetRegionTm`
(SantaLucia & Hicks 2004 table, `0.368 (N-1) ln[Na+]`, strand term ln(5 nM) where Primer3 uses
ln(5 nM / 4)), so the 66-76 C window applies as published; it reads only `na_mM` and `probe_nM`, and the
buffer's 30 % formamide does not enter it. Under Primer3 at 5 nM the same window rejected valid probes at
the low end and admitted ones up to 1.4 C too hot (17.7 % vs 9.1 % of random 30-mers outside 66-76). On
GAPDH, ACTB and EEF1A1 (hg38, Ensembl 113, no expression table) the formula gave 7, 11 and 38 probes
against 7, 11 and 33 with Primer3; 36 windows are shared, and each Primer3 set had one probe above 76 C
by the code (`references/uprobe/template-qa.md`).

What to edit: `encoding` (per gene, `readouts`: the 4 readout probe sequences of its code word, one
per on-bit, in any order), the two primer placeholders (20 N each), `equal_space` (92, 64, 48; 16 with bDNA),
`location_interval` (-20 for short genes), `gene_expression` (your tissue). `genome` needs `bowtie` in
`align_index`.

Readout rotation: `probe_index` numbers a gene's final probes 0, 1, 2, ... along the RNA, after
`post_process` has chosen them; probe k leaves out readout k % 4 and carries the other three, so the 4
combinations of 3 rotate along the RNA and each readout is on 3 of every 4 probes (69 of 92). The
final TSV has a `probe_index` column; check the rotation there (`encoding_probe.readout_1..3`).

Template departures:
- Readouts: rotated deterministically by `probe_index` (3 of 4 per probe, balanced along the RNA)
  instead of chosen at random, and always in configuration A (1 readout 5', 2 readouts 3'); Moffitt
  2016 also varied 1 or 2 readouts 5' at random.
- Specificity: Bowtie 1 (<= 2 substitutions, any gene) instead of 17-mer specificity indices; the
  15-nt rRNA/tRNA screen and probe-probe / probe-primer cross-homology are not modelled.
- Spacing: `equal_space` spreads 92 probes along the RNA; Moffitt 2016 took regions from the 5' end
  and then picked 92 at random.

## 9. Pitfalls

- Strand: `binding` = rc(RNA); readout sites = rc(readout probes). Inserting readout probe sequences
  directly makes the probe complementary to nothing.
- Readout orthogonality: all readouts of the panel must be screened against each other, the target
  regions, primers and the transcriptome (<= 11 nt contiguous homology).
- Primers: screen against all probe sequences (15-mers) and the T7 promoter; drop primers with
  3'-end homology to other primers.
- Highly expressed genes or paralog families: they crowd every image and remove probes of their
  paralogs; route to sequential rounds.
- Multi-colour rounds: Moffitt 2016 imaged two readouts per round (Cy5 and Alexa 750 probes), so the
  codebook's bit order must match the readout-to-dye assignment of the instrument.

## 10. Checklist

- [ ] 30-nt target regions, GC 43-63 % and Tm 66-76 C under the source formula (`tm_method: merfish_design`,
      300 mM, 5 nM; state them), no overlap except for transcripts too short.
- [ ] Every gene has its 4 readouts, each on ~75 % (3-readout design) or ~50 % (2-readout) of probes
      (template: count `encoding_probe.readout_1..3` per gene in the TSV).
- [ ] Probe = rc(readouts) + rc(RNA) + primers; A spacers present; length as designed (e.g. 132 nt).
- [ ] Probes per gene reached (92/64/48) or the gene moved to sequential rounds.
- [ ] Panel abundance: no gene above the per-gene cap; total within budget; blanks in the codebook.
- [ ] `target_mapped_genes` = 1 or explained (`_unexpressed`, `_copies`); no rRNA/tRNA 15-mers.
- [ ] Run last, with the published values (`binding_tm` is the design code's Tm):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --summary results/run1/summary.json \
    --transcripts <genome dir>/transcript.fa --gc 0.43 0.63 --tm 66 76 --tm-column binding_tm --min-gap 0 --min-probes 48
```

## 11. References

- Chen KH et al. 2015. Spatially resolved, highly multiplexed RNA profiling in single cells. Science 348:aaa6090. doi:10.1126/science.aaa6090
- Moffitt JR et al. 2016. High-throughput single-cell gene-expression profiling with multiplexed error-robust FISH. PNAS 113:11046-11051. doi:10.1073/pnas.1612826113
- Moffitt JR et al. 2016. High-performance multiplexed FISH in culture and tissue with matrix imprinting and clearing. PNAS 113:14456-14461. doi:10.1073/pnas.1617699113
- Moffitt JR, Zhuang X. 2016. RNA imaging with MERFISH. Methods Enzymol 572:1-49. doi:10.1016/bs.mie.2016.03.020
- Wang G, Moffitt JR, Zhuang X. 2018. Multiplexed imaging of high-density libraries of RNAs with MERFISH and expansion microscopy. Sci Rep 8:4847. doi:10.1038/s41598-018-22297-7. Author Correction (encoding probes ~50 uM, not ~300 uM): Sci Rep 8:6487. doi:10.1038/s41598-018-24844-8
- Moffitt JR et al. 2018. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 362:eaau5324. doi:10.1126/science.aau5324
- Xia C et al. 2019. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression. PNAS 116:19490-19499. doi:10.1073/pnas.1912459116
- Xia C, Babcock HP, Moffitt JR, Zhuang X. 2019. Multiplexed detection of RNA using MERFISH and branched DNA amplification. Sci Rep 9:7721. doi:10.1038/s41598-019-43943-8 (Table S1: amplifier sequences)
- Zhang M et al. 2021. Spatially resolved cell atlas of the mouse primary motor cortex by MERFISH. Nature 598:137-143. doi:10.1038/s41586-021-03705-x
- Blosser TR, Zhuang X. 2021. MERFISH probe design v1. protocols.io. doi:10.17504/protocols.io.bujfnujn
- Xu Q et al. 2009. Design of 240,000 orthogonal 25mer DNA barcode probes. PNAS 106:2289-2294. doi:10.1073/pnas.0812506106
- ZhuangLab MERFISH_analysis: example_scripts/library_design_example.m, probe_construction/TRDesigner.m (github.com/ZhuangLab/MERFISH_analysis)
- Vizgen. MERSCOPE Platform Grant Application Resource, Doc 90000111 Rev B (gene panel abundance limits).
