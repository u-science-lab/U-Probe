# TDDN-FISH primary probes (Wang 2025)

> Read this when: a user asks for TDDN-FISH (tetrahedral DNA dendritic nanostructure-enhanced FISH)
> probes, or for a few 35-nt primary probes that each carry one 30-nt readout for an enzyme-free DNA
> nanostructure amplifier. U-Probe designs the primary probes only. Related: `smifish.md` (one
> readout tail), `merfish.md` (source of the readouts), `hcr-v3.md`, `smfish.md` (template).

SD1 = Supplementary Data 1 of Wang 2025 (all probe and TDDN sequences). "Computed here" = parsed
from SD1 or RefSeq for this card.

## 1. What it is and when to choose it

- Unlabelled primary probes are hybridized once (8-12 h). Each round, a pre-assembled TDDN binds the
  readouts of that round (45 min-1 h, RT), is imaged, then cut off with TCEP (15 min). A TDDN is
  three layers of DNA tetrahedra (T0, T1, T2; 17-bp edges): strand T0a's 30-nt sticky end binds the
  readout, T2 binds dye oligos, 27 dyes per TDDN (authors, peer review file).
- One probe + TDDN is a spot: 3 probes per mRNA in cells (smFISH in the same paper: 48), 1 for
  miR-21. ~1 h per round after primary hybridization vs >= 8 h for HCR (Results).
- Codes: F colours x N rounds = F^N (3 colours, 2 rounds: 9 human mRNAs, 8 influenza vRNAs). Tissue
  cell typing gave the 10 marker genes of a cell type one code (60 and 90 genes).
- Choose it for short RNAs, RNAs with few accessible sites, fast low-plex rounds. Not for dense
  high-plex panels: no error-correcting code; ACTB-level RNAs overlap (Discussion).

## 2. Probe architecture (5' -> 3')

```
primary (65 nt): 5'-[ readout 30 = rc(T0a sticky end) ][ rc(RNA) 35 ]-3'      unmodified oligo
T0a:             5'-[ sticky end 30 ]-SS-[ tetrahedron strand 55 ]-3'          TCEP cuts at -SS-
dye oligo:       5'-Cy5/Cy3/FAM-CTCATCCATCGGTC-3'   one sequence for all colours; binds T2b-d
```

- Readout at the 5' end in all 185 SD1 primary probes, never 3' (computed here). The binding part
  is rc(RefSeq mRNA): ACTB 3/3, TGFB1, TNF and CCL2 6/6 each, PPIA 41/44, all 1,439 tissue windows.
- The 16 T0a sticky ends (T0a-hyb1-Cy5 ... T0a-hyb13) are exactly the 16 MERFISH readout probes of
  Moffitt & Zhuang 2016 Table 1 (computed here). Colour = the dye pre-loaded on that TDDN's T2.
- 2 SD1 probes are 70 nt: R1-FAM readout + TTGCG (the R1-Cy5 readout's 3' end) + a 35-nt Satb2
  site used with R1-Cy5 elsewhere; do not copy them as a 40-nt design.

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Binding length | 35 nt (183/185 SD1 probes) | | Methods "Primer probe design"; SD1 |
| Readout | 30 nt, 5' end, one per probe | | Methods; SD1 |
| Design tool | OligoArray 2.0 | | Methods (Rouillard 2003) |
| Tm | > 70 C | OligoArray conditions, not stated | Methods |
| Cross-hyb. | "minimum Tm of 72 C for hybridization to potential off-target sequences" | | Methods |
| Structure | "no internal secondary structures with a Tm lower than 76 C" (as printed; same thresholds as Chen 2015, `merfish.md`) | | Methods |
| Homopolymer | no run > 6 identical nt; 2 of 182 SD1 cell targets break it (A8, A15) | | Methods; computed here |
| GC | not stated; published 0.23-0.69, median 0.46 | 1,574 distinct targets | computed here |
| Tm, buffer A | published 42.9-65.7 C, median 52.3 | 2x SSC (390 mM Na+), 30 % FA, 16.7 nM, primer3 | computed here |
| Probes per RNA | 3 per round (cells); 1 (miR-21); "3-5 primary probes per target ... each targeting distinct transcript regions" | 2-round code: 6 per gene | Results; peer review file; SD1 |
| Spacing | not stated; ACTB sites 300-460 nt apart; PPIA 44 sites, 20-nt gaps | | computed here |
| Tissue panels | 10 windows per gene; 75 % overlap their neighbour (steps 10 or 20 nt); at the 3' end | | SD1; computed here |
| Primary hyb. | buffer A: 10 % dextran sulfate, 2x SSC, 30 % formamide, 20 mM RVC, 1 mg/mL yeast tRNA, 0.1 % BSA, probes "50 nM for each RNA" | 37 or 42 C, 8-12 h | Methods |
| Optimization | signal 37 > 42 > 47 C; 30 > 20 > 10 % formamide (ACTB) | | Suppl. Fig. 1d-g |
| Primary wash | 2x SSC, 30 % formamide, 20 mM RVC | 37 C, 3 x 10 min | Methods |
| TDDN hyb. | buffer B = A with 5 % formamide; wash 2x SSC, 5 % formamide | RT | Methods |
| TCEP | 0.25 uM (as printed) in 2x SSC, pH 7.0 | 15 min | Methods |

## 4. Target region

- Not stated; SD1 sites lie in CDS and UTRs; tissue windows end a median 42 nt from the transcript
  end (computed here). U-Probe: `source: mRNA`, `transcripts: mane`.
- miRNA: "miR-21 (72 nucleotides)" matches the MIR21 stem-loop (Ensembl ENST00000362134). Mature
  miR-21-5p is 22 nt (RNAcentral URS000039ED8D), so a 35-nt site covers precursor sequence
  (inference; the miR-21 probe is not in SD1); windows starting at stem-loop positions 1-8 contain
  all of miR-21-5p (computed here). The stem-loop lies in the 3' UTR of VMP1's MANE transcript
  (ENST00000262291, Ensembl 116): the probe also binds VMP1 mRNA.

## 5. Specificity

- Per-probe rule: one bound probe carries a 27-dye TDDN, enough for a spot (authors: > 6 dyes
  suffice); Reviewer 2 raised single-probe false positives. Keep every probe specific, as for HCR
  pairs and padlocks; no `off_target_pool`.
- U-Probe has no OligoArray cross-hybridization Tm. Proxy (U-Probe departs): Bowtie 1 within 3
  substitutions (exhaustive) plus Bowtie2 local (partial matches, indels), both `<= 1`, with the
  sample's `gene_expression` (`references/concepts/specificity.md`).
- Use the published readouts (MERFISH set "validated for human samples"); a new readout needs a
  transcriptome and cross-readout screen.

## 6. U-Probe: adapt `smfish_stellaris`

`uprobe templates show smfish_stellaris -o protocol.yaml`, keep `genome`, then replace these blocks
(two genes, 2-round codes as SD1's TNF "12" and IL6 "23"; passes `uprobe validate`):

```yaml
name: tddn_fish
targets: [TNF, IL6]
encoding:   # readout PROBE = T0a sticky end (MERFISH readout); one per round, from your codebook
  TNF: {readouts: [CGCAACGCTTGGGACGGTTCCAATCGGATC, GCGGGAAGCACGTGGATTAGGGCATCGACC]}   # hyb1-Cy5, hyb2-Cy3
  IL6: {readouts: [CGAATGCTCTGGCCTCGAACGAACGATAGC, AAGTCGTACGCCGATGCGCAGCAATTCACT]}   # hyb1-Cy3, hyb2-FAM
extracts:
  target_region: {source: mRNA, transcripts: mane, length: 35, step: 1}
probes:
  primary_probe:
    template: "{readout}{binding}"   # readout at the 5' end, as in all 185 published probes
    parts:
      readout: {expr: "rc(encoding[target]['readouts'][probe_index % 2])"}   # % number of rounds
      binding: {expr: "rc(target_region)"}
hybridization: {temperature_c: 37, na_mM: 390, formamide_pct: 30, probe_nM: 16.7, dntp_mM: 0}
attributes:
  target_gc: {target: target_region, type: gc_content}
  binding_tm: {target: primary_probe.binding, type: annealing_temperature}
  target_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie, max_mismatches: 3,
                        strand: sense, exempt_copies: {kmer: 30, min_shared: 0.5}}
  target_mapped_genes_bt2: {target: target_region, type: mapped_genes, aligner: bowtie2, strand: sense,
                            exempt_copies: {kmer: 30, min_shared: 0.5}}
post_process:
  filters:
    binding_tm: {condition: binding_tm >= 47}
    homopolymer: {type: sequence_pattern, target: target_region,
                  exclude_patterns: [AAAAAAA, CCCCCCC, GGGGGGG, TTTTTTT]}
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    target_mapped_genes_bt2: {condition: target_mapped_genes_bt2 <= 1}
  remove_overlap: {location_interval: 2}
  equal_space: {number_desired: 6}   # 3 per round x 2 rounds
  sorts: {is_ascending: [target_mapped_genes]}
summary: {report_name: rna_report, attributes: [target_gc, binding_tm, target_mapped_genes]}
```

- `probe_index` numbers the final probes along the RNA, so rounds alternate along it. R rounds:
  R readouts per gene, `% R`, `number_desired: 3R`; one round: a one-item list, `% 1`, 3.
- Genome: `align_index: [bowtie, bowtie2]`. Add `gene_expression` for the sample.
- miRNA (validated): `targets: [MIR21]`, one readout, `% 1`, `number_desired: 1`. `mane` falls back
  to the Ensembl canonical stem-loop (38 windows); host-gene transcripts that contain it (VMP1)
  appear in `target_mapped_genes_copies`. Unannotated RNA: a literal target
  `- ENSG...: <sequence>`, named by its gene id (`references/uprobe/yaml-reference.md`).
- Test run (Ensembl 113, Bowtie 1 attribute only): TNF 6/6 and IL6 6/6 with the readouts
  alternating as coded; MIR21 1/1 (window 1, VMP1 in `_copies`).

Departures (not published rules): Tm floor 47 C under buffer A (rejects 29 of 1,574 published
targets, 1.8 %); `probe_nM` 16.7 reads "50 nM for each RNA" as the 3-probe set; Bowtie proxies for
OligoArray's Tm tests; no structure filter (`fold_score` can rank); non-overlapping, evenly spaced
probes, not the overlapping tissue windows.

U-Probe cannot do: TDDN monomer, sticky-end or dye-oligo sequences, the -SS- modification of T0a,
assembly; the colour/round codebook or cell-type codes; OligoArray's cross-hybridization and
structure Tm tests; readout design and orthogonality; mature vs precursor miRNA discrimination.

## 7. Pitfalls

- Readout end: 5'. A 3' readout is not the published design.
- Strand: `encoding` holds the T0a sticky end; the probe gets its rc. SD1's "readout sequence" list
  (R1-Cy5 GATCCGATTGG...) is already the probe part: rc it before putting it in `encoding`.
- Never reuse a readout after TCEP: the cut leaves the 30-nt sticky end on the readout (duplex Tm
  71-75 C in 2x SSC, 5 % formamide, 50 nM; computed here, inference). Every published round has its
  own readout (16 T0a variants; PPIA: 15 readouts on 15 sets of 3 sites).
- A readout belongs to a (round, colour) and is shared by all genes with that colour in that round
  (R1-Cy5 on TGFB1, TNF, CCL2). Genes with the same full code cannot be told apart.
- Sites per gene grow with rounds (3 x R); a short RNA runs out of sites.

## 8. Checklist

- [ ] Each probe = 5'-rc(T0a sticky end) + 35-nt rc(RNA)-3', 65 nt.
- [ ] 3 probes per gene per round (1 for a miRNA), non-overlapping; each round's readout as coded.
- [ ] Every probe `target_mapped_genes` and `_bt2` = 1 (copies listed and explained).
- [ ] No run of > 6 identical nt; Tm stated with conditions (2x SSC, 30 % formamide, 37 C).
- [ ] No readout used in two rounds; T0a variant, dye and round recorded per readout.

## 9. References

- Wang YF, Chen HJ, He ZD, Wang ZG, Pang DW, Liu SL. 2025. Tetrahedral DNA dendritic nanostructure-enhanced FISH for high-speed, sensitive spatial transcriptomics. Nat Commun 16:9251. doi:10.1038/s41467-025-64294-1 (open: full text PMC12537822, Supplementary Information, SD1, peer review file)
- Wang YF, Liu SL, Chen HJ. 2025. Image-processing code (MATLAB; no probe design). Zenodo. doi:10.5281/zenodo.16993469
- Moffitt JR, Zhuang X. 2016. RNA imaging with multiplexed error-robust fluorescence in situ hybridization (MERFISH). Methods Enzymol 572:1-49. doi:10.1016/bs.mie.2016.03.020
- Chen KH et al. 2015. Spatially resolved, highly multiplexed RNA profiling in single cells. Science 348:aaa6090. doi:10.1126/science.aaa6090
- Rouillard JM, Zuker M, Gulari E. 2003. OligoArray 2.0: design of oligonucleotide probes for DNA microarrays using a thermodynamic approach. Nucleic Acids Res 31:3057-3062. doi:10.1093/nar/gkg426
