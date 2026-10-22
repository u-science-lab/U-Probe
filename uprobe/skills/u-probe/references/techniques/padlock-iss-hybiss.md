# Padlock probes for in situ sequencing (ISS, HybISS)

> Read this when: designing single padlock probes (PLPs) that are circularized by ligation and read out by
> rolling-circle amplification (Ke 2013 ISS, HybISS, direct-RNA padlocks), choosing cDNA vs RNA as the ligation
> template, or editing the `padlock_iss` template. For padlock + primer pairs see `starmap.md`, `mip-seq.md`, `double-hyb-rca.md`;
> general rules in `references/concepts/probe-architecture.md` (padlocks) and `references/concepts/specificity.md`.

## 1. What it is and when to choose it

A padlock is one linear oligo whose two ends (arms) hybridize side by side on the target. A ligase joins the
probe's 3'-OH to its own 5'-phosphate only when both termini are base-paired, so the probe becomes a circle
locked on the target. Phi29 polymerase copies the circle by RCA into a rolling-circle product (RCP) with
hundreds of copies of the backbone (Banér 1998), visible at 20x. One ligated probe gives one spot.

| Generation | Template | Readout of the backbone | Source |
| --- | --- | --- | --- |
| In situ padlock genotyping | cDNA (LNA RT primer), target-primed RCA | fluorescent detection oligo | Larsson 2010 |
| ISS, barcoded padlock | cDNA | 4-nt barcode by sequencing-by-ligation (SBL), 4^4 = 256 codes | Ke 2013 |
| ISS, gap-fill padlock | cDNA | target bases copied into the gap, read by SBL | Ke 2013; Mignardi 2015 (oligo gap-fill) |
| HybISS | cDNA (random decamers) | 20-nt ID read by bridge probes + 4 dye-labelled readout oligos | Gyllborg 2020 |
| Direct-RNA (dRNA) ISS | mRNA itself (PBCV-1 / SplintR ligase, or a commercial kit) | as HybISS | Schneider 2017; Lee 2022 |

Choose padlocks for single-nucleotide discrimination (SNVs, paralogs, splice junctions), large tissue areas at low
magnification, and panels of ~100-200 genes decoded combinatorially. The weakness is detection efficiency (RT,
hybridization and ligation each lose molecules; Lee 2022 gained fivefold by skipping RT); prefer smFISH/MERFISH-type probes when counting nearly every molecule matters more.

## 2. Probe architecture

HybISS PLP, 70 nt (Gyllborg 2020; all 1,195 published PLPs in Supp. Table 1 fit this layout exactly):

```text
5'-P-[arm_5: mRNA nt 16-30, 15][ID, 20][anchor, 20][arm_3: mRNA nt 1-15, 15]-OH-3'
Pvalb:  CTCGCGGCCCTCCTG  TAGACGTAGTGAGCATGACT  TCTACGAGTTTGCAGTCACG  CCCATCTCAGTCCTT
        (target window, mRNA sense: CCCATCTCAGTCCTT|CTCGCGGCCCTCCTG = nt 1-15|16-30)
```

On the cDNA template the arms carry the mRNA (sense) sequence. Ends of the bound PLP (backbone loops away):

```text
               arm_3 = mRNA 1-15  arm_5 = mRNA 16-30
PLP ends  5'-..CCCATCTCAGTCCTT|CTCGCGGCCCTCCTG..-3'   (3'-OH of arm_3 | 5'-P of arm_5)
cDNA      3'-..GGGTAGAGTCAGGAA GAGCGCCGGGAGGAC..-5'
                              ^ ligation junction = mRNA nt 15|16
```

- The 5' arm (5'-P donor) binds the upstream half of the cDNA window (complement of mRNA 16-30), the 3' arm (3'-OH
  acceptor) the downstream half (complement of mRNA 1-15). Check: arm_3 + arm_5 = the mRNA window (cDNA design).
- Direct-RNA PLP (template = mRNA): arms are reverse complements, 5'-P-[rc(mRNA 1-15)][backbone][rc(mRNA 16-30)]-OH-3';
  arm_3 + arm_5 = rc(window); junction still at mRNA nt 15|16.

```text
               arm_3 = rc(16-30)  arm_5 = rc(1-15)
PLP ends  5'-..CAGGAGGGCCGCGAG|AAGGACTGAGATGGG..-3'   (3'-OH of arm_3 | 5'-P of arm_5)
mRNA      3'-..GTCCTCCCGGCGCTC TTCCTGACTCTACCC..-5'   (nt 30..16 | 15..1)
                              ^ ligation junction = mRNA nt 15|16
```

Ke 2013 barcoded PLP (~89 nt; arm lengths obtained by mapping Supp. Table 5 probes onto RefSeq mRNAs):

```text
5'-P-[arm_5 ~18][TCCTCTATGATTACTGAC 18][TGCGTCTATTTAGTGGAGCC 20 = anchor/SBL primer site][NNNN 4 = barcode][CTATCTTCTTT 11][arm_3 ~18]-3'
```

Ke 2013 gap-fill PLP (ACTB, NM_001101.3): 3' arm on mRNA nt 117-136 (20), 4-nt gap at nt 137-140, 5' arm from nt 141
(~20). A polymerase extends the 3' arm across the gap and the ligase closes the circle; the anchor primer equals
the 3'-arm end (the beta-actin anchor primer of Supp. Table 4), so SBL reads into the copied gap bases. Oligo gap-fill (Mignardi 2015) instead
ligates a 5'-phosphorylated hexamer into a 6-nt gap (two ligations).

Readout oligos. Bridge probes, anchor oligos and SBL primers have the PLP's sense because they pair with the
RCP = rc(PLP); dye-labelled readout oligos pair with the bridge probe.

| Oligo | Layout | Source |
| --- | --- | --- |
| HybISS bridge probe (39 nt) | [ID nt 2-18, 17][2-nt linker][rc(readout oligo), 20] | Gyllborg 2020 + Supp. Table 1 |
| HybISS readout (detection) oligo | 19-20 nt, 5' dye; 4 sequences, one per colour | Supp. Table 1 |
| HybISS anchor oligo | 17-nt 3' end of the anchor, 3'-AF750, stains all RCPs | Supp. Table 1 |
| Ke 2013 anchor primer + interrogation nonamers | anchor primer's 3' end abuts the barcode; T4 ligase joins p-N(9)-dye nonamers, one queried position per cycle | Ke 2013 Supp. Table 4; protocols.io s8vehw6 |

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Arm length, HybISS | 15 + 15 (30-nt site) | nt | Gyllborg 2020, doi:10.1093/nar/gkaa792 |
| Arm length, Ke 2013 barcoded | 15-20 each, mostly 18 + 18 | nt | Ke 2013 Supp. Table 5 (mapped here), doi:10.1038/nmeth.2563 |
| Arm length, pciSeq | two arms totalling 40 | nt | Qian 2020, doi:10.1038/s41592-019-0631-4 |
| Gap, gap-fill ISS / oligo gap-fill | 4 / 6 | nt | Ke 2013 Supp. Table 4 (mapped here); Mignardi 2015, doi:10.1093/nar/gkv772 |
| Target Tm window (30-nt arm pair) | 65-75 | °C | Gyllborg 2020 |
| Tm conditions | 0.1 µM probe, 75 mM monovalent, 10 mM Mg2+, 20% formamide | - | multi_padlock_design readme.txt |
| Tm formula (repo) | Sugimoto 1996 NN; +16.6 log10([Na+] + 3.795 sqrt[Mg2+]); -0.72 °C per % formamide; in U-Probe `tm_method: multi_padlock_design` (same values to 1e-12 °C) | - | multi_padlock_design lib/screenseq.py |
| Space between targets | 15 | nt | Gyllborg 2020 |
| PLPs per gene | 5, chosen at random among passing targets | - | Gyllborg 2020 |
| ID / anchor | 20 (Affymetrix GenFlex Tag 16K) / 20 (6 shared variants, linker only) | nt | Gyllborg 2020 |
| Ke 2013 barcode / anchor site | 4 / 20 | nt | Ke 2013; Ke 2016 review, doi:10.1002/humu.23051 |
| HybISS decoding | 5 cycles, 4 colours, min. Hamming distance 2 | - | Gyllborg 2020 |
| GC | none stated; published HybISS targets span 43-83% (computed here) | % | Gyllborg 2020 Supp. Table 1 |
| Homopolymers | repo prefers PLPs without AAAA/CCCC/GGGG/TTTT, drops 6-nt runs and low-complexity targets | - | multi_padlock_design lib/finalizeprobes.py |
| PLP concentration | 10 each | nM | Gyllborg 2020; protocols.io xy4fpyw |
| cDNA ligation | Tth ligase 0.5 U/µl + RNase H 0.4 U/µl, 1x Tth buffer + 50 mM KCl (75 mM K+, 10 mM Mg2+), 20% formamide; 37 °C 30 min, then 45 °C 1.5 h | - | protocols.io xy4fpyw |
| Ke-era ligation | Ampligase 0.5 U/µl, 45 °C | - | Mignardi 2015 |
| RT primers | LNA-modified gene-specific (Ke 2013) / random decamers 5 µM (HybISS) | - | Ke 2013 Supp. Table 4; protocols.io xy4fpyw |
| RCA | phi29 1 U/µl + Exonuclease I 0.4 U/µl; 37 °C up to 4 h, then 30 °C overnight | - | protocols.io xy4fpyw |
| Direct-RNA ligation | PBCV-1 ligase 200 nM, PLP 1 µM + RCA primer 1 µM, 37 °C 30 min | - | Schneider 2017, doi:10.1261/rna.057836.116 |
| Bridge / readout hybridization | 10 nM bridge, 1 h RT; 100 nM readout, 2 h RT; 2x SSC, 20% formamide | - | Gyllborg 2020 |

**Protocols disagree:** arm length (15 + 15 in HybISS vs ~18 + 18 in Ke 2013 vs 40 nt total in pciSeq); RT priming
(gene-specific LNA primers vs random decamers); ligase (Ampligase vs Tth on cDNA; PBCV-1 on RNA).

## 4. Target region and template

| | cDNA template (Ke 2013, HybISS) | Direct RNA (SplintR, dRNA kits) |
| --- | --- | --- |
| Arms | mRNA sense | reverse complement of mRNA |
| Ligase, fidelity | thermostable NAD ligase (Tth, Ampligase) at 45 °C on DNA: high | PBCV-1 (SplintR): efficient on RNA (Lohman 2014), limited fidelity at both termini (Krzywkowski 2017) |
| RCA priming | target-primed: phi29's 3'->5' exonuclease trims the cDNA until it primes on the circle (Larsson 2010; Schneider 2017) | needs an RCA primer: Schneider 2017 used 44 nt = 20 nt on the RNA next to the PLP + 4-nt linker + 20 nt on the PLP backbone at its 3' end |
| Efficiency | limited by in situ RT | about 5x more RCPs than cDNA-HybISS (Lee 2022, commercial kit) |

- Region: HybISS took sites anywhere on the mRNA; multi_padlock_design designs on RefSeq NM/NR mRNAs, in the
  sequence common to all transcript variants (ClustalW alignment, lib/parmsa.py). Published sites include 3'-UTR
  sites, e.g. the two Pvalb sites shown here lie at nt 449 and 522 of NM_013645.4 (CDS 72-404).
- With gene-specific RT primers the PLP site must lie 5' of the primer on the mRNA, inside the cDNA it makes
  (Ke 2013 ACTB: primer on nt 207-231, gap-fill site on nt 117-161). Random priming removes this constraint.
- Isoforms: for gene-level counts put sites in exons shared by all variants; to tell isoforms apart, place
  the site, ideally the junction, on an isoform-specific exon or exon-exon junction.
- Direct RNA: site accessibility matters. Moving a PBCV-1 PLP 500 nt changed RCP counts threefold, and Lohman 2014
  reports faster RNA-templated ligation with A/T than with G/C at the donor (Schneider 2017). Test several sites.

## 5. Specificity

- One ligated PLP gives one RCP and one decoded spot. Nothing requires several probes to coincide, so every
  off-target ligation is a false spot carrying the gene's code. Filter per probe: `target_mapped_genes <= 1`.
  Do not use the smFISH-style `off_target_pool`.
- An off-target site needs (1) both arms to hybridize at the ligation temperature (a 15-nt arm alone is weak:
  the Pvalb arms compute to 35.5 and 48.1 °C under the HybISS ligation buffer with U-Probe's primer3), and
  (2) base-paired termini at the junction. Ligases discriminate most at the 3'-OH acceptor base
  (Schneider 2017: a 3'-terminal mismatch cut RCPs from 94 to 3 per cell). RNA-templated PBCV-1 ligation tolerates
  more mismatches (Krzywkowski 2017), so direct-RNA designs depend more on hybridization specificity (iLock improves it).
- multi_padlock_design rule: BLAST (word size 7, plus strand) of the 30-mer against RefSeq mRNAs, predicted
  XM/XR ignored. A hit in another gene rejects the target when it is longer than one arm (>50% coverage),
  more than 80% identical, and covers the ligation site +-5 nt. The target must also align perfectly to its own gene.
  In BLAST's output columns (lib/readblast.py): length > 15, pident > 80, qstart <= 10 and qend >= 21, so the hit
  spans nt 10-21 around the junction at nt 15|16. BLAST reports only hits below its default e-value of 10.
- The template approximates the rule with two Bowtie 1 searches, both filtered `<= 1`: `target_mapped_genes`
  (the whole 30-mer, within 2 substitutions) and `junction_mapped_genes` (nt 6-25, the 20 nt around the junction,
  a `layout` part, so no helper oligo; within 1 substitution). Both kinds of hit meet the repo rule, and by an
  estimate with blastn's default 1/-2 scoring on a RefSeq-sized database BLAST reports them (E < 10). Not
  searched: hits of 16-19 nt, gapped hits, and more divergent hits that BLAST still reports. With 2 substitutions
  the 20-nt search would also flag hits BLAST does not report (it keeps 863 of the 1,195 published targets, 1,143
  with 1). On the Ensembl 113 transcriptomes 1,155 published targets pass the 30-mer search and 1,138 pass both;
  the 17 more removed hit genes such as Dtnbp1, Lrrc58, NFYC, WDR82, CNTNAP3 and DUSP14, or lncRNAs, retained
  introns and processed pseudogenes, many of which RefSeq (the repo's database) lacks.
- Expression: an off-target RNA that is absent cannot ligate. Add `gene_expression` (`min_value`, e.g. TPM >= 1).
  A `min_fraction_of_target` threshold accepts false spots at roughly that fraction. `exempt_copies` means that
  pseudogene or retrocopy signal counts as the target's.
- Not seen by `mapped_genes`: intronic or unannotated RNA (random-primed RT copies them too). Genomic DNA is a
  minor concern because it is double-stranded; in situ DNA detection needs restriction digestion and exonucleolysis (Mignardi 2015).

## 6. U-Probe template

`uprobe templates show padlock_iss -o protocol.yaml` (HybISS design, cDNA template). Key lines:

```yaml
extracts:
  target_region:
    source: mRNA          # whole mRNA, UTRs included
    transcripts: mane     # MANE Select, else Ensembl canonical
    length: 30            # 2 x 15-nt arms
    layout:               # only for the junction search: nt 6-25 of the window
      template: "{arm_3_end}{junction}{arm_5_end}"
      parts: {arm_3_end: {length: 5}, junction: {length: 20}, arm_5_end: {length: 5}}
probes:
  padlock:
    template: "{arm_5}{id}TCTACGAGTTTGCAGTCACG{arm_3}"   # ID + anchor
    parts:
      arm_5: {expr: "target_region[15:30]"}   # mRNA sense (binds cDNA)
      id:    {expr: "encoding[target]['ID']"}
      arm_3: {expr: "target_region[0:15]"}
attributes:
  target_tm: {target: target_region, type: annealing_temperature,     # 30-mer, repo formula and conditions
              tm_method: multi_padlock_design,
              hybridization: {na_mM: 75, mg_mM: 10, probe_nM: 100, formamide_pct: 20}}
  target_mapped_genes: {target: target_region, type: mapped_genes,
                        aligner: bowtie, max_mismatches: 2, strand: sense}
  junction_mapped_genes: {target: target_parts.junction, type: mapped_genes,
                          aligner: bowtie, max_mismatches: 1, strand: sense}
post_process:
  filters:
    target_tm: {condition: target_tm >= 65 & target_tm <= 75}
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    junction_mapped_genes: {condition: junction_mapped_genes <= 1}
  remove_overlap: {location_interval: 15}
  equal_space: {number_desired: 5}
```

1. Replace `meta.replace` items: `genome` (its `align_index` must include `bowtie`, i.e. Bowtie 1), `targets`,
   and one unique 20-nt `encoding[gene].ID` per gene, taken from an orthogonal tag set. The anchor is a literal
   in the template; choose one of the six published anchors.
2. Optionally add a `gene_expression` table (commented block in the template).
3. `uprobe validate -p protocol.yaml -g genomes.yaml --check-targets`, then `uprobe run -p protocol.yaml -g genomes.yaml -o results`.
4. Order the PLPs 5'-phosphorylated, or treat them with T4 PNK. U-Probe does not add modifications.

Variants (each validated with `uprobe validate`):

| Variant | Edit |
| --- | --- |
| Direct RNA (SplintR) | `arm_5: rc(target_region[0:15])`, `arm_3: rc(target_region[15:30])`; set `hybridization` to the 37 °C ligation; an RCA primer next to the PLP is needed (see `mip-seq.md` for a padlock + primer layout) |
| Ke 2013 barcoded | `length: 36`; template `"{arm_5}TCCTCTATGATTACTGACTGCGTCTATTTAGTGGAGCC{bc}CTATCTTCTTT{arm_3}"`, `arm_5: target_region[18:36]`, `arm_3: target_region[0:18]`, 4-nt `bc` in `encoding`; layout parts `arm_3_end` and `arm_5_end` 8 nt each, so the junction part stays centred (18\|19) |
| Gap-fill | `length: 44` + `layout: {template: "{arm3_site}{gap}{arm5_site}", parts: {arm3_site: {length: 20}, gap: {length: 4}, arm5_site: {length: 20}}}`; `arm_5: target_parts['arm5_site']`, `arm_3: target_parts['arm3_site']`; this layout replaces the junction one, so drop `junction_mapped_genes` and its filter |

**U-Probe departs:**
- Tm model: none left. `target_tm` uses the repo formula (`tm_method: multi_padlock_design`), so all 1,195
  published HybISS targets pass 65-75 °C; with primer3, the default, 327 (27%) did. Primer3 under the same
  conditions is 3.4 +- 1.5 °C lower (-3.2 to +7.5 °C, larger for AT-rich targets; correlation with GC -0.71),
  so no primer3 window is equivalent: the best 10-°C one (59.1-69.1 °C) keeps 1,110, and 60-73 °C keeps 1,121.
  Mouse Pvalb: 70.5 °C (primer3 67.2). All filters together (Tm, both specificity searches on Ensembl 113) keep
  1,138 of the 1,195 published targets; the previous template kept 314.
- Selection. HybISS picked 5 sites at random; the template spreads 5 sites evenly along the mRNA with `equal_space`.
- Region. The repo uses the region common to all RefSeq variants; the template uses one reference transcript.
- Specificity. Two Bowtie 1 searches (30-mer, <= 2 substitutions; junction 20-mer, <= 1) approximate the BLAST
  junction rule (section 5), against the transcriptome of your genome instead of RefSeq mRNAs.
- Sequence filters. No homopolymer or low-complexity preference; add a `sequence_pattern` filter excluding
  AAAA/CCCC/GGGG/TTTT if wanted.
- Template. The arms suit a cDNA template only. `arm_5_tm`/`arm_3_tm` are reported, not filtered.

## 7. Pitfalls

- Arm orientation: sense arms for a cDNA template, rc arms for an RNA template; a swapped design binds nothing.
  In cDNA designs the 5' arm is the downstream half of the mRNA window. Without a 5' phosphate nothing ligates.
- Tm windows: apply 65-75 °C to the 30-nt arm pair under 75 mM K+ / 10 mM Mg2+ / 20% formamide, never to one arm
  or under 2x SSC, and with `tm_method: multi_padlock_design`: under primer3 the window rejects most published
  targets (section 6).
- Gene-specific RT: a PLP site 3' of the RT primer on the mRNA lies outside the cDNA and is never detected.
- Direct RNA without an RCA primer: no target-primed RCA, no signal (Schneider 2017).
- Reusing IDs or anchors that pair with transcripts or with each other gives background RCPs and false codes.
  Check the backbone with `mapped_genes`/BLAST and the bridge set for cross-hybridization.
- `source: CDS` misses the 3'-UTR sites HybISS used. Keep `mRNA` unless the CDS is required.
- Fewer than 5 passing sites gives a `partial` run (`summary.json`); acceptable, HybISS also kept fewer.

## 8. Verification checklist

- [ ] Each PLP = arm_5 (15) + ID (20) + anchor (20) + arm_3 (15) = 70 nt (HybISS), with a 5' phosphate.
- [ ] arm_3 + arm_5 equals the target window in the reference mRNA (cDNA design), or its reverse complement (RNA design).
- [ ] The window lies on the intended transcript, in exons shared by the isoforms you want to count.
- [ ] `target_tm` in the chosen window under the stated conditions; the model (primer3 vs repo formula) is recorded.
- [ ] `target_mapped_genes == 1` with Bowtie 1 (<= 2 substitutions) and `junction_mapped_genes == 1` (<= 1 in the
      20 nt around the junction), with empty `_off_targets` columns; the off-targets of any exempted copies or
      unexpressed genes were reviewed.
- [ ] Sites at least 15 nt apart, up to 5 per gene; with gene-specific RT, every site within its primer's cDNA.
- [ ] IDs unique per gene; bridge probe = ID nt 2-18 + 2-nt linker + rc(readout); codes at Hamming distance >= 2.
- [ ] Backbone (ID + anchor) has no transcriptome hit and no strong hairpin with the arms.
- [ ] Arms are 15 + 15 as published; no Tm is derived per arm (`arm_5_tm`/`arm_3_tm` are reported only) and no extra Tm band is added: the window is the 30-nt `target_tm` under the repo's formula and conditions (75 mM monovalent, 10 mM Mg2+, 20 % formamide, 0.1 uM).
- [ ] Run last, with the file's own values (Gyllborg 2020: Tm 65-75 C, space between targets 15 nt, 5 per gene; no GC is published):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --summary results/run1/summary.json \
    --transcripts <genome dir>/transcript.fa \
    --sense padlock.arm_5 --sense padlock.arm_3 --tm 65 75 --tm-column target_tm --min-gap 15 --min-probes 5
```

## 9. References

- Nilsson 1994, Science. Padlock probes: circularizing oligonucleotides for localized DNA detection. doi:10.1126/science.7522346
- Banér 1998, Nucleic Acids Res. Signal amplification of padlock probes by rolling circle replication. doi:10.1093/nar/26.22.5073
- Larsson 2004, Nat Methods. In situ genotyping individual DNA molecules by target-primed rolling-circle amplification of padlock probes. doi:10.1038/nmeth723
- Larsson 2010, Nat Methods. In situ detection and genotyping of individual mRNA molecules. doi:10.1038/nmeth.1448
- Ke 2013, Nat Methods. In situ sequencing for RNA analysis in preserved tissue and cells. doi:10.1038/nmeth.2563 (Supp. Tables 4-6 used; Online Methods paywalled)
- Lohman 2014, Nucleic Acids Res. Efficient DNA ligation in DNA-RNA hybrid helices by Chlorella virus DNA ligase. doi:10.1093/nar/gkt1032
- Mignardi 2015, Nucleic Acids Res. Oligonucleotide gap-fill ligation for mutation detection and sequencing in situ. doi:10.1093/nar/gkv772
- Ke 2016, Hum Mutat. Fourth generation of next-generation sequencing technologies: promise and consequences. doi:10.1002/humu.23051
- Schneider 2017, RNA. Efficient in situ detection of mRNAs using the Chlorella virus DNA ligase for padlock probe ligation. doi:10.1261/rna.057836.116
- Krzywkowski 2017, Nucleic Acids Res. Fidelity of RNA templated end-joining by chlorella virus DNA ligase and a novel iLock assay with improved direct RNA detection accuracy. doi:10.1093/nar/gkx708
- Qian 2020, Nat Methods. Probabilistic cell typing enables fine mapping of closely related cell types in situ. doi:10.1038/s41592-019-0631-4
- Gyllborg 2020, Nucleic Acids Res. Hybridization-based in situ sequencing (HybISS) for spatially resolved transcriptomics in human and mouse brain tissue. doi:10.1093/nar/gkaa792
- Gyllborg & Nilsson, protocols.io. HybISS: Hybridization-based In Situ Sequencing. doi:10.17504/protocols.io.xy4fpyw
- Gyllborg, Yokota & Nilsson, protocols.io. In situ sequencing for RNA analysis in tissue sections. doi:10.17504/protocols.io.s8vehw6
- Lee 2022, Sci Rep. Direct RNA targeted in situ sequencing for transcriptomic profiling in tissue. doi:10.1038/s41598-022-11534-9
- Moldia/multi_padlock_design (readme.txt; lib/screenseq.py, readblast.py, distributeprobes.py, finalizeprobes.py, parmsa.py). https://github.com/Moldia/multi_padlock_design
