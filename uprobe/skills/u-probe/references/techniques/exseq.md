# Targeted ExSeq: barcoded padlocks ligated on RNA in expanded tissue

> Read this when: designing or checking targeted expansion sequencing (ExSeq) padlocks, 2 x 16-nt
> arms joined by SplintR on the RNA, gene barcodes read by in situ sequencing in an expanded gel.
> Related: `padlock-iss-hybiss.md` (direct-RNA padlocks, section 6), `mip-seq.md` (SplintR padlock
> plus an RNA-bound primer); rules in `references/concepts/probe-architecture.md`, `specificity.md`.

## 1. What it is and when to choose it

Alon 2021: tissue RNA is anchored to a swellable gel with LabelX, the sample is expanded,
re-embedded in a non-expanding gel (3.3x linear, Tissue Preparation Step 17) and its carboxylates
passivated so enzymes work in the gel. Targeted ExSeq hybridizes barcoded padlocks to the
anchored RNA, ligates them on the RNA with SplintR (PBCV-1 DNA ligase; no RT), primes Phi29 RCA
with a universal primer on the padlock backbone, and reads the gene barcode by in situ sequencing
(Illumina SBS or SOLiD). Yield: ~62% of HCRv3.0-ExFISH for the same gene in the same cells
(slope 0.62, r = 0.991, 60 HeLa cells, 8 padlocks/gene; main text, Fig. 4B, SM).

Untargeted ExSeq (random-primed RT, cDNA circularization, reads matched to ex situ sequencing)
designs no probes; U-Probe has nothing to do there.

Choose it for subcellular maps of tens to ~300 genes (panels of 42, 35 and 297 genes) in intact
tissue. Without expansion, see `padlock_iss` (HybISS) or `mipseq`.

## 2. Probe architecture

```text
RNA window 5'->3':  [nt 1-16][nt 17-32]           ligation junction between nt 16 and 17
Padlock, 74 nt, 5' phosphate (/5Phos/):
5'P-[arm_5 16][AAA][SOLiD BC 7][TCTCGGGAACGCTGAAGACGGC 22][Illumina BC 7][AAA][arm_3 16]-OH 3'
    rc(nt 1-16)                 backbone: RCA primer site                       rc(nt 17-32)
RCA primer, 18 nt:   TCTTCAGCGTTCCCGA*G*A   = rc(backbone nt 1-18), * phosphorothioate
SBS primer, 22 nt:   TCTCGGGAACGCTGAAGACGGC = backbone; binds the amplicon, reads the Illumina BC
QC detection oligo:  /5Cy3/TCTCGGGAACGCTGAAGA
AKT1, NM_001014432 nt 264-295 (repo worked example, 5' UTR site):
5'P-tcggccttccctaagc AAA AGAACCA TCTCGGGAACGCTGAAGACGGC TGCGACA AAA tttgacccaggctggc-3'
```

- Arms are rc(RNA): rc(arm_5) + rc(arm_3) is the window, contiguous; the 5'-P pairs nt 16, the
  3'-OH nt 17. Checked here: all 5,259 panel padlocks of Table S10 map so to their RefSeq mRNA.
- Arms are target-dependent, barcodes gene-dependent; spacers and backbone are constant. The RCA
  primer binds the padlock only (unlike MiP-seq's RNA-bound primer). The ligated circle is wound
  around the gel-anchored RNA (32 bp; Library Prep Step 9), so no post-fixation (SM).

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Homology region | 32 = 2 x 16, no gap | nt | SM "Probe Design - Homology"; Fig. S12 |
| Probe length | 74 = 2n + 10 (n = 32); AAA spacers fill it | nt | Probe Generation Step 13.2 |
| Arm GC | 40-65 each | % | SM; Step 9.2 `GC_RANGE` |
| Arm Tm | both arms above a per-gene floor: cortex 55, lowered to >= 52 for genes with < 6 probes; hippocampus 50; repo example 55/53/52 | C | SM; Step 7; example_genelist.xlsx |
| Tm model | MATLAB oligoprop, mean of Tm methods 2-6, DNA/DNA, 0.3 M salt, 100 nM | - | screen_homology_physical.m l.34-35, 59-62 |
| Arm Tm difference | <= 8 (repo only) | C | Step 9.2 `MAX_TM_SEPARATION` |
| Sequence | no run of >= 5 identical bases; repo: all 4 bases present | - | SM; screen_homology_complexity.m |
| Structure | no hairpin (stem >= 7, loop 6) or dimer >= 10 in the 32-mer | - | SM; screen_homology_physical.m l.71 |
| Gap between windows | brain >= 5; repo and cancer panel >= 0 | nt | SM; Homology.m `SPACER_LENGTH`; cancer gaps from 0 (computed here) |
| Probes per gene | cortex <= 8, hippocampus <= 16 (random picks); cancer 9-16 stated (Table S10: 7-16, 212/297 genes at 16); repo <= 16, first passing from the 5' end | - | SM; Assembly.m l.461 |
| Hybridization | 2x SSC, 20% formamide, 37 C overnight; 100 nM/probe suggested (paper: ~162, 100, ~32) | - | SM; Library Prep Step 8 |
| Ligation | SplintR 1.25 U/uL, 1x SplintR buffer, 4 C 6 h, then fresh mix 37 C overnight (protocol: 6 h) | - | SM; Step 9 |
| RCA | primer 500 nM, wash buffer, 37 C 2 h; Phi29 1 U/uL, 4 C 6 h, then 30 C overnight + 40 uM aminoallyl-dUTP; BS(PEG)9 cross-link | - | SM; Steps 10-12 |

The SM lists "melting temperature for both 16-mer halves > threshold" among exclusion criteria;
the code and the published arms (primer3 medians 55-57 C at 0.3 M Na+, 100 nM; computed here)
show it is a floor.

## 4. Barcodes and readout

- One logical code (digits 0-3) per gene, shared by all its padlocks.
- Brain (Illumina SBS, 4 rounds): R4 = (R1 + R2 + R3) mod 4; only perfect 4-base matches kept.
  Illumina BC = code mapped 0/1/2/3 -> A/C/G/T (padlock sense) + 3 unread bases that vary between
  a gene's probes (Table S10). Code 2222 (Rgs5) was dropped in silico as a simple repeat.
- Cancer (SOLiD, 5 primers): DNABarcodes, length 7, Hamming distance 3 (one error corrected);
  the same logical code in colour space (SOLiD BC) and base space (Illumina BC).
- Repo: `utilities/prioritizedBarcodeLibrary.csv`, 326 length-7 codes, distance 3; its 58
  priority-1/2 codes keep distance 2 in positions 1-4 (4-round runs). Full code: 7 SBS rounds.
  Assign codes so each round's colours are balanced, weighted by expression (Step 6).
- SOLiD BC depends on the dibase specification (default `embedded_len4`; the cancer panel used the
  sequential one): compute it with the repo, never by hand.

## 5. Target region and specificity

- Window search on the whole RefSeq mRNA, UTRs included (Homology.m l.154-159). For several
  variants take the one sharing most exons (Step 2); > 1.5 kb usually yields 16 probes.
- SM rule: reject a window whose BLAST hit in another transcript is > 12 nt and spans the
  junction by >= 3 nt each side. Repo: blastn plus strand, E <= 0.05, species RefSeq mRNAs; reject
  on an HSP with >= 13 identities covering nt 13-20 (>= 4 nt each side); hits naming "(GENE" or
  the accession (own pseudogenes, e.g. KRT17P1) pass; Step 10 re-BLASTs for manual review.
- One ligated padlock gives one read: filter each probe (`<= 1`), never pool off-targets.

## 6. U-Probe

Nearest template: `padlock_iss` in its direct-RNA form (`padlock-iss-hybiss.md` section 6). Not
`mipseq`: no primer binds the RNA. Copy `padlock_iss`, set `genome` and `targets`, then replace
these blocks (the copy validated with 0 errors, 0 warnings):

```yaml
encoding:
  AKT1: {SOLiD: AGAACCA, Illumina: TGCGACA}   # code 3212010, repo worked example
extracts:
  target_region:
    source: mRNA
    transcripts: mane
    length: 32
    layout:                       # junction search only: nt 7-26
      template: "{flank_5}{junction}{flank_3}"
      parts: {flank_5: {length: 6}, junction: {length: 20}, flank_3: {length: 6}}
probes:
  padlock:
    template: "{arm_5}AAA{solid}TCTCGGGAACGCTGAAGACGGC{illumina}AAA{arm_3}"
    parts:
      arm_5: {expr: "rc(target_region[0:16])"}
      solid: {expr: "encoding[target]['SOLiD']"}
      illumina: {expr: "encoding[target]['Illumina']"}
      arm_3: {expr: "rc(target_region[16:32])"}
hybridization: {temperature_c: 37, na_mM: 390, formamide_pct: 20, probe_nM: 100,
                mg_mM: 0, dntp_mM: 0}
attributes:     # replace all but target_mapped_genes, junction_mapped_genes
  arm_5_gc: {target: padlock.arm_5, type: gc_content}            # arm_3_gc alike
  arm_5_tm: {target: padlock.arm_5, type: annealing_temperature,  # arm_3_tm alike
             hybridization: {na_mM: 300, probe_nM: 100, formamide_pct: 0}}
  arm_tm_diff: {type: difference, of: [arm_5_tm, arm_3_tm]}
summary: {report_name: rna_report, attributes: [arm_5_tm, arm_3_tm, arm_tm_diff,
          target_mapped_genes, junction_mapped_genes]}   # no target_gc, target_tm left
post_process:
  filters:      # drop target_tm; keep the two mapped_genes filters (<= 1)
    arm_gc: {condition: [arm_5_gc >= 0.40, arm_5_gc <= 0.65, arm_3_gc >= 0.40, arm_3_gc <= 0.65]}
    arm_tm: {condition: arm_5_tm >= 50 & arm_3_tm >= 50}
    arm_tm_diff: {condition: arm_tm_diff <= 8}
    homopolymers: {type: sequence_pattern, target: target_region,
                   exclude_patterns: [AAAAA, CCCCC, GGGGG, TTTTT]}
  remove_overlap: {location_interval: 0}   # 5 for the brain panels
  equal_space: {number_desired: 16}        # 8 for the cortex panel
```

Run on a toy genome of the repo's three example mRNAs (Bowtie2 there), it rebuilt the published
AKT1 padlocks at nt 264, 629 and 980 base for base. Use Bowtie 1 (`bowtie` in `align_index`).

**U-Probe departs:** Tm: primer3 floor 50 C at 300 mM Na+/100 nM stands in for oligoprop's 52-55;
it keeps all 334 cortex, 4,348/4,353 cancer and 486/540 hippocampus padlocks (48 keeps all;
computed here). Specificity: Bowtie 1 on the transcriptome (32-mer <= 2, junction 20-mer <= 1
substitutions), not BLAST on RefSeq; add `exempt_copies` to tolerate own pseudogenes. Selection:
`equal_space` spreads sites; the repo takes the first 16 from the 5' end, the brain panels random
ones. No hairpin/dimer screen, no all-four-bases test.

**U-Probe cannot:** generate or assign codes (checksum, DNABarcodes, colour balance); convert
codes to SOLiD dibase barcodes; relax the Tm floor per gene (run short genes in a second protocol);
add the 5' phosphate (order /5Phos/); anything of expansion (LabelX, gels, passivation), sequencing
or decoding (ExSeqProcessing).

## 7. Pitfalls

- Sense arms (HybISS cDNA layout) bind nothing on RNA; arm_5 is rc of the upstream half.
- No 5' phosphate, no circle. Missing RCA primer, no amplicon.
- Codes: low-complexity codes look like artefacts (2222); 4-round codes need the checksum or the
  58-code subset (distance 2); the SOLiD BC must match the specification the sequencing uses.
- Short mRNAs give few probes (3 of 300 cancer genes dropped as too short).

## 8. Verification checklist

- [ ] Each padlock = 16 + AAA + 7 + backbone 22 + 7 + AAA + 16 = 74 nt, ordered /5Phos/.
- [ ] rc(arm_5) + rc(arm_3) equals the 32-nt window on the intended transcript.
- [ ] Arm GC 40-65%, arm Tm floor met with its model and conditions recorded, |dTm| <= 8.
- [ ] `target_mapped_genes` and `junction_mapped_genes` = 1 with Bowtie 1; off-targets reviewed.
- [ ] Codes unique, distance >= 3 (7 rounds) or checksum/distance 2 (4 rounds); Illumina BC =
      code mapped A/C/G/T; SOLiD BC from the repo; all padlocks of a gene share the code.
- [ ] Probes per gene (8-16) and spacing recorded; shortfalls named.

## 9. References

- Alon S, Goodwin DR, Sinha A, et al. 2021, Science 371:eaax2656. Expansion sequencing: spatially
  precise in situ transcriptomics in intact biological systems. doi:10.1126/science.aax2656 (PMC
  author manuscript PMC7900882 with Supplementary Materials (SM) text and Tables S9/S10 used).
- Sinha A, Chen F, Alon S, Boyden E 2021, protocols.io. Targeted ExSeq -- Probe Generation V.2.
  doi:10.17504/protocols.io.byrbpv2n
- Sinha A, Wassie AT, Chen F, Cui Y, Boyden E 2021, protocols.io. Targeted ExSeq -- Sequencing
  Library Preparation. doi:10.17504/protocols.io.bgp7jvrn
- Cui Y, Sinha A, Boyden E 2021, protocols.io. Targeted ExSeq -- In Situ Sequencing (Illumina
  Chemistry). doi:10.17504/protocols.io.bgs4jwgw; collection doi:10.17504/protocols.io.b2e4qbgw
- Sinha A, Cui Y, Alon S, Wassie AT, Chen F, Boyden E 2021, protocols.io. Targeted ExSeq --
  Tissue Preparation V.2. doi:10.17504/protocols.io.b2emqbc6
- anusinha/TargetedExSeqProbeGen, commit 2fe7748 (Homology.m, Assembly.m, screen_homology_*.m,
  utilities/). https://github.com/anusinha/TargetedExSeqProbeGen
