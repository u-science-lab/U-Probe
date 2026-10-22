# RAEFISH: reverse padlock + priming probe + encoding probe

> Read this when: designing or checking RAEFISH (Reverse-padlock Amplicon Encoding FISH) sets: a
> reverse padlock and a priming probe on adjacent 30-nt RNA sites, plus an encoding probe on the
> RCA amplicon. Related: `starmap.md` (padlock + primer `layout`), `mip-seq.md` (primer tail on
> the padlock), `merfish.md` (MHD4 codes, readouts, pool amplification, crowding).

## 1. What it is and when to choose it

Cheng 2025. A reverse padlock (80 nt) and a priming probe (60 nt) bind adjacent 30-nt sites of one
RNA. The padlock's invariant ends face away from the RNA: a common splint lets T4 DNA ligase close
the circle, a toehold oligo removes the splint, and the priming probe's 10-nt 3' tail, paired with
the padlock, primes phi29 RCA; "non-specific binding of the reverse-padlock probe or the linear
priming probe alone cannot initiate RCA" (Methods). The amplicon copies the RNA site; an encoding
probe binds it with four readout sites, the on-bits of a 94-bit weight-4 MHD4 code read in 47
two-colour rounds (adapter oligos + two common dye oligos).

Choose it for transcriptome-scale targeted panels (23,312 human, 21,955 mouse genes) and short RNAs
(shortest target 78 nt): all oligos have invariant ends, so the three sub-libraries are amplifiable
pools. Costs: detection 8.3% of MERFISH's (A549, Fig. S1D), 4.2 days of imaging for 6.9 mm2.
Smaller panels: `starmap.md`, `mip-seq.md`; long genes, higher efficiency: `merfish.md`.

## 2. Probe architecture

```text
RNA 5'-[ padlock site 30 ][gap 0-1][ priming site 30 ]-3'                  window 60-61 nt
Padlock 80 nt:  5'P-ACTCTCACTTGGGTTCCAAT-[rc(padlock site)]-TTCCTCTAACTCTCTAATCCACTAATACCC-3'
                    nt 1-6: splint; 11-20: priming tail              nt 75-80: splint
Priming 60 nt:  5'-GTGGTAAAGCTCCGCGGCTT-[rc(priming site)]-ATTGGAACCC-3'  (tail = rc nt 11-20)
Splint   GAGAGTGGGTATA*T*C/3InvdT/  nt 1-6 pair padlock 1-6, 7-12 pair 75-80, ATC overhang
Toehold  GA+TAT+ACC+CAC+TC+TC (+ = LNA)  full complement of the splint's 15 nt
Encoding 150 nt as hybridized (a < b < c < d = the gene's on-bits; Ad = adapter nt 1-20):
5'-TTGTCGCACGTTCGGTGTCG-rc(Ad_d)-rc(Ad_c)-[rc(padlock site)]-rc(Ad_b)-rc(Ad_a)-
   ACGTTGCACCATAACGACGC-3'
```

- Target-dependent: the three binding parts and the four readout sites. The encoding binding part
  equals the padlock's ("a 30-nt sequence same as the RNA-binding region on the corresponding
  reverse padlock probe"): RCA copies the circle, so the amplicon holds the RNA-sense site. Adapter
  nt 21-60 carry two sites for the 647 or the 750 dye oligo. All else is constant in every Table S1
  row; the pairings above are computed here.
- Padlock site 5' of the priming site, gap 0 (82% human, 63% mouse) or 1, in 59,700 of the 59,721
  sets located on Ensembl 113 transcripts (computed here).
- Pool templates (Table S1): padlock as above (the phosphorylated RT primer adds the 5'-P); priming
  + TATAGTGAGTCGTATTATCCGTACGGCGATTGGTCAATGCGGAACGA (rc T7 + PCR site) = 107 nt; encoding =
  rc(probe), "because our probe amplification procedure reverse-complements this sequence".

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Padlock / priming site | 30 / 30, "without overlap ... no or a 1-nt gap" | nt | Probe design |
| Code | 94 bits, weight 4, MHD4 (checked in Table S1); 47 rounds x 2 colours | - | Probe design |
| Sets per gene | human 1; mouse 2, one set twice if only one passed (4,959 genes) | sets | Probe design |
| Panel | GENCODE v46 / vM35 coding + lncRNA; 40 / 26 top-expressed genes left out | - | Probe design |
| GC, Tm, off-target, isoform rules | not stated | - | Methods |
| Published parts, GC | 53-90%; >= 57% in all but 28 of 124,196 | % | computed here |
| Published parts, Tm (U-Probe; 2x SSC, 50% FA, 3 nM) | >= 41.6 (medians 43-44); random 30-mers: 76% < 41.5 | C | computed here |
| Homopolymer | no run >= 6 nt (2.0% of GC-matched random 30-mers have one) | nt | computed here |
| Primary hybridization | 3 nM each oligo, 2x SSC, 50% formamide, 37 C, 16 h | - | protocols |
| Ligation; splint removal | 1 uM splint, T4 ligase, RT 2 h (tissue 16 h); 800 nM toehold, 37 C | - | protocols |
| RCA | phi29, 30 C, 2 h (tissue 3 h) | - | protocols |
| Encoding hybridization | 3 nM each, 2x SSC, 40% FA, 10% dextran sulfate, 37 C | - | protocols |

## 4. Target region

- RNA itself, no RT. Region unstated. Of the 19,733 human windows on the MANE/canonical Ensembl 113
  transcript, 9,037 lie in CDS, 3,455 in 5' UTR, 1,968 across start/stop, 797 in 3' UTR, 4,476 on
  non-coding RNA; 18% span an exon junction; 2,903 more are not on it (computed here).
- 5'-biased: starts at 3/14/36% of the reference transcript (quartiles). Mouse second sets never
  overlap the first (median 106 nt apart).
- Rows `<gene>NegativeControl` (40 human, 26 mouse: the excluded genes) carry codes and RNA-SENSE
  binding parts (all 40 human found verbatim in their mRNA), i.e. non-binding controls (computed
  here; not described in the text).
- Perturb-RAEFISH (gRNAs): 40-nt padlock site (20-nt spacer + 20 nt of the U6 promoter), one
  shared priming probe on the scaffold (45 nt), 14-bit HD2 codes, 7 rounds; nothing for U-Probe.

## 5. Specificity

- False spot = padlock and priming probe adjacent on one molecule: map the 60-61-nt window (sense),
  keep `target_mapped_genes <= 1`.
- Wrong identity: encoding probes meet every amplicon, each a copy of a padlock site. A padlock site
  within 2 substitutions of another gene's RNA could light that gene's amplicons with this gene's
  bits; `padlock_site_genes <= 1` excludes it. Both filters are U-Probe's; the paper states none.
- Under-detected genes had higher Tm and GC (Fig. S1H-I); structure or cross-hybridization is the
  authors' hypothesis (Limitations).
- MHD4 corrects 1-bit errors; foci on unassigned codes (of 25,281) estimate the error rate.

## 6. U-Probe: start from `starmap`

`starmap` pairs a padlock site and a primer site via `layout`, both `rc()`, and maps the whole
window (`mipseq` ligates on the RNA instead). Copy it and replace these blocks:

```yaml
genome: GRCh38_e113                 # bowtie in align_index
targets: [DPM1]
encoding:                           # adapter nt 1-20 of the 4 on-bits, ascending bit order
  DPM1:
    A1: GACTTTGTGGGTACGTTGAT
    A2: GACACGTCCCTCAATCGAAC
    A3: TTGGACCGGCATGTCAAGAG
    A4: TGGACTTCCGGTCGATATCC
extracts:
  target_region:
    source: mRNA
    transcripts: mane               # canonical outside human
    length: [60, 61]
    step: 1
    layout:
      template: "{padlock_site}{gap}{priming_site}"
      parts: {padlock_site: {length: 30}, gap: {length: [0, 1]}, priming_site: {length: 30}}
probes:
  padlock:                          # order with 5'-phosphate
    template: "ACTCTCACTTGGGTTCCAAT{binding}TTCCTCTAACTCTCTAATCCACTAATACCC"
    parts: {binding: {expr: "rc(target_parts['padlock_site'])"}}
  priming:
    template: "GTGGTAAAGCTCCGCGGCTT{binding}ATTGGAACCC"
    parts: {binding: {expr: "rc(target_parts['priming_site'])"}}
  encoding_probe:
    template: "TTGTCGCACGTTCGGTGTCG{r4}{r3}{binding}{r2}{r1}ACGTTGCACCATAACGACGC"
    parts:
      binding: {expr: "rc(target_parts['padlock_site'])"}
      r1: {expr: "rc(encoding[target]['A1'])"}
      r2: {expr: "rc(encoding[target]['A2'])"}
      r3: {expr: "rc(encoding[target]['A3'])"}
      r4: {expr: "rc(encoding[target]['A4'])"}
  priming_pool:                     # Table S1 "Priming", 107 nt
    template: "{p}TATAGTGAGTCGTATTATCCGTACGGCGATTGGTCAATGCGGAACGA"
    parts: {p: {expr: "priming"}}
  encoding_pool:                    # Table S1 "Encoding"
    template: "{e}"
    parts: {e: {expr: "rc(encoding_probe)"}}
hybridization: {temperature_c: 37, na_mM: 390, probe_nM: 3, formamide_pct: 50, dntp_mM: 0}
attributes:
  padlock_tm: {target: padlock.binding, type: annealing_temperature}
  priming_tm: {target: priming.binding, type: annealing_temperature}
  target_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie,
                        max_mismatches: 2, strand: sense}
  padlock_site_genes: {target: target_parts.padlock_site, type: mapped_genes, aligner: bowtie,
                       max_mismatches: 2, strand: sense}
post_process:
  filters:
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    padlock_site_genes: {condition: padlock_site_genes <= 1}
    padlock_tm: {condition: padlock_tm >= 41.5}
    priming_tm: {condition: priming_tm >= 41.5}
    homopolymers: {type: sequence_pattern, target: target_region,
                   exclude_patterns: [AAAAAA, CCCCCC, GGGGGG, TTTTTT]}
  remove_overlap: {location_interval: 0}
  equal_space: {number_desired: 1}  # 2 for a mouse-style library
  sorts: {is_ascending: [target_mapped_genes]}
```

- On the published windows of DPM1, FGR (gap 0), HS3ST1 and FKBP4 (gap 1) these templates rebuild
  Table S1's three strings exactly (computed here).
- `encoding`: your codebook's on-bits mapped to adapters (Table S1, sheet "94-bit RAEFISH readout
  oligos"; DPM1 = bits 1, 92, 93, 94). Pools: order `padlock`, `priming_pool`, `encoding_pool`;
  single oligos: `padlock` (5'-P), `priming`, `encoding_probe`.
- `equal_space: 1` keeps the 5'-most passing window, like the published bias; 2 keeps the 5'-most
  and 3'-most. A mouse gene with one window was ordered twice in the paper.
- **U-Probe departs:** MANE/canonical mRNA (isoform rule unpublished); the Tm floor (recompute it if
  `probe_nM` or formamide change) and the homopolymer filter reproduce the published sets, not a
  stated rule; both specificity filters.

U-Probe cannot do:
- codes: no 94-bit weight-4 MHD4 codebook generator; use published codes or build one outside
  (`merfish.md`) and keep unassigned codes as blanks;
- panel choice and exclusion of crowding genes;
- common oligos (splint, toehold, adapters, dye and blocking oligos, pool primers): copy Table S1;
- modifications (5'-P, phosphorothioate, LNA, inverted dT) and pool amplification;
- a direct check of encoding probes against the panel's amplicons (`padlock_site_genes` stands in).
  Sense negative controls: a second run with `expr: "target_parts['padlock_site']"` (no rc).

## 7. Pitfalls

- Table S1 "Encoding" is rc of the binding probe: synthesized as is, it misses the amplicon.
  "Priming" ends in the T7/PCR tail: as a single oligo its 3' end is off the padlock, no priming.
- Encoding binding part = rc(padlock site), never the priming site: only the padlock is copied.
- Keep `{padlock_site}{gap}{priming_site}`; swapped, the priming tail and padlock nt 11-20 sit at
  opposite ends of the window.
- Highly expressed genes crowd their rounds and merge into single foci (undercount, Limitations);
  drop them or image them separately (`merfish.md` section 5).
- Bowtie 1 must be installed and in `align_index`; Bowtie2 misses sites of a 30-nt query.

## 8. Verification checklist

- [ ] Window = padlock site + 0-1 nt + priming site on the intended transcript, padlock site 5'.
- [ ] Padlock 80 nt with both constant ends and 5'-P; priming 60 nt ending ATTGGAACCC.
- [ ] `encoding_probe.binding == padlock.binding`; readout sites = rc(adapter nt 1-20) of the
      on-bits; `encoding_pool == rc(encoding_probe)`.
- [ ] Codes unique, weight 4, no two genes sharing 3 on-bits; blanks kept.
- [ ] `target_mapped_genes == 1`, `padlock_site_genes == 1`; Tm given at 2x SSC, 50% FA, 3 nM.
- [ ] 1 set per gene (2 mouse-style); shortfalls named.
- [ ] Run last (the paper states no GC, Tm, off-target or isoform rule: the Tm floor 41.5 C at 2x SSC, 50 % FA, 3 nM is U-Probe's calibration on the published parts; published parts have GC 53-90 %):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --summary results/run1/summary.json \
    --transcripts <genome dir>/transcript.fa \
    --binding padlock.binding --binding priming.binding --tm 41.5 100 --tm-column padlock_tm --tm-column priming_tm \
    --min-gap 0 --min-probes 1
```

## 9. References

- Cheng Y, Dang S, Zhang Y, et al. 2025, Cell 188:6953-6970.e12. Sequencing-free whole-genome
  spatial transcriptomics at single-molecule resolution. doi:10.1016/j.cell.2025.09.006. Read: PMC
  author manuscript PMC12662569 (Results; STAR Methods: Probe design, Probe synthesis, RAEFISH
  protocols, Image analysis; Fig. S1; Key resources table) and Table S1 "All oligos" (mmc1.xlsx).
- Preprint bioRxiv doi:10.1101/2025.03.06.641951 and Mendeley Data doi:10.17632/8kbv637pxh (code,
  codebooks): not used.
