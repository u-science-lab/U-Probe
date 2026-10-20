# MiP-seq: padlock + RCA initiator primer, pairwise sequencing

> Read this when: designing or checking MiP-seq probes (a padlock ligated directly on RNA plus an
> RCA initiator primer that bridges RNA and circle), editing the `mipseq` template, or comparing
> it with `double_hyb_rca` (double-hyb-rca.md) and STARmap SNAIL (starmap.md). General rules:
> `references/concepts/probe-architecture.md` (padlock + primer units), `references/concepts/specificity.md`.

## 1. What it is and when to choose it

Multi-omics in situ pairwise sequencing (Wu et al. 2024). Per target site, two oligos bind one ~40-nt
RNA window:
- a **padlock** whose two 13-nt arms meet on the RNA and are joined by **SplintR** (PBCV-1 DNA ligase,
  which ligates DNA ends on an RNA splint; no reverse transcription);
- an **RCA initiator primer** whose 5' part binds the RNA just downstream of the padlock and whose
  3' end binds the padlock backbone, so Phi29 amplifies only a circle that sits next to its primer.

The padlock carries two 3-nt barcodes (B1, B2) read at the same time by sequencing by ligation (SBL)
from two anchor sites: 10 colour combinations per round (4 single + 6 two-colour), so N rounds
decode 10^N codes instead of 4^N (Wu 2024, peer-review file). Proteins (antibody-conjugated oligos),
DNA loci and small molecules (aptamers) are read with the same padlock/RCA chemistry.

Choose it for: targeted panels of tens to ~1,000 genes with few sequencing rounds; direct RNA
detection without RT; SNV / allele / m6A readouts at the ligation junction.
Choose instead: `padlock_iss` (padlock-iss-hybiss.md) for cDNA-based HybISS; `starmap`
(starmap.md) when ligation should be templated by DNA (the primer), not by RNA;
`double_hyb_rca` (double-hyb-rca.md) only for U-Probe's 17-nt-barcode RCA layout.

## 2. Probe architecture

```text
RNA window (sense), nt 1..40, 5'->3':
  [arm5 site nt 1-13][arm3 site nt 14-26][gap nt 27][primer site nt 28-40]

(a) Oligos, 5'->3'
Padlock, 58 nt, 5'-phosphate (T4 PNK):
5'P-[arm_5 13]-CC-[B1 3]-TGCGTCTATTT-GT-[B2 3]-TAGTGGAGCCT-[arm_3 13]-OH 3'
     rc(nt 1-13)  2  code  anchor-L 11  2  code  anchor-R 11   rc(nt 14-26)
                                              (= primer tail site)
RCA initiator primer, 25 nt in the template (published: 25 in 247 of 432, 26 in 185):
5'-[RNA part 13]-[s 1]-AGGCTCCACTA-3'      (26-nt form: one C between s and the tail)
    rc(nt 28-40)  nt 27  rc(TAGTGGAGCCT), 11 nt: pairs with the padlock backbone
                  (sense base, unpaired spacer)

(b) Bound (RNA written 3'<-5' so the oligos read 5'->3' above it)
      primer                                padlock 3' arm         padlock 5' arm
 5'-[rc nt 28-40]-s-AGGCTCCACTA-3'       5'-[rc nt 14-26]-OH3' ^ 5'P-[rc nt 1-13]-3'
    |||||||||||||    (to backbone)           |||||||||||||   junction  |||||||||||||
 3'-[nt 40 .. 28]-[27]------------------- [nt 26 ....... 14]---------[nt 13 ..... 1]-5'  RNA
```

- Padlock 5' end (phosphate) sits on nt 13, its 3'-OH on nt 14: the **ligation junction is between
  RNA nt 13 and 14**, inside the 26-nt arm site. arm_5 (5' end of the oligo) binds the 5'-most RNA
  segment; arm_3 binds the segment downstream of it.
- The backbone loop leaves the RNA at nt 26 through TAGTGGAGCCT, right where the primer's RNA part
  ends (nt 28): the primer's 3' tail AGGCTCCACTA pairs with it, and the **primer 3' end sits on the
  circle**, which Phi29 copies by strand displacement.
- The rolling-circle product (RCP) is rc(circle) in tandem. SBL: 5'-phosphorylated anchors
  TGCGTCTATTT (anchor-L) and TAGTGGAGCCT (anchor-R) hybridize to the RCP; fluorescent query 6-mers
  5'-dye-NNNNN[b]-3' are ligated to the anchors' 5' phosphate, so the query's 3' base reads the
  barcode base next to the anchor (B1 and B2 base 3). Anchors shifted by N / NN (N-TGCGTCTATT,
  NN-TGCGTCTAT, N-TAGTGGAGCC, NN-TAGTGGAGC) read bases 2 and 1. Dyes: A Alexa488, G Alexa546,
  C Alexa594, T Alexa647 (Supplementary Table, sheet "In Situ Sequencing probe").
- Backbones are modular: the published GAPDH/MTOR sets also use hybridization-readout backbones
  (e.g. ...ACTAGCTCGATTCAA-GT-ACTAGCTCGATTCAA... with a Cy3 detection oligo) and primers whose
  tail pairs that backbone instead of TAGTGGAGCCT.

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Padlock arm length | 12-16 (Methods); 13 + 13 in all 433 58-nt padlocks of the 217-gene set | nt | Wu 2023 preprint Methods, doi:10.1101/2023.01.07.523058; Wu 2024 Supp. Table (sheet "Fig 3"), doi:10.1038/s41551-024-01205-7 |
| Gap between arms | 0 (arms abut; no gap fill) | nt | Th-1/Th-2 padlocks mapped to mouse Th NM_009377.2 |
| Padlock length | 58 | nt | Wu 2024 Supp. Table |
| Backbone | CC-B1-TGCGTCTATTT-GT-B2-TAGTGGAGCCT (32) | nt | Wu 2024 Supp. Table, e.g. Th-1-PLP-1 |
| Barcodes B1, B2 | 3 each; 217 genes, 217 distinct per-position colour codes | nt | Wu 2024 Supp. Table |
| Coding capacity | 10^N for N rounds (3 rounds = 1,000) | codes | Wu 2024 peer-review file |
| Primer RNA part | 12-16 (Methods); 13 in 417 of the 430 mapped pairs, 14 in 13 | nt | Wu 2023 Methods; this KB, 430 units mapped to mouse Ensembl 113 |
| Gap padlock site -> primer site | 0-1 (1 in 417 of 430) | nt | mapping above; template: 1 |
| Primer spacer before tail | the sense base of nt 27 alone in 234 pairs (25-nt primers), that base + C in 163 (26-nt), C alone in 14, another 2-nt spacer in 13 | nt | mapping above; template: the sense base alone |
| Primer tail | AGGCTCCACTA (all 432 primers) | nt | Wu 2024 Supp. Table |
| Primer length | 25 (247 primers) or 26 (185); the template builds 25 | nt | Wu 2024 Supp. Table |
| Sites per gene | 2, "in an exon close to the 5' end"; 8 pairs raised detection (Supp. Fig. 7) | pairs | Wu 2023 Methods; Wu 2024 SI |
| Site position (no number published) | median start nt 885 (quartiles 455, 1,650); 54% within the first 1,000 nt, 71% within 1,500; template cap 1,000 | nt | this KB, 367 sites on the mouse canonical transcript |
| Tm / GC criteria | none stated ("self-developed" program) | - | Wu 2023 Methods |
| Arm Tm, estimate | Th-1 arms 37-39 C; primer RNA part 47 C (2x SSC, 10% FA, 30 nM) | C | U-Probe/primer3 estimate, not a published criterion |
| Hybridization | 30 nM each padlock + equimolar primer, 2x SSC, 10% formamide, 20 mM RVC, 37 C overnight | - | Wu 2023 Methods |
| Ligation (RNA) | SplintR 1 U/uL, 0.2 U/uL RNase inhibitor, 37 C, 2 h | - | Wu 2023 Methods |
| RCA | Phi29 1 U/uL, 30 C, 2 h | - | Wu 2023 Methods |
| SBL | T4 DNA ligase 0.2 U/uL, 1 uM anchor + 1 uM query, 25 C, 2 h | - | Wu 2023 Methods |
| DNA targets | denature in 70% formamide/2x SSC at 85 C; ligate with T4 DNA ligase | - | Wu 2023 Results/Methods |

**Protocols disagree:** SplintR fidelity. Wu 2024 (SI Fig. 1) shows SplintR rejecting a mismatch next
to the ligation point for ACTB; Krzywkowski & Nilsson 2017 found "an overall limited end-joining
fidelity" of PBCV-1 ligase on RNA, with mismatch-specific variation. Do not rely on SplintR alone to
discriminate a SNV; validate with wild-type/mutant controls.

## 4. Target region and template

- Template: the **RNA itself** (no RT). Arms carry rc(RNA), i.e. antisense sequence.
- Region: Methods choose "the exon close to the 5' end"; CDS and UTR are not distinguished. Of the
  published padlock sites that map to mouse RefSeq mRNAs, 290/380 lie in the CDS, 64 in a UTR, and
  343 lie in every RefSeq variant of their gene (template meta). "Close to the 5' end" is not a
  number: on the mouse canonical transcripts the 367 mapped sites start at a median of nt 885, 54%
  within the first 1,000 nt. The template turns it into `target_start <= 1000` (U-Probe's choice);
  raise the cap for long transcripts, and check `summary.json` for `partial`.
- Isoforms: choose sites shared by the isoforms you want to count; U-Probe's `transcripts: mane`
  designs on one reference transcript, and with one transcript windows can span exon junctions
  (they then bind only spliced RNA).
- SNV / allele designs: put the variant at the ligation junction, i.e. under the 3'-terminal base
  of arm_3 (pairs nt 14) or the 5'-terminal base of arm_5 (pairs nt 13); Wu 2024 tested "a single
  unmatched base adjacent to the ligation point". U-Probe has no variant-aware option; build such
  probes by hand or from a direct-sequence target with a fixed window.

## 5. Specificity

- One padlock-primer unit gives one RCP, which is detected: there is no multi-probe averaging as in
  smFISH. Keep every unit clean (`target_mapped_genes <= 1`).
- A signal needs (1) both arms hybridized and ligated on the same molecule, and (2) the primer bound
  next to the circle. An unligated padlock is not amplified; a primer alone has no circle; the 11-nt
  tail alone (Tm about 30 C under the hybridization conditions, primer3 estimate) does not hold the
  primer at 37 C without its RNA part. An off-target is therefore a transcript holding the 26-nt arm
  site and the primer site together. The template maps the whole 40-nt `target_region` (the union).
- Short arms help: a mismatch in a 13-nt arm with Tm near the 37 C hybridization strongly reduces
  binding; junction mismatches also reduce ligation (limited with SplintR, see section 3).
- Bowtie2 on 40-nt windows (template default): flagged 97% of random 40-nt mRNA windows that have a
  site within 2 substitutions in another gene and listed 87% of those genes; it also reports partial
  matches (about 25 nt or more, e.g. the arm site alone), which is conservative. For an exhaustive
  search within 2 substitutions use `aligner: bowtie` (Bowtie 1, needs `bowtie` in `align_index`).
- Expression-aware counting: add `gene_expression` so unexpressed off-target genes are not counted;
  add `exempt_copies` for pseudogene-rich targets. Repeats: a window in a repeat hits many genes and
  is removed by `<= 1`.

## 6. U-Probe template

```bash
uprobe templates show mipseq -o protocol.yaml
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets
uprobe run -p protocol.yaml -g genomes.yaml -o results --threads 8
```

Key lines (uprobe/templates/mipseq.yaml):

```yaml
encoding:
  Th: {B1: GAG, B2: CAT}            # replace per gene; published Th-1-PLP-1 values
extracts:
  target_region:
    source: mRNA                     # whole mRNA; target_start counts from its 5' end
    transcripts: mane
    length: 40
    layout:
      template: "{arm5_site}{arm3_site}{gap}{primer_site}"
      parts: {arm5_site: {length: 13}, arm3_site: {length: 13}, gap: {length: 1}, primer_site: {length: 13}}
probes:
  padlock:
    template: "{arm_5}CC{b1}TGCGTCTATTTGT{b2}TAGTGGAGCCT{arm_3}"
    parts: {arm_5: {expr: "rc(target_parts['arm5_site'])"}, arm_3: {expr: "rc(target_parts['arm3_site'])"}, ...}
  rca_primer:
    template: "{binding}{gap_base}AGGCTCCACTA"     # 25 nt; the 26-nt form adds a C before the tail
    parts: {binding: {expr: "rc(target_parts['primer_site'])"}, gap_base: {expr: "target_parts['gap']"}}
hybridization: {temperature_c: 37, na_mM: 390, probe_nM: 30, formamide_pct: 10, dntp_mM: 0}
post_process:
  filters:
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    near_5p: {condition: target_start <= 1000}     # "an exon close to the 5' end"
  remove_overlap: {location_interval: 0}
  equal_space: {number_desired: 2}                 # 2 sites per gene
```

What to edit:
- `genome`, `targets`, and `encoding` (B1/B2 per gene; see pitfalls for code collisions).
- Sites per gene and region: `near_5p` plus `equal_space` give the paper's 2 sites near the 5' end.
  Raise the 1,000-nt cap (U-Probe's number, not the paper's) for long transcripts, and raise
  `number_desired` for more pairs (8 raised detection in SI Fig. 7).
- Primer length: 25 nt here; for the 26-nt form write `"{binding}{gap_base}CAGGCTCCACTA"`.
- Optional Tm/GC filters on `arm_5_tm`, `arm_3_tm`, `primer_tm`, `target_gc` (reported, not
  filtered): none are published, so choose them yourself and say so.
- `aligner: bowtie` + `max_mismatches: 2` for an exhaustive off-target search.

**U-Probe departs:** (1) `source: mRNA` on the MANE/canonical transcript with a 1,000-nt cap, not
"an exon near the 5' end" (no distance is published); (2) `equal_space` spreads the 2 sites evenly
over the passing windows below the cap, while the paper's two sites are a median 223 nt apart;
(3) the specificity filter is U-Probe's approximation, the paper states no criterion.

## 7. Pitfalls

- Orientation: arms are rc(RNA). Do not swap arms: arm_5 must bind nt 1-13 (5' side of the window),
  or the ends do not meet. The primer's RNA part must be at its 5' end and the tail at its 3' end.
- Barcode collisions: decoding sees an unordered base pair per position. B1=GAG/B2=CAT and
  B1=CAG/B2=GAT give the same colours ({C,G},{A},{G,T}); check uniqueness of the per-position sets.
  B1 = B2 gives a single colour every round (12 published genes); avoid codes that look like a
  missing barcode in your decoder.
- Order the padlock 5'-phosphorylated (or phosphorylate with T4 PNK); SplintR needs it.
- SplintR ligation fidelity is limited (section 3); T4 DNA ligase ligates poorly on RNA splints
  (Lohman 2014), and the authors chose SplintR after comparing ligases (SI Fig. 1).
- SplintR is "strongly inhibited" above 100 mM NaCl and slower when the 5'-phosphorylated base
  (padlock 5' end, pairing nt 13) is dC or dG (Lohman 2014): ligate in ligase buffer after washing
  out 2x SSC; prefer windows with A/T at the padlock 5' end if you rank by hand (not a MiP-seq rule).
- 13-nt arms have Tm near 37 C: high-AT windows may hybridize poorly; consider filtering
  `arm_5_tm`/`arm_3_tm` against the 37 C hybridization temperature.
- Bowtie2 reports at most 100 alignments; windows in repeats or gene families may be under-counted.
- The primer's spacer base is the *sense* nt 27, so it cannot pair with nt 27; it is a spacer, not
  an extra binding base. If you write primers by hand, keep this (25 nt, the published majority),
  add a C after it (26 nt) or use C alone (published variants).
- A gene whose first 1,000 nt hold no passing window gets no probes with the `near_5p` filter: the
  run is `partial` (`summary.json`), so raise the cap for that gene's panel rather than ignoring it.

## 8. Verification checklist

- [ ] Padlock = 58 nt: 13 + CC + 3 + TGCGTCTATTT + GT + 3 + TAGTGGAGCCT + 13; 5' phosphate noted.
- [ ] rc(arm_5) + rc(arm_3) equals RNA nt 1-26 of the window, contiguous (junction 13|14).
- [ ] Primer 5' part = rc(RNA nt 28-40) (or 27-40); primer ends with AGGCTCCACTA (3' end on circle)
      and is 25 nt (or 26 with the extra C).
- [ ] Window lies in the intended transcript region/isoforms; MANE/canonical choice logged.
- [ ] `target_mapped_genes` = 1 for every unit (with the aligner and expression table you intended).
- [ ] Per-position B1/B2 colour codes are unique across genes; dyes match the SBL scheme.
- [ ] Number of sites per gene matches the plan (2 near the 5' end, or more); sites do not overlap.
- [ ] Conditions recorded: 2x SSC, 10% FA, 37 C hybridization; SplintR; Phi29 at 30 C.

## 9. References

- Wu X, Xu W, Deng L, et al. 2024, Nature Biomedical Engineering. Spatial multi-omics at subcellular
  resolution via high-throughput in situ pairwise sequencing. doi:10.1038/s41551-024-01205-7
  (Supplementary Information, Supplementary Table of probe sequences, peer-review file).
- Wu X, Xu W, Deng L, et al. 2023, bioRxiv. Development of Multiomics in situ Pairwise Sequencing
  (MiP-Seq) for Single-cell Resolution Multidimensional Spatial Omics. doi:10.1101/2023.01.07.523058
- Lohman GJS, Zhang Y, Zhelkovsky AM, Cantor EJ, et al. 2014, Nucleic Acids Research. Efficient DNA ligation in
  DNA-RNA hybrid helices by Chlorella virus DNA ligase. doi:10.1093/nar/gkt1032
- Krzywkowski T, Nilsson M. 2017, Nucleic Acids Research. Fidelity of RNA templated end-joining by
  chlorella virus DNA ligase and a novel iLock assay with improved direct RNA detection accuracy.
  doi:10.1093/nar/gkx708
- Wang X, Allen WE, Wright MA, et al. 2018, Science. Three-dimensional intact-tissue sequencing of
  single-cell transcriptional states (SNAIL). doi:10.1126/science.aat5691
- Zhang Q, Cai H, Zhang J, et al. 2026, openRxiv (bioRxiv). Universal agentic probe design for
  imaging-based spatial-omics (U-Probe). doi:10.64898/2026.04.12.717982
