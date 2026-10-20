# Double-hybridization RCA (circle probe + amplification primer)

> Read this when: using or adapting U-Probe's `double_hyb_rca` template, deciding between it,
> `mipseq` (mip-seq.md) and `starmap` (starmap.md), or checking which way round the amplification
> primer must be built (its 3' end must sit on the circle, not on the RNA). General rules:
> `references/concepts/probe-architecture.md` (several probes per window), `references/concepts/specificity.md`.

## 1. What it is and when to choose it

`double_hyb_rca` is the RCA probe configuration of the U-Probe preprint (Zhang et al. 2026, which
shares authors with MiP-seq), shipped as `tests/data/double_hyb_rca.yaml` and
ported to a template. The preprint describes it as "a circle probe with target-complementary arms
flanking barcodes, plus an amplification probe referencing both the circle probe and target region".
Both oligos hybridize to the RNA ("double hybridization"):
- the **circle probe** (padlock) has two 13-nt arms that meet on the RNA and two 17-nt barcodes;
- the **amplification probe** (primer) binds the RNA next to the circle with its 5' part and binds
  barcode 2 of the circle with its 3' part, so it can prime rolling circle amplification (RCA).

The two-oligo geometry (13 + 13 arms, one skipped nt, 13-nt primer site in a 40-nt window) is that of
MiP-seq; the backbone is not. No published experiment, readout scheme or protocol is attached to
this exact layout (as far as could be checked), so treat it as a **design scaffold**, not a
validated protocol.

| Design | Ligation template / ligase | What the primer's 3' end binds | RNA sites | Readout |
| --- | --- | --- | --- | --- |
| STARmap SNAIL | primer tail (DNA splint, 6 + 6 nt) / T4 DNA ligase | padlock ends it just splinted | padlock 19-25 + primer 19-25 nt, gap 0-2 | SEDAL sequencing of a gene code |
| MiP-seq | RNA (13 + 13 arms) / SplintR | constant 11-nt backbone site TAGTGGAGCCT | 26 + 13 nt, gap 0-1 | pairwise SBL of B1/B2 (3 nt) |
| double_hyb_rca | RNA (13 + 13 arms) / not stated (needs an RNA-templated ligase, e.g. SplintR) | gene-specific barcode 2 (17 nt) | 26 + 13 nt, gap 1 | not stated (17-nt barcodes) |
| PLAYR | two probe tails template a backbone + insert / T4 DNA ligase | circle, from the 3'-OH of one of the two probes | two adjacent probes | insert-specific detection oligo |

Decision rule:
- Validated chemistry and SBL decoding with few rounds -> `mipseq`.
- Ligation templated by DNA (T4 ligase), longer 19-25-nt sites, STARmap workflow -> `starmap`.
- Single padlock on cDNA (RT first), HybISS readout -> `padlock_iss` (padlock-iss-hybiss.md).
- `double_hyb_rca` only when you want this 17-nt-barcode layout (e.g. your own hybridization
  readout sequences) or a starting point for a custom two-oligo RCA design; then fix the items in
  sections 6-7 first.

## 2. Probe architecture

```text
RNA window (sense), nt 1..40, 5'->3':
  [arm site A: nt 1-13][arm site B: nt 14-26][nt 27, unbound][primer site: nt 28-40]

(a) Oligos, 5'->3'
circle_probe, 61 nt (needs a 5' phosphate for ligation; the template does not mark it):
5'-[part1 13]-[barcode1 17]-N-[barcode2 17]-[part3 13]-3'
   rc(nt 1-13)                                rc(nt 14-26)
amp_probe (template order "{part2}N{part1}"), 31 nt:
5'-[part2 13]-N-[part1 17]-3'
   rc(nt 28-40)  rc(barcode2): pairs with the circle

(b) Bound (RNA written 3'<-5' so the oligos read 5'->3' above it)
      amp_probe                               circle 3' arm          circle 5' arm
 5'-[rc nt 28-40]-N-[rc(barcode2)]-3'    5'-[rc nt 14-26]-OH3' ^ 5'P-[rc nt 1-13]-3'
    |||||||||||||    (to barcode2)           |||||||||||||   junction  |||||||||||||
 3'-[nt 40 .. 28]-[27]------------------- [nt 26 ....... 14]---------[nt 13 ..... 1]-5'  RNA

(c) Old order in tests/data/double_hyb_rca.yaml ("{part1}N{part2}") -- do not use:
 5'-[rc(barcode2) 17]-N-[rc(nt 28-40) 13]-3'      3' end on the RNA, circle part at the 5' end
```

- part1 (the circle's 5' arm) binds the 5' segment nt 1-13; part3 (3' arm) binds nt 14-26. The 5'
  end sits on nt 13 and the 3'-OH on nt 14: **ligation junction between nt 13 and 14**.
- barcode2 lies directly 5' of part3, so it leaves the RNA at nt 26, next to the primer's RNA part
  (nt 28-40). The primer's 3' part pairs with it; Phi29 extends the **primer's 3' end around the
  circle**.
- Why the old order fails: a polymerase extends only a 3'-OH. With rc(barcode2) at the 5' end, the
  part paired with the circle cannot be extended, and the 3' end is paired with RNA nt 28, where any
  extension would copy RNA (and run into the circle's arm), never the circle. MiP-seq states the
  rule: "The 5' end of the RCA initiator primer is complementary to the 3' end of the targeting
  sequence (12-16 bp), and the 3' end is complementary to the padlock probe" (Wu 2023), and its
  published primers (e.g. Th-1-RCA-2 GCCTCTGGGCAGG-C-AGGCTCCACTA) and U-Probe's
  `tests/data/RNA_format.yaml` put the RNA part at the 5' end. The template uses that order.

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Window | 40, step 20 | nt | U-Probe test protocol (overlap 20); U-Probe preprint, doi:10.64898/2026.04.12.717982 |
| Circle arms | 13 + 13, adjacent | nt | test protocol; same as MiP-seq (Wu 2024, doi:10.1038/s41551-024-01205-7) |
| Barcodes | 17 + N + 17 (placeholders) | nt | test protocol barcode lengths; sequences (unverified) |
| Primer | 13 (RNA) + N + 17 (rc barcode2) = 31 | nt | test protocol; orientation per Wu 2023, doi:10.1101/2023.01.07.523058 |
| Unbound nt between sites | 1 (nt 27) | nt | test protocol (primer site = window[-13:]) |
| Arm / primer-site Tm | 35-45 under 50 mM Na+, 0.6 mM dNTP, 50 nM, 0% formamide | C | test protocol (legacy U-Probe Tm conditions); criterion (unverified) |
| Same window, MiP-seq conditions | about +1.5 to +2.3 C vs legacy, i.e. ~37-47 | C | primer3 estimate (2x SSC, 10% FA, 30 nM) |
| Off-target genes | <= 1 on the 40-nt window (`n_mapped_genes`, Bowtie2) and <= 1 on the 26-nt arm site (`arm_mapped_genes`, Bowtie 1, 2 substitutions); both include the target | genes | template (U-Probe choice; the test protocol allowed <= 4) |
| Circle fold filter | none; `circle_fold_score` only sorts (ascending) | kcal/mol | template (U-Probe choice; the test protocol's `<= 80` removed nothing, section 7) |
| Hybridization | 2x SSC, 10% formamide, 37 C | - | borrowed from MiP-seq (Wu 2023); (unverified) for this design |
| Ligation / RCA | SplintR 37 C; Phi29 30 C | - | borrowed from MiP-seq; (unverified) for this design |
| Region | CDS of MANE/canonical transcript | - | U-Probe choice; no published source |

## 4. Target region and template

- Template molecule: the RNA itself (both oligos are rc(RNA)); ligation on an RNA splint needs
  SplintR or another RNA-templated ligase (T4 DNA ligase is inefficient there, Lohman 2014).
- Region: `source: CDS` (U-Probe choice; the test protocol used `exon`, which in U-Probe means CDS).
  Switch to `source: mRNA` to include UTRs, as the `mipseq` template does.
- Isoforms: `transcripts: mane` designs on one reference transcript; windows may span its exon
  junctions and then bind only spliced RNA.

## 5. Specificity

- One ligated circle plus a primer gives one RCP, which is counted as a molecule: each unit must be
  free of expressed off-targets. The template filters `n_mapped_genes <= 1` and `arm_mapped_genes <= 1`
  (the test protocol's `<= 4` left up to 3 off-target genes per unit).
- An off-target needs the circle's 26-nt arm site ligated on one molecule; the primer must then
  reach the circle. In MiP-seq the 11-nt tail alone (Tm ~30 C under the hybridization conditions,
  primer3 estimate) needs the primer's RNA part, so the AND gate covers the whole 40-nt window. Here
  the 17-nt barcode2 duplex (Tm ~44-50 C for typical 17-mers, same estimate) may hold the primer on
  any ligated circle by itself, so specificity rests mostly on the **26-nt arm site and the junction**.
  The template therefore searches the arm site separately: a `layout` splits `target_region` into
  `arm_site` (26 nt), `gap` (1) and `primer_site` (13), and `arm_mapped_genes` searches `arm_site`
  with Bowtie 1 (exhaustive within 2 substitutions; no helper oligo, so nothing extra is ordered).
- Bowtie2 on 40-nt windows flagged 97% of random windows with a site within 2 substitutions in
  another gene and listed 87% of those genes; its local alignment also reports partial hits of about
  25 nt or more, which covers an arm-site-only match. For an exhaustive search of the whole window
  set `n_mapped_genes` to `aligner: bowtie`, `max_mismatches: 2` (needs `bowtie` in `align_index`,
  which `arm_mapped_genes` needs anyway).
- Add `gene_expression` (unexpressed off-targets not counted) and `exempt_copies` if needed.

## 6. U-Probe template

```bash
uprobe templates show double_hyb_rca -o protocol.yaml
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets
uprobe run -p protocol.yaml -g genomes.yaml -o results --threads 8
```

Key lines (uprobe/templates/double_hyb_rca.yaml):

```yaml
encoding:
  GAPDH: {barcode1: NNNNNNNNNNNNNNNNN, barcode2: NNNNNNNNNNNNNNNNN}   # 17-nt placeholders
extracts:
  target_region:
    {source: CDS, transcripts: mane, length: 40, step: 20}
    layout:                                 # only for the arm-site off-target search
      template: "{arm_site}{gap}{primer_site}"
      parts: {arm_site: {length: 26}, gap: {length: 1}, primer_site: {length: 13}}
probes:
  circle_probe:
    template: "{part1}{part2}{part3}"
    parts:
      part1: {length: 13, expr: "rc(target_region[0:13])"}
      part2: {template: "{barcode1}N{barcode2}", parts: {...}}
      part3: {length: 13, expr: "rc(target_region[13:26])"}
  amp_probe:
    template: "{part2}N{part1}"            # was "{part1}N{part2}" in tests/data
    parts:
      part1: {expr: "rc(circle_probe.part2.barcode2)"}
      part2: {expr: "rc(target_region[-13:])"}
hybridization: {temperature_c: 37, na_mM: 390, formamide_pct: 10}   # from MiP-seq
attributes:
  tm1: {target: "circle_probe:part1", type: annealing_temperature,
        hybridization: {na_mM: 50, mg_mM: 0, dntp_mM: 0.6, probe_nM: 50, formamide_pct: 0}}
  n_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie2, strand: sense}
  arm_mapped_genes: {target: target_parts.arm_site, type: mapped_genes,
                     aligner: bowtie, max_mismatches: 2, strand: sense}
post_process:
  filters:
    tm1: {condition: tm1 >= 35 & tm1 <= 45}          # tm2 (part3), tm3 (amp part2) alike
    n_mapped_genes: {condition: n_mapped_genes <= 1}
    arm_mapped_genes: {condition: arm_mapped_genes <= 1}
  sorts: {is_ascending: [circle_fold_score, amp_fold_score, tm1]}   # least structured first
  remove_overlap: {location_interval: 0}
```

What to edit:
1. `encoding`: real barcode1/barcode2 per gene (17 nt each, or your readout lengths).
2. Replace the literal `N` spacers in `part2` and `amp_probe` with defined bases (a synthesized N is
   a degenerate mix).
3. Optional `gene_expression`; `aligner: bowtie`, `max_mismatches: 2` on `n_mapped_genes` too for an
   exhaustive search of the whole window.
4. `step: 1` for more candidates (the `mipseq` template uses 1); keep `remove_overlap`.
5. Decide sites per gene (`equal_space: {number_desired: N}`) and region (`source: mRNA` if UTRs
   are allowed).
6. Once the barcodes are final, add a `circle_fold_score` threshold if you want one (section 7).

**U-Probe departs** (template vs its own test protocol): primer order reversed (the fix above);
`step: 20` instead of `overlap: 20` (same windows); `source: CDS` with `transcripts: mane` instead
of `exon` on all transcripts; Bowtie2 with `strand: sense` (the test used
`similarity_threshold_ratio: 0.75`) plus the Bowtie 1 search of the arm site; `n_mapped_genes <= 1`
instead of `<= 4`; no fold filter; top-level hybridization conditions borrowed from MiP-seq while
tm1-tm3 keep the legacy conditions under which the 35-45 C window was set.
**U-Probe departs** (template vs MiP-seq): 17-nt gene-specific primer tail on barcode2 instead of
the constant 11-nt TAGTGGAGCCT site; a fixed 13-nt primer site; Tm and fold filters MiP-seq does not
state; `step: 20` (the `mipseq` template scans every position).

## 7. Pitfalls

- `fold_score` in U-Probe is **-MFE** (ViennaRNA, RNA parameters, kcal/mol, at `temperature_c`), so
  it is >= 0 and higher means more structure, and the sorts list it ascending: least structured first.
  The old `circle_fold_score <= 80` removed nothing: on 4,375 circle probes (875 human CDS windows x 5
  random 17-nt barcode pairs) the score was 14.3 (median), 21.0 (p90) and 37.2 at most. The template
  now filters nothing on it, because a threshold depends on the barcodes you choose; `<= 21` would
  keep 90% of those probes. Sorting runs last, so it changes only the output order, not which units
  are kept.
- Fold scores replace non-ACGU characters with A: with placeholder barcodes (all N) or literal N
  spacers, `circle_fold_score` and `amp_fold_score` describe poly-A, not your probe.
- `fold_score` (unsuffixed) targets `circle_probe:part1`, a 13-nt arm; it is computed but unused.
- Tm windows are on the legacy conditions; do not compare them with Tm computed under the
  top-level `hybridization` block (they differ by ~2 C here, more for other buffers).
- The circle needs a 5'-phosphate (order it or use T4 PNK); U-Probe does not annotate it.
- Do not loosen `n_mapped_genes` or `arm_mapped_genes` beyond 1: one unit makes one spot, so every
  off-target unit is a false molecule (the test protocol's `<= 4` kept units with 3 off-target genes).
- A gene-specific primer tail means one primer sequence per gene and barcode; barcode2 must also be
  unique enough that one gene's primer does not prime another gene's circle.

## 8. Verification checklist

- [ ] amp_probe = rc(RNA nt 28-40) + spacer + rc(barcode2), 5'->3' (RNA part at the 5' end).
- [ ] rc(part1) + rc(part3) = RNA nt 1-26, contiguous; junction 13|14; circle 5'-phosphorylated.
- [ ] barcode2 in the circle is immediately 5' of part3; primer tail = its exact reverse complement.
- [ ] No literal N left; barcodes are the intended readout sequences, not placeholders.
- [ ] `n_mapped_genes` and `arm_mapped_genes` = 1 per unit with the intended aligner/expression table.
- [ ] Fold/Tm filters recomputed after barcodes are final; arm Tm near the hybridization temperature.
- [ ] Ligase (RNA-templated), hybridization and RCA conditions stated, and marked as borrowed.

## 9. References

- Zhang Q, Cai H, Zhang J, et al. 2026, openRxiv (bioRxiv). Universal agentic probe design for
  imaging-based spatial-omics (U-Probe). doi:10.64898/2026.04.12.717982
- Wu X, Xu W, Deng L, et al. 2024, Nature Biomedical Engineering. Spatial multi-omics at subcellular
  resolution via high-throughput in situ pairwise sequencing. doi:10.1038/s41551-024-01205-7
- Wu X, Xu W, Deng L, et al. 2023, bioRxiv. Development of Multiomics in situ Pairwise Sequencing
  (MiP-Seq) for Single-cell Resolution Multidimensional Spatial Omics. doi:10.1101/2023.01.07.523058
- Wang X, Allen WE, Wright MA, et al. 2018, Science. Three-dimensional intact-tissue sequencing of
  single-cell transcriptional states (STARmap, SNAIL). doi:10.1126/science.aat5691
- Frei AP, Bava F-A, Zunder ER, et al. 2016, Nature Methods. Highly multiplexed simultaneous
  detection of RNAs and proteins in single cells (PLAYR). doi:10.1038/nmeth.3742
- Lohman GJS, Zhang Y, Zhelkovsky AM, Cantor EJ, et al. 2014, Nucleic Acids Research. Efficient DNA
  ligation in DNA-RNA hybrid helices by Chlorella virus DNA ligase. doi:10.1093/nar/gkt1032
