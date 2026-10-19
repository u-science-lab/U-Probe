# HCR v3.0 split-initiator probe pairs

Read this when a user asks for HCR / HCR v3 / HCR RNA-FISH probes, needs enzyme-free signal
amplification (embryos, thick or autofluorescent tissue, flow cytometry), or must assign B1-B5
amplifiers. Covers pair geometry, initiator halves and spacers, pair counts and buffers.

## 1. Method and when to choose it

- Each probe pair binds a 52-nt RNA window: two 25-nt sites with a 2-nt gap. Each probe carries half
  of a 36-nt HCR initiator; only when both bind side by side is the full initiator formed, which
  triggers polymerization of two dye-labelled 72-nt hairpins (H1, H2; 60 nM each).
- Split initiators give "automatic background suppression": a probe bound alone cannot amplify.
  Choi 2018 measured ~50-fold suppression in situ (~60-fold in gels).
- Choose it for: whole-mount embryos, thick or autofluorescent tissue, low-abundance RNA, 1-5 targets
  at once (one amplifier each; 10-plex with spectral imaging, Schulte 2024), quantitative relative
  expression (qHCR, overnight amplification) or single-molecule dots (dHCR, 45-90 min amplification).
- Do not choose it for >~10 genes in the same sample (see `merfish.md`, `seqfish-plus.md`).

## 2. Probe architecture (5' -> 3')

```
odd probe  (45 nt): 5'-[initiator half 1, 18][spacer, 2][rc(window nt 1-25), 25]-3'
even probe (45 nt): 5'-[rc(window nt 28-52), 25][spacer, 2][initiator half 2, 18]-3'

bound on the RNA (initiator halves meet over the 2-nt gap):
mRNA      5'-[ nt 1 .............. 25 ][26-27][ nt 28 ............. 52 ]-3'
odd        3'-[ binding 25            ]-sp-+  +-sp-[ binding 25       ]-5'   even
                            init half 1 (18)    init half 2 (18)
```

- The odd probe binds the 5' (upstream) segment of the window and carries the FIRST 18 nt of the
  initiator at its 5' end; the even probe binds the 3' segment and carries the LAST 18 nt at its
  3' end (Choi 2018 SI S4 probe-table headers; confirmed on all published human pairs).
- Pairs are numbered along the mRNA (1/2, 3/4, ...). Amplifier halves and spacers (Choi 2014 SI S6
  initiators, split 18/18; spacers as in every Choi 2018 probe set):

| Amplifier | Half 1 (odd, 5' end) | Odd spacer | Even spacer | Half 2 (even, 3' end) |
| --- | --- | --- | --- | --- |
| B1 | GAGGAGGGCAGCAAACGG | AA | TA | GAAGAGTCTTCCTTTACG |
| B2 | CCTCGTAAATCCTCATCA | AA | AA | ATCATCCAGTAAACCGCC |
| B3 | GTCCCTGCCTCTATATCT | TT | TT | CCACTCAACTTTAACCCG |
| B4 | CCTCAACCTACCTCCAAC | AA | AT | TCTCACCATATTCGCTTC |
| B5 | CTCACTCCCAATCTCTAT | AA | AA | CTACCCTACAAATCCAAT |

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Binding length | 25 nt per probe | | Choi 2018, doi:10.1242/dev.165753 |
| Gap within a pair | 2 nt (pair footprint 52 nt) | Table S1 model target; all published pairs | Choi 2018 SI |
| Initiator | 36 nt, split 18 + 18 | B1-B5 orthogonal | Choi 2014, doi:10.1021/nn405717p |
| Spacer | 2 nt, amplifier-specific (table above) | | Choi 2018 SI S4 |
| GC / Tm / homopolymer | none stated | published sets: GC 8-96 % (median 52 %), runs up to 13 nt | Choi 2018 SI (parsed) |
| Spacing between pairs | 2 nt in the published sets (54-nt grid) | not stated as a rule | Choi 2018 SI (parsed) |
| Pairs per target, qHCR | 20 ("we recommend using 20 split-initiator probe pairs") | Fig. 3 tested 5, 10, 20 | Choi 2018 Discussion |
| Pairs per target, dHCR | >= 25 preferred (Choi 2018); >= 30 (Molecular Technologies FAQ) | maximize for single molecules | Choi 2018; moleculartechnologies.org/info/faq |
| Probe hyb. buffer | 30 % formamide, 5x SSC, 9 mM citric acid pH 6.0, 0.1 % Tween 20, 50 ug/mL heparin, 1x Denhardt's, 10 % dextran sulfate | | Choi 2018 SI S2.1.4 |
| Probe concentration | 4 nM each probe (qHCR); 10 nM (dHCR, Choi 2018) vs 16 nM (MI RNA-FISH Rev11) | e.g. 1.2 pmol / 300 uL, 0.4 pmol / 100 uL, 4 pmol / 1 mL | Choi 2018; MI protocols |
| Hybridization | 37 C, overnight (12-16 h); prehyb 30 min 37 C | | Choi 2018 |
| Probe wash | buffer above minus Denhardt's and dextran; 37 C, 4 x 5 min (cells) / 4 x 15 min (embryos); then 5x SSCT RT | FFPE (MI): graded 75/50/25/0 % wash buffer in 5x SSCT, 15 min each | Choi 2018; MI FFPE Rev1 |
| Amplification | 5x SSC, 0.1 % Tween 20, 10 % dextran sulfate; hairpins 60 nM each, snap-cooled (95 C 90 s, 30 min RT, separate tubes) | RT, overnight (qHCR) or 45-90 min (dHCR) | Choi 2018; MI |

Signal versus probe-set size (Choi 2018 Tables S10-S11, chick Sox10): v3 background stays flat
(29, 26, 28 for 5, 10, 20 pairs) while signal-to-background rises 28 -> 42 -> 65; v2 standard probes
(full initiator on each probe) lose S/B as probes are added (7.8 -> 2.0).

Community tools disagree on filters Choi 2018 never set; none is canonical:

| Tool | GC | Spacing | Pairs | Specificity |
| --- | --- | --- | --- | --- |
| insitu_probe_generator (Kuehn 2022) | none | >= 3 nt in code ("skips >= 2 bp" in paper) | all, or cap 33 | optional blastn-short of the 52-mer |
| HCRProbeMakerCL | optional, default 0-100 % | >= 3 nt | cap 36 | BLAST flags only |
| gofflab HCRProbeDesign | 45-55 % | 0 nt | 20 | Bowtie2 genome, 1 hit |
| JeffLeeLab HCRv3_probe_design | 40-60 %, Tm 60-90 C | 3 nt | max 33 | blastn, 1 gene |
| Easy_HCR (Elagoz 2022) | ends not GC/CG | as insitu_probe_generator | >= 20, cap 33 | |

## 4. Target region

- Whole mRNA, UTRs included: Choi 2018 designed on NCBI mRNAs; 9 of the 10 human ACTB pairs lie in
  the 5' or 3' UTR (likely to avoid the other actin isoforms), GAPDH/PGK1/BRAF sets reach the 3' UTR.
- Use the reference transcript (`transcripts: mane`); list transcript ids for an isoform-specific set.
- Short RNAs: a 20-pair set needs >= 20 x 54 = ~1,080 nt of usable sequence; 5-10 pairs still work
  but with lower S/B.
- Introns / nascent RNA: same pair design on `source: genome` intervals (not covered by the template).

## 5. Specificity

- Unit of specificity = the pair. One probe alone does not amplify, so a 25-nt off-target match of
  one probe matters little. A pair that binds an off-target RNA side by side DOES amplify, and every
  amplified pair is a bright dot in dHCR, so each 52-nt window must be unique among expressed RNAs.
  This is stricter than smFISH, where one off-target probe is invisible.
- Paralogs and pseudogenes: put pairs in UTRs or divergent regions; check `target_mapped_genes_off_targets`.
- Cross-talk between targets: only through amplifier choice; B1-B5 are orthogonal.

## 6. U-Probe template: `hcr_v3`

```bash
uprobe templates show hcr_v3 -o protocol.yaml
```

```yaml
encoding:
  GAPDH: {amplifier: B1, init_odd: GAGGAGGGCAGCAAACGG, spacer_odd: AA, spacer_even: TA, init_even: GAAGAGTCTTCCTTTACG}
extracts:
  target_region:
    source: mRNA
    transcripts: mane
    length: 52
    layout: {template: "{left}{gap}{right}", parts: {left: {length: 25}, gap: {length: 2}, right: {length: 25}}}
probes:
  probe_odd:  {template: "{initiator}{spacer}{binding}", ...}   # binding = rc(target_parts['left'])
  probe_even: {template: "{binding}{spacer}{initiator}", ...}   # binding = rc(target_parts['right'])
hybridization: {temperature_c: 37, na_mM: 975, probe_nM: 4, formamide_pct: 30}   # 5x SSC
attributes:
  target_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie2, strand: sense}
post_process:
  filters: {target_mapped_genes: {condition: target_mapped_genes <= 1}}
  remove_overlap: {location_interval: 3}
  equal_space: {number_desired: 20}
```

What to edit:
- `targets` and `encoding`: one amplifier per target imaged together; copy halves and spacers from
  the table in section 2 (the template ships B1 and B2 only).
- `equal_space.number_desired`: 20 for qHCR; 25-30+ for dHCR (all passing pairs if fewer).
- Add `gene_expression` so only expressed off-target genes remove pairs.
- Optional exhaustive search: add `bowtie` to `align_index` and set `aligner: bowtie`, `max_mismatches: 2`.

Template departures:
- Spacing: `location_interval: 3` (insitu_probe_generator code); Choi 2018 sets sit 2 nt apart.
- Off-targets: Bowtie2 local on the 52-nt window (heuristic, <= 100 alignments; flagged 99 % of windows
  with a <= 2-substitution site in another gene). It also counts partial matches of ~25 nt, i.e. one
  half: conservative, it drops pairs that could not amplify on that gene.
- No GC/Tm/homopolymer filter (as Choi 2018). `odd_tm`/`even_tm` and GC are reported only.

## 7. Pitfalls

- Swapped halves: half 1 always on the odd (5'-window) probe's 5' end, half 2 on the even probe's
  3' end. Mixing halves of different amplifiers in one pair gives no signal.
- Spacers are amplifier-specific (B1 even = TA, B4 even = AT); do not reuse B1 spacers for B4.
- Amplifier names beyond B5 (B7, B9-B17 from Wang et al. 2018 bioRxiv; B6-B10 from Molecular
  Technologies) may not denote the same sequences; get sequences from the supplier's documentation.
- Molecular Instruments' current "HCR Gold" amplifiers (X1-X10) are proprietary; custom pairs need
  matching public initiators (B1-B5) and hairpins.
- Ordering: tool caps (33/36 pairs) follow IDT oPool price tiers, not biology.
- Repeats in 3' UTRs (Alu, poly(A)-rich stretches) give windows with many hits; they are filtered,
  which can leave short UTR-rich genes with few pairs.

## 8. Checklist

- [ ] Odd = 45 nt `half1 + spacer + rc(left 25)`; even = 45 nt `rc(right 25) + spacer + half2`.
- [ ] Windows 52 nt, pairs non-overlapping, gap 2 nt inside each pair.
- [ ] Each target has its own amplifier; halves and spacers match that amplifier's row.
- [ ] 20 pairs (qHCR) or >= 25-30 (dHCR) per target, or the shortfall reported.
- [ ] `target_mapped_genes` = 1 for every pair (or off-targets unexpressed / copies of the target).
- [ ] Buffer plan: 30 % formamide, 5x SSC, 4 nM per probe, 37 C overnight; hairpins 60 nM snap-cooled.

## 9. References

- Choi HMT et al. 2018. Third-generation in situ hybridization chain reaction: multiplexed, quantitative, sensitive, versatile, robust. Development 145:dev165753. doi:10.1242/dev.165753
- Choi HMT, Beck VA, Pierce NA. 2014. Next-generation in situ hybridization chain reaction: higher gain, lower cost, greater durability. ACS Nano 8:4284-4294. doi:10.1021/nn405717p
- Schwarzkopf M et al. 2021. Hybridization chain reaction enables a unified approach to multiplexed, quantitative, high-resolution immunohistochemistry and in situ hybridization. Development 148:dev199847. doi:10.1242/dev.199847
- Schulte SJ et al. 2024. HCR spectral imaging: 10-plex, quantitative, high-resolution RNA and protein imaging in highly autofluorescent samples. Development 151:dev202307. doi:10.1242/dev.202307
- Kuehn E et al. 2022. Segment number threshold determines juvenile onset of germline cluster expansion in Platynereis dumerilii (insitu_probe_generator). J Exp Zool B 338:225-240. doi:10.1002/jez.b.23100
- Elagoz AM et al. 2022. Optimization of whole mount RNA multiplexed in situ HCR ... octopus embryonic neurogenesis (Easy_HCR). Front Physiol 13:882413. doi:10.3389/fphys.2022.882413
- Molecular Instruments protocols (mammalian cells on slide Rev6/Rev11, zebrafish Rev7, FFPE Rev1): files.molecularinstruments.com
