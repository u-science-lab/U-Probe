# Hybridization, Tm and secondary structure

Read this when you set window lengths, a Tm or GC window, the `hybridization:` block, or
structure/homopolymer filters, or when a Tm filter removes everything. It explains what sets duplex
stability, how U-Probe computes Tm, and how to place a Tm window relative to the hybridization temperature.

## What sets duplex stability

| Factor | Effect on Tm | Source |
| --- | --- | --- |
| Length | longer -> higher, with diminishing gain per nt; one mismatch costs less in a longer duplex | NN model, SantaLucia 1998 |
| Sequence / GC | GC stacks are more stable; Tm depends on nearest-neighbour stacks, not only GC fraction | SantaLucia 1998 |
| Monovalent cations (Na+) | higher [Na+] -> higher Tm (log-linear); 2x SSC = 0.39 M Na+ | Owczarzy et al. 2004 |
| Mg2+ | strongly stabilizing; dNTPs chelate it | Owczarzy et al. 2008; von Ahsen et al. 2001 |
| Formamide | lowers Tm linearly: 2.4-2.9 °C per mol/L depending on GC, ~0.6-0.7 °C per % (v/v) | Blake & Delcourt 1996; Hutton 1977 |
| Strand concentration | higher probe concentration -> higher Tm (bimolecular) | SantaLucia 1998 |
| Duplex type | RNA:DNA hybrid stability depends on sequence (most stable with a pyrimidine-rich DNA strand, i.e. purine-rich RNA); in high formamide RNA:DNA hybrids are much more stable than DNA:DNA | Sugimoto et al. 1995; Lesnik & Freier 1995; Casey & Davidson 1977 |
| Temperature | binding is a two-state equilibrium; fraction bound falls steeply around Tm | SantaLucia 1998 |

- SantaLucia J Jr 1998 PNAS 95:1460, doi:10.1073/pnas.95.4.1460.
- Owczarzy R et al. 2004 Biochemistry 43:3537, doi:10.1021/bi034621r; 2008 Biochemistry 47:5336, doi:10.1021/bi702363u.
- von Ahsen N, Wittwer CT, Schütz E 2001 Clin Chem 47:1956, doi:10.1093/clinchem/47.11.1956
  ([Na+eq] = [monovalent] + 120 sqrt([Mg2+] − [dNTP]), mM).
- Blake RD, Delcourt SG 1996 Nucleic Acids Res 24:2095, doi:10.1093/nar/24.11.2095.
- Hutton JR 1977 Nucleic Acids Res 4:3537, doi:10.1093/nar/4.10.3537 (0.60 °C per % formamide, 0.035-0.88 M NaCl).
- Sugimoto N et al. 1995 Biochemistry 34:11211, doi:10.1021/bi00035a029.
- Lesnik EA, Freier SM 1995 Biochemistry 34:10807, doi:10.1021/bi00034a013.
- Casey J, Davidson N 1977 Nucleic Acids Res 4:1539, doi:10.1093/nar/4.5.1539 (80% formamide, 0.3 M
  NaCl: rRNA:DNA hybrid melted at 57 °C vs 38 °C for DNA:DNA).

## How U-Probe computes Tm (`annealing_temperature`)

- Primer3 `calc_tm` (Untergasser A et al. 2012 Nucleic Acids Res 40:e115, doi:10.1093/nar/gks596):
  SantaLucia 1998 nearest-neighbour **DNA/DNA** parameters, SantaLucia salt correction, Mg2+ converted
  to a Na+ equivalent (von Ahsen 2001), conditions from `hybridization:` (`na_mM`, `mg_mM`, `dntp_mM`,
  `probe_nM`, `formamide_pct`). Unset fields keep PCR-like defaults (50 mM Na+, 0.6 mM dNTP, 50 nM,
  no formamide): **always declare the real conditions**.
- RNA targets: the value is a DNA/DNA proxy for a DNA:RNA hybrid. Use it to rank windows and to
  reproduce windows that papers computed the same way; do not read it as the true hybrid Tm (at high
  formamide the true hybrid is likely more stable, Casey & Davidson 1977).
- Strand concentration: Primer3 assumes equal concentrations of both strands (the C_T/4 term of
  SantaLucia 1998). In FISH the probe is in large excess over its target, where the term is about the
  probe concentration; by the same equations this raises Tm by about 1.4-2 °C (computed with U-Probe's Primer3 call: +2.0 °C for 20-mers at 50 nM,
  +1.4 °C for 30-mers at 1 nM). Uniform offset; it rarely changes which windows pass.

### Formamide: Blake & Delcourt vs the 0.65 °C/% rule

| Model | Tm drop | Used by |
| --- | --- | --- |
| Blake & Delcourt 1996 (in Primer3): (2.88 − 0.453 × f_GC) °C per mol/L; 1% (v/v) = 0.2516 mol/L | 0.725 °C/% at 0% GC, 0.668 at 50%, 0.611 at 100% | U-Probe |
| Fixed linear rule, 0.65 °C/% | 0.65 °C/% regardless of GC | Biopython `chem_correction` default (its docstring: "factors between 0.6 and 0.72"), hence tools built on it such as OligoMiner's `blockParse` |
| Measured linear coefficient | 0.60 °C/% (Hutton 1977) | - |

At 50% formamide and 30-70% GC, U-Probe's Tm is 2 °C lower to 0.3 °C higher than under the 0.65 rule
(lower for AT-rich windows: 0.691 vs 0.65 °C/% at 30% GC). A window published under the 0.65 rule,
applied to U-Probe Tm, is off by up to ~2 °C; widen it by that much or note it. OligoMiner: Beliveau BJ et al. 2018 PNAS 115:E2183,
doi:10.1073/pnas.1714530115 (`blockParse` defaults 36-41 nt, Tm 42-47 °C, 390 mM Na+, 50% formamide).

Typical values (U-Probe/Primer3, median of random 45%-GC windows, 50 nM probe):

| Length | 50 mM Na+ (default) | 2x SSC | 2x SSC + 10% FA | + 30% FA | + 50% FA |
| --- | --- | --- | --- | --- | --- |
| 20 nt | 50 | 60 | 54 | 40 | 27 |
| 25 nt | 56 | 66 | 59 | 46 | 32 |
| 30 nt | 60 | 71 | 64 | 51 | 37 |
| 40 nt | 66 | 77 | 70 | 56 | 43 |

## Choosing the Tm window

1. Declare the hybridization conditions of the protocol (`temperature_c`, `na_mM`, `mg_mM`,
   `formamide_pct`, `probe_nM`). SSC: 1x = 150 mM NaCl + 15 mM trisodium citrate = 195 mM Na+.
2. If the method's paper gives a Tm window **with** the conditions it was computed under, reproduce it:
   put those conditions in a per-attribute `hybridization:` override and use the published window
   (MERFISH: 66-76 °C at 300 mM Na+, 5 nM probe, no formamide, 30-nt target regions, Moffitt et al.
   2016 PNAS 113:11046, doi:10.1073/pnas.1612826113; OligoMiner: 42-47 °C at 390 mM Na+ and 50%
   formamide, Beliveau et al. 2018).
3. Otherwise set the window **relative to the hybridization temperature under the real conditions**.
   Published protocols, recomputed with U-Probe under their own conditions (median of random windows
   in their GC range), sit at these margins above the hybridization temperature:

| Design (source) | Window | Conditions | T_hyb | Median Tm − T_hyb |
| --- | --- | --- | --- | --- |
| smFISH (Raj et al. 2008) | 20 nt | 2x SSC, 10% FA | 37 | +18 °C |
| MERFISH (Moffitt et al. 2016) | 30 nt | 2x SSC, 30% FA | 37 | +16 °C |
| HCR v3 half (Choi et al. 2018) | 25 nt | 5x SSC, 30% FA, 4 nM | 37 | +13 °C |
| HybISS arm pair (Gyllborg et al. 2020) | 2 x 15 nt | 75 mM monovalent, 10 mM Mg2+, 20% FA (design assumptions; ligation 37 then 45 °C) | 45 | +14 °C (pair); −9 °C (one arm) |
| seqFISH+ (Eng et al. 2019) | 28 nt | 2x SSC, 40% FA, 1 nM | 37 | +5 °C |
| Oligopaint (Beliveau et al. 2018) | 35-41 nt | 2x SSC, 50% FA | 42 | +3 °C |

   So: RNA oligo sets run ~5-20 °C above T_hyb under the real conditions; DNA Oligopaints run close to
   T_hyb at high formamide. Computed by this skill's authors from the templates' conditions; treat as
   orientation, not as a rule from one paper.
4. Keep the window **narrow** (MERFISH 10 °C, OligoMiner 5 °C) so all probes of a set bind similarly
   under one stringency; a variable `length: [min, max]` lets each site reach the window.
5. Stringency is the margin: the closer T_hyb is to Tm, the more a mismatch matters. U-Probe/Primer3
   thermodynamic alignment, 2x SSC: one central mismatch lowers a 20-mer's Tm by ~6 °C (median; 3-9),
   a 30-mer's by ~3 °C; two central mismatches ~13 °C (20-mer) and ~7 °C (30-mer); a central 1-nt
   bulge (an indel) ~9-11 °C (20-mer) and ~6-8 °C (30-mer); a terminal mismatch < 1 °C (estimates from
   Primer3's parameters, computed for this skill; a 1-nt bulge costs +4.0 kcal/mol at 1 M NaCl in
   SantaLucia & Hicks 2004 Annu Rev Biophys Biomol Struct 33:415, doi:10.1146/annurev.biophys.32.110601.141800).
   With a 15-20 °C margin, near-matches still bind: specificity must come from the off-target search
   (`references/concepts/specificity.md`), not from the Tm window.
6. Auto-generated Tm filters (no `post_process`) are [T_hyb, T_hyb + 40] °C: a sanity bound, not a design window.

Classic kinetics: for long DNA the renaturation rate is optimal at Tm − 25 ± 5 °C, with a flat
temperature response (Hutton 1977). Oligo probes in large excess bind fast; for them the equilibrium
margin and the mismatch penalties above decide what stays bound.

## GC content

GC fraction is a coarse proxy for stability; the Tm window does the real work. Use GC filters to
reproduce a protocol (MERFISH 43-63%, Moffitt et al. 2016; seqFISH+ 45-65%, Eng et al. 2019; OligoMiner
20-80%) or to drop extremes. Raj et al. 2008 aimed for GC "as close to 45% as possible" (Supplement).

## Secondary structure (`fold_score`, `self_match`)

- Target structure competes with probe binding; probe hairpins and self-dimers reduce the free probe
  and can make probes stick to each other. Formamide weakens all base pairing, structure included;
  many-probe sets tolerate a few poorly binding probes (the spot is the sum of many).
- `fold_score` on `target_region` folds the window alone: it sees local hairpins, not long-range
  structure of the whole RNA or proteins bound to it.
- `fold_score` = −MFE (kcal/mol, >= 0) of the sequence folded alone by ViennaRNA (Lorenz R et al. 2011
  Algorithms Mol Biol 6:26, doi:10.1186/1748-7188-6-26) with RNA parameters (Turner 2004; Mathews DH
  et al. 2004 PNAS 101:7287, doi:10.1073/pnas.0401799101) at `temperature_c`; formamide and salt are not
  modelled. **Higher = more structured, so lower is better.** Filter with `fold <= X` and sort it in
  `is_ascending` (least structured first), never `is_descending`.
- For a DNA probe, RNA parameters overstate hairpin stability; compare windows with each other rather
  than reading absolute values. None of the papers cited here sets a universal threshold. Heuristic
  (uncertain, not from a paper): drop only the most structured tail of a target's windows (e.g. the top
  5-10% of `fold_score`), and only when the target has candidates to spare.
- `self_match` counts 4-mer complementarities within the sequence; lower is better. Constant flanks
  (readouts, initiators, backbones) are fixed: check that they do not pair with the binding region of
  each probe (fold the whole probe: `fold_score` on `probe`, not only on the binding part).

## Homopolymers and low-complexity sequence

- Exclude runs of >= 5 identical bases with a `sequence_pattern` filter: OligoMiner's default prohibited
  sequences are `AAAAA,TTTTT,CCCCC,GGGGG` (Beliveau et al. 2018); seqFISH+ dropped "probe sequences that
  contained five or more consecutive bases of the same kind" (Eng CHL et al. 2019 Nature 568:235,
  doi:10.1038/s41586-019-1049-y, Methods).
- Mechanisms: G runs can form G-quadruplexes (Burge S et al. 2006 Nucleic Acids Res 34:5402,
  doi:10.1093/nar/gkl655); a T run in the probe (A run in `target_region`) can bind the poly(A) tail of
  every mRNA; homopolymers make register-shifted (slipped) partial duplexes possible. Filter
  `target_region` for A runs when the probe is `rc(target_region)`.
- Dinucleotide repeats and other low-complexity windows hit many loci; they are handled by the
  off-target search or `kmer_count` (`references/concepts/specificity.md`).
