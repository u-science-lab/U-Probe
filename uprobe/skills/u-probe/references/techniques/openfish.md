# OpenFISH split probe pairs + reusable padlock

> Read this when: designing OpenFISH (OpenFISH-MALDI) gene probes, or any split-probe pair that
> brings a generic padlock to the RNA. Same pair geometry as `hcr-v3.md`; padlock/RCA background
> in `padlock-iss-hybiss.md` and `starmap.md`. Sources: the bioRxiv preprint (Methods,
> Supplementary Note, Tables S1-S3) and the authors' code; the Neuron paper was not accessible.

## 1. What it is and when to choose it

Two unmodified 46-nt gene probes (P1 = probe A, P2 = probe B) bind adjacent 30-nt RNA sites.
Their 14-nt overhangs are two halves of a 30-nt bridge, the docking site of a gene-independent
72-nt padlock with two 20-nt readout (RO) sites. T4 DNA ligase closes the padlock on the
overhangs and phi29 RCA starts at P2's 3' end: no signal without P1, P2 and the padlock (preprint
Supp. Note, Supp. Fig. 1a). Readout: 5 dyes; each gene = two different colours in one round, 10
genes per round (109 genes + 1 false-positive control in 11 rounds, 20x objective). Choose it for
panels of tens to ~100 genes in tissue on a standard epifluorescence microscope, when only the
gene probes should change between panels. One pair + padlock = one spot.

## 2. Probe architecture (5'->3')

```text
mRNA    5'-[ nt 1 ......... 30 ][31-32][ nt 33 ......... 62 ]-3'
probe A 5'-[bridge nt 1-14][TA][rc(nt 1-30)]-3'            46 nt, unmodified
probe B 5'-[rc(nt 33-62)][TA][bridge nt 17-30]-3'          46 nt, unmodified
padlock 5'-P-[rc(bridge 1-7)][RO x, 20][2 nt][RO y, 20][rc(bridge 8-30), 23]-3'   72 nt
```

- A binds the upstream site with its overhang at its 5' end, B the downstream site with its
  overhang at its 3' end; both point into the 2-nt gap (utils.py:12-86). All 225 published pairs
  found on one Ensembl GRCm39 transcript (of 231) have A upstream and a 2-nt gap, none reversed
  (computed here from Table S1).
- The padlock's 3' arm pairs with bridge nt 8-14 (on A) and 17-30 (on B), its 5' end with nt 1-7
  (on A): the nick lies mid-way in A's 14-nt overhang, where ligation "peaked" (Supp. Fig. 1b).
  Bridge nt 15-16 are on neither probe. The RCA product is rc(padlock), so each RO probe has the
  sequence of its padlock RO site (Table S1).
- Padlocks are 5'-phosphorylated; a "placeholder" padlock without phosphate dilutes abundant
  genes (Supp. Fig. 4; ratios in Table S3). Library: 110 padlocks, 60 ROs of A/C/T only (Table S1).

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Arm / gap / overhang / spacer | 30 / 2 / 14 / TA (written `ta`) | nt | utils.py:28-35,80-86; Methods; Table S1 |
| Arm GC | 30-70 (Para default); notebook 20-80; published arms 33-70 | % | generator.py:48; 05_Probe_design cell 10 |
| Arm sequence | no N, no GGGGG/CCCCC, Shannon entropy >= 1 | - | filter.py:464-491 |
| Arm Tm floor | >= 47 (sCelsius) under the author formula | C | filter.py:494-496; nb 05 cell 10 |
| Author Tm | Biopython Tm_NN, RNA/DNA table (Sugimoto 1995), 20 nM, 390 mM Na+, 0 Mg; + (0.453 fGC - 2.88) x 6.67 M (30 % formamide) | C | filter.py:445-456; generator.py:76 |
| Published arms, author formula | 45.3-65.6 (median 50.0); 98 of 462 below 47 | C | computed here, Table S1 |
| Published arms, U-Probe Tm | 42.2-64.2 (median 51.2) at 390 mM Na+, 30 % FA, 20 nM | C | computed here |
| Off-target | blastn-short per arm vs RefSeq RNA (taxid); drop the pair if a non-target hit has Tm > 37 C | - | cmdtools.py:4-22; filter.py:544-587 |
| Pairs per gene | 2 (98 of 109 genes; 1-4 after a pilot) | pairs | Methods "Panel Design"; Table S3 |
| Spacing | none; 10 of 96 neighbouring published pairs overlap | nt | computed here |
| Gene probe hyb. | 2x SSC, 30 % formamide, 1 % Tween 20, 0.02 % BSA, yeast tRNA; 20 nM/oligo; 42 C, 90 min; wash 30 % FA 2 x 15 min 42 C | - | Methods "OpenFISH protocol" |
| Padlock, ligation, RCA | 20 nM, no formamide, 37 C 60 min; T4 DNA ligase 7 U/uL RT 60 min; phi29 0.5 U/uL, 250 uM dNTP, 30 C 3.5 h | - | same |
| Readout | 2 nM RO + 20 nM blockers, 2x SSC 20 % FA, RT 15 min; strip 80 % FA 37 C | - | Methods "OpenFISH Imaging" |

Versions disagree: BLAST e-value 100/5 (notebook), 100/50 (Para), 10 (Methods text); the Methods
also say `selfcomp=True` and `fmd=30`, the code `False` and 6.67 M. Spacer: the Methods say "AT",
the code and every Table S1 probe `ta`.

## 4. Target region

- Authors: NCBI RefSeq records (NM_/XM_/NR_) of the gene, whole record by default
  (`gtype='gene'`); windows shared by several records are collapsed and favoured by the score
  (generator_multi.py:86-129; utils.py:26-56,184-196).
- 209 of the 225 located published pairs lie on the Ensembl canonical transcript (computed here).
- Gjb2 (low abundance): 4 pairs with 30-nt arms beat 10 pairs with 20-nt arms (Supp. Fig. 2).

## 5. Specificity

- The authors screen each 30-nt arm alone (a non-target hit with Tm > 37 C under their formula
  fails the pair), then BLAST full probes against a blacklist of failed probes + mm39 tRNAs
  (>= 16 matched nt fails) and against the padlocks in use (only A may span its own padlock's
  ligation junction) (generator.py:155-253; filter.py:624-682). Expression only annotates.
- U-Probe: per-arm `mapped_genes` with Bowtie 1, <= 2 substitutions, `<= 1` gene each (rule 4:
  30-nt arms need an exhaustive search). Not the authors' Tm rule; a `gene_expression` table
  relaxes it to expressed genes (U-Probe departs).

## 6. U-Probe: start from `hcr_v3`

Same geometry as HCR v3 (overhang at A's 5' end and B's 3' end, 2-nt gap). Replace from `name:` on:

```yaml
genome: mm39                     # align_index must include bowtie
targets: [Gad2]
encoding:
  Gad2: {bridge: TTTCGAGAGCCTTTAAGATAATGTGCCACA}   # bridge of the padlock assigned to Gad2 (PadLock_001)
extracts:
  target_region:
    source: mRNA
    transcripts: canonical       # mane in human
    length: 62
    step: 1
    layout:
      template: "{left}{gap}{right}"
      parts: {left: {length: 30}, gap: {length: 2}, right: {length: 30}}
probes:
  probe_a:
    template: "{overhang}TA{binding}"
    parts:
      overhang: {expr: "encoding[target]['bridge'][0:14]"}
      binding: {expr: "rc(target_parts['left'])"}
  probe_b:
    template: "{binding}TA{overhang}"
    parts:
      binding: {expr: "rc(target_parts['right'])"}
      overhang: {expr: "encoding[target]['bridge'][-14:]"}
hybridization: {temperature_c: 42, na_mM: 390, probe_nM: 20, formamide_pct: 30, mg_mM: 0, dntp_mM: 0}
attributes:
  a_gc: {target: probe_a.binding, type: gc_content}
  b_gc: {target: probe_b.binding, type: gc_content}
  a_tm: {target: probe_a.binding, type: annealing_temperature}
  b_tm: {target: probe_b.binding, type: annealing_temperature}
  a_mapped_genes: {target: target_parts.left, type: mapped_genes, aligner: bowtie, max_mismatches: 2, strand: sense}
  b_mapped_genes: {target: target_parts.right, type: mapped_genes, aligner: bowtie, max_mismatches: 2, strand: sense}
post_process:
  filters:
    a_gc: {condition: a_gc >= 0.3 & a_gc <= 0.7}
    b_gc: {condition: b_gc >= 0.3 & b_gc <= 0.7}
    a_tm: {condition: a_tm >= 42}
    b_tm: {condition: b_tm >= 42}
    a_mapped_genes: {condition: a_mapped_genes <= 1}
    b_mapped_genes: {condition: b_mapped_genes <= 1}
    no_g5_c5: {type: sequence_pattern, target: target_region, exclude_patterns: [GGGGG, CCCCC]}
  sorts: {is_ascending: [a_mapped_genes, b_mapped_genes]}
  remove_overlap: {location_interval: 2}
  equal_space: {number_desired: 2}
summary: {report_name: rna_report, attributes: [a_gc, b_gc, a_tm, b_tm, a_mapped_genes, b_mapped_genes]}
```

- One bridge per gene: it selects the padlock, hence the gene's round and colour pair. Take it
  from the authors' `SEQUENCES/Bridge_sequences.csv` (ids = padlock ids) and their codebook.
- U-Probe departs: Tm floor 42 C (the hybridization temperature) with primer3, not 47 C under the
  authors' formula (the two differ per arm by -4.4 to +4.9 C, 5th-95th percentile: do not port 47);
  it keeps all 231 published pairs. GGGGG/CCCCC is checked on the window (gap included). Pairs
  never overlap (`location_interval: 2`, not a published rule).
- `equal_space` with 2 keeps the 5'-most and 3'-most passing pairs; the authors sort by a score
  and pick by hand. To choose by quality, request more and pick rows (never edit sequences).

U-Probe cannot do: RO and bridge library design, padlock assembly and its NUPACK checks, panel
selection and round/colour assignment (scGIST, MAGIC, genetic algorithm), placeholder ratios,
decoding; nor the authors' Tm formula, entropy filter, BLAST Tm > 37 C rule, blacklist, padlock
cross-check and NUPACK binding-efficiency score.

## 7. Pitfalls

- Keep bridge nt 1-14 at A's 5' end and nt 17-30 at B's 3' end (`[0:14]`, `[-14:]`): the nick must
  fall on A's overhang and RCA starts at B's 3' end; swapped halves do not match the padlock.
- Two genes sharing a bridge share a code. Padlocks need a 5'-phosphate (placeholders lack it).
- With 2 pairs requested, a gene with one passing pair keeps one: report it (Table S3 has 7).
- Do not shorten arms to 20 nt for more pairs (Supp. Fig. 2).

## 8. Verification checklist

- [ ] A = bridge[0:14] + TA + rc(window 1-30); B = rc(window 33-62) + TA + bridge[-14:]; 46 nt.
- [ ] Windows 62 nt on the reference mRNA; every pair of a gene carries that gene's bridge.
- [ ] Each arm: GC 30-70 %, Tm >= 42 C at 390 mM Na+/30 % FA/20 nM, `*_mapped_genes` = 1.
- [ ] `check_probes.py ... --gc 0.3 0.7 --tm 42 100 --min-gap 2 --min-probes 2` passes.
- [ ] 2 pairs per gene, or the shortfall named.

## 9. References

- Li X et al. 2025. OpenFISH enables integrated high-resolution spatial transcriptomics and
  metabolomics on a single tissue section. bioRxiv. doi:10.1101/2025.08.19.671030 (v1; full text,
  Supplementary Note, Tables S1-S3 read; source of every wet-lab and table number above)
- Li X et al. 2026. OpenFISH enables same-section spatial transcriptomics and MALDI-MSI
  integration. Neuron. doi:10.1016/j.neuron.2026.09.007 (closed; abstract only)
- Lixinyoung/OpenFISH_MALDI_Protocol, GitHub, commit be490eb (openDesign/*.py, notebooks 01-05);
  Lixinyoung/openfish-maldi, commit 01ea661 (identical openDesign code)
