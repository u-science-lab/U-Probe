# PRISM: colour-intensity barcoded padlocks, one imaging round

> Read this when: a user names PRISM (Chang 2025, "profiling of RNA in situ through single-round
> imaging"), wants up to 64 genes in one staining round on a widefield microscope, or asks for
> RNA-templated padlocks with an external RCA primer. Padlock basics and the direct-RNA variant:
> `padlock-iss-hybiss.md`; RNA-bound primer instead: `mip-seq.md`; rules:
> `references/concepts/specificity.md`.

## 1. What it is and when to choose it

A 5'-phosphorylated padlock binds a 40-nt mRNA window with two 20-nt arms, SplintR ligates it on
the RNA, an external primer starts Phi29 RCA, and all imaging probes are applied in one pooled
staining round. Each padlock carries four 20-nt segments, one per channel; each segment sequence
stands for one intensity level of its channel.

- Identity is read as intensities that the imaging-probe mix assigns to the segments, not by
  sequencing: every segment sequence has a dye-labelled and an unlabelled twin, premixed at the
  level's ratio; competitive binding sets the spot's intensity per channel, and the intensity
  ratios are decoded in colour space.
- Capacity: Texas Red + AF488 + Cy5 levels in quarters summing to 4 give 15 codes, times Cy3 0/1
  = 30, plus Cy3 alone = 31; in fifths summing to 5 and Cy3 0/1/2: 21 x 3 + 1 = 64 (v3 Fig. 1e and
  Methods; article Extended Data Fig. 7). The fixed sum keeps codes apart when RCP size scales all
  channels together ("radius vector filtering").

Choose it for <= 64 genes, conventional microscopes, no fluidics, thick tissue. Choose MERFISH,
seqFISH+ or ISS for more genes; smFISH when counting matters (PRISM sensitivity "about 10%", v3
Discussion).

## 2. Probe architecture

```text
Padlock, 122 nt (31-panel; 120 nt without TC in the 64-panel), 5'->3', ordered with "5' P":
5'P-[arm_5 20][TC][G 20][R 20][B 20][Y 20][arm_3 20]-OH 3'
arm_5 = rc(RNA nt 1-20); G = Cy3, R = Cy5, B = AF488, Y = Texas Red segment; arm_3 = rc(RNA 21-40)

Bound (backbone loops away):
     5'-[arm_3 .............]-OH 3' ^ 5' P-[arm_5 ..............]-3'
     3'-[RNA nt 40 ........ 21]        [nt 20 ................ 1]-5'    junction = RNA nt 20|21
RCA primer, 20 nt: rc(R segment), two 3'-terminal phosphorothioates; binds the circle, not the RNA
Imaging probe, 20 nt: = segment sequence (binds the RCP); 3' Cy5 (R), 3' TxRed (Y), 3' AF488 (B),
                      5' Cy3 (G; 64-panel 2G: both ends)
Example Slc17a7: GGACGTAAAGAAGCGCCTCC TC ACACGACGCTCTTCCGATCT CGCAGCAGATAAATGAGCCA
                 TCTCGCTAGATCTTGTCAGT GGTCTTAGCACCTCGTGATC ATGATGGCATAGACGGGCAT (barcode "0040")
```

- Only the arms depend on the target; the 80-82-nt backbone is the gene's barcode (ST2 column I).
- Codeword digits read Ch1-Ch4 = Y, B, R, G (TxRed, AF488, Cy5, Cy3): "0040" = R level 4. The
  physical order is G-R-B-Y. One primer per R sequence is pooled (5 in the 31-panel, 12 in the 64).

## 3. Design parameters

| Parameter | Value | Source |
| --- | --- | --- |
| Window, arms | 40 nt split 20 + 20; 214 of 215 published padlocks; HBX (viral) 14 + 14 | v3 Methods; ST2 (parsed here) |
| Orientation | arms = rc(RNA): 182 of 184 distinct windows match their gene's mRNA exactly (Ensembl GRCm39, GRCh38) | computed here from ST2 |
| Tm | 65-75 C on the 40-nt sequence; model and conditions not stated | v3 Methods |
| Sequence | no 5 consecutive G; no significant DNA secondary structure (dG > -10) | v3 Methods |
| Isoforms | cover all transcript variants if possible | v3 Methods |
| Specificity | "no blast results matching other genes in the NCBI database" | v3 Methods |
| GC | none stated; published windows 37.5-67.5 % | computed here |
| Padlocks per gene | text: 2 (mouse), 3 (HCC); ST2: 1 per mouse gene, 3 per HCC gene; "30-plex PRISM (one RCA probe per gene)" | v3 Methods, validation; ST2 |
| Hybridization | 200 nM each padlock, 50 mM KCl, 20 % formamide in Ampligase buffer (1x: 25 mM KCl, 10 mM MgCl2; so 75 mM K+), 55 C 15 min then 45 C 2 h | v3; Biosearch A1905B |
| Ligation, RCA | SplintR 2.5 U/uL, 37 C 2 h; Phi29 0.25 U/uL, 300 nM RCA primer, 30 C overnight | v3 Methods |
| Imaging probes | 31-panel: 17 labelled/unlabelled pairs, 120 nM per pair, 20 % formamide, 2x SSC, 50 C 9 min then 37 C 21 min; 64-panel: 39 sequences | v3 Methods; ST2 |

ST2 = Supplementary Table 2; v3 = bioRxiv v3 (2025-10-05), the open text closest to the article.
The text's "2 per mouse gene" contradicts ST2 and the validation section: follow ST2, use more for
sensitivity ("may be alleviated by increasing the number of probes per target gene", v3 Discussion).

## 4. Target region

mRNA, direct RNA (no RT); CDS and UTR not distinguished. No positional rule is published. Abundant
genes were diluted with "irreplicable" padlocks (1/2 to 1/36 functional, ST1) to limit crowding.

## 5. Specificity

One ligated padlock gives one spot: keep each unit strict, `target_mapped_genes <= 1` with Bowtie 1
on the whole 40-nt window (175 of 184 published windows pass on Ensembl; the others hit Gapdh
pseudogenes, Mylf-ps, Gm28309, FCGR3B (three FCGR3A sites) or unnamed genes sharing C1QA (two) or
EPCAM sequence; computed here). SplintR end-joining fidelity on RNA is limited (Krzywkowski 2017),
so specificity rests on hybridization.

The authors' probe-designer (PRISM-Code README.md:40; commit d7326be, publication-time version not
identified) filters differently: per-arm Tm 50-65 C by Biopython R_DNA_NN1 at 50 mM Na+
(configs/config_single_sequence.yaml:22-24; filtering/__init__.py:146-154, 207-212), G fraction
0.3-0.7 per arm (:92-104), no G run >= 5 in arms or window (:106-129), RNA MFE >= -10 (:156-165),
megablast against refseq_rna (config :30-34), any hit in another gene rejects (:1027-1061).
Only 78 and 47 of the 184 published windows pass its Tm and G-fraction rules (computed here): these
are not PRISM rules.

## 6. U-Probe

Nearest template: `padlock_iss` with its direct-RNA arms (`padlock-iss-hybiss.md` section 6);
`mipseq` does not fit (its primer binds the RNA). Edits to
`uprobe templates show padlock_iss -o protocol.yaml`:

```yaml
encoding:
  Slc17a7: {G: ACACGACGCTCTTCCGATCT, R: CGCAGCAGATAAATGAGCCA,
            B: TCTCGCTAGATCTTGTCAGT, Y: GGTCTTAGCACCTCGTGATC}   # ST2 barcode 1, "0040"
extracts:
  target_region: {source: mRNA, transcripts: mane, length: 40, step: 1}   # layout removed
probes:
  padlock:
    template: "{arm_5}TC{g}{r}{b}{y}{arm_3}"   # 64-panel backbones: no TC
    parts:
      arm_5: {expr: "rc(target_region[0:20])"}
      g: {expr: "encoding[target]['G']"}
      r: {expr: "encoding[target]['R']"}
      b: {expr: "encoding[target]['B']"}
      y: {expr: "encoding[target]['Y']"}
      arm_3: {expr: "rc(target_region[20:40])"}
hybridization: {temperature_c: 45, na_mM: 75, mg_mM: 10, probe_nM: 200,
                formamide_pct: 20, dntp_mM: 0}
attributes:   # target_gc, arm_*_tm, target_mapped_genes kept; junction_mapped_genes removed
  target_tm: {target: target_region, type: annealing_temperature, tm_method: multi_padlock_design,
              hybridization: {na_mM: 75, mg_mM: 10, probe_nM: 100, formamide_pct: 20, dntp_mM: 0}}
post_process:
  filters:    # target_tm and target_mapped_genes (<= 1) kept; junction filter removed
    g_run: {type: sequence_pattern, target: target_region, exclude_patterns: [GGGGG]}
  sorts: {is_ascending: [target_mapped_genes]}
  remove_overlap: {location_interval: 0}
  equal_space: {number_desired: 1}             # 3 as in the HCC panel
summary: {attributes: [target_gc, target_tm, arm_5_tm, arm_3_tm, target_mapped_genes]}
```

Put the gene's four segments (ST2 column I, or new orthogonal ones) in `encoding`, one barcode per
gene. On a toy genome holding Slc17a7 this rebuilds the published padlock exactly from its window.
Check with `check_probes.py --binding padlock.arm_5 --binding padlock.arm_3 --tm 65 75
--tm-column target_tm --min-gap 0`.
Pseudogene-rich genes (Gapdh): add `exempt_copies: true` to `target_mapped_genes`.

**U-Probe departs** (computed here on the 184 published windows; together they keep 164):
- Tm model unstated: `target_tm` keeps the template's HybISS formula and design conditions, which
  match the PRISM buffer but for 100 vs 200 nM probe: 174 in 65-75 C (none below, 10 at
  75-77.7 C). Primer3 in the PRISM buffer (200 nM) keeps 148 (34 below 65 C).
- G run checked on the window only (3 published windows hold GGGGG; 4 hold CCCCC, a G run in
  the arm).
- No structure filter: `fold_score` uses RNA parameters (24 windows have -MFE > 10); with DNA
  parameters (ViennaRNA 2.7, Mathews 2004, 37 C) all 184 have dG >= -10.
- One reference transcript, not "all variants"; sites spread by `equal_space`; Bowtie 1 on Ensembl
  instead of BLAST on NCBI. The template's junction search would remove 9 more published windows
  (e.g. Krt5 vs Krt4/6a/6b/7, Cdh1 vs Cdh3). For paralog-rich panels keep it with layout
  parts 10, 20 and 10 nt (junction = nt 11-30).

U-Probe cannot do:
- choose intensity levels, codewords (radius-vector filtering) or the gene-to-barcode assignment;
- design imaging-probe mixes (labelled:unlabelled ratios) or tune them to a microscope;
- design new orthogonal 20-nt segments, the RCA primers or the imaging probes (fixed ST2 sets);
- selective amplification: irreplicable padlocks and dilution ratios (construction not given);
- add the 5' phosphate, phosphorothioates or dyes; fold with DNA parameters; decode spots.

## 7. Pitfalls

- Swapped arms or sense arms (cDNA style) bind nothing; no 5' phosphate, no ligation.
- Codeword digits (Y, B, R, G) are not the segment order (G, R, B, Y). The 64-panel workbook writes
  Cy3 levels 0/2/5 where Extended Data Fig. 7 says 0/1/2.
- Use the empirical pooling, not the theoretical one (3R: 9.6:2.4 uL vs 3:1, ST2); ratios "can be
  finetuned based on actual optical setups" (v3).
- High-Cy3 64-panel codes need their own R/B/Y families ("FRET-like effect", ST2 64-panel note).
- Single-colour codes (0040, 0052) had more false calls (5.85 %, 7.20 %); avoid them for key
  genes (v3).
- Three 64-panel imaging cells (5B, 1G, 2G) start with U+FEFF: strip it before use.

## 8. Verification checklist

- [ ] Padlock = rc(RNA 1-20) + [TC] + G + R + B + Y + rc(RNA 21-40): 122 nt (120 without TC).
- [ ] rc(arm_5) + rc(arm_3) is the 40-nt window on the intended transcript.
- [ ] `target_tm` 65-75 C (model and conditions recorded); no GGGGG; `target_mapped_genes` = 1.
- [ ] One unique barcode per gene; its R segment has a primer, rc(R), in the RCA pool.
- [ ] Order notes: 5' phosphate on padlocks, phosphorothioates on primers, dyes on imaging probes.

## 9. References

- Chang T, Zhao S, Deng K, et al. 2025, Nat Biotechnol. High-plex spatial RNA imaging in one round
  with conventional microscopes using color-intensity barcodes. doi:10.1038/s41587-025-02883-7
  (paywalled: abstract, Extended Data legends, Supplementary Information and Tables 1-2 used).
- Chang T et al. 2025, bioRxiv v3 (2025-10-05), same title. doi:10.1101/2024.06.29.601330
  (Methods used; v1 and v2 not retrieved).
- HuangLab-PKU/probe-designer, commit d7326be; HuangLab-PKU/PRISM-Code, commit ab7b72d (GitHub).
- Krzywkowski T, Nilsson M 2017, Nucleic Acids Res. doi:10.1093/nar/gkx708
- Biosearch Technologies, Ampligase 10X Reaction Buffer (A1905B), product page.
