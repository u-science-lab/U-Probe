# smFISH with singly labelled oligos (Raj 2008 / Stellaris)

Read this when a user asks for smFISH, Stellaris or "Raj-style" probes: 20-nt oligos, each with one
3' fluorophore, 25-48 per RNA. Covers design rules, buffers, where Raj, Stellaris and tissue
protocols disagree, and the `smfish_stellaris` template.

## 1. Method and when to choose it

- 25-48 short DNA oligos tile one RNA; each carries one dye. One probe is invisible; dozens on one
  molecule make a diffraction-limited spot, so RNAs are counted one by one without amplification.
- Choose it for 1-3 RNAs per sample (one dye each), absolute counts per cell, nascent transcription
  sites (with intron probes), cultured cells, yeast, worms, flies, thin tissue sections.
- Choose smiFISH (`smifish.md`) for the same readout at lower cost or with switchable dyes; HCR v3
  (`hcr-v3.md`) for thick or autofluorescent samples or RNAs too short for ~25 probes;
  MERFISH / seqFISH+ for many genes.

## 2. Probe architecture (5' -> 3')

```
5'-[ rc(RNA), 20 nt (17-22) ]-3'-amine -> NHS-ester dye (TMR, Alexa 594, Cy5, Quasar, CAL Fluor)
RNA 3'-.......[ site 1 ]--2-3 nt--[ site 2 ]--2-3 nt--[ site 3 ]....... (x 25-48) -5'
```

- No flanks, no readouts: the probe is only the binding sequence; the dye is a synthesis
  modification (3' amine, coupled after synthesis, HPLC-purified), not part of the designed sequence.

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Length | 17-22 nt ("roughly 17-22 bases"); sets used 20 nt (Actb, 48) and 17 nt (d2EGFP, 48; Suppl. Fig. 4 labels these "20 bp") | | Raj 2008 Suppl. Methods and probe list, doi:10.1038/nmeth.1253 |
| Length (Stellaris) | 18-22 nt, typically 20; 18-19 nt in GC-rich regions, 21-22 nt in AT-rich regions | | Biosearch designer FAQ / blog |
| GC | "as close to 45 % as possible" (dynamic programming); the Actb set actually spans 35-80 % (mean 56 %) | no hard window | Raj 2008 |
| GC (Stellaris) | no number; thermodynamic placement, skewed-GC stretches skipped | | Biosearch FAQ |
| Tm | none published (Raj, Stellaris) | | |
| Spacing | >= 3 nt (Raj 2008); default 2 nt, may be 1 (Stellaris; CAL Fluor Red 635 needs 2) | | Raj 2008; Biosearch blog |
| Probes per RNA | 48 (Raj 2008); Stellaris max 48 per set, >= 25 recommended; tissue: > 20, 48 preferred, up to 96 | 12 or 24 probes: spots not reliably identifiable (Raj Suppl. Fig. 3) | Raj 2008; Stellaris guide; Lyubimova 2013 |
| Target length | >= ~1 kb for 48 x 20-mers with 2-nt gaps (1,054 nt); ~600 nt for a reduced set | | Lyubimova 2013; Stellaris FAQ |
| Labelling | 3' amine + NHS-ester dye, one per oligo | | Raj 2008 |
| Hyb. buffer (Raj) | 10 % dextran sulfate, 2x SSC, 10 % formamide, 2 mM VRC, 0.02 % BSA, 1 mg/mL E. coli tRNA; probe 5-50 ng per 50-100 uL | 37 C, 3-24 h (worms 30 C) | Raj 2008 Suppl. |
| Hyb. (Stellaris) | 10 % formamide in Stellaris buffer (older recipe: 100 mg/mL dextran sulfate, 10 % formamide, 2x SSC); 125 nM total probe (~2.6 nM per oligo for 48) | 37 C, 4-16 h | Stellaris adherent-cell protocol |
| Wash (Raj) | 10-20 % formamide, 2x SSC, 2 x 30 min | 30 C | Raj 2008 Suppl. |
| Wash (Stellaris) | Wash A (10 % formamide, 2x SSC) 30 min, then DAPI in Wash A 30 min; Wash B (2x SSC) 2-5 min | 37 C; Wash B RT | Stellaris protocol |
| Tissue (Lyubimova) | hybridization overnight and washes 30 min, 10 % formamide; raise to 25 % formamide if background is high | 30 C | Lyubimova 2013, doi:10.1038/nprot.2013.109 |

Disagreements to state to the user: spacing (Raj 3 nt vs Stellaris/Lyubimova 2 nt); wash temperature
(Raj and Lyubimova 30 C vs Stellaris 37 C); formamide in washes (10 % default, 15-25 % to cut
background); region (CDS only vs UTRs, below); probe dose (Stellaris 125 nM total vs Raj ~7-140 nM).

## 4. Target region

- Raj 2008 tiled the coding region of each endogenous gene (ATG to stop, one RefSeq/GenBank mRNA);
  Stellaris prefers "unique coding sequences" and suggests adding 5'/3' UTR to reach more probes, or
  using UTRs to tell paralogs apart; Lyubimova 2013 avoids UTRs where possible.
- U-Probe template: whole MANE Select mRNA (`source: mRNA`), so paralog-rich genes still reach 48;
  `source: CDS` reproduces Raj 2008.
- Isoforms: probes on constitutive exons count all isoforms; to target one isoform, list its transcript
  id and keep probes on its unique exons or junctions (check `exon_name` / `transcript_names`).
- Introns (transcription sites, iceFISH): up to 8 kb of introns shared by all isoforms, avoiding the
  last intron (Stellaris); use `source: genome` intervals.

## 5. Specificity

- One off-target probe gives no spot (a "single off-target smFISH probe is below detection limits",
  Calvo 2021); a false spot needs many probes of the set on the same off-target RNA. Abundant
  off-target RNAs bound by a few probes raise diffuse background.
- Stellaris rule: exclude 20-mers with > 16 nt matching other RNAs; if one off-target gene matches
  >= 5 probes with >= 16 nt, remove those probes. "The binding of a few probes to other targets is
  expected, and does not interfere."
- Highly expressed paralogs and processed pseudogenes (ACTB/ACTG1, GAPDH, ribosomal proteins) are the
  main problem: they remove most probes. In U-Probe use `exempt_copies` (copies of the target) and a
  `gene_expression` table with `min_fraction_of_target`, not a looser per-probe filter.

## 6. U-Probe template: `smfish_stellaris`

```bash
uprobe templates show smfish_stellaris -o protocol.yaml
```

```yaml
extracts:
  target_region: {source: mRNA, transcripts: mane, length: 20, step: 1}
probes:
  smfish_probe: {template: "{binding}", parts: {binding: {expr: "rc(target_region)"}}}
hybridization: {temperature_c: 37, na_mM: 390, formamide_pct: 10, probe_nM: 2.6, dntp_mM: 0}   # 2x SSC, 10 % FA; 125 nM / 48 oligos
attributes:
  target_mapped_genes:
    {target: target_region, type: mapped_genes, aligner: bowtie, max_mismatches: 2, strand: sense,
     exempt_copies: {kmer: 30, min_shared: 0.5}}
post_process:
  filters: {target_mapped_genes: {condition: target_mapped_genes <= 1}}
  remove_overlap: {location_interval: 2}
  equal_space: {number_desired: 48}
# gene_expression: {table: expression.tsv, value_column: TPM, min_value: 1, min_fraction_of_target: 0.05}
```

What to edit: `targets`; `genome` with `bowtie` in `align_index` (Bowtie2 cannot align 20-mers and
would report no off-targets); uncomment `gene_expression` with the sample's table;
`location_interval: 3` for Raj spacing; `source: CDS` for a CDS-only set. `length: [18, 22]` makes
windows of every length; only together with a GC or Tm filter does it let GC-skewed sites take
another length (`remove_overlap` keeps the shortest passing window at each start).

Template departures:
- Off-targets: Bowtie 1 within 2 substitutions (>= 18 of 20 nt). The Stellaris "> 16 nt" rule also
  catches 17-nt matches; `max_mismatches: 3` comes closer. Off-targets are counted per probe, not
  per off-target gene across the set (`off_target_pool` is the opt-in for that).
- GC: reported, not filtered or optimized toward 45 %; `equal_space` picks 48 by position.
- Region: whole mRNA (Raj used the CDS).
- Length: 20 nt, as Raj's Actb set and the Stellaris default; Raj's GFP (d2EGFP) set was 48 17-mers.

## 7. Pitfalls

- Strand: probes are the reverse complement of the RNA; never order the mRNA sequence itself
  (`strand: sense` for `target_region`).
- Short RNAs (< ~600 nt) cannot hold 25 probes; use HCR v3 (amplified) or accept dimmer spots and
  validate detection against a reference set.
- Genes with expressed copies (pseudogenes, retrocopies, readthrough genes) or close paralogs: check
  `target_mapped_genes_copies` / `_off_targets`; a set that detects a family must be reported as such.
- Masking: Stellaris masking level 5 (species-specific) excludes abundant-RNA cross-hybridization;
  U-Probe's equivalent is the off-target search with an expression table.
- "125 nM" in Stellaris protocols is total oligo, not per probe.
- Multi-gene experiments: one dye per gene; check that no probe of gene A lists gene B in
  `target_mapped_genes_off_targets`.

## 8. Checklist

- [ ] All probes 17-22 nt (20 default), reverse complement of the mRNA.
- [ ] >= 25 probes per RNA (48 target); `summary.json` status `ok` or shortfall explained.
- [ ] Gaps >= 2 nt (3 nt for Raj rules); no overlapping probes.
- [ ] GC near 45 % on average; no probe with extreme GC (< 30 % or > 70 %) if avoidable.
- [ ] `target_mapped_genes` = 1, or off-targets listed as unexpressed / copies; no off-target gene hit by
      >= 5 probes.
- [ ] Conditions stated: 10 % formamide, 2x SSC, 37 C hybridization; wash temperature chosen.

## 9. References

- Raj A, van den Bogaard P, Rifkin SA, van Oudenaarden A, Tyagi S. 2008. Imaging individual mRNA molecules using multiple singly labeled probes. Nat Methods 5:877-879. doi:10.1038/nmeth.1253
- Raj A, Tyagi S. 2010. Detection of individual endogenous RNA transcripts in situ using multiple singly labeled probes. Methods Enzymol 472:365-386. doi:10.1016/S0076-6879(10)72004-8
- Lyubimova A et al. 2013. Single-molecule mRNA detection and counting in mammalian tissue. Nat Protoc 8:1743-1758. doi:10.1038/nprot.2013.109
- Calvo L, Ronshaugen M, Pettini T. 2021. smiFISH and embryo segmentation for single-cell multi-gene RNA quantification in arthropods. Commun Biol 4:352. doi:10.1038/s42003-021-01803-0
- LGC Biosearch Technologies. Stellaris RNA FISH protocol for adherent cells; Stellaris Getting Started Guide; Stellaris Probe Designer FAQ (biosearchtech.com).
