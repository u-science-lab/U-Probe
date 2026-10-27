# smiFISH: gene-specific probes with a FLAP readout (Tsanov 2016)

Read this when a user asks for smiFISH, "FLAP" probes, or cheap unlabelled smFISH probes read by one
dye-labelled secondary oligo. Covers primary-probe layout, FLAP-X/Y/Z sequences, Oligostan rules,
pre-annealing and buffers. U-Probe has no smiFISH template; a validated protocol is given below.

## 1. Method and when to choose it

- 12-48 (typically 24) unlabelled primary probes tile one RNA. Each has a gene-specific part and a
  shared 28-nt FLAP-binding site. A 28-nt secondary "FLAP" oligo, labelled at both ends, is
  pre-annealed to the primaries in vitro; the duplexes are hybridized like smFISH probes.
- Same single-molecule readout as smFISH (`smfish.md`), at a fraction of the cost: primaries are
  ordinary unmodified oligos, and one labelled FLAP serves every gene; changing dye = changing FLAP.
- Choose it for many genes tested one or two at a time, dye flexibility, super-resolution (Tsanov
  2016 used it with STED and expansion). For thick tissue / amplification see `hcr-v3.md`.

## 2. Probe architecture (5' -> 3')

```
primary (<= 60 nt): 5'-[ rc(RNA), 26-32 nt ][ FLAP-binding site, 28 nt ]-3'
FLAP (secondary)  : 5'-dye-[ FLAP, 28 nt ]-dye-3'        (pairs with the FLAP-binding site)

RNA  3'-......[ 26-32 nt site ].....-5'
               | | | | | | | |
primary     5'-[ binding     ][ FLAP-binding ]-3'
                              | | | | | | | |
FLAP                     3'-dye[ FLAP 28   ]dye-5'
```

| FLAP | Secondary oligo (labelled), 5'->3' | Site on the primary = rc(FLAP), 5'->3' |
| --- | --- | --- |
| X | CACTGAGTCCAGCTCGAAACTTAGGAGG | CCTCCTAAGTTTCGAGCTGGACTCAGTG |
| Y | AATGCATGTCGACGAGGTCCGAGTGTAA | TTACACTCGGACCTCGTCGACATGCATT |
| Z | CTTATAGGGCATGGATGCTAGAAGCTGG | CCAGCTTCTAGCATCCATGCCCTATAAG |

- FLAP end on the primary: 3' in Oligostan (`paste(probe, FLAP)`), the Tsanov repository probe table
  ("added at the 3' end of target sequence") and Querido 2020; 5' in Calvo 2021 and Bouchet 2023.
  The paper text does not say. Both work; keep one convention per order.
- Naming trap: some protocols call the primary's site "FLAP X" (CCTCC...); others the labelled oligo
  (CACTG...). Always check that site = rc(labelled FLAP).
- Two or three genes at once: give each gene's primaries a different FLAP (X, Y, Z) and a different dye.

## 3. Design parameters

| Parameter | Value | Conditions / note | Source |
| --- | --- | --- | --- |
| Gene-specific length | 26-32 nt | chosen per site to hit the dG target | Tsanov 2016 Fig. 1A, doi:10.1093/nar/gkw784; Oligostan |
| Affinity | dG37 of the RNA/DNA duplex in [-36, -28] kcal/mol (Fig. 1A). Set-level rule: Oligostan.r tries each target dG from -36 to -28 in 0.5 steps, keeps the probes within +-1 kcal/mol of it (score 1 - 0.1 x \|dG - target\| >= 0.9) and takes the target that leaves the most probes (or a fixed one, e.g. -32) | dG = sum of 16 dinucleotide values over the probe 5'->3' minus a 0.115 M salt term, no initiation or formamide term; U-Probe's `rna_dna_dg37` (section 6) | Tsanov 2016 Fig. 1A; Oligostan.r |
| dG vs Tm design | dG-based sets beat Tm-based sets; 3 of 7 Tm-designed sets gave no signal | | Tsanov 2016 Suppl. Table 1 |
| GC | 40-60 % | | Oligostan |
| Sequence rules | "PNAS rules" (Xu 2009) on the probe: 1 A < 28 %; 2 no AAAA; 3 C 22-28 %; 4 no CCCC; 5 no 4 C in any 6 consecutive nt within the first 12; Oligostan script default 1, 2, 4 (its documentation says all 5) | default 1, 2, 4; all five when the user follows the documentation. Hard filters: never drop one to fill the count (SKILL.md rule 10) | Oligostan; Xu 2009, doi:10.1073/pnas.0812506106 |
| Repeats | <= 10 % RepeatMasker-masked | | Oligostan |
| Spacing | >= 2 nt | | Oligostan |
| Probes per RNA | 12-48; 24 "a good compromise between cost and signal quality" | tested 12, 24, 45 | Tsanov 2016 |
| FLAP | 28 nt, dG37 -35 to -31 kcal/mol, dye at both ends | | Tsanov 2016 |
| Pre-annealing | 40 pmol primary set (24 probes, equimolar) + 50 pmol FLAP in 10 uL 1x NEBuffer 3; 85 C 3 min, 65 C 3 min, 25 C 5 min | duplexes survive 2-5 freeze-thaws | Tsanov 2016 Suppl. Protocol |
| Hyb. mix (100 uL) | 1x SSC, 15 % formamide, ~10.6 % dextran sulfate, 0.34 mg/mL tRNA, 0.2 mg/mL BSA, 2 mM VRC, 2 uL duplexes (~80 nM total primary, ~100 nM FLAP) | 37 C overnight, after 15 min in 15 % formamide / 1x SSC at RT | Tsanov 2016 |
| Washes | 15 % formamide, 1x SSC, 2 x 30 min | 37 C | Tsanov 2016 |
| Alternative buffer | 10 % formamide, 2x SSC, 100 mg/mL dextran sulfate (Stellaris buffers also work) | 37 C | Tsanov 2016; Calvo 2021 |

Later protocols changed the probes:

| Protocol | Gene-specific part | FLAP | Buffer | Probe dose |
| --- | --- | --- | --- | --- |
| Tsanov 2016 (cells) | 26-32 nt, dG-matched (Oligostan) | X/Y/Z, 3' | 15 % FA, 1x SSC | ~80 nM total |
| Querido 2020 (hTR, 451 nt) | >= 23 nt, GC up to 85 %, 15 probes (for that RNA) | Y, 3' | 15 % FA, 1x SSC | 80 nM |
| Calvo 2021 (arthropod embryos) | 20 nt, Stellaris designer, <= 48 | X, 5' | 10 % FA, 2x SSC | 80 nM, 37 C 14 h |
| Bouchet 2023 (Drosophila) | 20 nt, GC > 40 %, <= 48 | X, 5' | 10 % FA, 2x SSC, 5 % DS | 0.4-2 uM |

## 4. Target region

- Oligostan: "We usually use only the coding sequence and exclude the UTR"; add UTRs when the CDS is
  short (< 500 bp). Querido 2020's length and GC (table above) belong to that protocol: use them only
  when the user chooses it, not to fill a short Tsanov-style set.
- In U-Probe: `source: CDS`, `transcripts: mane`; switch to `source: mRNA` for short-CDS genes or
  when fewer than 24 probes pass.
- Isoforms and introns: as in `smfish.md` section 4.

## 5. Specificity

- As in smFISH, a spot needs many primaries on one RNA; a single off-target primary is not detected.
  The protocol keeps the templates' off-target rule (Bowtie 1 within 2 substitutions, plus an expression
  table; SKILL.md rules 4-5); a stricter search is not a published smiFISH rule (SKILL.md rule 10).
- Abundant paralogs/pseudogenes remove probes, as in smFISH; use `exempt_copies` and `gene_expression`.
- The FLAP must not bind the transcriptome or the gene-specific parts; X, Y, Z were designed for that.
  Do not invent a new FLAP without screening it.

## 6. U-Probe protocol (no built-in template; adapted from `smfish_stellaris`, passes `uprobe validate`)

`rna_dna_dg37` computes Oligostan's dG37 (the script's 16 dinucleotide values summed over the probe read
5'->3', plus its 0.115 M salt term; it reproduces the dG37 column of the Oligostan tutorial output to
1e-9). Put it on the probe's binding part, as below; on `target_region` (the RNA sense sequence) it needs
`sequence: target`, otherwise it gives the other strand's value and `validate` warns.

```yaml
name: smifish
genome: hg38                 # needs bowtie in align_index
targets: [GAPDH]
encoding:
  GAPDH: {flap: CACTGAGTCCAGCTCGAAACTTAGGAGG}   # FLAP-X, the dye-labelled secondary
extracts:
  target_region: {source: CDS, transcripts: mane, length: [26, 32], step: 1}
probes:
  primary:
    template: "{binding}{flap_site}"            # FLAP-binding site at the 3' end (Oligostan)
    parts:
      binding: {expr: "rc(target_region)"}
      flap_site: {expr: "rc(encoding[target]['flap'])"}   # CCTCCTAAGTTTCGAGCTGGACTCAGTG
hybridization: {temperature_c: 37, na_mM: 195, probe_nM: 3.3, formamide_pct: 15}   # 1x SSC
attributes:
  target_gc: {target: target_region, type: gc_content}
  dg37: {target: primary.binding, type: rna_dna_dg37}   # the probe 5'->3' (sequence: probe, the default)
  binding_tm: {target: primary.binding, type: annealing_temperature}   # reported, not filtered
  target_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie, max_mismatches: 2,
                        strand: sense, exempt_copies: {kmer: 30, min_shared: 0.5}}
post_process:
  filters:
    target_gc: {condition: target_gc >= 0.40 & target_gc <= 0.60}
    dg37: {condition: dg37 >= -33 & dg37 <= -31}   # one target dG (-32) +- 1, inside [-36, -28]
    pnas_rules_2_4: {type: sequence_pattern, target: primary.binding, exclude_patterns: [AAAA, CCCC]}
    target_mapped_genes: {condition: target_mapped_genes <= 1}
  remove_overlap: {location_interval: 2}
  equal_space: {number_desired: 24}
  sorts: {is_ascending: [target_mapped_genes]}
summary: {report_name: rna_report, attributes: [target_gc, dg37, binding_tm, target_mapped_genes]}
```

How it departs from Oligostan (say so to the user):
- Target dG: the filter applies Oligostan's set-level rule (every probe within +-1 kcal/mol of one target,
  score >= 0.9) for a fixed target, -32, the example of the script's `dG37Desiree` setting. Oligostan's
  default instead tries every target from -36 to -28 in 0.5 steps and keeps the one that leaves the most
  probes: run the protocol with a few targets (band `T-1` to `T+1`), keep the run with the most probes and
  report T. `length: [26, 32]` makes every length, and `remove_overlap` keeps the shortest passing window
  at each start (Oligostan takes the length nearest the target); both stay inside the band.
- PNAS rules 2 and 4 are filters; rule 1 (A < 28 %, part of the default) and rules 3 and 5 (when all five
  apply) need a post-filter on the TSV before `equal_space` picks the set, e.g. on a `--raw` run
  `b = df['primary.binding']; df = df[b.str.count('A') / b.str.len() < 0.28]`.
- The script's dG table is not cited to Wu 2002 in the code and was not checked against the published
  tables; compare with Oligostan output for critical genes.
- RepeatMasker masking is not applied; the off-target search removes most repeat windows.
- Optional: `gene_expression` as in `smfish_stellaris`; FLAP at the 5' end: `"{flap_site}{binding}"`.

## 7. Pitfalls

- Wrong FLAP strand (most common error): the primary must carry rc(labelled FLAP).
- Probe length > 60 nt raises synthesis cost (Tsanov's limit: 32 + 28).
- Forgetting pre-annealing, or using too little FLAP: unlabelled primaries compete with labelled
  duplexes. Use a slight FLAP excess (50 vs 40 pmol).
- Mixing conventions: 20-nt Stellaris-style primaries + FLAP (Calvo, Bouchet) work but are a different
  design from Tsanov's dG-matched 26-32-mers; state which one you designed.
- Two genes with the same FLAP in one sample cannot be told apart.

## 8. Checklist

- [ ] Gene-specific parts 26-32 nt (or the user's chosen 20-nt variant), rc of the RNA.
- [ ] dG37 (`rna_dna_dg37` on the binding part) of every gene-specific part in [-36, -28] and within +-1 kcal/mol of one target dG T (section 6); T reported.
- [ ] FLAP-binding site = rc(chosen FLAP), at the agreed end; total length <= 60 nt.
- [ ] 24 probes per RNA (12-48), gaps >= 2 nt; GC 40-60 %; PNAS rules 1, 2, 4 (A < 28 %, no AAAA / CCCC; also 3 and 5 if all five apply).
- [ ] `target_mapped_genes` = 1 or explained; one FLAP (and dye) per gene imaged together.
- [ ] Pre-annealing and buffer plan stated (15 % formamide, 1x SSC, 37 C).
- [ ] Run last, with the file's own values (the dG window of Fig. 1A; with T = -32, `--range dg37 -33 -31` checks the band):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --summary results/run1/summary.json \
    --transcripts <genome dir>/transcript.fa --binding primary.binding --gc 0.40 0.60 --range dg37 -36 -28 \
    --min-gap 2 --min-probes 12
```

## 9. References

- Tsanov N et al. 2016. smiFISH and FISH-quant - a flexible single RNA detection approach with super-resolution capability. Nucleic Acids Res 44:e165. doi:10.1093/nar/gkw784
- Xu Q, Schlabach MR, Hannon GJ, Elledge SJ. 2009. Design of 240,000 orthogonal 25mer DNA barcode probes. PNAS 106:2289-2294. doi:10.1073/pnas.0812506106
- Wu P, Nakano S, Sugimoto N. 2002. Temperature dependence of thermodynamic properties for DNA/DNA and RNA/DNA duplex formation. Eur J Biochem 269:2821-2830. doi:10.1046/j.1432-1033.2002.02970.x
- Querido E, Sfeir A, Chartrand P. 2020. Imaging of telomerase RNA by single-molecule inexpensive FISH combined with immunofluorescence. STAR Protoc 1:100104. doi:10.1016/j.xpro.2020.100104
- Calvo L, Ronshaugen M, Pettini T. 2021. smiFISH and embryo segmentation for single-cell multi-gene RNA quantification in arthropods. Commun Biol 4:352. doi:10.1038/s42003-021-01803-0
- Bouchet M et al. 2023. A simple smiFISH pipeline to quantify mRNA at the single-cell level in 3D. STAR Protoc 4:102316. doi:10.1016/j.xpro.2023.102316
- Oligostan (R script and documentation), FISH-quant repository: bitbucket.org/muellerflorian/fish_quant (Oligostan/)
