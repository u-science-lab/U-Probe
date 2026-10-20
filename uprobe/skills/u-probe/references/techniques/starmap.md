# STARmap SNAIL probes (padlock + primer pairs)

> Read this when: designing STARmap / STARmap PLUS / RIBOmap-style SNAIL probe pairs, editing the
> `starmap` template, choosing gene codes for SEDAL readout, or judging the specificity of a
> two-oligo padlock+primer unit on RNA. General rules: `references/concepts/probe-architecture.md`
> (padlock + primer units) and `references/concepts/specificity.md`.

## 1. What it is and when to choose it

STARmap (Wang 2018) detects RNA without reverse transcription. A SNAIL pair (specific
amplification of nucleic acids via intramolecular ligation) is two DNA oligos on adjacent sites of
one RNA: a padlock and a primer. The primer's 3' tail is the splint that circularizes the padlock
(T4 DNA ligase), and the same 3' end primes phi29 rolling-circle amplification (RCA). The amplicon
(with aminoallyl-dUTP) is copolymerized into a hydrogel and the gene code in it is read by SEDAL
sequencing by ligation. Only when both oligos sit on the same RNA does an amplicon form, which
"excludes noise" from single-probe mis-hybridization (Wang 2018).

Choose it for highly multiplexed (10^2-10^3 genes) single-cell spatial transcriptomics in tissue,
when bright spots (RCA) matter, transcripts are short, or 3D/thick samples are planned
(STARmap 3D, Deep-STARmap). Variants: STARmap PLUS (Zeng 2023; protein co-staining, larger code
space), RIBOmap (tri-probe: a third splint oligo on 18S rRNA circularizes the padlock, so only
ribosome-bound mRNA is amplified), Deep-STARmap (Sui 2025; a 24-nt 5' linker on the primer pairs
with an acrydite + CNVK adaptor that anchors it in the gel). Prefer a padlock ligated on the target
itself (cDNA in HybISS, RNA in MiP-seq; see `padlock-iss-hybiss.md`, `mip-seq.md`) when
single-nucleotide discrimination at a junction is needed: in SNAIL the ligation template is the
primer (DNA), not the RNA.

## 2. Probe architecture

Target window: RNA (sense) nt 1..L, L = 40-46, read 5'->3'. Padlock site is UPSTREAM, primer site
DOWNSTREAM (checked: 4,178 of 4,178 located published pairs, 2018 / 64-gene / 2,766-gene sets and the
4,086 pairs of Shi 2023, mapped to Ensembl 113 mouse cDNAs).

```text
RNA 5'-[ padlock_site 19-25 ][gap 1-2][ primer_site 19-25 ]-3'        L = 40-46   (0 allowed in the text, never used)
         ^nt 1                                     nt L^

Padlock (U-Probe template = STARmap PLUS 64-gene backbone; 5'-phosphate), 52-58 nt:
5'-P- A c1c2c3 TA | rc(padlock_site) | AATTATTAC | c1c2c3 | CATACACTAAAGATA -3'
      [ flap 6  ]   [ 19-25, on RNA ]   [  9    ]  [code 3]  [15: reading-probe site
       pairs primer                                           CATACACTAAAG; last 6 AAGATA
       tail nt 1-6                                            pair primer tail nt 7-12]
Primer (no phosphate), 31-37 nt:
5'- rc(primer_site) | TA c3'c2'c1' T | TATCTT -3'      c' = complement; TA rc(code) T
    [ 19-25, on RNA ] [6: pairs padlock 5' end] [6: pairs padlock 3' end] -> RCA starts here
```

Bound state. The padlock binding part's 3' end sits on nt 1 and its 5' end at the gap; the
primer binding part's 3' end sits at the gap, so its 12-nt tail and the padlock's 5' flap meet
there. The padlock's 27-29-nt 3' arm loops back so its last 6 nt pair with the end of the tail:

```text
primer tail   5'-  T  A  c3' c2' c1' T  T  A  T  C  T  T  -3'  (primer 3'-OH primes RCA)
                   |  |  |   |   |   |  |  |  |  |  |  |
padlock ends  3'-  A  T  c3  c2  c1  A ^A  T  A  G  A  A  ...   (written 3'<-5')
                  [padlock 5' end, 5'-P][padlock 3' end, 3'-OH]
                                       ^ ligation junction: on the primer (DNA), not on the RNA
```

Each 6-bp half has a Tm below room temperature, so primer and padlock do not pair during the
40 C hybridization but are ligated afterwards at room temperature (Wang 2018, Fig. S1A). Because
ligation is templated by the primer, the 1-2 nt gap on the RNA needs no gap fill.

Backbones of the published sets (derived from the supplementary sequence tables; "bc" is the
SEDAL-read code, "s" the splint code):

| Set | Padlock 5' -> 3' | Primer tail | Code | SEDAL |
| --- | --- | --- | --- | --- |
| Wang 2018 (160, 1,020 genes) | P-ACATTA-[site]-AATTATTAC-bc5-CATACACTAAAGATA | TAATGT TATCTT (constant 6+6) | bc 5 nt (1,024 codes) | 6 rounds, 2-base |
| Shi 2023 (1,022 genes, STARmap PLUS) | P-A s3 TA-[site]-AATTATTAC-bc5-CATACACTAAAGATA | TA rc(s3) T TATCTT | bc 5 nt + splint s3 | 6 rounds, 2-base |
| Zeng 2023 PLUS 64 genes (= U-Probe template) | P-A c3 TA-[site]-AATTATTAC-c3-CATACACTAAAGATA | TA rc(c3) T TATCTT | one 3-nt code, used twice (64 codes) | 4 rounds |
| Zeng 2023 PLUS 2,766 genes | P-GTAGACTA-[site]-AAC n3 CATACACTAAAGATT n3' TC | TAGTCTAC GA rc(n3') A (8+6) | two 3-nt blocks; paper text: "10 nt" | 8 rounds |

- The splint code (s3) blocks cross-gene ligation: it "reduce[s] ... false positives" when the
  primer of transcript A circularizes the padlock on transcript B (Shi 2023 Methods).
- SEDAL (2018 Table S2): reading probes R1-R6 are 11-nt, 5'-phosphate, 3'-inverted-dT copies of
  the padlock sequence `CATACACTAAAG`, shifted by one degenerate N per cycle; 8-nt decoding
  probes `NNNNNNXY` (16 dinucleotides, 4 dyes) ligate to the reading probe's 5' end. From these
  sequences, cycles 1-6 read the padlock-sense dinucleotides (b5,C) (b4,b5) ... (C,b1): the code
  between two known C (G in the amplicon), so a single color error corrupts the final known base
  and the read is rejected (Fig. S3H). Products are stripped with 60% formamide between cycles.
  The 64-gene set's 4 rounds fit the same scheme for a 3-nt code (inferred, not stated).
- The 28-gene 3D set (2018) swapped the code for per-gene orthogonal 11-nt reading sites
  (OR1-OR7) read with 1-base probes, 4 genes per cycle, and 5'-amino primers.

## 3. Design parameters

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Window (padlock + gap + primer site) | 40-46 | nt | Wang 2018 SM, doi:10.1126/science.aat5691 |
| Each half | 19-25 (text says 20-25; Fig. S1A and all published sets 19-25) | nt | Wang 2018 SM; Zeng 2023 Supp Table 1 |
| Gap between halves | 0-2 stated (2018), 1-2 (Ren 2025); 1 (1,286) or 2 (2,892) in every located published pair, never 0 | nt | Wang 2018 SM; Ren 2025; this KB, 4,178 pairs |
| Half Tm | 60 designed, halves Tm-matched (Picky 2.2; salt/oligo conditions not stated). Recomputed at 2x SSC, 10% FA, 100 nM: median 58.0, 5th-95th percentile 53.4-63.9 (8,444 published halves); the two halves of a pair differ by a median 1.0 C (p90 3.1) | C | Wang 2018 Fig. S1A; this KB |
| Window GC, Tm (2025) | GC 30-70 %, Tm 56-70 C at 5 nM oligo and 300 mM salt (Picky's salt term: [Na+]), set in Picky on the whole 40-46-nt sequence before the split (section 6) | %, C | Ren 2025, "Padlock and primer probe design"; Chou 2010 |
| Splint overlap | 6 + 6 (2018, 1,022, 64-gene sets); 8 + 6 (2,766-gene set) | nt | Fig. S1A; supp tables |
| Splint Tm | below room temperature | C | Wang 2018 Fig. S1A |
| Pairs per gene | 4 (2018, 1,022, 64-gene); 2 (2,766-gene); 3-6 kept, > 6 designed as filtering removes 10-20 % (Ren 2025) | pairs | supp tables; Ren 2025 |
| Oligo length | padlock 54-60, primer 31-37 (2018); "under 60 nt" | nt | Wang 2018 SM, Table S1 |
| Code length | 5 (2018, 1,022); 3 (64-gene); 2 x 3 (2,766) | nt | tables; Zeng 2023; Shi 2023 doi:10.1038/s41586-023-06569-5 |
| Hybridization | 2x SSC, 10% formamide, 1% Tween-20, 20 mM RVC, 40 C | - | Wang 2018 SM |
| Probe concentration / time | 100 nM per oligo, overnight (2018); 1 nM, 36 h (PLUS); 5 nM, 24-36 h (Shi 2023); 100 nM, or 1-5 nM above 100 genes, 12-36 h (Ren 2025 Steps 16-18) | nM | papers' Methods |
| Ligation | T4 DNA ligase, room temperature, 2 h | - | Wang 2018; Zeng 2023 doi:10.1038/s41593-022-01251-x |
| RCA | phi29, 250 uM dNTP + 20 uM aminoallyl-dUTP, 30 C, 2 h | - | Wang 2018; Zeng 2023 |
| SEDAL | reading probe 11 nt (10 uM) + decoding 8 nt, T4 ligase, RT 3 h; strip 60% formamide | - | Wang 2018 SM |

**Protocols disagree:** probe concentration (100 vs 1-5 nM), code length and splint design
differ between versions; half length is 20-25 in the 2018 text but 19-25 in Fig. S1A and in
the published sequences.

## 4. Target region and template

- Template: RNA directly (no RT, no cDNA). Wang 2018: "only the shortest isoforms were considered
  and the coding regions were used except for non-coding RNAs". All 36 sampled 2018 pairs lie in
  the CDS of the canonical Ensembl transcript; in later sets 11-17% of sampled pairs lie outside
  the CDS (64-gene: 11/64; 2,766-gene: 4/36), so PLUS is not strictly CDS-only.
- Both halves must be on the same molecule, so the 40-46 nt window must be contiguous on the
  mature RNA. With one selected transcript, U-Probe also tiles across exon junctions.
- Isoforms: a pair detects every isoform containing the whole window. Pick a transcript shared by
  the isoforms you want to count; of the 136 pairs of the 2018, 64-gene and 2,766-gene samples, 28 lay
  on all Ensembl isoforms of their gene and 103 more on at least half.
- Spacing: published pairs of one gene are usually hundreds of nt apart but sometimes overlap
  (by up to 11 nt in the 2018 sample); no minimum spacing is stated. Ren 2025 removes "overlaps and
  containments" (Dedupe, BBMap), yet its preselected windows still overlap by up to 19 nt (this KB).
- Ren 2025 Supp Tables 1-2: preselected pairs for 25,916 human / 23,519 mouse genes; columns Primer
  and Padlock are the binding parts, 5'->3' (2,766 sampled mouse rows located: padlock site upstream,
  gap 1-2 nt; this KB). Cross-check designs against them; their "Tm" columns (51-90 C per half) are
  on an unstated scale.

**U-Probe departs:** `source: CDS` of the MANE Select (else Ensembl canonical) transcript instead
of the shortest isoform; `remove_overlap: 0` forbids overlapping pairs (stricter than published).

## 5. Specificity

- One SNAIL pair makes one amplicon, i.e. one countable spot. A false spot needs the padlock
  and primer to bind ADJACENTLY on the SAME RNA molecule (both halves of a 40-46 nt region on one
  transcript, at the right spacing). A padlock alone is never circularized (no splint) and a
  primer alone has no circle to prime, so single-oligo off-target binding gives no signal
  (Wang 2018, Fig. 1C-D: mismatched pairs give no labeling).
- Hence the unit of off-target search is the whole window: the template maps `target_region`
  (padlock site + gap + primer site) with `strand: sense` and keeps `target_mapped_genes <= 1`
  (no other gene carries the window). Because one unit is enough for a spot, do not tolerate
  off-target genes per pair as smFISH does; use `gene_expression` to ignore unexpressed genes and
  `exempt_copies` for pseudogenes/readthrough copies (`references/concepts/specificity.md`).
- Second failure mode: intermolecular cross-pairing (primer of gene A splinting the padlock of
  gene B on a neighboring molecule) gives a wrong code. Gene-specific splint codes (1,022-,
  64- and 2,766-gene sets) reduce it; the 2018 constant splint (ACATTA/TAATGT) does not.
- Neither Wang 2018 nor PLUS states a specificity criterion beyond Picky's design.
  **U-Probe departs:** specificity is the U-Probe approximation `target_mapped_genes <= 1`
  (Bowtie2 local search by default, which can miss some sites; see section 7).

## 6. U-Probe template

```bash
uprobe templates show starmap -o protocol.yaml        # copy, then edit meta.replace items
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets
uprobe run -p protocol.yaml -g genomes.yaml -o results --threads 10
```

Key lines of uprobe/templates/starmap.yaml (condensed excerpt; `...` marks omitted lines):

```yaml
encoding:
  Aldoc:
    code: AAA                       # 3-nt code (64 possible), unique per gene
extracts:
  target_region:
    source: CDS
    transcripts: mane
    length: [40, 46]
    layout:
      template: "{padlock_site}{gap}{primer_site}"     # padlock upstream, primer downstream
      parts: {padlock_site: {length: [19, 25]}, gap: {length: [1, 2]}, primer_site: {length: [19, 25]}}
probes:
  padlock:
    template: "A{code_5}TA{binding}AATTATTAC{code}CATACACTAAAGATA"
    ...  binding: {expr: "rc(target_parts['padlock_site'])"}; code_5 and code: encoding code
  primer:
    template: "{binding}TA{code_rc}TTATCTT"
    ...  binding: {expr: "rc(target_parts['primer_site'])"}; code_rc: rc(encoding code)
hybridization: {temperature_c: 40, na_mM: 390, probe_nM: 100, formamide_pct: 10, dntp_mM: 0}
attributes:
  padlock_tm: {target: padlock.binding, type: annealing_temperature}   # primer_tm alike
  tm_balance: {type: difference, of: [padlock_tm, primer_tm]}          # |padlock_tm - primer_tm|
  target_mapped_genes: {target: target_region, type: mapped_genes, aligner: bowtie2, strand: sense}
post_process:
  filters:
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    padlock_tm: {condition: padlock_tm >= 50 & padlock_tm <= 65}      # primer_tm alike
  remove_overlap: {location_interval: 0, prefer: tm_balance}           # best Tm-matched split per start
  equal_space: {number_desired: 4}
```

What to edit:
- `genome`, `targets`, and `encoding`: one code per gene, the same code feeds `code_5`, `code`
  and `code_rc`. Order the padlock with a 5'-phosphate (`/5Phos/`) - the template string has none.
- More than 64 genes: switch to the 1,022-gene backbone, e.g. padlock
  `A{s}TA{binding}AATTATTAC{bc}CATACACTAAAGATA`, primer `{binding}TA{s_rc}TTATCTT` with a 5-nt
  `bc` (6 SEDAL rounds) and a separate 3-nt splint code `s`; codes can come from
  `uprobe generate-barcodes`.
- Exhaustive off-target search: add `bowtie` to the genome's `align_index`, then set
  `aligner: bowtie` and `max_mismatches: 2` (Bowtie2 flagged 97-99% of 40-46-nt windows that
  have a <=2-substitution site in another gene; README).
- Tm-balanced halves: the template keeps halves of 50-65 C (both halves of 96.8% of the 9,910
  published pairs lie there under its conditions) and, at each start, the split with the smallest
  `tm_balance` (`remove_overlap.prefer`), as Picky picked the best Tm match. Narrow the window to your
  conditions if you change `probe_nM`.
- Probe concentration: `probe_nM: 100` is Wang 2018's; use 1 or 5 nM for a PLUS-style hybridization,
  which lowers every Tm by 6.4 or 4.2 C.
- Ren 2025 windows (checked here on the 2,766 located preselected mouse pairs): GC applies to the
  whole window (all 32.6-70.7 %; 12.5 % of pairs have a half outside 30-70 %). Do not filter a
  whole-window Tm at 56-70: with an attribute-level `hybridization: {na_mM: 300, probe_nM: 5,
  formamide_pct: 0}` U-Probe puts these windows at 72-87 C (1st-99th percentile), so Picky's Tm is on
  another scale. Their halves run hotter: under the template's conditions 50-65 C keeps both halves
  of 80.5 % of the pairs, 50-70 C 97.3 % (not a published rule):

```yaml
attributes:
  target_gc: {target: target_region, type: gc_content}     # the whole 40-46-nt window
post_process:
  filters:
    target_gc: {condition: target_gc >= 0.30 & target_gc <= 0.70}
    padlock_tm: {condition: padlock_tm >= 50 & padlock_tm <= 70}   # primer_tm alike
```

## 7. Pitfalls

- Orientation: swapping the halves (padlock downstream) puts the primer tail away from the
  padlock 5' flap; keep `{padlock_site}{gap}{primer_site}` and `rc()` on both binding parts.
- The 3-nt code is used three times (padlock 5' flap, padlock code, primer tail); a typo in one
  breaks either ligation (splint) or decoding (code). Check that `code_rc == rc(code)`.
- A 3-nt code caps the panel at 64 genes (the published 64-gene set uses every code, AAA
  included); more genes need the 5-nt (1,022-gene) or two-block (2,766-gene) backbone.
- Without `prefer`, `remove_overlap` keeps the first candidate at each start: the shortest split,
  19 + 1 + 20 nt, for every unit. The published pairs vary: each site 19-25 nt (19-21 in 69%), the
  window 40-46 nt and the gap 2 nt in 72%; the template's `prefer: tm_balance` lets them vary too.
- The designed 60 C Tm was computed by Picky under unstated conditions, so the template records it and
  does not aim at it. Under the hybridization block (2x SSC, 10% formamide, 100 nM) the published halves
  sit at 58.0 C (median; 5th-95th percentile 53.5-63.3), which is close but not the same quantity; the
  50-65 C filters keep both halves of 96.8% of the pairs. Do not filter at exactly 60.
- Probe concentration differs between the protocols: 100 nM per oligo overnight (Wang 2018, the
  template's `probe_nM`), 1 nM for 36 h (STARmap PLUS, whose 64-gene backbone the template uses) and
  5 nM for 24-36 h (Shi 2023). Tm follows: the same halves are 6.4 C lower at 1 nM and 4.2 C lower at
  5 nM, so set `probe_nM` to what you will hybridize at before judging any Tm window.
- The gap is `[1, 2]`: abutting sites (the 0 the text allows) leave the padlock's ligated ends and
  the primer's 3' end no room between the two duplexes, and no published pair uses them.
- Bowtie2 is a heuristic with at most 100 alignments: windows in repeats can pass with
  `target_mapped_genes <= 1`; prefer `aligner: bowtie` when the index exists.
- No fold or self-complementarity filter is applied; a padlock whose binding part pairs with its
  own backbone or with the primer tail will not hybridize/ligate as intended. Ren 2025 drops strong
  NUPACK secondary structures and BLASTs each probe against the minus strand of the whole pool
  (heterodimers); U-Probe has `fold_score` (ViennaRNA MFE as RNA; no published cutoff) and
  `self_match`, but no probe-vs-pool screen: check the final pool externally.

## 8. Verification checklist

- [ ] Each pair: rc(padlock binding) + 1-2 nt + rc(primer binding) is a contiguous 40-46 nt
      substring of the target mRNA, padlock site first (5'); the sites never abut.
- [ ] Each binding part 19-25 nt; padlock 52-60 nt with 5'-phosphate; primer 31-37 nt, no phosphate.
- [ ] Primer 3' tail = rc(padlock 3'-terminal 6 nt + padlock 5'-terminal 6 nt) (8 nt at the 5'
      end for the 2,766-gene backbone): a splint for its own padlock and no other gene's.
- [ ] Codes unique per gene; Hamming distance as large as the code space allows; read positions
      flanked by the known bases of the reading-probe design.
- [ ] `target_mapped_genes == 1` (after `gene_expression` / `exempt_copies`) for every pair;
      off-target lists empty.
- [ ] 4 pairs per gene (3-6 in Ren 2025; 2 for very large panels), non-overlapping, on the transcript
      you mean to count; `summary.json` status `ok` and no `fewer_probes_than_requested` warning.
- [ ] Half Tms reported at the probe concentration you will use and reasonably balanced (published
      pairs differ by ~1 C); no strong hairpins in either oligo, no heterodimers within the pool.

## 9. References

- Wang X, et al. 2018, Science. Three-dimensional intact-tissue sequencing of single-cell
  transcriptional states. doi:10.1126/science.aat5691 (SM: probe design, Fig. S1A, Fig. S3;
  Tables S1-S2 sequences)
- Zeng H, et al. 2023, Nature Neuroscience. Integrative in situ mapping of single-cell
  transcriptional states and tissue histopathology in a mouse model of Alzheimer's disease
  (STARmap PLUS). doi:10.1038/s41593-022-01251-x (Supp Table 1: 2,766- and 64-gene sets)
- Shi H, et al. 2023, Nature. Spatial atlas of the mouse central nervous system at molecular
  resolution. doi:10.1038/s41586-023-06569-5 (Methods: splint barcode; Supp Tables 1-2)
- Zeng H, et al. 2023, Science. Spatially resolved single-cell translatomics at molecular
  resolution (RIBOmap). doi:10.1126/science.add3067
- Sui X, et al. 2025, Nature Methods. Scalable spatial single-cell transcriptomics and
  translatomics in 3D thick tissue blocks (Deep-STARmap). doi:10.1038/s41592-025-02867-0
- Frei AP, et al. 2016, Nature Methods. Highly multiplexed simultaneous detection of RNAs and
  proteins in single cells (PLAYR; an earlier two-probe ligation assay that shares authors with STARmap).
  doi:10.1038/nmeth.3742
- Ren J, Zeng H, Huang J, et al. 2025, Nature Protocols. Spatially resolved in situ profiling of
  mRNA life cycle at transcriptome scale in intact cells and tissues using STARmap PLUS, RIBOmap and
  TEMPOmap. doi:10.1038/s41596-025-01248-3 (PMC12648639 author manuscript: "Padlock and primer
  probe design", Steps 16-18, Supp Tables 1-2)
- Chou HH, et al. 2004, Bioinformatics. Picky: oligo microarray design for large genomes.
  doi:10.1093/bioinformatics/bth347
- Chou HH. 2010, BMC Bioinformatics. Shared probe design and existing microarray reanalysis using
  PICKY. doi:10.1186/1471-2105-11-196 (salt correction on [Na+])
- Nilsson M, et al. 1994, Science. Padlock probes: circularizing oligonucleotides for localized
  DNA detection. doi:10.1126/science.7522346
