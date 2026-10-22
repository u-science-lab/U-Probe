# Specificity: off-targets, copies, paralogs and repeats

Read this when you choose or tune the off-target attribute (`mapped_genes`, `mapped_sites`,
`kmer_count`), its threshold, an expression table, `exempt_copies` or `off_target_pool`, or when
specificity filters leave too few probes. Each rule below rests on a cited result or a stated mechanism.

## Why search near-matches exhaustively

Mechanism: FISH oligos hybridize well below their Tm, so a site with one or two mismatches is still
largely bound. Under 2x SSC, one central mismatch lowers a 20-mer's Tm by ~6 °C and a 30-mer's by
~3 °C, two central mismatches by ~13 °C and ~7 °C (Primer3 thermodynamic alignment, computed for this
skill; terminal mismatches cost < 1 °C), while published RNA-FISH designs run 5-20 °C above the
hybridization temperature (`references/concepts/hybridization-and-tm.md`). Mismatch energetics: Allawi HT & SantaLucia J
Jr 1997 Biochemistry 36:10581, doi:10.1021/bi962590c; Peyret N et al. 1999 Biochemistry 38:3468,
doi:10.1021/bi9825091; SantaLucia J Jr & Hicks D 2004 Annu Rev Biophys Biomol Struct 33:415,
doi:10.1146/annurev.biophys.32.110601.141800.

So every transcript site within a few substitutions must be found, not just the best hits:
- **Bowtie 1** (`aligner: bowtie`, `bowtie -v N -a`): reports every end-to-end alignment with <= N
  substitutions (N = 0-3), i.e. an exhaustive search within that radius (Langmead B et al. 2009 Genome
  Biol 10:R25, doi:10.1186/gb-2009-10-3-r25).
- **Bowtie2** (`aligner: bowtie2`, `--very-sensitive-local -k 100`): a seed-and-extend heuristic with
  local alignment (Langmead B & Salzberg SL 2012 Nat Methods 9:357, doi:10.1038/nmeth.1923). Its
  default local minimum score 20 + 8 ln(L) exceeds the perfect-match score 2L for L < 22, so **reads
  shorter than 22 nt never align** and count 0 off-targets (verified for this skill: 21-mers 0/5 aligned,
  22-mers aligned).

U-Probe README measurements (Ensembl 113 human transcriptome):
- All 15,295 18-22-mers of the ACTB and GAPDH MANE mRNAs: Bowtie 1 (`max_mismatches: 2`) reported
  exactly the sites a brute-force scan found (sample of 1,000; 147,426 sites); Bowtie2 found 1.5% of the
  off-target genes with a site within 2 substitutions (none for 18-21-mers; 56% for 30-mers, 82% for 40-mers).
- Per window (1,500 random windows per length, 400 MANE mRNAs, `strand: sense`): Bowtie2 flagged 0% of
  20-mers, 69% of 25-mers, 88% of 28-mers, 87% of 30-mers and 97-99% of 40-52-mers that have a site
  within 2 substitutions in another gene.
- Cost: Bowtie 1 searched 15,295 sequences in 4.6 s (8 threads); its index took ~5 min and 650 MB.

Rule: windows <= 30 nt -> `aligner: bowtie`, `max_mismatches: 2` (the templates' choice). Windows
>= 40 nt -> Bowtie2 also finds partial matches and indels but is a heuristic (it misses 1-3 % of the
windows with a <= 2-substitution site in another gene); where the unit of signal is the whole window
(HCR pairs, SNAIL and padlock sites) the exhaustive substitution search is the decisive check, so Bowtie 1 alone
is a sound default there (it also spares the Bowtie2 transcriptome index, the slowest build of an RNA run),
and Bowtie2 is an optional second attribute for indels and one-half matches. For 28-30-nt
windows, extrapolating the numbers above, a 3-substitution site can still have its Tm above the
hybridization temperature, so `max_mismatches: 3` is defensible (derived, not from a paper).

## Substitutions, indels and partial matches

- Bowtie 1 sees substitutions only. A 1-nt bulge (an indel between probe and site) costs +4.0 kcal/mol
  at 1 M NaCl (SantaLucia & Hicks 2004, Table 4), with strong sequence dependence (Tanaka F et al. 2004
  Biochemistry 43:7143, doi:10.1021/bi036188r). Primer3's parameters put a central 1-nt bulge at a
  ~9-11 °C Tm loss for a 20-mer and ~6-8 °C for a 30-mer, about two substitutions' worth (computed for
  this skill). So a 1-indel site can bind about as well as a 2-substitution site, and Bowtie 1 misses it.
  Bowtie2's local alignment reports indels and partial matches (from 40 nt on, README).
- Partial matches: a long probe can bind through a contiguous sub-stretch. For 50-mer microarray probes,
  non-targets > 75% similar over the 50 nt may cross-hybridize, and the probe must not contain a
  complementary stretch > 15 contiguous bases (Kane MD et al. 2000 Nucleic Acids Res 28:4552,
  doi:10.1093/nar/28.22.4552); He Z et al. 2005 Appl Environ Microbiol 71:3753,
  doi:10.1128/AEM.71.7.3753-3760.2005 set <= 85% identity, <= 15-nt stretch and >= −30 kcal/mol for
  50-mers (<= 20 nt, >= −40 kcal/mol for 70-mers). By the Tm table in `references/concepts/hybridization-and-tm.md`, a perfect
  20-nt duplex has Tm ~40 °C in 2x SSC + 30% formamide: a 20-nt identity inside a 30-nt MERFISH window
  is partly bound at 37 °C, and Bowtie 1's end-to-end search over the whole window cannot see it.
  For long windows (Oligopaint, 40-52-nt units) apply the stretch rule with Bowtie2 or `kmer_count`.
- Stellaris guidance excludes 20-mers with > 16 nt matching other RNAs (vendor rule, quoted in the
  smFISH template); U-Probe approximates it with <= 2 substitutions over the window.

## Strand of the search

`target_region` is the RNA sense sequence. Search it with `strand: sense`: hits are transcripts that
contain the window, i.e. RNAs the probe (`rc(target_region)`) can bind. `strand: both` (the default)
also counts transcripts antisense to the window, which the probe cannot bind (it has their sequence).
For DNA targets both strands are present; `kmer_count` counts canonical k-mers (both strands) by
default. `strand: sense` counts the k-mers as written on the FASTA's forward strand, as OligoMiner's and
PaintSHOP's dictionaries do; the oligopaint template uses it to keep OligoMiner's threshold.

## Expression-aware off-target counting

Mechanism: a probe can only bind RNAs present in the sample, and the false signal an off-target makes
scales with its abundance. 1 TPM is a small amount: cells hold ~50,000-300,000 mRNAs (Marinov GK et
al. 2014 Genome Res 24:496, doi:10.1101/gr.161034.113), so 1 TPM ~ 0.05-0.3 molecules per cell on average.

Precedent: MERFISH removed probes with > 14 nt of homology to rRNA/tRNA or > 17 nt of homology to
highly expressed genes (FPKM > 10,000) (Chen KH et al. 2015 Science 348:aaa6090,
doi:10.1126/science.aaa6090, Supplement), and later penalized off-target 17-mers by the summed RNA-seq
abundance of the isoforms containing them, keeping probes with gene/isoform specificity indices of
0.75-1 and no > 15-nt homology to rRNA/tRNA (Moffitt JR et al. 2016 PNAS 113:11046,
doi:10.1073/pnas.1612826113, SI).

U-Probe (`gene_expression:` block, applies to all `mapped_genes` attributes):
- `min_value` (absolute): off-target genes below it are not counted (listed in `<attr>_unexpressed`).
  TPM >= 1 is the template's choice; it is a convention, not a measured cutoff.
- `min_fraction_of_target` (relative): an off-target gene counts only above max(`min_value`, f x
  target's value). Rationale: an RNA at a few percent of the target's level adds at most a few percent
  false spots. README: 0.05 raised human GAPDH from 7 to 14 probes, ACTB 38 -> 44, mouse Eef1a1 17 -> 48.
- Default configuration. With the sample's expression table, start with all of: `level: transcript` (or
  gene), `min_value: 1`, `min_fraction_of_target: 0.05` and `exempt_copies` on the off-target attribute. They are
  the specificity model, not relaxations: without them the pool shrinks for no imaging reason and leaves less
  room for the method's other windows. An absolute 1 TPM floor alone (`absent: count`) is a stricter, not a
  "conservative", model: a negligible RNA next to a highly expressed target still removes every probe on it.
- `absent: count` (default) counts genes missing from the table as expressed; `ignore` trusts the table
  to list every expressed gene (right only for a table on the same annotation).
- Transcript vs gene level. Only molecules that carry the site can make false signal. A gene-level
  table values an off-target gene by all its isoforms, so a site on a rarely used isoform of a highly
  expressed gene still counts (conservative). With `level: transcript` (salmon/kallisto/RSEM transcript
  tables), U-Probe sums only the transcripts the sequence hits; prefer it when you have one. README
  (smFISH template, salmon tables of HeLa-S3 and NIH/3T3 from one quantification, TPM >= 1, fraction
  0.05): the transcript table counted 61% fewer off-target genes for six housekeeping genes and
  10-11% fewer for random genes; mouse Gapdh went from 26 to 48 probes (its windows hit two lncRNA genes
  mostly on isoforms the cells hardly express); most other genes changed little.
- Isoforms of the target gene are never off-targets at either level (see
  `references/concepts/targets-and-isoforms.md`).
- The table must come from the same species, a matching annotation and the sample's cell types. Bulk
  tables under-represent rare cell types: an off-target expressed only in a minority population can
  be "below threshold" in bulk; prefer the maximum over cell-type pseudobulks when the tissue is mixed.
- README (smFISH template, HeLa TPM >= 1): 54% instead of 39% of 20-mers of 30 random genes had no
  off-target gene; six housekeeping genes kept 180 instead of 119 probes (193 with fraction 0.05).

## Near-identical copies: unavoidable

Processed pseudogenes and retrocopies are reverse-transcribed copies of mRNAs; human has thousands:
2,090 processed ribosomal-protein pseudogenes (Zhang Z, Harrison PM, Gerstein M 2002 Genome Res
12:1466, doi:10.1101/gr.331902); ~8,000 processed pseudogenes (Zhang Z et al. 2003 Genome Res 13:2541,
doi:10.1101/gr.1429003); of 11,216 GENCODE pseudogenes, 8,248 processed, and at least 9% (876)
transcribed (Pei B et al. 2012 Genome Biol 13:R51, doi:10.1186/gb-2012-13-9-r51). Readthrough genes
and duplicated loci also contain the target.

Mechanism: if a gene contains most of the target's sequence, every window of the target also lies in
it; no probe choice can avoid it, and counting it removes all probes. Decide instead:
- `exempt_copies` (k = 30, `min_shared` 0.5): genes with a transcript holding >= half of the target's
  30-mers are listed in `<attr>_copies` and not counted. It exempts only such near-identical copies (every
  probe of the set hits them anyway); it does not exempt paralogs or genes that share some windows, and
  it does not say that a copy is silent. Do not drop it "to tell the target from its paralogs": the paralogs
  stay counted with it. README shares: EEF1A1P5 holds 81% of the EEF1A1
  merged CDS 30-mers; mouse Gapdh retrocopy Gm13394 76% (668 TPM in NIH/3T3); paralog ACTG1 only 5-7% of
  ACTB (stays counted); SPDYE8-15 hold 84-89% of SPDYE17 (the probes detect the family).
- Then report: if a copy is expressed, the probe set measures target + copy. Say so to the user.
- Partial escape: old copies diverge; ligation-based single-probe methods can discriminate a single
  base at the ligation junction (Larsson C et al. 2004 Nat Methods 1:227, doi:10.1038/nmeth723), so a
  padlock junction placed on a diagnostic base can separate near-identical sequences. Many-probe sets
  cannot.

## Paralogs: avoidable

Paralogs share stretches, not the whole sequence (ACTG1: 4.8% of ACTB mRNA 30-mers). Windows in the
divergent parts are specific, and UTRs usually diverge more than CDS between paralogs (actin 3' UTRs
"are isotype specific but are conserved in evolution" and served as isotype-specific probes: Ponte P,
Gunning P, Blau H, Kedes L 1983 Mol Cell Biol 3:1783, doi:10.1128/mcb.3.10.1783-1791.1983; U-Probe
README: paralogs and pseudogenes share coding sequence more than UTRs). So: `source: mRNA`, Bowtie 1
within 2 substitutions, `target_mapped_genes <= 1`; never exempt paralogs.

## Pool rule vs per-probe rule

| Method class | Unit that makes a spot | Rule |
| --- | --- | --- |
| Many singly labelled oligos (smFISH, smiFISH), MERFISH/seqFISH primary sets | many probes on one molecule | pool rule allowed |
| Padlock / RCA (ISS, HybISS, STARmap, MiP-seq) | one ligated padlock -> one RCA product | strict per probe (unit) |
| HCR v3 split-initiator pairs | one pair -> one amplification polymer | strict per pair |
| Oligopaint DNA FISH | many probes in a genomic interval | pool rule; avoid clustered off-target hits |

- Pool rule (mechanism): a spot is called when many probes of the set sit on one molecule; an off-target
  RNA hit by k of N probes carries k/N of a spot's fluorescence. Raj A et al. 2008 Nat Methods 5:877,
  doi:10.1038/nmeth.1253 (Supp. Fig. 3a): with the first 12 or 24 of 48 probes, particles were dimmer and
  "not readily identifiable"; heavily labelled single probes instead "suffer from false positives ...
  owing to individual probe misbinding". smFISH protocols ask for > 20 probes (Lyubimova A et al. 2013
  Nat Protoc 8:1743, doi:10.1038/nprot.2013.109); smiFISH separated true from false spots with 24-45
  probes (Tsanov N et al. 2016 Nucleic Acids Res 44:e165, doi:10.1093/nar/gkw784). So a few probes on
  one off-target RNA stay below detection. Hence a set may tolerate a few probes per off-target gene:
  `off_target_pool: {max_probes_per_gene: 2}` with a per-probe cap such as `<= 2`. README: this added
  only 9-18 probes per gene set while 34-47% of final probes came to carry an off-target, because
  `equal_space` picks by position; an expression table helped more. Default to `<= 1`.
- Per-probe rule, padlock/RCA: a single circularized padlock is amplified into one bright, countable
  product (Larsson et al. 2004; Ke R et al. 2013 Nat Methods 10:857, doi:10.1038/nmeth.2563), so one
  off-target ligation is one false spot. Keep `target_mapped_genes <= 1` (with expression) on the whole
  site both arms (and the primer) bind; no pool. A second `mapped_genes` attribute on part of the site
  catches what the whole-window search cannot: the padlock templates give it a `layout` part and filter
  it `<= 1` too (`padlock_iss`: the 20 nt around the ligation junction, within 1 substitution, after
  the HybISS design code's BLAST rule; `double_hyb_rca`: the 26-nt arm site, within 2, because Bowtie2
  searches the 40-nt window only locally).
- Per-pair rule, HCR v3: the initiator is split between two probes and amplification needs both bound
  next to each other (Choi HMT et al. 2018 Development 145:dev165753, doi:10.1242/dev.165753). A lone
  probe is suppressed, but a pair that binds an off-target site grows an amplified polymer. Search the
  whole pair window (e.g. 52 nt) and keep it free of off-targets.
- DNA FISH: scattered single off-target probes give diffuse background; many probes landing in one
  off-target region form a false focus. `avoid_otp` (DNA mode, `mapped_sites`) removes probes whose
  off-target genome hits cluster.

## Repeats

Interspersed repeats make up ~45% of the human genome; Alu alone has > 1 million copies (~10%) (Lander
ES et al. 2001 Nature 409:860, doi:10.1038/35057062), and repeats also occur in UTRs and introns.
- Repeats within the target RNA itself (tandem copies) are not off-targets but make one oligo bind several
  sites: count distinct oligos and use all sites for spacing (`references/concepts/targets-and-isoforms.md`).
- RNA designs: repeat windows hit hundreds of genes. Bowtie2 stops at 100 alignments (counts saturate);
  Bowtie 1 lists them all (slow; `_off_targets` cells beyond Excel's 32,767 characters are cut in the
  `.xlsx`, full in the TSV). A per-probe cap (`<= 1` or
  `<= 2`) removes them; `off_target_pool` alone does not (each gene is hit by one probe only).
- DNA designs: use `kmer_count` (Jellyfish genome k-mer counts) with a low maximum, and
  `sequence_pattern` for simple repeats. OligoMiner's `kmerFilter` drops a probe if any of its 18-mers
  occurs >= k times in the genome; the paper ran k = 5 (Beliveau BJ et al. 2018 PNAS 115:E2183,
  doi:10.1073/pnas.1714530115). U-Probe's `kmer_count` reports the highest 18-mer count, so the
  equivalent filter is `<= 4`. Clone-based probes needed preannealing with competitor genomic DNA to
  suppress repeats (Lichter P et al. 1988 Hum Genet 80:224, doi:10.1007/BF01790090); oligo libraries
  avoid repeats by design.

## Recipe by method

| Method | Attribute on | Aligner, radius | Filter | Expression | Copies | Pool |
| --- | --- | --- | --- | --- | --- | --- |
| smFISH 18-22 nt | `target_region`, `sense` | bowtie, 2 | `<= 1` | recommended (+ fraction 0.05) | `exempt_copies` | optional |
| MERFISH / seqFISH+ 28-30 nt | `target_region`, `sense` | bowtie, 2 (3 defensible) | `<= 1` | recommended | `exempt_copies` | optional |
| HCR v3 pair 52 nt | whole pair window, `sense` | Bowtie 1 alone, 2 (exhaustive within the radius; Bowtie2 only as an optional second attribute for indels and partial matches) | `<= 1` | recommended | report copies | no |
| Padlock / SNAIL 30-46 nt | whole site, `sense`; plus the junction (`padlock_iss`) or the arm site (`double_hyb_rca`) | bowtie (<= 30 nt) / Bowtie2 | `<= 1` on each | recommended | report copies | no |
| Oligopaint 35-41 nt | homology, windows from a soft-masked genome (`masked: drop`) | `kmer_count` (k 18, `strand: sense`) + `mapped_sites` | `kmer_count <= 4`, `<attr>_num == 1` | n/a | n/a | `avoid_otp` |

"Report copies": with single-unit methods every probe of a copied target also detects an expressed
copy, so each copy molecule makes a real false spot. Exempting it is the only way to get probes at all;
do it deliberately, name the copies (`<attr>_copies`) and their expression in the report, or place a
padlock junction on a diagnostic base. The technique files in `references/techniques/` take precedence
for method-specific numbers.
