# Oligopaint / OligoMiner DNA FISH

> Read this when: designing oligo DNA FISH probes for a genomic interval (a locus, a chromosome region, a whole chromosome or chromatin-tracing segments), adapting the `oligopaint` template, or judging the specificity of genomic-homology oligos. Also covers repeat targets (centromeres, telomeres) briefly. General specificity rules (repeats, `avoid_otp`): `references/concepts/specificity.md`.

## 1. What it is and when to choose it

Oligopaints are pools of synthetic oligos, each a short genomic homology region (~30-50 nt) flanked by non-genomic sequence (PCR primers, readout sites). The pool is array-synthesized, amplified and converted to single-stranded probes; the genome is denatured and thousands of probes tile the target (Beliveau 2012, 2018). Like smFISH, a spot is the sum of many independent oligos; there is no ligation and no pairing between probes.

| Use | Typical target | Source |
| --- | --- | --- |
| Locus FISH | 10 kb to a few Mb, 5-20 probes/kb | Beliveau 2012; Nguyen & Joyce 2019 |
| Chromosome paint | multi-Mb to whole chromosome, ~0.5-2 probes/kb | Nguyen & Joyce 2019 |
| Chromatin tracing | segments of 2-100 kb, each with its own readout sequence | Wang 2016; Bintu 2018; Mateo 2019 |
| Superresolution, homolog-specific (SNP) imaging | 20 kb to Mb | Beliveau 2015, 2018 |

Choose something else for RNA (`smfish.md`, `merfish.md` and the other RNA files in this folder; see `README.md`) and for tandem repeats, where one or a few short oligos or PNA probes suffice (section 4).

## 2. Probe architecture

There is no ligation junction and no arm order: each oligo hybridizes on its own, and neighbours need not be adjacent.

```text
Library oligo (U-Probe template layout), 5'->3':
5'-[primer_1, 21][homology, 35-41][primer_2 site = rc(reverse primer), 21]-3'    77-83 nt
                 |<- copied 5'->3' from the + strand of chrom:start-end ->|
Beliveau 2012 format: 21 + 32 + 21 = 74 nt.

Final probe, one route (Beliveau 2012): PCR with a dye-labelled primer_1 and a reverse
primer carrying a nicking-endonuclease site; the nick lies 3' of the homology on the
labelled strand; a denaturing gel keeps the labelled 53-mer:
5'-Dye-[primer_1, 21][homology, 32]-3'                                           53 nt

After denaturation (heat + 50% formamide) both genomic strands are single-stranded:
probe (+ strand copy)  5'-Dye-[primer_1]-ACCTGATTC...GGTCA-3'
                                         |||||||||...|||||
genomic - strand                      3'-TGGACTAAG...CCAGT-5'   bound
genomic + strand                      5'-ACCTGATTC...GGTCA-3'   same sequence, not bound
```

Indirect labelling: the non-genomic 5' tail ("MainStreet") stays single-stranded on the hybridized probe and carries binding sites for labelled secondary/readout oligos.

```text
Beliveau 2015:  5'-[secondary site(s), 32 each][primer][homology]-3'  (sites added 5' of the primer by PCR)
Wang 2016:      MainStreet holds one readout sequence shared by the 1000 probes of a TAD
Bintu 2018:     20-nt readout unique to each 30-kb segment, plus primers
ORCA (Mateo 2019): 120-nt template = 40-nt homology + 20-nt segment barcode
                   + 20-nt fiducial/RT site + forward and reverse index primers
```

Which strand survives depends on the ssDNA route (labelled strand after nicking, RT product after in vitro transcription, or the strand lambda exonuclease does not digest); check it when strand matters (section 4).

## 3. Design parameters

Homology selection. Tm values compare only under the same conditions.

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Homology length, Oligopaints 2012 | 32 | nt | Beliveau 2012, doi:10.1073/pnas.1213818110 |
| OligoMiner Coverage / Balance / Stringent | 26-32 / 36-41 / 40-46 | nt | Beliveau 2018 Fig. 2A, doi:10.1073/pnas.1714530115 |
| **Protocols disagree:** Balance length | 36-41 (Fig. 2A; Table S1, the imaged sets; blockParse default `-l 36`) vs 35-41 (main text; the genome-wide hg38 Balance set OligoMiner hosts, `hg38b.tar.gz`, whose probes are 35-41 nt with 35 the most common length; PaintSHOP Table 1). The template follows the hosted set | nt | Beliveau 2018; Hershberg 2021, doi:10.1038/s41592-021-01187-3 |
| PaintSHOP newBalance | 30-37, Tm 42-47 C, may include repeat-masked bases | nt | Hershberg 2021 |
| ORCA / chromatin tracing | 40 | nt | Mateo 2019, doi:10.1038/s41586-019-1035-4; Mateo 2021, doi:10.1038/s41596-020-00478-x |
| Tm window, Coverage / Balance / Stringent | 37-42 / 42-47 / 47-52 (390 mM Na+, 50% formamide, 25 nM each strand, Biopython DNA_NN3, 0.65 C per % formamide; U-Probe `tm_method: oligominer`) | C | Beliveau 2018 Fig. 2A; blockParse.py, probeTm.py |
| Hybridization temperature per preset | 37 / 42 / 47 (the Tm window is [T_hyb, T_hyb + 5]) | C | Beliveau 2018 Fig. 2A |
| Tm, other designs | 75-90 (OligoArray, 32-mers; 34.2-49.7 when recomputed for PaintSHOP); >65 (ORCA, conditions not stated) | C | Beliveau 2012; Hershberg 2021 Table 1; Mateo 2021 |
| GC | 20-80 (blockParse default, ORCA); 35-80 (2012) | % | Beliveau 2012, 2018; Mateo 2021 |
| Prohibited runs | AAAAA, TTTTT, CCCCC, GGGGG and any N (blockParse); 2012: >=5 A/T or >=4 G/C | - | blockParse.py; Beliveau 2012 |
| Spacing between probes | 0 (blockParse -S); >=10 (2012) | nt | blockParse.py; Beliveau 2012 |
| Repeat masking | candidates are mined from the repeat-masked assembly (UCSC, RepeatMasker; blockParse rejects every window with an N); the Bowtie2 index and Jellyfish dictionary are built from the unmasked one | - | Beliveau 2018 Methods; OligoMiner README; blockParse.py |
| k-mer filter | 18-mers (Balance, Stringent), 16-mers (Coverage); drop a probe if any k-mer occurs >= 5 times. The dictionary counts each k-mer on the forward strand only (`jellyfish count -s 3300M -m 18 --out-counter-len 1 -L 2`, no `-C`); PaintSHOP builds it the same way | - | Beliveau 2018 Methods; kmerFilter.py; PaintSHOP_pipeline Snakefile |
| Alignment filter | Bowtie2 `--very-sensitive-local -k 2` (Unique Mode, UM; the README uses `-k 100`) or `--local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -k 2` (LDA Mode, LDM). UM keeps probes that aligned with no `XS:i` tag, i.e. without a second alignment; LDM (models at 32-57 C, 5-C steps; class "binds" = duplex probability >= 0.2; reject at predicted probability >= 0.5) | - | Beliveau 2018 Methods; outputClean.py |
| Secondary structure (optional) | NUPACK structureCheck (dna1998, at -T + 0.65 C per % formamide); keep probes whose probability of being linear is above -t (default 0.1) at -T (default 47 C). Not applied to any Beliveau 2018 set (Table S1) | - | structureCheck.py |
| Repeat homology (ORCA) | <14 nt to RepBase repeats; <12 nt to any barcode | nt | Mateo 2019, 2021 |

Density, library and wet lab.

| Parameter | Value | Unit | Source |
| --- | --- | --- | --- |
| Density, OligoArray sets | ~10 (9.7-12.8) | probes/kb | Beliveau 2012 |
| Density, OligoMiner hg38 after all filters, Coverage / Balance / Stringent | 1.6-3.1 / 4.6-4.8 / 2.6-2.9 | probes/kb | Beliveau 2018 |
| Density of imaged sets | 3.2-6.1 (Xq28, 19p13.2; 20 kb with 104 oligos) | probes/kb | Beliveau 2018 Table S1 |
| Typical locus / whole chromosome | 5-20 / ~0.5-2 | probes/kb | Nguyen & Joyce 2019, doi:10.1007/978-1-4939-9674-2_11 |
| Minimum set | >200 probes over >10 kb (10 kb was imaged with 200 oligos in 2012); ORCA >=20 probes per 2-10-kb segment | probes | Hershberg 2021; Beliveau 2012; Mateo 2019 |
| Chromatin tracing segments | 1000 probes per 100-kb TAD; ~300 per 30-kb segment | probes | Wang 2016, doi:10.1126/science.aaf8084; Bintu 2018, doi:10.1126/science.aau1783 |
| Primer flanks | 21 each side (2012); 20-22 (2018 Table S2); T7 promoter TAATACGACTCACTATAGGG prepended for IVT/RT | nt | Beliveau 2012, 2018 |
| Readout / secondary site | 32 (2015 secondaries); 20 (Bintu 2018, ORCA); 42 bridge (PaintSHOP chromosome patterns) | nt | Beliveau 2015, doi:10.1038/ncomms8147; Bintu 2018; Hershberg 2021 |
| Hybridization buffer | 2x SSCT, 50% formamide, 10% dextran sulfate (+ RNase A for DNA-only FISH) | - | Beliveau 2012, 2018 |
| Probe amount | 1-20 pmol per 25 uL (2012); 1.6-2.5 uM (2018); 50-200 pmol, signal lost above ~600 pmol (Nguyen & Joyce) | - | as cited |
| **Protocols disagree:** denaturation | 92 C 2.5 min (2012; Nguyen & Joyce, or 80 C 5 min); 78 C 2.5-3 min (2015 coverslips; 2018 cells); metaphase spreads 70% formamide 70 C 90 s (2018); 0.1 M HCl 5 min then 90 C 10 min (ORCA) | - | as cited |
| **Protocols disagree:** hybridization | overnight 37 or 42 C (2012); 42 or 47 C (2015); >24 h at 47-52 C (2018); 42 C overnight (ORCA) | C | as cited |
| Post-hybridization wash | 2x SSCT at 60 C (15 min in 2012; 4 x 5 min in 2018) | - | Beliveau 2012, 2018 |

## 4. Target region and template

- **Source.** A genomic interval `chrom:start-end`; there is no mRNA/CDS/isoform choice. Probes cover introns, intergenic DNA and both alleles. U-Probe reads the window from the + strand of the FASTA, so the homology is a + strand copy that binds the - strand.
- **One strand per position.** After denaturation either strand is a target, so probes from one strand suffice. Probes from both strands at the same position are complementary and anneal to each other in solution. OligoMiner's `probeRC` flips a set for strand-specific DNA or RNA FISH.
- **RNA cross-talk.** Over transcribed regions a DNA probe antisense to the RNA also labels it. Add RNase A (Beliveau 2018 buffer), or design the probe sense to the transcript (ORCA). ORCA images RNA before DNA when doing both.
- **Masked vs unmasked genome.** OligoMiner mines repeat-masked FASTA (UCSC, RepeatMasker) and builds its Bowtie2 and Jellyfish indexes from the unmasked assembly (Beliveau 2018 Methods; README). U-Probe does both with one soft-masked FASTA, repeats in lowercase (`uprobe genome fetch --soft-masked`, Ensembl `dna_sm`; or UCSC `hg38.fa`): Bowtie2 and Jellyfish ignore case, and `extracts.target_region.masked: drop` skips every window with a lowercase or N base. The oligopaint template sets it. `masked: keep`, the default, upper-cases the windows and ignores masking. With an unmasked FASTA (Ensembl `dna`, a lab's own assembly) give the repeats as a file instead of making a masked copy of the genome: `mask: <file>` next to `masked: drop`, a BED of RepeatMasker intervals of the same assembly (UCSC `rmsk`; 0-based, may be `.gz`) or a FASTA of soft-masked loci named by region (`>chr19:10000001-10050000`); its chromosome names must match the targets' (`validate` warns otherwise). Without either, `drop` skips only N windows; the run log and `uprobe validate --check-targets` warn when a target has no masked base. Do not use a hard-masked FASTA (`dna_rm`) as the genome: the indexes would no longer see the repeats. PaintSHOP's newBalance sets keep soft-masked windows and flag them instead (Hershberg 2021).
- **Locus vs whole chromosome.** Locus: dense tiling (Balance or Coverage), every passing window, and at least ~200 probes. Chromosome paint: subsample to ~0.5-2 probes/kb (`equal_space`). Use universal primers on all oligos plus nested region-specific primers so one pool yields sub-paints (Nguyen & Joyce 2019). Chromatin tracing: split the region into equal segments and give each segment its own readout and probe count (Wang 2016; Bintu 2018; Mateo 2019). Repeats and assembly gaps leave holes, so check density per bin, not only the mean.
- **Coordinates.** Regions are 1-based and inclusive, as genome browsers and Ensembl write them: `chr19:10000001-10050000` is the 50,000 bases from 10,000,001 to 10,050,000. (Older U-Probe versions read regions 0-based, one base off.) Chromosome names must match the FASTA: Ensembl FASTAs from `uprobe genome fetch` use `19`, UCSC FASTAs `chr19`.

Repeat targets need a different design: thousands of tandem copies give the signal, so one or a few short oligos (or a PNA) against the repeat consensus are enough. Specificity is the choice of repeat family or chromosome-specific variant, not genomic uniqueness.

| Target | Repeat | Probe | Source |
| --- | --- | --- | --- |
| Telomeres (vertebrates) | (TTAGGG)n | labelled (CCCTAA)3 PNA, quantitative FISH | Moyzis 1988, doi:10.1073/pnas.85.18.6622; Lansdorp 1996, doi:10.1093/hmg/5.5.685 |
| Human centromeres | alpha satellite, 170-bp units; 17-bp CENP-B box (absent from Y) | alphoid oligos: pan-centromeric or chromosome-specific (2, 12); 30-min hybridizations | Masumoto 1989, doi:10.1083/jcb.109.5.1963; Matera & Ward 1992, doi:10.1093/hmg/1.7.535 |
| Mouse centromeres | major satellite (pericentric, ~6 Mb of 234-bp units); minor satellite (centric, ~600 kb of 120-bp units) | satellite probes | Guenatri 2004, doi:10.1083/jcb.200403109 |
| Drosophila X | 359-bp satellite | a single synthetic oligo carrying a 32-nt MainStreet | Beliveau 2015 |

In U-Probe a repeat consensus can be given as a direct-sequence target (`targets: [{majsat: <consensus>}]`). Drop the k-mer filter for it, because it rejects every repeat oligo by design. No template covers this case.

## 5. Specificity

- **Signal logic.** A spot needs many probes on one locus. An off-target site bound by one oligo adds one dye: diffuse background, not a spot. False foci come from: (1) **repeats**, where one oligo has hundreds to millions of sites (Alu, L1, satellites); (2) **off-target peaks**, where many different oligos of the set hit one other region (segmental duplications, paralog clusters, retrocopy loci); (3) oligos that pair with each other or with the readouts.
- **k-mer counting** (Jellyfish) catches high-copy repeats, including matches too short for an aligner seed. OligoMiner uses 18-mers for its 42 and 47 C presets and 16-mers for the 37 C preset, chosen from simulated k-mer duplexing (Beliveau 2018 Fig. S6); PaintSHOP's "Max K-mer Count" uses 18-mers. It gives no positions, so it cannot see off-target peaks. Copies at 2-4 sites pass a `<= 4` threshold.
- **Strand of the k-mer count.** A probe binds its k-mers on either strand of the denatured genome. U-Probe's `kmer_count` counts canonical k-mers (both strands) by default; OligoMiner and PaintSHOP count the probe's k-mers on the forward strand only, so a repeat present in both orientations counts about half as often, and their `-k 5` threshold refers to these counts. `strand: sense` reproduces them (index `jf_genome/<stem>.k<k>.stranded.jf`); the template uses it. No canonical cut-off is equivalent: a canonical count is the forward count of the k-mer plus that of its reverse complement.
- **Alignment** gives positions. Bowtie2 local alignment flags partial matches (OligoMiner Unique Mode keeps probes that align once); LDA Mode keeps probes whose second-best hit is not thermodynamically relevant at the hybridization temperature. PaintSHOP replaces the classifier with an XGBoost duplex-probability score summed over off-target sites.
- **`mapped_genes` is the wrong measure for DNA FISH.** It searches spliced transcripts and counts genes, so it misses introns, intergenic DNA and most repeats. Use `kmer_count` (genome k-mers) and `mapped_sites` (Bowtie2 `--very-sensitive-local -k 100` against the genome, both strands; `<attr>_num` counts the alignments). `<attr>_num == 1` is OligoMiner's Unique Mode: Bowtie2 found the probe's own site and no second one (outputClean keeps a probe without an `XS:i` tag). It is local, not end-to-end, alignment, so partial matches count.
- **`avoid_otp`** (DNA-only `post_process` step; needs a `mapped_sites` attribute) visits probes in order and adds their alignments to an interval tree. It drops a probe when one of its off-target hits lands where the hits already kept reach `density_thresh` per bp within +-1 Mb (`search_range`). It also drops on-target overlaps (`avoid_target_overlap`). It writes `on_target`, `off_target` and `target_ratio`. The thresholds have no published calibration. `tests/data/DNA_format.yaml` uses `1e-5`, i.e. 20 hits per 2 Mb; the code default is `1e-4`.

## 6. U-Probe template

```bash
uprobe templates show oligopaint -o protocol.yaml
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets   # also warns if the FASTA is not soft-masked
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --threads 10
```

Key lines (OligoMiner Balance, as in its genome-wide hg38 set):

```yaml
targets: [chr19:10000001-10050000]          # chrom:start-end, 1-based, in your FASTA's naming
extracts: {target_region: {source: genome, length: [35, 41], step: 1,
                           masked: drop}}   # skip soft-masked (lowercase) or N windows; or add mask: rmsk.bed
probes:
  oligopaint_oligo:
    template: "{primer_1}{homology}{primer_2_rc}"   # primers: 21-N placeholders to replace
    parts: {homology: {expr: "target_region"}}
attributes:
  homology_tm: {target: oligopaint_oligo.homology, type: annealing_temperature, tm_method: oligominer,
                hybridization: {na_mM: 390, formamide_pct: 50, probe_nM: 25}}   # OligoMiner's formula
  homology_kmer_count: {target: oligopaint_oligo.homology, type: kmer_count, kmer_len: 18, aligner: jellyfish,
                        strand: sense,        # forward-strand counts, as OligoMiner's dictionary
                        locus_only: true}     # only the candidates' k-mers, no genome-wide table
  homology_sites: {target: oligopaint_oligo.homology, type: mapped_sites, aligner: bowtie2}
post_process:
  filters:
    homology_gc: {condition: homology_gc >= 0.2 & homology_gc <= 0.8}
    homology_tm: {condition: homology_tm >= 42 & homology_tm <= 47}
    homopolymers: {type: sequence_pattern, target: oligopaint_oligo.homology,
                   exclude_patterns: [AAAAA, TTTTT, CCCCC, GGGGG]}
    n_bases: {type: sequence_pattern, target: oligopaint_oligo.homology, exclude_patterns: [N]}
    homology_kmer_count: {condition: homology_kmer_count <= 4}
    homology_sites_num: {condition: homology_sites_num == 1}    # outputClean Unique Mode
  remove_overlap: {location_interval: 0}
```

The genome entry needs a soft-masked FASTA (or a `mask`) and `bowtie2` in `align_index` (a genome-wide Bowtie2 index; `uprobe build-index`). `kmer_count` with `locus_only: true` needs Jellyfish 2 on `PATH`, but no `jellyfish: true` and no genome-wide table: it lists the candidates' 18-mers and counts them in one pass over the genome FASTA (`jellyfish count --if`), with the counts of the whole-genome table. On the 50-kb template locus the final table was identical, row for row and count for count, to one made with a table of the locus' k-mers; the pass added 37-47 s per run (human, 8 threads), and the temporary table is a few MB. A repeat `mask` works as well as soft-masking: on that locus an uppercase chr19 with a BED of UCSC `rmsk` intervals kept 150,896 windows and 447 probes (GC, Tm and homopolymer filters only), the soft-masked chr19 150,392 and 446 (it masks 38 bases of the locus that `rmsk` does not), the uppercase chr19 without a mask 349,741 and 983.

Measured on chr19:10,000,001-10,200,000 (200 kb of 19p13.2; hg38 primary assembly with UCSC soft-masking; RepeatMasker covers 54% of it). Repeat fraction = probes overlapping a UCSC RepeatMasker element.

| Set | Probes | Per kb | In repeats | Largest gap |
| --- | --- | --- | --- | --- |
| Template before (36-41 nt, unmasked windows, canonical 18-mers, no alignment filter) | 1,881 | 9.4 | 29% | 1.8 kb |
| Template with Primer3 Tm (U-Probe 0.3.0) | 1,486 | 7.4 | 0% | 5.3 kb |
| Current template (`tm_method: oligominer`) | 1,479 | 7.4 | 0% | 5.3 kb |
| OligoMiner hosted hg38 Balance set (`hg38b`) | 1,313 | 6.6 | 0% | 5.3 kb |

Masking makes the difference (1,881 to 1,388 probes); forward-strand k-mers add 33, the unique-alignment filter removes 10, and 35-mers add 75. OligoMiner's own Tm formula hardly changes the count (1,486 to 1,479; 432 to 432 on the 50-kb template locus) but changes the set: under OligoMiner's formula 473 of the 1,486 probes selected with Primer3 (32%) were above 47 C, and only 232 windows are shared. With it, every hosted OligoMiner probe overlaps a template probe (97.3% before), and 89.7% of the template probes overlap a hosted one (87.0%); both sets are mostly 35-mers.

What to edit:
- `genome`, `targets` and `encoding` keys. Replace both primer placeholders: `primer_1` is the forward PCR primer as ordered, `primer_2_rc` the reverse complement of the reverse primer, each 20-22 nt, absent from the genome and unpaired with each other and the homologies, one pair per pool or sub-pool, e.g. a published set (PaintSHOP; ORCA-public; Beliveau 2018 Table S2). The T7 promoter (for IVT/RT), the dye or a nicking site go on the ordered primers; add readout sites as parts if you image indirectly.
- Presets: Coverage `length: [26, 32]`, Tm 37-42, `kmer_len: 16`; Stringent `[40, 46]`, Tm 47-52, `kmer_len: 18`. Set `hybridization.temperature_c` to the preset's temperature.
- Optional additions, all supported by the code:

```yaml
post_process:
  avoid_otp: {chr19:10000001-10050000: {target_regions: [chr19:10000001-10050000], density_thresh: 0.00001}}
  equal_space: {number_desired: 250}   # cap the count, e.g. for paints; runs after remove_overlap
```

**U-Probe departs** from OligoMiner in these ways:
1. Tm: no departure. `tm_method: oligominer` is OligoMiner's own formula (Biopython `Tm_NN`, DNA_NN3, salt correction 5, strand term `dnac1 - dnac2/2`, -0.65 C per % formamide; it reproduces the Tm column of OligoMiner's example probes to 0.01 C). Keep it for every preset window (37-42, 42-47, 47-52 C): Primer3 reads about 1.5 C lower for these windows (SantaLucia salt correction, Blake-Delcourt formamide model), and a Primer3 42-47 C filter left 23-32 % of its probes above 47 C for OligoMiner (random windows; the locus above).
2. Only outputClean's Unique Mode is reproduced, not its LDA Mode, which the 19p13.2 sets of Table S1 used: Unique Mode drops a probe with any second alignment, however weak. There is no structureCheck, which no Beliveau 2018 set used either.
3. OligoMiner's dictionary records only k-mers seen at least twice and caps counts at 255 (`-L 2 --out-counter-len 1`); neither changes the drop-at->=5 decision.
4. U-Probe filters all windows first and then tiles greedily, shortest window first. OligoMiner tiles first (blockParse) and filters afterwards, which leaves gaps. U-Probe therefore returns denser sets for the same parameters (1,486 vs 1,313 probes in the locus above).
5. The 21-nt primer flanks come from the 2012 format, which paired them with a 32-nt homology. `temperature_c: 42` only affects `fold_score`, which the template does not compute.

## 7. Pitfalls

- **Locus-sized targets do not need genome-scale k-mer tables.** Without `locus_only`, `kmer_count` builds a whole-genome Jellyfish table at its first use (GBs for a mammalian genome). With `locus_only: true`, as in the template, it counts only the candidates' k-mers, in one pass over the genome FASTA, with the same counts; check that a known repeat of the locus gets a high count. A genome-wide table pays off only when many runs or a whole-chromosome paint reuse it (`locus_only: false`, `jellyfish: true`). `mapped_sites` still needs the whole-genome Bowtie2 index: share one built index through `indexes:` in `genomes.yaml` rather than building it again. Never write a masked copy of the whole genome: take a soft-masked FASTA as published (UCSC `hg38.fa`, Ensembl `dna_sm`), or give an unmasked one a repeat `mask` (section 4).
- **Yield before `masked: drop`.** In a repeat-dense locus `masked: drop` can skip most windows (58 % in the locus above) and leave fewer non-overlapping probes than the request. Count the windows it keeps (the unmasked fraction of the locus) against the probes needed before you commit. OligoMiner mines the masked assembly only (`blockParse` upper-cases the sequence and rejects windows with an N, so its masking is by Ns); PaintSHOP's newBalance "can include repeat-masked bases" and lets the user choose whether to allow them, with its k-mer (18-mer count) and Bowtie2 alignment filters doing the work (Hershberg 2021). If the masked set is too sparse, `masked: keep` with the k-mer and alignment filters is that published alternative; report which you used.
- **Jellyfish indexes** (without `locus_only`). Each `kmer_len` has its own index, built by the first `kmer_count` that needs it: `jf_genome/<fasta stem>.k<k>.jf` (canonical, `strand: both`) or `<stem>.k<k>.stranded.jf` (`strand: sense` or `antisense`), next to the FASTA or in the entry's `index_dir`. An index left by older versions, `jf_genome/<stem>.jf`, is reused only for the k in its header. Either way, check that a known repeat gets a high count: the 40-mer of an AluJo at chr19:10,000,001-10,000,040 has a maximal 18-mer count of 706 on both strands and 357 on the forward strand.
- Every start position yields one window per length (7 per bp for 35-41 nt at `step: 1`), and every attribute is computed on every window before the filters. `masked: drop` skips the repeat windows first (58% of them in the locus above), which also spares Bowtie2 its slowest queries. For multi-Mb paints, split the region into targets of a few Mb, or start from a precomputed OligoMiner/PaintSHOP set.
- An unmasked FASTA (Ensembl `dna`) with `masked: drop` and no `mask` skips no repeat: only the k-mer and alignment filters remove them. Use a soft-masked one, or add a RepeatMasker `mask`.
- Templates ship with N primers; ordering them as-is gives an unamplifiable library.
- Tm windows from OligoArray (75-90 C) or ORCA (>65 C) do not transfer to 50% formamide conditions; recompute the Tm under the declared conditions. OligoMiner's preset windows hold under its own formula (`tm_method: oligominer`), not under Primer3's.
- Too few probes give dim or missing spots. A 5-kb target at 5 probes/kb has only 25 probes; use Coverage settings, denser tiling or signal amplification, or enlarge the target.
- RNA labelling by antisense probes (section 4) and segmental duplications (k-mer counts of 2-4; `homology_sites_num` > 1 removes those Bowtie2 finds) are the usual sources of a second, false spot.
- `avoid_otp` without a `mapped_sites` attribute is skipped (with an error in the log); `uprobe validate` also flags it.

## 8. Verification checklist

- [ ] Targets resolve in the FASTA (`--check-targets`); chromosome names and the assembly match the coordinates.
- [ ] Every homology: length within the preset, GC 20-80%, Tm in the window under OligoMiner's formula and conditions (`tm_method: oligominer`, 390 mM, 25 nM, 50 %), no homopolymer runs, no N.
- [ ] The genome FASTA is soft-masked (or a RepeatMasker `mask` is given) and `masked: drop` is set; a target without masked bases means neither is.
- [ ] Max 18-mer count <= 4 on the forward strand (`strand: sense`), from `locus_only` counts or a table built at that k; an Alu control gives a high count.
- [ ] `homology_sites_num == 1` (one Bowtie2 alignment); no cluster of off-target hits (`avoid_otp` or a BLAT/BLAST spot check).
- [ ] Probes do not overlap. Probes per kb meet the goal (locus >= ~5, paint ~0.5-2); the largest gap is known; total >= ~200 for a >= 10-kb locus.
- [ ] The homology strand is consistent across the set; over transcribed regions the probe is sense to the RNA, or RNase A is in the protocol.
- [ ] Primers and readouts are real sequences (no N), unique per sub-pool or segment, absent from the genome, and not complementary to the homologies; T7 promoter present if needed.
- [ ] Library oligo length is within the synthesis limit (74-120 nt in the cited designs).
- [ ] Run last, with the file's own values (Balance preset; Coverage 37-42 C and Stringent 47-52 C change `--tm`; `homology_tm` is OligoMiner's Tm; no `--transcripts`):

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv --sense oligopaint_oligo.homology \
    --gc 0.2 0.8 --tm 42 47 --tm-column homology_tm --min-gap 0
```

## 9. References

- Beliveau 2012, PNAS. Versatile design and synthesis platform for visualizing genomes with Oligopaint FISH probes. doi:10.1073/pnas.1213818110
- Beliveau 2014, Curr Protoc Mol Biol. Visualizing genomes with Oligopaint FISH probes. doi:10.1002/0471142727.mb1423s105
- Beliveau 2015, Nat Commun. Single-molecule super-resolution imaging of chromosomes and in situ haplotype visualization using Oligopaint FISH probes. doi:10.1038/ncomms8147
- Beliveau 2018, PNAS. OligoMiner provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. doi:10.1073/pnas.1714530115 (code: github.com/beliveau-lab/OligoMiner, blockParse.py, outputClean.py, kmerFilter.py, structureCheck.py; hosted hg38 Balance set: yin.hms.harvard.edu/oligoMiner/probe_seqs/hg38/hg38b.tar.gz)
- Hershberg 2021, Nat Methods. PaintSHOP enables the interactive design of transcriptome- and genome-scale oligonucleotide FISH experiments. doi:10.1038/s41592-021-01187-3 (pipeline: github.com/beliveau-lab/PaintSHOP_pipeline)
- Nguyen & Joyce 2019, Methods Mol Biol. Programmable chromosome painting with Oligopaints. doi:10.1007/978-1-4939-9674-2_11
- Wang 2016, Science. Spatial organization of chromatin domains and compartments in single chromosomes. doi:10.1126/science.aaf8084
- Bintu 2018, Science. Super-resolution chromatin tracing reveals domains and cooperative interactions in single cells. doi:10.1126/science.aau1783
- Mateo 2019, Nature. Visualizing DNA folding and RNA in embryos at single-cell resolution. doi:10.1038/s41586-019-1035-4
- Mateo 2021, Nat Protoc. Tracing DNA paths and RNA profiles in cultured cells and tissues with ORCA. doi:10.1038/s41596-020-00478-x
- Boyle 2011, Chromosome Res. Fluorescence in situ hybridization with high-complexity repeat-free oligonucleotide probes generated by massively parallel synthesis. doi:10.1007/s10577-011-9245-0
- Marcais 2011, Bioinformatics. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers (Jellyfish). doi:10.1093/bioinformatics/btr011
- Blake 1996, Nucleic Acids Res. Thermodynamic effects of formamide on DNA stability. doi:10.1093/nar/24.11.2095
- Moyzis 1988, PNAS. A highly conserved repetitive DNA sequence, (TTAGGG)n, present at the telomeres of human chromosomes. doi:10.1073/pnas.85.18.6622
- Lansdorp 1996, Hum Mol Genet. Heterogeneity in telomere length of human chromosomes. doi:10.1093/hmg/5.5.685
- Masumoto 1989, J Cell Biol. A human centromere antigen (CENP-B) interacts with a short specific sequence in alphoid DNA, a human centromeric satellite. doi:10.1083/jcb.109.5.1963
- Matera 1992, Hum Mol Genet. Oligonucleotide probes for the analysis of specific repetitive DNA sequences by fluorescence in situ hybridization. doi:10.1093/hmg/1.7.535
- Guenatri 2004, J Cell Biol. Mouse centric and pericentric satellite repeats form distinct functional heterochromatin. doi:10.1083/jcb.200403109
