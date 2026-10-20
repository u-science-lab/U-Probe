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
5'-[primer_1, 21][homology, 36-41][primer_2 site = rc(reverse primer), 21]-3'    78-83 nt
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
| **Protocols disagree:** Balance length | 36-41 (Fig. 2A, Table S1, blockParse default) vs 35-41 (main text; hosted hg38 set in PaintSHOP Table 1) | nt | Beliveau 2018; Hershberg 2021, doi:10.1038/s41592-021-01187-3 |
| PaintSHOP newBalance | 30-37, Tm 42-47 C, may include repeat-masked bases | nt | Hershberg 2021 |
| ORCA / chromatin tracing | 40 | nt | Mateo 2019, doi:10.1038/s41586-019-1035-4; Mateo 2021, doi:10.1038/s41596-020-00478-x |
| Tm window, Coverage / Balance / Stringent | 37-42 / 42-47 / 47-52 (390 mM Na+, 50% formamide, 25 nM each strand, Biopython DNA_NN3, 0.65 C per % formamide) | C | Beliveau 2018 Fig. 2A; blockParse.py |
| Hybridization temperature per preset | 37 / 42 / 47 (the Tm window is [T_hyb, T_hyb + 5]) | C | Beliveau 2018 Fig. 2A |
| Tm, other designs | 75-90 (OligoArray, 32-mers; 34.2-49.7 when recomputed for PaintSHOP); >65 (ORCA, conditions not stated) | C | Beliveau 2012; Hershberg 2021 Table 1; Mateo 2021 |
| GC | 20-80 (blockParse default, ORCA); 35-80 (2012) | % | Beliveau 2012, 2018; Mateo 2021 |
| Prohibited runs | AAAAA, TTTTT, CCCCC, GGGGG and any N (blockParse); 2012: >=5 A/T or >=4 G/C | - | blockParse.py; Beliveau 2012 |
| Spacing between probes | 0 (blockParse -S); >=10 (2012) | nt | blockParse.py; Beliveau 2012 |
| k-mer filter | 18-mers (Balance, Stringent), 16-mers (Coverage); drop a probe if any k-mer occurs >= 5 times | - | Beliveau 2018 Methods; kmerFilter.py |
| Alignment filter | Bowtie2 local; Unique Mode keeps probes with one alignment; LDA Mode (models at 32-57 C, 5-C steps; class "binds" = duplex probability >= 0.2; reject at predicted probability >= 0.5) | - | Beliveau 2018; outputClean.py |
| Secondary structure (optional) | NUPACK structureCheck; keep probes whose probability of being linear is above -t (default 0.1) at -T (default 47 C) | - | structureCheck.py |
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
- **Masked vs unmasked genome.** OligoMiner mines repeat-masked FASTA (UCSC, RepeatMasker) and builds its Bowtie2 and Jellyfish indexes from the unmasked assembly (Beliveau 2018 Methods; README). U-Probe uses one FASTA for both. `uprobe genome fetch` takes the unmasked primary assembly, and `extract_fasta` upper-cases windows, so soft-masking is lost. Repeat removal therefore rests on the k-mer filter (section 5). Do not hard-mask the U-Probe FASTA: the indexes would no longer see those repeats.
- **Locus vs whole chromosome.** Locus: dense tiling (Balance or Coverage), every passing window, and at least ~200 probes. Chromosome paint: subsample to ~0.5-2 probes/kb (`equal_space`). Use universal primers on all oligos plus nested region-specific primers so one pool yields sub-paints (Nguyen & Joyce 2019). Chromatin tracing: split the region into equal segments and give each segment its own readout and probe count (Wang 2016; Bintu 2018; Mateo 2019). Repeats and assembly gaps leave holes, so check density per bin, not only the mean.
- **Coordinates.** U-Probe slices `fa[chrom][start:end]` (pyfaidx, 0-based), so `chr19:10000000-10050000` covers 1-based bases 10,000,001-10,050,000. Chromosome names must match the FASTA: Ensembl FASTAs from `uprobe genome fetch` use `19`, not `chr19`.

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
- **Alignment** gives positions. Bowtie2 local alignment flags partial matches (OligoMiner Unique Mode keeps probes that align once); LDA Mode keeps probes whose second-best hit is not thermodynamically relevant at the hybridization temperature. PaintSHOP replaces the classifier with an XGBoost duplex-probability score summed over off-target sites.
- **`mapped_genes` is the wrong measure for DNA FISH.** It searches spliced transcripts and counts genes, so it misses introns, intergenic DNA and most repeats. Use `kmer_count` (genome k-mers) and `mapped_sites` (Bowtie2 `--very-sensitive-local -k 100` against the genome; `<attr>_num` counts the alignments).
- **`avoid_otp`** (DNA-only `post_process` step; needs a `mapped_sites` attribute) visits probes in order and adds their alignments to an interval tree. It drops a probe when one of its off-target hits lands where the hits already kept reach `density_thresh` per bp within +-1 Mb (`search_range`). It also drops on-target overlaps (`avoid_target_overlap`). It writes `on_target`, `off_target` and `target_ratio`. The thresholds have no published calibration. `tests/data/DNA_format.yaml` uses `1e-5`, i.e. 20 hits per 2 Mb; the code default is `1e-4`.

## 6. U-Probe template

```bash
uprobe templates show oligopaint -o protocol.yaml
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets   # needs jellyfish: true in genomes.yaml
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --threads 10
```

Key lines (OligoMiner Balance defaults):

```yaml
targets: [chr19:10000000-10050000]          # chrom:start-end in your FASTA's naming
extracts: {target_region: {source: genome, length: [36, 41], step: 1}}
probes:
  oligopaint_oligo:
    template: "{primer_1}{homology}{primer_2_rc}"   # primers: 21-nt N placeholders
    parts: {homology: {expr: "target_region"}}
attributes:
  homology_tm: {target: oligopaint_oligo.homology, type: annealing_temperature,
                hybridization: {na_mM: 390, formamide_pct: 50, probe_nM: 25}}
  homology_kmer_count: {target: oligopaint_oligo.homology, type: kmer_count, kmer_len: 18, aligner: jellyfish}
post_process:
  filters:
    homology_gc: {condition: homology_gc >= 0.2 & homology_gc <= 0.8}
    homology_tm: {condition: homology_tm >= 42 & homology_tm <= 47}
    homopolymers: {type: sequence_pattern, target: oligopaint_oligo.homology,
                   exclude_patterns: [AAAAA, TTTTT, CCCCC, GGGGG]}
    homology_kmer_count: {condition: homology_kmer_count <= 4}
  remove_overlap: {location_interval: 0}
```

What to edit:
- `genome`, `targets` and `encoding` keys. Replace both primer placeholders with real primers that are orthogonal to the genome and to each other, e.g. a published set (PaintSHOP offers 13 sets; ORCA-public). Add the T7 promoter if you make probes by IVT/RT, and readout sites if you image indirectly.
- Presets: Coverage `length: [26, 32]`, Tm 37-42, `kmer_len: 16`; Stringent `[40, 46]`, Tm 47-52, `kmer_len: 18`. Set `hybridization.temperature_c` to the preset's temperature.
- Add `N` to `exclude_patterns`: blockParse rejects N, the template does not.
- Optional additions, all supported by the code:

```yaml
attributes:
  homology_sites: {target: oligopaint_oligo.homology, type: mapped_sites, aligner: bowtie2}  # bowtie2 in align_index
post_process:
  filters:   # add to the template's filters
    homology_sites_num: {condition: homology_sites_num <= 1}   # Unique Mode analogue
  avoid_otp: {chr19:10000000-10050000: {target_regions: [chr19:10000000-10050000], density_thresh: 0.00001}}
  equal_space: {number_desired: 250}   # cap the count, e.g. for paints; runs after remove_overlap
```

**U-Probe departs** from OligoMiner in these ways:
1. The template takes 36-41 nt, but the paper text and the hosted hg38 Balance set use 35-41.
2. Tm comes from primer3 (SantaLucia; Blake-Delcourt formamide correction, 0.61-0.72 C per %), not Biopython Tm_NN DNA_NN3 at 0.65 C per %. Oligos near 42 or 47 C can fall on the other side of the window.
3. Jellyfish counts canonical k-mers (`-C`, both strands). OligoMiner's dictionary counts forward k-mers only (`jellyfish count -m 18 ... -L 2`, no `-C`). The same threshold (drop at >= 5) is therefore stricter in U-Probe for repeats found in both orientations.
4. There is no outputClean (Unique or LDA Mode) and no structureCheck. The k-mer filter is the only specificity filter.
5. Windows are cut from the unmasked FASTA, not a repeat-masked one (section 4).
6. U-Probe filters all windows first and then tiles greedily, shortest window first. OligoMiner tiles first (blockParse) and filters afterwards, which leaves gaps. U-Probe can therefore return denser sets for the same parameters (inferred from the code, not measured).
7. The 21-nt primer flanks come from the 2012 format, which paired them with a 32-nt homology. `temperature_c: 42` only affects `fold_score`, which the template does not compute.

## 7. Pitfalls

- **Stale Jellyfish index (confirmed in code).** The index is `<fasta dir>/jf_genome/<fasta stem>.jf`, with no k in the name, and is built by the first `kmer_count` run. Auto-generated DNA protocols and `tests/data/DNA_format.yaml` use `kmer_len: 35`. A later run with `kmer_len: 18` queries the 35-mer index, while `cal_kmer_count` splits the output assuming 18-mers, so counts shift between probes and most probes get 0. The filter then passes nearly everything without an error. Delete `jf_genome/` whenever you change `kmer_len`, and check that a known repeat (an Alu) gets a high count.
- Every start position yields one window per length (6 per bp for 36-41 nt at `step: 1`). For multi-Mb paints, split the region into targets of a few Mb, or start from a precomputed OligoMiner/PaintSHOP set.
- Templates ship with N primers; ordering them as-is gives an unamplifiable library.
- Tm windows from OligoArray (75-90 C) or ORCA (>65 C) do not transfer to 50% formamide conditions; recompute the Tm under the declared conditions.
- Too few probes give dim or missing spots. A 5-kb target at 5 probes/kb has only 25 probes; use Coverage settings, denser tiling or signal amplification, or enlarge the target.
- RNA labelling by antisense probes (section 4) and segmental duplications (k-mer counts of 2-4) are the usual sources of a second, false spot.
- `avoid_otp` without a `mapped_sites` attribute is skipped (with an error in the log); `uprobe validate` also flags it.

## 8. Verification checklist

- [ ] Targets resolve in the FASTA (`--check-targets`); chromosome names and the assembly match the coordinates.
- [ ] Every homology: length within the preset, GC 20-80%, Tm in the window under the declared conditions, no homopolymer runs, no N.
- [ ] Max 18-mer count <= threshold, with the Jellyfish index built at that k; an Alu control gives a high count.
- [ ] Optional: `homology_sites_num == 1` for most probes; no cluster of off-target hits (`avoid_otp` or a BLAT/BLAST spot check).
- [ ] Probes do not overlap. Probes per kb meet the goal (locus >= ~5, paint ~0.5-2); the largest gap is known; total >= ~200 for a >= 10-kb locus.
- [ ] The homology strand is consistent across the set; over transcribed regions the probe is sense to the RNA, or RNase A is in the protocol.
- [ ] Primers and readouts are real sequences (no N), unique per sub-pool or segment, absent from the genome, and not complementary to the homologies; T7 promoter present if needed.
- [ ] Library oligo length is within the synthesis limit (74-120 nt in the cited designs).

## 9. References

- Beliveau 2012, PNAS. Versatile design and synthesis platform for visualizing genomes with Oligopaint FISH probes. doi:10.1073/pnas.1213818110
- Beliveau 2014, Curr Protoc Mol Biol. Visualizing genomes with Oligopaint FISH probes. doi:10.1002/0471142727.mb1423s105
- Beliveau 2015, Nat Commun. Single-molecule super-resolution imaging of chromosomes and in situ haplotype visualization using Oligopaint FISH probes. doi:10.1038/ncomms8147
- Beliveau 2018, PNAS. OligoMiner provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. doi:10.1073/pnas.1714530115 (code: github.com/beliveau-lab/OligoMiner, blockParse.py, outputClean.py, kmerFilter.py, structureCheck.py)
- Hershberg 2021, Nat Methods. PaintSHOP enables the interactive design of transcriptome- and genome-scale oligonucleotide FISH experiments. doi:10.1038/s41592-021-01187-3
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
