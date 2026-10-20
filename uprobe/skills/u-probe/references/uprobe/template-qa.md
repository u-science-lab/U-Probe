# Templates on real genomes: runtimes, yields and checks

Read this before running a template, to know what a run costs and what it usually yields, and after
it, to tell an expected shortfall from a broken design. The numbers come from one end-to-end test of
every template, run as `references/uprobe/workflow.md` describes; your genes will differ.

## Setup of the test

- U-Probe 0.2.0, 2026-10-19; macOS, 18 cores, `uprobe run -t 8`, indexes built beforehand.
- Human: UCSC hg38 primary assembly (soft-masked) with the Ensembl 113 GTF. Mouse: Ensembl 113
  GRCm39 primary assembly and GTF.
- RNA templates on 3 genes each, `transcripts: mane` (human) or `canonical` (mouse). Expression,
  where the method's sample has a public table: salmon transcript TPM of HeLa-S3 (ENCODE
  ENCSR000CPR) or NIH/3T3 (4DN 4DNESQAI59VQ), `level: transcript`, `min_value: 1`,
  `min_fraction_of_target: 0.05`. Brain targets of the padlock methods had no table.
- Placeholders were replaced by random sequences (MERFISH readouts and primers, seqFISH+ readouts,
  HybISS IDs, Oligopaint primers, double_hyb_rca barcodes and spacers); HCR used B1-B3.

## Index builds (once per genome)

| Index | Time | Disk |
| --- | --- | --- |
| Bowtie 1 transcriptome (`bowtie`), human | about 5 min | 650 MB |
| Bowtie2 transcriptome, human | 6 min (8 threads) | 855 MB |
| Bowtie2 genome + transcriptome (`bowtie2`), mouse | 16 min (9 threads) | 3.5 GB + 570 MB |

A run builds only the indexes its attributes search, so the RNA templates that use Bowtie2
(`hcr_v3`, `starmap`, `mipseq`, `double_hyb_rca`) build its transcriptome index and not the genome
index, which only `mapped_sites` (Oligopaint) needs; `uprobe build-index --all` builds both. The
Oligopaint run used Jellyfish 18-mer indexes restricted to the k-mers of its locus
(`jellyfish count --if`); a genome-wide one was not built in this test.

## Results

Wall time of `uprobe run` with the indexes in place. About 30 s of a human RNA run (20 s of a mouse
one) is parsing the GTF, once for target validation and extraction; the design itself takes 5-25 s.

| Template | Targets | Probes / requested | Time | Stage that removed most | Notes |
| --- | --- | --- | --- | --- | --- |
| `smfish_stellaris` | GAPDH, TFRC, EEF2 (HeLa) | 48/48 each | 55-75 s | Remove overlap | without a table GAPDH keeps 6 (its pseudogenes); TFRC keeps a 20-mer of 5% GC, Tm 34 C at a 37 C hybridization |
| `hcr_v3` | GAPDH, TFRC, EEF2 (HeLa) | 20/20 each | 60 s | Remove overlap | same pairs with `aligner: bowtie` |
| `merfish` | GAPDH, TFRC, EEF2 (HeLa) | 37, 92, 76 of 92 | 55-75 s | Filters (GC, Tm) | GAPDH's 1,285-nt mRNA holds at most 42 non-overlapping 30-mers |
| `seqfish_plus` | Gapdh, Tfrc, Eef2 (NIH/3T3, mouse) | 2, 24, 24 of 24 | 36 s | Filters (off-target) | Gapdh copies (Gm13394 and others) count; `exempt_copies` gives 24/24 |
| `padlock_iss` | Pvalb, Gad2, Aldoc (mouse) | 5/5 each | 37 s | Filters | target Tm 65-73 C (multi_padlock_design) |
| `starmap` | Aldoc, Pvalb, Gad2 (mouse) | 4/4 each | 56 s | Remove overlap | every unit is the shortest layout, 19 + 1 + 20 nt (below) |
| `mipseq` | Pvalb, Gad2, Aldoc (mouse) | 2/2 each | 35 s | Filters (`near_5p`) | no Tm filter: one Gad2 unit has an arm of Tm 24 C and a primer of 20 C at 37 C |
| `oligopaint` | chr19:10,000,001-10,200,000 as two 100-kb targets | 749 + 737 (7.4 per kb) | 55 s | Filters | largest gap 5.3 kb; OligoMiner's hosted hg38 set has 6.6 per kb here |
| `double_hyb_rca` | GAPDH, TFRC, EEF2 (HeLa) | 11, 5, 21 (no request) | 57 s | Filters (Tm 35-45 C) | `step: 20` (a 40-nt window every 20 nt) leaves few candidates; the technique file suggests `step: 1` |

- No final probe of any template listed an off-target gene (`<attr>_off_targets` empty); pseudogenes
  below the expression thresholds were listed in `_unexpressed` (GAPDHP1, GAPDHP65 and others).
- Options on `smfish_stellaris`: a gene-level table gave the same final counts as the transcript
  table but removed more candidates (TFRC 4,427 vs 4,803 passing); `off_target_pool` with `<= 2`
  without a table lifted GAPDH from 6 to 7 probes while 43 of the 103 probes carried an off-target
  gene; `-t 1` and `-t 8` gave identical tables.
- A misspelled gene stops `run` with exit 1 in about 25 s (the GTF read); `--continue-invalid`
  designs the others (`partial`, exit 0) or ends `no_valid_targets` (exit 3). A filter nothing
  passes ends `no_probes` (exit 3) and writes the raw table.

## Expected limits of the templates

- `smfish_stellaris`, `hcr_v3` and `mipseq` filter no GC and only a Tm floor taken from the published
  probe sets: smFISH and each HCR half at least the 37 C hybridization temperature (4% of the Raj
  2008 probes and 8% of the Choi 2018 pairs are below it), MiP-seq arms and primer site at least 28 C
  (2% of the published units are below it). AT-rich regions therefore give fewer probes.
- `starmap`: `remove_overlap` keeps the shortest window at each start, so the 40-46-nt layout always
  ends as 19 + 1 + 20 nt (padlock 52 nt, primer 32 nt); the published halves are 19-25 nt, balanced
  by Tm.
- `merfish` and `seqfish_plus` request more probes than short or copy-rich genes hold; report the
  shortfall (`fewer_probes_than_requested`) or use the iteration order of SKILL.md step 9.

## check_probes.py options per template

`check_probes.py` finds the parts named `<probe>.binding` by itself; name the others.

| Template | Options |
| --- | --- |
| `smfish_stellaris` | `--min-gap 2 --min-probes 25` |
| `hcr_v3` | `--min-gap 2` |
| `merfish` | `--gc 0.43 0.63 --tm 66 76 --tm-column binding_tm --min-gap 0` |
| `seqfish_plus` | `--gc 0.45 0.65 --min-gap 0` |
| `padlock_iss` | `--sense padlock.arm_5 --sense padlock.arm_3 --tm 65 75 --tm-column target_tm --min-gap 15` |
| `starmap` | `--min-gap 0` |
| `mipseq` | `--binding padlock.arm_5 --binding padlock.arm_3 --binding rca_primer.binding --min-gap 0` |
| `oligopaint` | `--sense oligopaint_oligo.homology --gc 0.2 0.8 --tm 42 47 --tm-column homology_tm --min-gap 0` (no `--transcripts`) |
| `double_hyb_rca` | `--binding circle_probe.part1 --binding circle_probe.part3 --binding amp_probe.part2 --tm 35 45 --min-gap 0` |
