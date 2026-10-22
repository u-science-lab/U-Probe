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
Oligopaint template counts its 18-mers with `locus_only: true`: no genome-wide Jellyfish table, one
pass over the genome FASTA per run (37-47 s for human with 8 threads; Jellyfish 2.3.0).

## Results

The `smfish_stellaris` rows below were measured before the template filtered GC to 40-60 % and kept 3 nt
between probes (it then used 2 nt and no GC filter); expect fewer probes per gene with them, and report the
shortfall (`references/uprobe/troubleshooting.md`). The other rows are unchanged.

Wall time of `uprobe run` with the indexes in place. About 30 s of a human RNA run (20 s of a mouse
one) is parsing the GTF, once for target validation and extraction; the design itself takes 5-25 s.

| Template | Targets | Probes / requested | Time | Stage that removed most | Notes |
| --- | --- | --- | --- | --- | --- |
| `smfish_stellaris` | GAPDH, TFRC, EEF2 (HeLa) | 48/48 each | 34 s | Remove overlap | Tm 37.4-65 C (floor 37); without a table GAPDH keeps 6 (its pseudogenes) |
| `hcr_v3` | GAPDH, TFRC, EEF2 (HeLa) | 20/20 each | 33 s | Remove overlap | halves 37.8-63.8 C (floor 37); same pairs with `aligner: bowtie` |
| `merfish` | GAPDH, TFRC, EEF2 (HeLa) | 37, 92, 76 of 92 | 32 s | Filters (GC, Tm) | GAPDH's 1,285-nt mRNA holds at most 42 non-overlapping 30-mers; Tm then by Primer3 (below for the design code's formula) |
| `seqfish_plus` | Gapdh, Tfrc, Eef2 (NIH/3T3, mouse) | 24/24 each | 21 s | Filters (off-target) | Gapdh's copies (Gm13394 and 43 others) exempt; without `exempt_copies` Gapdh keeps 2 |
| `padlock_iss` | Pvalb, Gad2, Aldoc (mouse) | 5/5 each | 21 s | Filters | target Tm 65-73 C (multi_padlock_design) |
| `starmap` | Aldoc, Pvalb, Gad2 (mouse) | 4/4 each | 24 s | Remove overlap | sites of 19-25 nt, gaps of 1 and 2 nt; halves 54.7-64.9 C, differing by 0.1 C (median) |
| `mipseq` | Pvalb, Gad2, Aldoc (mouse) | 2/2 each | 20 s | Filters (`near_5p`) | arms and primer site 29.1-50.1 C (floor 28) |
| `oligopaint` | chr19:10,000,001-10,200,000 as two 100-kb targets | 749 + 737 (7.4 per kb) | 32 s | Filters | largest gap 5.3 kb; OligoMiner's hosted hg38 set has 6.6 per kb here; report 1.2 MB; Tm then by Primer3 (below for OligoMiner's formula) |
| `double_hyb_rca` | GAPDH, TFRC, EEF2 (HeLa) | 11, 5, 21 (no request) | 32 s | Filters (Tm 35-45 C) | `step: 20` (a 40-nt window every 20 nt) leaves few candidates; the technique file suggests `step: 1` |

- No final probe of any template listed an off-target gene (`<attr>_off_targets` empty); pseudogenes
  below the expression thresholds were listed in `_unexpressed` (GAPDHP1, GAPDHP65 and others).
- Options on `smfish_stellaris`: a gene-level table gave the same final counts as the transcript
  table but removed more candidates (TFRC 4,427 vs 4,803 passing); `off_target_pool` with `<= 2`
  without a table lifted GAPDH from 6 to 7 probes while 43 of the 103 probes carried an off-target
  gene; `-t 1` and `-t 8` gave identical tables.
- A misspelled gene stops `run` with exit 1 in about 25 s (the GTF read); `--continue-invalid`
  designs the others (`partial`, exit 0) or ends `no_valid_targets` (exit 3). A filter nothing
  passes ends `no_probes` (exit 3) and writes the raw table.

## Design-code Tm formulas and locus-only k-mers (2026-10-22)

`merfish` now computes `binding_tm` with `tm_method: merfish_design` and `oligopaint` `homology_tm` with
`tm_method: oligominer`, the formulas their published windows were defined under (66-76 C at 300 mM and
5 nM; 42-47 C at 390 mM, 25 nM and 50 % formamide); `oligopaint` counts k-mers with `locus_only: true`.
Measured before (main at 1279145: Primer3 windows, a k-mer table of the locus' 18-mers) and after, U-Probe
at 1279145 with these template changes, `uprobe run -t 8`, indexes in place; UCSC hg38 primary assembly
(soft-masked) with the Ensembl 113 GTF, placeholders replaced by random sequences, and no expression
table (every annotated gene counts as an off-target, so the MERFISH counts are below the HeLa ones
above).

| Template | Target | Before | After | Notes |
| --- | --- | --- | --- | --- |
| `merfish` | GAPDH | 7 of 92 (105 windows pass the filters) | 7 (106) | off-target filter removes most |
| `merfish` | ACTB | 11 (147) | 11 (147) | |
| `merfish` | EEF1A1 | 33 (515) | 38 (585) | 36 of the 51 and 56 windows of the three genes are shared; each Primer3 set held one probe above 76 C by the code's formula |
| `oligopaint` | chr19:10,000,001-10,050,000 (50 kb, the template's) | 432 (45,402 pass) | 432 (47,358) | 72 windows shared; 117 of the 432 Primer3 probes (27 %) above 47 C by OligoMiner's formula |
| `oligopaint` | chr19:10,000,001-10,200,000 (two 100-kb targets) | 1,486 | 1,479 | 232 shared; 473 (32 %) above 47 C by OligoMiner's formula; hosted hg38b probes overlapping a template probe 97.3 % -> 100 % |

- Run times: `merfish` 34 s either way; `oligopaint` 10 s -> 56 s (50 kb) and 33 s -> 80 s (200 kb), the
  difference being the `locus_only` pass over the genome, which replaces building a genome-wide 18-mer
  table (the before runs used a table of the locus' k-mers built ahead, 49 s).
- `locus_only` alone (the old template with `locus_only: true`) gave the same 432 rows with the same k-mer
  counts as the table of the locus' k-mers.
- Primer3 at `probe_nM: 20` (4 x 5 nM, the earlier way to match the MERFISH code) gave the same 7, 11 and 38
  probes as `merfish_design`.

## Expected limits of the templates

- `smfish_stellaris` filters GC to the Stellaris designer's 40-60 % (Xiao 2023) and, like `hcr_v3` and
  `mipseq`, a Tm floor taken from the published probe sets: smFISH and each HCR half at least the 37 C
  hybridization temperature (4% of the Raj 2008 probes and 8% of the Choi 2018 pairs are below it),
  MiP-seq arms and primer site at least 28 C (2% of the published units are below it). `hcr_v3` and
  `mipseq` filter no GC. AT-rich regions therefore give fewer probes.
- `starmap` keeps halves of 50-65 C and, at each start, the split whose halves' Tm match best
  (`remove_overlap.prefer: tm_balance`); it does not aim at Picky's designed 60 C.
- `merfish` requests more probes than short genes hold; report the shortfall
  (`fewer_probes_than_requested`) or use the iteration order of SKILL.md step 9. `merfish`,
  `seqfish_plus` and `smfish_stellaris` exempt copies of the target: name the genes listed in
  `target_mapped_genes_copies` in the report.

## check_probes.py options per template

`check_probes.py` finds the parts named `<probe>.binding` by itself; name the others. Use the method's own
values from its technique file (the same commands end each technique checklist), never the values of the
protocol you just wrote: a filter written too loosely then verifies itself. `--min-probes` counts distinct
oligos; `--part-gc COLUMN MIN MAX` and `--part-tm COLUMN MIN MAX` check one part (a half of a split unit)
with its own window, `--range COLUMN MIN MAX` any numeric column. The Tm columns of `merfish`,
`oligopaint` and `padlock_iss` are computed with their design code's formula (`tm_method`), so their
published windows are checked as they stand.

| Template | Options |
| --- | --- |
| `smfish_stellaris` | `--gc 0.40 0.60 --min-gap 3 --min-probes 25 --tm 37 100` (`--min-gap 2` only for the Stellaris designer's spacing) |
| `hcr_v3` | `--min-gap 2 --min-probes 20 --tm 37 100` |
| `merfish` | `--gc 0.43 0.63 --tm 66 76 --tm-column binding_tm --min-gap 0 --min-probes 48` |
| `seqfish_plus` | `--gc 0.45 0.65 --min-gap 0 --min-probes 24` |
| `padlock_iss` | `--sense padlock.arm_5 --sense padlock.arm_3 --tm 65 75 --tm-column target_tm --min-gap 15 --min-probes 5` |
| `starmap` | `--min-gap 0 --min-probes 3 --part-gc target_region 0.30 0.70 --part-tm padlock_tm 50 65 --part-tm primer_tm 50 65 --part-gc padlock.binding 0.30 0.70 --part-gc primer.binding 0.30 0.70` |
| `mipseq` | `--binding padlock.arm_5 --binding padlock.arm_3 --binding rca_primer.binding --tm 28 100 --min-gap 0 --min-probes 2` |
| `oligopaint` | `--sense oligopaint_oligo.homology --gc 0.2 0.8 --tm 42 47 --tm-column homology_tm --min-gap 0`; no `--transcripts` |
| `double_hyb_rca` | `--binding circle_probe.part1 --binding circle_probe.part3 --binding amp_probe.part2 --tm 35 45 --min-gap 0` |
