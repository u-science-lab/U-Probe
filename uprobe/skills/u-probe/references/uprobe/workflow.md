# U-Probe workflow: from a request to delivered probes

Read this before running U-Probe for a design request, and again before you hand results back.
It covers genome setup, templates, editing the protocol, `validate`, `run`, reading the outputs,
iterating, and the checks to make before delivery. Field details: `references/uprobe/yaml-reference.md`;
failures: `references/uprobe/troubleshooting.md`. Paths are relative to the skill folder.

## 0. Environment

- `uprobe --version` (or `python -m uprobe`). Python 3.10+.
- External tools, called by name from `PATH`: `bowtie` and `bowtie-build` (Bowtie 1, bioconda package `bowtie`),
  `bowtie2`/`bowtie2-build`, `jellyfish` (only for `kmer_count`), `makeblastdb`/`mmseqs` (only if listed in `align_index`).
  Install missing ones with `conda install -c bioconda bowtie bowtie2 jellyfish`. U-Probe never installs them itself.
  When conda is slow or lacks a package, Bowtie 1 (`bowtie`, `bowtie-build`, `bowtie-inspect`) also comes as
  prebuilt binaries in the release assets of github.com/BenLangmead/bowtie (zip per platform, e.g.
  `bowtie-1.3.1-linux-x86_64.zip`); unzip and put the folder on `PATH`. Never replace it with your own scanner.
- Use one Python environment for `uprobe`, its dependencies and the commands you run (mixing interpreters
  gives import errors such as a numpy mismatch).
- Runs write a scratch directory `.uprobe-tmp-*` in the current directory, and resolve relative paths
  (`genomes.yaml` entries, `gene_expression.table`) from it. Run from a writable directory.

## 1. Get a genome

```bash
uprobe genome search mouse                     # Ensembl species names, assemblies, taxids
uprobe genome fetch --species mus_musculus --release latest --out ./genomes --genomes genomes.yaml
uprobe genome fetch --species homo_sapiens --release 113 --dry-run   # URLs and entry, no download
uprobe genome fetch --species homo_sapiens --soft-masked   # repeats in lowercase (oligopaint: masked: drop)
uprobe genome list -g genomes.yaml [--json]    # files present? indexes built?
```

- Ensembl and Ensembl Genomes (`--division plants|fungi|metazoa|protists`); not Ensembl Bacteria.
  Primary-assembly FASTA when available, else toplevel; checksums verified; downloads resume.
- The entry name defaults to `<species>_<assembly>_e<release>`; the command prints it. Put it in `genome:`.
- `align_index` defaults to `[bowtie2]` plus `bowtie` when `bowtie-build` is on `PATH`. **Make sure `bowtie`
  is in it** for short probes (<= 30 nt): only Bowtie 1 searches off-targets exhaustively (`references/concepts/specificity.md`).
- Pin `--release` for reproducible work; `latest` is resolved and recorded in the entry's `source`.
- Your own FASTA + GTF: write the entry by hand (`fasta`, `gtf`, `align_index`, `jellyfish`). GTF `gene_name`
  attributes are what `targets` match; MANE/canonical selection needs `tag "MANE_Select"` /
  `tag "Ensembl_canonical"` attributes on the transcript's GTF lines (Ensembl or GENCODE write them).
- Only a transcriptome FASTA (spliced transcripts) plus a transcript table with MANE/canonical tags: treat each
  transcript as a one-exon contig and write the two files yourself, then use the entry like any other. The
  design then needs no genome, and only the index the decisive check needs (Bowtie 1 for windows of
  <= 30 nt; in HCR also, `references/concepts/specificity.md`):

```python
# rows: (transcript_id, gene_id, gene_name, tag) with tag "MANE_Select", "Ensembl_canonical" or ""
with open('transcripts.gtf', 'w') as gtf:
    for tid, gid, name, tag in rows:
        attrs = f'gene_id "{gid}"; gene_name "{name}"; transcript_id "{tid}";' + (f' tag "{tag}";' if tag else '')
        gtf.write(f'{tid}\ttable\texon\t1\t{lengths[tid]}\t.\t+\t.\t{attrs}\n')   # FASTA record = transcript id
```

```yaml
tx_hs:                                   # genomes.yaml
  fasta: /data/tx/transcripts.fa         # one record per transcript id, spliced, + strand
  gtf: /data/tx/transcripts.gtf
  align_index: [bowtie]                  # add bowtie2 only if an attribute needs it
  jellyfish: false
```

  Targets are the table's `gene_name`s (a different spelling than the user's symbols is the usual miss).
  U-Probe cuts `transcript.fa` (records `<gene_id>_<transcript_id>`) from these two files and indexes it:
  budget about three times the FASTA in disk (FASTA, `transcript.fa`, index; a Bowtie 1 index of the human
  Ensembl 113 transcriptome is 650 MB, below). An index of the same
  transcriptome built earlier is reused when `bowtie_transcript/` (files named `<fasta stem>.*`) and
  `transcript.fa` sit next to the FASTA (symlink them); the same holds for `jf_genome/` indexes.
- Before any genome-scale build (Bowtie2 genome index, genome-wide Jellyfish table, soft-masked copy of the
  FASTA) estimate the size against the free disk (`df -h .`): a Bowtie2 genome index of a mammalian genome takes
  a few GB (mouse GRCm39: 3.5 GB, 13.5 min, below), a genome-wide Jellyfish table grows with the number of
  distinct k-mers of the whole genome, and temporary files come on top. Prefer an existing or
  shared index, and for a locus-sized DNA FISH target count only the locus' k-mers
  (`references/techniques/oligopaint-dna-fish.md`, section 6). Never copy a genome to mask it.

### Indexes

`uprobe run` builds the missing indexes of `align_index` that the protocol's attributes search before it
designs anything; build them ahead with `uprobe build-index -p protocol.yaml -g genomes.yaml -t 16` (the
same indexes; `--all`: every index of `align_index`) or `genome fetch --build-index` (every one).

| `align_index` entry | Builds (next to the FASTA) | Needed by |
| --- | --- | --- |
| `bowtie` | `bowtie_transcript/` (transcriptome only) | `mapped_genes` with `aligner: bowtie` |
| `bowtie2` | `bowtie2_genome/` (for `mapped_sites`) and `bowtie2_transcript/` (for `mapped_genes`), each only when an attribute searches it | `mapped_sites`; `mapped_genes` with `aligner: bowtie2` |
| `blast`, `mmseqs` | genome (and BLAST transcript) databases | nothing in the current attributes |
| `jellyfish: true` | `jf_genome/<stem>.k<kmer_len>.jf`, one per k-mer length (`.stranded.jf` for `strand: sense`/`antisense`), at the first `kmer_count` with it | `kmer_count` |

- The transcriptome is `transcript.fa` (every GTF transcript, spliced, records `<gene_id>_<transcript_id>`),
  extracted next to the FASTA once.
- Bowtie 1 transcriptome index of human Ensembl 113: about 5 min and 650 MB (README). A Bowtie2 *genome*
  index of a mammalian genome takes much longer (mouse GRCm39: 13.5 min, 3.5 GB); RNA designs, which
  only search the transcriptome, never build it.
- `uprobe run --threads/-t` (default 10) sets the threads of index building and of every aligner and
  Jellyfish step of an attribute without its own `threads:` option; an attribute's `threads:` wins.

## 2. Start from a template

```bash
uprobe templates list
uprobe templates show smfish_stellaris -o protocol.yaml
```

Templates: `smfish_stellaris`, `merfish`, `seqfish_plus`, `hcr_v3`, `padlock_iss` (HybISS), `starmap`,
`mipseq`, `oligopaint` (DNA), `double_hyb_rca`. Each has a `meta` block (ignored by the run) with
citations, `parameter_sources` (value -> paper) and `replace` (what you must change). Read `meta.replace`
and `meta.limitations` before editing. The technology-specific choices are in `references/techniques/`.

## 3. Edit the protocol

Always change: `genome` (a key of `genomes.yaml`), `targets` (GTF `gene_name`s, case-sensitive; or
`chrom:start-end`, 1-based inclusive as in a genome browser, with `source: genome`; or `name: SEQUENCE`), `encoding` (one entry per target that the
expressions read; `{}` when none), placeholder `N` sequences (primers, readouts).

Then decide, with the concept files:
- target region and transcripts (`references/concepts/targets-and-isoforms.md`);
- window length, Tm window under the real hybridization conditions (`references/concepts/hybridization-and-tm.md`);
- off-target attribute, threshold, expression table, copies, pool rule (`references/concepts/specificity.md`);
- probes per target (`equal_space`), spacing (`remove_overlap`).

Write `attributes` and `post_process` explicitly. If either is missing or empty, U-Probe generates
defaults (GC 0.2-0.8, Tm [T, T+40] °C on **every** probe and part, `target_mapped_genes <= 4`), which
are permissive and can filter on constant parts; see `references/uprobe/yaml-reference.md`.

## 4. Validate

Write the protocol file first, then validate, then run, one command after the other. Parallel commands made
`validate` read a file that did not exist yet (exit 2, a confusing message).

```bash
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets [--json]
```

- Exit 0: no errors (warnings may remain); exit 1: errors; exit 2: bad command line.
- Checks the schema, expression grammar (and builds every probe on synthetic windows), encodings,
  filters/sorts columns, hybridization and gene_expression blocks, aligners on `PATH`, indexes the
  attributes need, targets in the GTF/FASTA, `transcripts` tags. Each finding has a YAML line and a fix.
- Fix every error. Read every warning: several (Bowtie2 on <= 30-nt windows, `source: exon` on coding
  genes, `equal_space` without `remove_overlap`, Tm conditions without `temperature_c`) change results silently.
- `validate` does not check that index files exist (the run builds them) nor expression-table contents
  beyond columns.

## 5. Run

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o results/run1 --json > results/run1.json 2> results/run1.log
echo $?
```

| Exit | Meaning | `status` in `summary.json` |
| --- | --- | --- |
| 0 | probes designed | `ok`, or `partial` (a target has no probes or fewer than `equal_space` asked, or invalid targets were skipped with `--continue-invalid`) |
| 1 | the run failed | `error` (`stopped_at`, `error` say where and why) |
| 2 | bad command-line usage | no summary written |
| 3 | finished without any probe | `no_probes`, `no_probes_constructed`, `no_target_sequences`, `no_valid_targets` |

- Order of steps: normalize protocol -> validate targets (RNA) -> build indexes -> extract windows ->
  construct probes -> attributes -> post_process (`filters`, `avoid_otp`, `remove_overlap`,
  `off_target_pool`, `equal_space`, `sorts`) -> report.
- A target missing from the GTF stops the run before indexes are built (status `error`, exit 1; the
  error lists the targets). `--continue-invalid` skips them instead and lists them in `summary.json`
  (`inputs.invalid_targets`, warning `invalid_targets`; status `partial`, or `no_valid_targets` with
  exit 3 when no target is left). Run `validate --check-targets` first to catch them.
- `--json` prints `summary.json` on stdout; logs and tool output go to stderr.
- `--raw` also writes all candidates with their attributes before filtering (`*_raw.xlsx`); a run that
  ends with no probes writes it anyway. Excel limits never fail a run: a table longer than a worksheet
  is written as `*_raw.tsv` / `*_probes.tsv` only, and cells over 32,767 characters are cut in the
  `.xlsx` with a marker naming the TSV that keeps them (warning `excel_limits`).
- Use a fresh `-o` directory per run: tables are timestamped and accumulate, `summary.json` is overwritten.

## 6. Read the outputs

Files in `-o` (names in `summary.json -> output_files`): `<name>_<time>_probes.tsv` (final table, read
this one), `<name>_<time>.xlsx` (same, formatted), `<name>_report_<time>.html` (final table, filter
counts, a 100-row preview of the raw candidates and a link to the `--raw` table), `summary.json`.

`summary.json` keys to read first:
- `status`, `message`, `warnings[].code` (`invalid_targets`, `no_candidate_windows`,
  `probe_construction_failed`, `stage_removed_all`, `targets_without_probes`,
  `fewer_probes_than_requested`, `off_target_genes_not_in_expression_table`,
  `targets_not_in_expression_table`, `copies_not_checked`).
- `targets.<name>`: `candidate_windows`, `probes_constructed`, `after_stage` (count after each stage),
  `final`, `requested_probes`, `top_removal_stage`, `transcripts` (`selection`, `chosen_by`, `ids`).
- `stages`: before/after/removed per stage, all targets together. `Filters` is one stage: to see which
  filter removed candidates, rerun with `--raw` and count passes per condition on the raw table.
- `attributes`: min/median/max of the `summary.attributes` on the final table.
- `off_targets.attributes.<attr>`: aligner, `max_mismatches`, `strand`, `off_target_genes`,
  `final_off_target_genes`, copies (`copies`, `final_copy_genes`), expression counts; `off_targets.pool`.
- `inputs.hybridization_conditions`: the conditions each Tm/fold attribute was computed under.

Final table columns (TSV): `probe_id`, `target`, `exon_name`, `transcript_names`, `start`/`end` (within the
exon), `target_start`/`target_end` (along the target), `target_region` (the target RNA/DNA sense
sequence of the window), one column per probe and part (`probe`, `probe.part`), one per attribute, and
for each `mapped_genes` attribute `<attr>_off_targets` (`GENE:mismatches, ...`) and, as configured,
`<attr>_unexpressed`, `<attr>_hit_transcripts`, `<attr>_copies`. Rows are in `sorts` order, not along the target: sort by `target`, `target_start` to
see the tiling.

## 7. Report, then relax

Change one thing at a time and rerun into a new `-o`. Typical loop:
1. Run with the method's published rules, the sample's expression table, `min_fraction_of_target` and
   `exempt_copies` already in place (they are the specificity model, not relaxations). If a target has fewer
   probes than requested, write the shortfall down first (probes, requested, `top_removal_stage`).
2. Look at `top_removal_stage`, then `references/uprobe/troubleshooting.md` ("Too few probes"). Relax in this
   order, one step per run, each reported with its yield: widen the region (`CDS` -> `mRNA`), then Tm/GC
   windows within published ranges, then `location_interval` down to the method's smallest published value.
   Add the expression settings of step 1 first if they were missing. Do not relax specificity
   (`max_mismatches`, `<= 1`) for single-probe-signal methods, and do not change the standard length or arm
   lengths of the set to fill the count.
3. Never loosen a filter you cannot justify from `references/concepts/` or the technique file, and do not add
   a hard filter the method does not publish; if you do, compare the yield with and without it.

## 8. Check before delivering

Do these on the final TSV and `summary.json`; report any failure instead of hiding it.

| Check | How | Expect |
| --- | --- | --- |
| Status | `status`, exit code | `ok`; `partial` only with an explanation per target |
| Probes per target | `targets.<t>.final` vs `requested_probes`; `check_probes.py` distinct oligos | the technique's minimum (technique file; e.g. Stellaris recommends at least 25 of up to 48 smFISH probes); count distinct oligos, not positions, for repeated targets |
| Right transcripts | `targets.<t>.transcripts` | the isoform the user asked for (`chosen_by`) |
| Strand | for a few rows: binding part == reverse complement of `target_region` (RNA probes; `references/concepts/probe-architecture.md`) | true; padlocks on cDNA are the exception |
| Own-gene hit | `target_mapped_genes` on `target_region` with `strand: sense` | >= 1 on every row (0 means the window was not found in its own gene: wrong `strand`, wrong attribute target, or Bowtie2 on < 22 nt); literal-sequence targets absent from the annotation are 0 by design |
| Off-targets | `<attr>_off_targets`, `off_targets.*.final_off_target_genes`, `_copies` | empty for `<= 1`; copies named and acceptable |
| Tm, GC | `attributes` stats; histogram from TSV | inside the windows you set; Tm computed under the declared conditions |
| Spacing | consecutive `target_start` differences per target; `check_probes.py --min-gap` (also over every site of a repeated window) | >= window length + `location_interval` |
| Method's own values | the last command of the technique checklist (`--gc`, `--tm`, `--min-gap`, per-half `--part-gc`/`--part-tm`) | the values of the technique file's table, not those of the protocol you wrote |
| Sequences | probe columns | only ACGT (no `N` placeholders left), expected lengths |
| Barcodes | encoding parts per target | the intended readout/initiator per gene; no two genes share a code unless intended |

Deliver: the TSV (or a derived order sheet), the protocol YAML used, `summary.json`, the U-Probe version
(`uprobe_version`), genome release (`inputs.genome_source`), and a short note of every non-default choice.

## Step-by-step and Python API (rarely needed)

`uprobe build-index`, `validate-targets`, `generate-targets`, `construct-probes`, `post-process`,
`generate-report` run single stages. `uprobe.UProbeAPI(protocol_config=Path(...), genomes_config=Path(...),
output_dir=Path(...)).run_workflow(threads=10)` returns the final DataFrame and writes the same files.
`uprobe generate-barcodes --strategy max_orthogonality --num-barcodes 16 --length 8 --alphabet ACT -o dir`
writes candidate barcodes (`barcodes.csv`, `barcodes.txt`).
