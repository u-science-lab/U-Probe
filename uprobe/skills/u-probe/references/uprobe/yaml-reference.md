# U-Probe protocol YAML reference

Read this when writing or editing `protocol.yaml` or `genomes.yaml`: every block, its options, its
defaults and a minimal working example. Derived from the U-Probe code (0.2.x, `uprobe/core`), not from
memory; when the README and the code differ, this file follows the code and says so.

## Top level

| Key | Required | Notes |
| --- | --- | --- |
| `genome` | yes | a key of `genomes.yaml` |
| `targets` | yes | list (or mapping) of targets, below |
| `extracts` | yes | only `extracts.target_region` is read |
| `encoding` | yes (may be `{}`) | per-target barcodes/readouts that expressions read |
| `probes` | yes | mapping of probe name -> `template` + `parts` |
| `name` | no | prefix of output files |
| `attributes`, `post_process`, `summary` | no | generated when missing or empty (see "Defaults") |
| `hybridization`, `gene_expression` | no | conditions for Tm/fold; expression-aware off-targets |
| `meta`, `description` | no | ignored by the run |

Unknown top-level keys are ignored (`validate` warns).

## genomes.yaml

```yaml
GRCh38_e113:
  fasta: /data/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa   # plain (not gzipped) FASTA
  gtf: /data/genomes/Homo_sapiens.GRCh38.113.gtf
  align_index: [bowtie, bowtie2]   # bowtie2 | bowtie | blast | mmseqs
  jellyfish: false                 # true enables kmer_count
```

Only a transcriptome FASTA (and a transcript table)? Write each transcript as a one-exon contig: a GTF
with one `exon` line per transcript (`1` to its length, `+`, `gene_id`, `gene_name`, `transcript_id`, and
`tag "MANE_Select";` or `tag "Ensembl_canonical";` where the table says so), `fasta:` the transcript FASTA
(records = transcript ids) and `align_index: [bowtie]`; `references/uprobe/workflow.md` (section 1) has the
script and the sizes. Indexes built earlier are found when `bowtie_transcript/<fasta stem>.*`,
`transcript.fa` and `jf_genome/` sit next to the FASTA (symlinks work).

`uprobe genome fetch` writes such entries plus a `source` record. With `--soft-masked` it fetches Ensembl's
soft-masked FASTA (`dna_sm`: repeats in lowercase; default name suffix `_sm`). Every template can use it: the
aligners and Jellyfish ignore case, and only `extracts.target_region.masked: drop` reads the masking. See
`references/uprobe/workflow.md` for what each index builds.

## targets

```yaml
targets:
  - GAPDH                       # GTF gene_name (sources mRNA, exon, CDS, UTR); case-sensitive
  - chr19:10000000-10050000     # region (source: genome)
  - my_reporter: ACGTTGCA...    # a literal sequence, 5'->3', RNA sense; its windows become target_region
```

- Gene names are matched against GTF `gene_name` (falls back to `gene`), not gene ids. All gene ids
  sharing that name are pooled.
- Regions are 1-based and inclusive, as genome browsers and Ensembl show them: `chr1:1000-2000` is
  bases 1000-2000 (1,001 nt); a start of 0 is an error. `avoid_otp` `target_regions` read the same way.
  Windows come from the reference (+) strand; `sub_region` `s-e` is the window along the region, i.e.
  genomic bases `start+s-1` to `start+e-1`. (U-Probe before this fix read targets 0-based: bases
  1001-2000; designs from old protocols shift by one base.)
- A literal sequence named after its Ensembl gene id (as written in the GTF) lets `mapped_genes`
  recognize that gene as the target rather than as an off-target.

## extracts.target_region

```yaml
extracts:
  target_region:
    source: mRNA          # genome | mRNA | exon | CDS | UTR   (required)
    transcripts: mane     # all (default) | canonical | mane | [ENST..., ...]   (RNA sources)
    length: 20            # int, or inclusive range [min, max] (every length per start)
    step: 1               # distance between window starts (default 1); legacy overlap: step = length - overlap
    masked: keep          # source genome: keep (default) | drop
```

`masked: drop` (source `genome` only) skips every window with a soft-masked (lowercase) or N base, as
OligoMiner mines a repeat-masked genome; the other windows keep their `probe_id`. It needs a soft-masked
FASTA (`uprobe genome fetch --soft-masked`, UCSC `hg38.fa`); the run log and `validate --check-targets`
warn when a target has no lowercase base. `keep` upper-cases the windows and ignores masking.

| `source` | Windows come from |
| --- | --- |
| `genome` | the region string, + strand |
| `mRNA` | all `exon` features: 5' UTR + CDS + 3' UTR (whole processed RNA of non-coding genes) |
| `exon` | **CDS features** of coding genes (no UTRs); all exons only for genes without CDS. Same as `CDS`; `validate` warns |
| `CDS` | CDS features (exons for genes without CDS) |
| `UTR` | `UTR`, `five_prime_utr`, `three_prime_utr` features |

`transcripts`: `all` = union of every transcript's features (merged exon model; windows never cross a
splice junction; segments shorter than `length` give no window, except with `mRNA` where they still
count in positions; a target with no long-enough segment gets no windows, the others are unaffected). `canonical` = GTF tag `Ensembl_canonical`; `mane` = `MANE_Select`, else
`Ensembl_canonical` (logged); a list = those ids (version optional), each target needs one. When one
transcript per target is selected, `mRNA`/`exon`/`CDS` tile its **spliced** sequence, so windows also
span exon-exon junctions (`exon_name` then lists the joined exons with `;`). `UTR` never splices.
Missing tags or ids stop the run; `validate --check-targets` reports them first.

`layout` (optional) names contiguous sub-parts of each window, e.g. split probes:

```yaml
    length: 52
    layout:
      template: "{left}{gap}{right}"
      parts: {left: {length: 25}, gap: {length: 2}, right: {length: 25}}   # lengths: int or [min, max]
```

Expressions read them as `target_parts['left']`; attributes can target `target_parts.left`; columns
`target_parts.<p>`, `.length`, `.start`, `.end` and `layout_id` are added; with ranges, every split that
sums to the window length becomes its own candidate (`probe_id` suffix `__layoutN`).

## encoding

```yaml
encoding:
  GAPDH: {R1: ATCCTCCTTCAATACATCCC, R2: ACACTACCACCATTTCCTAT}
  ACTB:  {R1: ..., R2: ...}
```

Free-form mapping read by `encoding[target]['R1']`. Every target an expression looks up needs every
key it reads (`validate` checks). Use `encoding: {}` (or `GENE: {}`) when no expression reads it.

## probes and the expression language

```yaml
probes:
  encoding_probe:                         # top level: template + parts (no bare expr)
    template: "{readout}A{binding}"        # {part} placeholders + literal bases
    parts:
      binding: {expr: "rc(target_region)"}
      readout: {expr: "rc(encoding[target]['R1'])"}
  amp_probe:
    template: "{part}TTT{fixed}"
    parts:
      part:  {expr: "rc(target_region[-13:])", length: 13}   # length: declared, checked by validate only
      fixed: {template: "CAGGCTCCACTA", parts: {}}            # a constant sequence
```

- Parts are `expr`, or `template` + `parts` (nested). Names must not contain `.` or `:`; a part is
  addressed as `probe.part` (or `probe:part`), nested `probe.part.sub`.
- Allowed in `expr`: string literals (`"'ACGT'"`, note the inner quotes) and integers; names
  `target_region`, `target` (the target's name), `encoding`, `target_parts`, `probe_index`, and any
  probe/part name (including dotted, e.g. `circle_probe.part2.barcode2`); indexing and slicing
  (`target_region[0:20]`, `[-13:]`, `[::-1]`, a list in `encoding`: `encoding[target]['readouts'][2]`);
  `+` (concatenate strings or add integers); `%` (remainder of two integers); unary `-` on integers;
  `rc(x)` (reverse complement, one argument). Nothing else: no other functions, methods, attributes,
  comparisons, conditionals, f-strings or comprehensions. Expressions are parsed, never `eval`ed.
- `probe_index`: the probe's number within its target, 0, 1, 2, ... along the target, for parts that
  vary from probe to probe, e.g. readout rotation: `rc(encoding[target]['readouts'][(probe_index + 1) % 4])`
  with `encoding: {GAPDH: {readouts: [R_a, R_b, R_c, R_d]}}`. The final probes are numbered after
  `post_process` (column `probe_index`) and the parts that read it are rebuilt, with their attributes;
  filters and sorts see the candidates numbered in extraction order (`validate` warns about attributes
  on such parts). The raw table keeps the candidates' numbering.
- An expression can read other probes/parts (a DAG); cycles are errors.
- A candidate whose expression fails (e.g. a missing encoding key) gets an empty probe; if all do,
  status `no_probes_constructed`.

## attributes

Each attribute adds a column named after it, computed on `target` (a sequence column:
`target_region`, `probe`, `probe.part`, `target_parts.<p>`). Name attributes with letters, digits and
`_` only (filter conditions parse `\w+`).

| `type` | Value | Options |
| --- | --- | --- |
| `gc_content` | GC fraction 0-1 | - |
| `length` | sequence length | - |
| `annealing_temperature` | Tm (°C), primer3 SantaLucia NN + salt correction | `hybridization` override; `tm_method: multi_padlock_design` = the HybISS design code's formula (Sugimoto 1996 NN, `16.6 log10([Na+] + 3.795 sqrt[Mg2+])`, -0.72 °C per % formamide; no dNTP term), for windows published under it |
| `fold_score` | **-MFE** (kcal/mol, >= 0) of the sequence folded as RNA by ViennaRNA at `temperature_c`; higher = more structure, so lower is better (filter `<=`, sort ascending) | `hybridization` (only `temperature_c` used) |
| `self_match` | number of 4-mer pairs between the sequence and its reverse complement (crude self-complementarity); lower is better | - |
| `mapped_genes` | distinct genes with a site in the transcriptome, the target included | `aligner` (required: `bowtie` or `bowtie2`), `max_mismatches` (bowtie: 0-3, default 2), `strand` (`both` default, `sense`, `antisense`), `exempt_copies`, `threads` |
| `mapped_sites` | list of genome sites in `<attr>`, their count in `<attr>_num` (Bowtie2 `--very-sensitive-local -k 100`, both strands: at most 100; `<attr>_num == 1` is OligoMiner's Unique Mode) | `aligner: bowtie2` (required), `threads` |
| `kmer_count` | highest genome count among the sequence's k-mers (Jellyfish); 0 if shorter than k; lower is better | `aligner: jellyfish` (required), `kmer_len` (35; each k has its own index `jf_genome/<stem>.k<k>.jf`), `strand` (`both`, the default: canonical k-mers, both strands; `sense`: the k-mers as written, counted on the FASTA's forward strand, as OligoMiner's dictionary; `antisense`: their reverse complements; `sense` and `antisense` read `jf_genome/<stem>.k<k>.stranded.jf`), `size` (hash, `1G`), `threads`; genome `jellyfish: true` |
| `difference` | `|a - b|` of two attributes defined before it, e.g. the Tm of two probe halves; no `target` | `of: [a, b]` (required) |
| `blocks` | legacy; do not use | - |

`threads` (aligner and Jellyfish attributes): default `uprobe run --threads` (10); set it on an attribute
only to override the run's value.

```yaml
attributes:
  target_gc:   {target: target_region, type: gc_content}
  probe_tm:    {target: smfish_probe.binding, type: annealing_temperature}
  target_fold: {target: target_region, type: fold_score}
  target_mapped_genes:
    target: target_region
    type: mapped_genes
    aligner: bowtie          # exhaustive within max_mismatches substitutions; bowtie2 = local heuristic
    max_mismatches: 2
    strand: sense            # target_region: sites where the window itself occurs
    exempt_copies: {kmer: 30, min_shared: 0.5}   # or true
```

`mapped_genes` details (see `references/concepts/specificity.md` for when to use what):
- `strand: sense` for `target_region` (RNA sense windows); `antisense` when `target` is a probe's
  binding part (the reverse complement); `both` also counts antisense transcripts the probe cannot bind.
- Bowtie2 runs `--very-sensitive-local -k 100`: sequences < 22 nt never align (count 0), 22-30 nt are
  searched incompletely, at most 100 hits. Bowtie 1 runs `-v N -a`: every end-to-end site with <= N
  substitutions, no indels.
- Extra columns: `<attr>_off_targets` (`GENE:mismatches`, names only with bowtie2), `<attr>_unexpressed`
  (with `gene_expression`), `<attr>_hit_transcripts` (with `level: transcript`), `<attr>_copies` (with
  `exempt_copies`). `<attr> <= 1` = no counted off-target gene.
- `exempt_copies`: a gene with a transcript containing >= `min_shared` of the target's distinct
  `kmer`-mers (same orientation, exact) is not counted (listed in `_copies`). `kmer` 12-32. The
  target's sequence is read from its windows: keep `step` <= window length.

## hybridization

```yaml
hybridization:          # all optional; units below
  temperature_c: 37     # °C: fold_score temperature; start of auto Tm filters
  na_mM: 390            # monovalent cations (2x SSC = 390 mM Na+)
  mg_mM: 0
  dntp_mM: 0
  probe_nM: 50
  formamide_pct: 30     # % v/v; primer3 Blake & Delcourt correction
```

Unset fields keep the legacy PCR-like defaults (`na_mM` 50, `mg_mM` 0, `dntp_mM` 0.6, `probe_nM` 50,
`formamide_pct` 0; fold at 37 °C). Every template sets `dntp_mM: 0`, since FISH buffers hold no dNTPs
(without Mg2+ they change no Tm), and `probe_nM` per oligo where its protocol states it (smFISH: 125 nM
for the 48-oligo set, so 2.6). An `annealing_temperature` or `fold_score` attribute can carry its own
`hybridization:` mapping, merged over the block (used to reproduce a published Tm window computed under
other conditions, or the Tm of one part of a split unit, which keeps a published window computed under
its own concentrations separate from the experiment's buffer: one attribute per part, each with its
override, filtered per part and checked with `check_probes.py --part-tm`). Values >= 0; `na_mM`, `probe_nM` > 0; `temperature_c`, `formamide_pct` <= 100.

## gene_expression

```yaml
gene_expression:               # a table of transcripts (preferred when you have one)
  table: quant.sf              # salmon quant.sf, kallisto abundance.tsv, RSEM isoforms.results; may be .gz
  level: transcript            # gene (default) | transcript
  transcript_column: Name      # level: transcript only; default: first column
  value_column: TPM            # default: the only numeric column
  min_value: 1                 # required: off-target genes below it are not counted
  min_fraction_of_target: 0.05 # optional 0-1: nor those below this fraction of the target's value
  absent: count                # count (default) | ignore: identifiers missing from the table
---
gene_expression:               # a table of genes
  table: expression.tsv        # TSV (CSV if the name contains .csv), may be .gz; GTEx .gct works
  gene_column: gene_id         # level: gene only; or [gene_id, gene_name]; default: first column
  value_column: TPM
  min_value: 1
```

- `level: gene`: rows are looked up by Ensembl gene id (version stripped), then by gene name; a gene
  listed twice takes its highest value. An off-target gene is valued by its whole-gene expression.
- `level: transcript`: rows are transcript ids (with or without version; a GENCODE FASTA header by its
  first `|` field); the genome's GTF assigns them to genes. An off-target gene is valued by the **sum
  of the transcripts the sequence hits**, so a gene expressed only on isoforms without the site does not
  count; the target's value is the sum of all its transcripts. A hit transcript missing from the table
  makes its gene count (`absent: count`) or adds nothing (`absent: ignore`). Extra column
  `<attr>_hit_transcripts`, e.g. `Mir100hg: ENSMUST00000212079=0.014, ENSMUST00000362798=0` (`NA` = not
  listed). Bowtie2 reports <= 100 alignments, so genes with many isoforms can be under-summed.
- Applies to every `mapped_genes` attribute; the target gene always counts. A table that matches no
  annotated gene stops the run; so does a table whose ids fit the other `level` (`validate` reports it).

## post_process

Steps always run in this order, whatever their order in the YAML:
`filters` -> `avoid_otp` (DNA only) -> `remove_overlap` -> `off_target_pool` -> `equal_space` -> `sorts`.

```yaml
post_process:
  filters:
    probe_tm:   {condition: probe_tm >= 47 & probe_tm <= 60}
    target_mapped_genes: {condition: target_mapped_genes <= 1}
    gc:         {condition: [target_gc >= 0.35, target_gc <= 0.6]}     # a list = all must hold
    homopolymers:
      type: sequence_pattern
      target: target_region                     # any sequence column
      exclude_patterns: [AAAAA, TTTTT, CCCCC, GGGGG]   # regexes, case-insensitive
  remove_overlap: {location_interval: 2}        # min gap (nt) between kept windows
  off_target_pool: {attribute: target_mapped_genes, max_probes_per_gene: 2}
  equal_space: {number_desired: 48}             # or per target: {GAPDH: {number_desired: 30}}
  sorts: {is_ascending: [target_mapped_genes], is_descending: []}
```

- `condition`: `<column> <op> <number>` joined by `&`; ops `<=`, `>=`, `<`, `>`, `==`; numbers in plain,
  non-negative decimal (no `-`, no `1e-5`, no `|`/or: use separate filters). A condition on a column that
  does not exist is **silently skipped** at run time (`validate` reports it).
- `remove_overlap`: per target, visits candidates by position (`target_start` in RNA mode) and keeps
  one when it starts more than `location_interval` nt after the end of the last kept one. It packs from
  the 5' end; of the candidates of one start (a `length` range or a `layout` with ranges) it keeps the
  shortest, or with `prefer: <attribute>` the one with the lowest value (the STARmap template: the split
  whose halves' Tm match best; keeping a standard length: `references/uprobe/troubleshooting.md`, "Length
  ranges are not a way to fill the count"). It is first-fit: it never looks back, so when the passing
  pool is barely enough it leaves no alternative for a weak window. It does not choose between different
  starts, so put quality into `filters`, and add only filters the method publishes (compare the yield with and
  without any other).
- `off_target_pool`: per target, walks probes in `sorts` order and drops a probe whose off-target gene
  (from `<attr>_off_targets`) is already hit by `max_probes_per_gene` kept probes. `attribute` may be a
  list. Combine with a per-probe filter (e.g. `<= 2`); `<= 1` makes it a no-op (`validate` warns).
- `equal_space`: keeps `number_desired` per target at evenly spaced ranks along the target (first and
  last kept); fewer candidates -> all kept, status `partial`. Use with `remove_overlap`.
- `sorts`: orders the final table (and the pool's ranking). Unknown columns are dropped. Lower is better
  for `fold_score`, `self_match`, `kmer_count` and `mapped_genes`: put them in `is_ascending` (`validate`
  warns when one is in `is_descending`).
- `avoid_otp` (DNA, needs a `mapped_sites` attribute): per target `{target_regions: [chr:s-e] (1-based),
  density_thresh: 1e-4, search_range: [-1000000, 1000000], avoid_target_overlap: true}`; removes probes
  whose off-target genome hits cluster (a density of other probes' hits within `search_range`), which
  would form a false focus. Adds `on_target`, `off_target`, `target_ratio`.

## summary

```yaml
summary:
  report_name: rna_report      # rna_report | dna_report; empty = no HTML report
  attributes: [target_gc, probe_tm, target_mapped_genes]   # columns summarized in summary.json and the report
  report:
    embed_raw: false           # true: embed every raw candidate as a workbook (large for big designs)
```

The HTML report embeds the final table and a 100-row preview of the raw candidates, and links the raw
table that `--raw` writes next to it. With `embed_raw: true` it also embeds the whole raw table for
download; for a 200-kb Oligopaint design (581,793 candidates) that is 63 MB of a 66 MB report.

## Defaults when `attributes` / `post_process` / `summary` are missing or empty

- RNA mode attributes: `target_gc`, `target_tm`, `target_fold`, `target_self_match`, `target_mapped_genes`
  (Bowtie 1, 2 mismatches, `sense` when `bowtie` is in `align_index` and installed; Bowtie2 otherwise), plus
  `_gc`, `_tm`, `_fold`, `_self_match` for every probe and part. DNA mode: `target_kmer_count` on the
  genomic window (jellyfish, k 35), and those four for every part plus `_mapped_sites` (bowtie2) and
  `_kmer_count` (jellyfish, k 35).
- Filters: Tm [T, T+40] °C at `hybridization.temperature_c` (37 °C if conditions are declared without
  it; 50-90 °C without a block) on **every** Tm attribute, constant parts included; GC 0.2-0.8 on
  `target_gc` and whole probes; RNA `target_mapped_genes <= 4` (up to 3 off-target genes); DNA
  `target_kmer_count <= 100` (no 35-mer of the window more than 100 times in the genome: a loose repeat
  filter, far looser than OligoMiner's 18-mers < 5; see `references/concepts/specificity.md`).
- Sorts: all ascending: GC, Tm, mapped genes, then `_fold`, `_self_match`, `_kmer_count` (least
  structured, least self-complementary, least repetitive first). No `remove_overlap` or `equal_space`.
Write your own blocks for any real design.
