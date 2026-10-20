---
name: u-probe
description: Design and check FISH and in situ hybridization probes with the U-Probe CLI (smFISH and Stellaris, smiFISH, HCR v3, MERFISH, seqFISH+, padlock probes for in situ sequencing such as ISS and HybISS, STARmap, MiP-seq, targeted ExSeq, PRISM, RAEFISH, OpenFISH and other RCA designs, TDDN-FISH, Oligopaint and DNA FISH, oligo pools and libraries for them). Use when a user wants probes against an RNA (gene, transcript, isoform, intron) or a genomic DNA locus or chromosome, asks which binding length, Tm, GC, formamide, spacing, probe count or off-target rule to use, names a FISH or imaging-based spatial transcriptomics method, or mentions uprobe or U-Probe. Not for PCR or qPCR primers or CRISPR guides.
---

# U-Probe: FISH and in situ probe design

U-Probe turns a declarative protocol (YAML) into probes: it cuts windows from the target, builds
each probe from parts (binding region, readouts, initiators, primers, padlock backbones), computes
attributes (GC, Tm under the stated conditions, off-target genes with Bowtie 1 or Bowtie2, genome
k-mer counts), then filters and tiles. Design with it instead of writing your own scripts or aligner.

Paths below are relative to this skill's folder; `uprobe skill path` prints the copy shipped with
the installed package. Read the references only when a step sends you there.

## When to use

- Probes against RNA or genomic DNA for imaging or in situ sequencing: one gene, a panel, an oligo
  pool, or a review of probes someone else designed.
- Questions on probe parameters: binding length, Tm and formamide, GC, spacing, probes per target,
  isoforms, off-target rules.

## Setup

```bash
pip install uprobe                                   # Python >= 3.10; or `pip install -e .` in a clone
conda install -c bioconda bowtie bowtie2 jellyfish   # external tools; U-Probe never installs them
uprobe --version
which bowtie bowtie-build bowtie2 bowtie2-build jellyfish
```

- Bowtie 1 (`bowtie`, `bowtie-build`) is required for probes of ~30 nt or less (rule 4). Bowtie2
  serves genome searches and long windows; Jellyfish only `kmer_count` (DNA FISH).
- A missing tool is installed or reported to the user, never replaced by your own code.
- Work in a writable directory: runs write `.uprobe-tmp-*` there and resolve relative paths from it.
- Details: `references/uprobe/workflow.md` (section 0).

## Workflow

Copy this checklist into your notes and tick it off:

```text
- [ ] 1. Clarify: organism, cell type, technology, targets, probes per target, conditions, deliverable
- [ ] 2. Choose the technology and template (references/techniques/README.md)
- [ ] 3. Fetch genome + annotation; Bowtie 1 in align_index; build indexes
- [ ] 4. Copy and edit the template (transcript, source, hybridization, expression, probe count)
- [ ] 5. uprobe validate --check-targets: no errors, every warning read
- [ ] 6. uprobe run --json into a fresh output directory
- [ ] 7. Inspect summary.json and the *_probes.tsv
- [ ] 8. Verify: strand, transcript, coverage, off-targets, counts, Tm/GC, spacing
- [ ] 9. Iterate, one change per run
- [ ] 10. Deliver in the requested format, with the run's provenance
```

### 1. Clarify the request

Establish, from the request or by asking:
- organism (and the genome assembly or annotation release, if the user has one);
- sample and cell type, and whether expression data exist for it;
- technology, or the experimental goal it must serve (step 2);
- targets: gene symbols, a named isoform, a region (mRNA, CDS, UTR, intron) or a genomic interval;
- probes per target; hybridization conditions (formamide %, SSC or Na+, temperature, probe concentration);
- the deliverable (step 10).

Ask when the answer changes the design and cannot be inferred (organism, technology, isoform).
Otherwise take the method's published default and say so.

### 2. Choose the technology and template

Read `references/techniques/README.md` (goal -> technology -> file -> template), then the
technique file it names. `uprobe templates list` shows the templates; smiFISH has none, its file
gives a protocol.

### 3. Genome, annotation and indexes

```bash
uprobe genome search human
uprobe genome fetch --species homo_sapiens --release <N> --out ./genomes --genomes genomes.yaml \
    --align-index bowtie,bowtie2
uprobe genome list -g genomes.yaml
uprobe build-index -p protocol.yaml -g genomes.yaml -t 16   # after step 4; `run` builds missing ones too
```

- Pin `--release`. MANE Select tags exist for human only; Ensembl canonical for every species.
- `align_index` must include `bowtie` for RNA probes of ~30 nt or less. `run` and `build-index` build
  only the indexes the attributes search: the slow Bowtie2 genome index only for `mapped_sites`
  (DNA FISH); `build-index --all` builds every one. DNA FISH needs `jellyfish: true`.
- Own FASTA + GTF: write the `genomes.yaml` entry by hand (`references/uprobe/yaml-reference.md`).

### 4. Copy and edit the template

```bash
uprobe templates show smfish_stellaris -o protocol.yaml
```

Read its `meta.replace` and `meta.limitations`, then set:
- `genome`, `targets` (GTF `gene_name`, case-sensitive: `GAPDH` in human, `Gapdh` in mouse) and
  `encoding` (real readout, initiator, barcode or primer sequences; no `N` left);
- the designated transcript, `extracts.target_region.transcripts`: `mane` (human) or `canonical`,
  or `[ENST...]` for an isoform the user names; `source` as the technique file says;
- `hybridization:` the user's conditions (`temperature_c`, `na_mM` with 2x SSC = 390,
  `formamide_pct`, `probe_nM` per oligo, `mg_mM`, `dntp_mM: 0`); keep per-attribute overrides that
  reproduce a published Tm window under its own conditions;
- `gene_expression:` when the cell type is known (see "Expression profile" below);
- `post_process.equal_space.number_desired`: the probe count, with `remove_overlap` at the
  method's spacing.

Keep `attributes` and `post_process` explicit: when they are empty U-Probe generates permissive
defaults. Every field: `references/uprobe/yaml-reference.md`.

### 5. Validate

```bash
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets --json
```

Exit 0 means no errors. Fix every error, then read every warning: Bowtie2 on windows of 30 nt or
less, `source: exon` on coding genes, `equal_space` without `remove_overlap` and conditions without
`temperature_c` change results silently.

### 6. Run

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o results/run1 -t 8 --json > results/run1.json 2> results/run1.log
```

Exit 0: probes designed (`ok`, or `partial` when a target is short); 1: failed (`stopped_at`,
`error`); 2: bad command line; 3: no probes. Use a new `-o` for every run.

### 7. Inspect

- `summary.json`: `status`, `warnings[].code`; per target `final` vs `requested_probes`,
  `top_removal_stage`, `transcripts` (`chosen_by`, `ids`); `attributes` (min/median/max);
  `off_targets`; `inputs.hybridization_conditions`.
- `<name>_<time>_probes.tsv` (`output_files.tsv`): one row per probe. Rows follow `sorts`; sort by
  `target`, `target_start` to see the tiling. Columns: `references/uprobe/workflow.md` (section 6).

### 8. Verify

```bash
python <skill>/scripts/check_probes.py results/run1/<name>_<time>_probes.tsv \
    --summary results/run1/summary.json --transcripts <genome dir>/transcript.fa
```

The script (standard library only; `--help`) fails when a binding part is not the reverse
complement of `target_region`, a window is not on the transcript U-Probe designed on, or a sequence
holds non-ACGT letters; it reports per target the count vs request, the span along the transcript,
gaps and overlaps, binding length, GC, Tm and off-target columns. Add `--gc`, `--tm`, `--min-gap`,
`--min-probes` with the method's values; `--sense` for padlock arms that bind cDNA (options per
template: `references/uprobe/template-qa.md`).
`transcript.fa` (records `<gene_id>_<transcript_id>`) appears next to the genome FASTA once a
transcriptome index is built. Then confirm:

| Check | Expect |
| --- | --- |
| Strand | binding part = rc(RNA); readout sites = rc(readout probe) |
| Transcript | `chosen_by` is `MANE_Select`, `Ensembl_canonical` or the user's isoform; probes spread over the mRNA |
| Off-targets | `<attr>_off_targets` empty; genes in `_copies` and `_unexpressed` named and acceptable |
| Counts | the technique's minimum per target, or the shortfall explained |
| Tm, GC | inside the windows; Tm computed under the declared conditions |
| Spacing | gaps at least the method's spacing; no overlaps unless intended |

Finish with the checklist at the end of the technique file and `references/uprobe/workflow.md`
(section 8).

### 9. Iterate

Change one thing per run. Too few probes: follow `top_removal_stage` into
`references/uprobe/troubleshooting.md` ("Too few probes"). Relax in this order and report each
step: region (`CDS` -> `mRNA`), expression table, `exempt_copies`, `min_fraction_of_target`, Tm/GC
within the method's published ranges, spacing. Never relax the off-target rule of a single-unit
method (HCR, padlock/RCA), and never loosen a filter the references do not justify.

### 10. Deliver

Produce the format the user asks for (TSV, CSV, FASTA, an order sheet) from the final TSV with code,
never by retyping. Include sequences 5'->3' with probe ids and targets; the protocol YAML and
`summary.json`; the U-Probe version (`uprobe_version`), genome release (`inputs.genome_source`) and
transcript ids; the hybridization conditions the Tm values refer to; modifications U-Probe does not
add (dyes, 5' phosphate on padlocks); every non-default choice, template departure and shortfall.

## Core rules

1. A probe binds its reverse complement: an RNA probe's binding part is `rc(target_region)` and
   readout sites are rc of the labelled oligo; only cDNA padlock arms carry the RNA sense
   (`references/concepts/probe-architecture.md`).
2. Design on the reference transcript (`transcripts: mane` in human, `canonical` elsewhere) unless
   the user names isoforms; not the longest isoform and not `all`, which add alternative exons, 3'
   ends and retained introns the reference mRNA lacks (`references/concepts/targets-and-isoforms.md`).
3. Use the chosen method's published windows (binding length, GC, Tm with its conditions, spacing,
   probes per target) from its technique file; do not widen them by default, and state every
   departure (`references/techniques/`, `references/concepts/hybridization-and-tm.md`).
4. Windows of ~30 nt or less need an exhaustive off-target search: `mapped_genes` with
   `aligner: bowtie`, `max_mismatches: 2`, `strand: sense`; Bowtie2 cannot align < 22 nt and misses
   sites up to 30 nt; never write your own aligner (`references/concepts/specificity.md`).
5. A near-match (<= 2 substitutions) in another expressed gene excludes the probe
   (`target_mapped_genes <= 1` as a filter, not a score); only many-probe methods (smFISH, smiFISH,
   MERFISH, seqFISH+, Oligopaint) may add the pool rule (`off_target_pool`); HCR pairs and
   padlock/RCA units stay strict per unit (`references/concepts/specificity.md`).
6. "Expressed" means expressed in the user's sample: give `gene_expression` a table for that cell
   type; without one every annotated gene counts (`references/concepts/specificity.md`).
7. Give every Tm with its conditions; a Tm window without them is meaningless
   (`references/concepts/hybridization-and-tm.md`).
8. U-Probe's output is the source of truth: never hand-edit or retype sequences; if you change one
   (trim, add a flank, swap a readout), re-run U-Probe or re-check it with `scripts/check_probes.py`
   and the off-target search (`references/uprobe/workflow.md`).
9. Report shortfalls: name each target with fewer probes than requested and why; never pad a set
   with probes that fail the method's rules (`references/uprobe/troubleshooting.md`).

## Where to read what

| Situation | Read |
| --- | --- |
| Choosing a technology and template | `references/techniques/README.md` |
| smFISH, Stellaris, Raj-style probes | `references/techniques/smfish.md` |
| smiFISH, FLAP readout (no template; protocol in the file) | `references/techniques/smifish.md` |
| HCR v3 split-initiator pairs, amplifiers B1-B5 | `references/techniques/hcr-v3.md` |
| MERFISH encoding probes, readouts, codebooks, abundance limits | `references/techniques/merfish.md` |
| seqFISH+ primary probes, readout rotation | `references/techniques/seqfish-plus.md` |
| Padlock probes, ISS, HybISS, direct-RNA padlocks | `references/techniques/padlock-iss-hybiss.md` |
| STARmap SNAIL padlock + primer pairs | `references/techniques/starmap.md` |
| MiP-seq padlock + RCA primer | `references/techniques/mip-seq.md` |
| `double_hyb_rca` circle + amplification probe | `references/techniques/double-hyb-rca.md` |
| Targeted ExSeq padlocks (expanded tissue, SplintR, SBS/SOLiD barcodes) | `references/techniques/exseq.md` |
| PRISM colour-intensity padlocks (one imaging round) | `references/techniques/prism.md` |
| RAEFISH reverse padlock + priming probe + encoding probe | `references/techniques/raefish.md` |
| OpenFISH split probe pairs + reusable padlock | `references/techniques/openfish.md` |
| TDDN-FISH primary probes (DNA-nanostructure amplifier) | `references/techniques/tddn-fish.md` |
| Oligopaint, DNA FISH, chromosome paints, chromatin tracing, repeats | `references/techniques/oligopaint-dna-fish.md` |
| Commands, outputs, pre-delivery checks | `references/uprobe/workflow.md` |
| Runtimes, probe yields and `check_probes.py` options of each template on real genes | `references/uprobe/template-qa.md` |
| Writing `protocol.yaml` or `genomes.yaml` | `references/uprobe/yaml-reference.md` |
| `validate` or `run` fails, too few probes, wrong strand, odd off-target counts | `references/uprobe/troubleshooting.md` |
| Probe parts, orientation, split pairs, padlock arms | `references/concepts/probe-architecture.md` |
| Tm window, salt, formamide, GC, structure, homopolymers | `references/concepts/hybridization-and-tm.md` |
| Aligners, mismatches, expression tables, pseudogenes, paralogs, pool rule, repeats | `references/concepts/specificity.md` |
| mRNA vs CDS vs UTR vs intron, MANE/canonical, isoforms, junctions, short targets | `references/concepts/targets-and-isoforms.md` |

## Expression profile for a named cell type

When the user names a cell type or tissue but gives no expression table:

1. Ask for their own RNA-seq of the same cells (bulk or single-cell): the best source.
2. Otherwise find a public dataset that matches the species, the cell line, primary cell type or
   tissue (and its state), and the identifiers of the genome's annotation (Ensembl/GENCODE ids,
   ideally the same release), with TPM or another per-gene abundance. Look in processed RNA-seq of
   consortia and repositories (ENCODE, GEO, ArrayExpress), the cell-line and tissue panels of
   expression atlases (Human Protein Atlas, GTEx), and single-cell atlases aggregated to pseudobulk
   per cell type. No dataset is the default: choose one that fits and record its accession.
3. Mixed tissue: take per gene the maximum over the cell types present, so an off-target expressed in
   a minority population still counts (`references/concepts/specificity.md`).
4. Transcript level, when reads are available: quantify them with salmon (Patro R et al. 2017 Nat
   Methods 14:417, doi:10.1038/nmeth.4197) or kallisto (Bray NL et al. 2016 Nat Biotechnol 34:525,
   doi:10.1038/nbt.3519) against the transcript FASTA of the genome's annotation release (headers =
   transcript ids), and set `level: transcript`: only the isoforms that carry a site then count.
5. Configure `gene_expression` as the template comments show (`min_value: 1` TPM, optionally
   `min_fraction_of_target: 0.05`; `absent: count` unless the table covers the whole annotation), and
   read the `off_target_genes_not_in_expression_table` warning (`references/uprobe/yaml-reference.md`).
6. With no fitting dataset, run without a table (every annotated gene counts) and say so. Never
   estimate expression from memory.
