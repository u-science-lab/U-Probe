<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/uprobe-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="./assets/uprobe.svg">
    <img src="./assets/uprobe.svg" alt="uprobe logo" height="80" align="absmiddle">
  </picture>
  Universal Agentic Probe Design Platform
</h1>
<br>

[![PyPI](https://img.shields.io/pypi/v/uprobe)](https://pypi.org/project/uprobe/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/UFISH-Team/U-Probe/blob/main/LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.10+-blue)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://ufish-team.github.io/uprobe-official/)
[![Website](https://img.shields.io/badge/Website-Official-green.svg)](https://ufish-team.github.io/uprobe-official/)

**U-Probe** is a universal and agentic probe design platform tailored for imaging-based spatial-omics. It overcomes the architectural limitations of existing tools and lowers the expertise barrier by combining an innovative declarative configuration system with LLM-based AI agents.

Whether you are designing standard probes for established protocols or developing entirely novel architectures, U-Probe provides a comprehensive, automated workflow from target sequence extraction to rigorous thermodynamic filtering.

## Features

- **Universal Probe Architecture**: Features a declarative configuration system and a directed acyclic graph (DAG)-based assembly engine. Easily design arbitrary and complex multi-part probe structures (e.g., DNA-FISH, smFISH, MERFISH, seqFISH, $\pi$-FISH, MiP-seq, RCA) without modifying any source code.
- **Agentic AI Workflow**: An Agent Skill (`uprobe/skills/u-probe/`) lets AI coding agents such as Claude Code or Codex design probes with U-Probe, and the web server offers conversational design through hierarchical LLM-based AI agents (via PantheonOS, an optional extra). Users can describe experimental goals in plain language or provide scRNA-seq data, and the agents will autonomously construct configurations, select parameters, and execute the design pipeline.
- **Comprehensive Quality Filtering**: Automatically computes and filters candidates based on GC content, melting temperature, secondary structure stability (via ViennaRNA), off-target mapping (via Bowtie2, or exhaustively within a number of mismatches via Bowtie 1), and k-mer frequency (via Jellyfish).
- **End-to-End Automation**: Streamlines the entire process from target sequence extraction to advanced post-processing, including overlap removal and equal spacing for tiling designs.
- **Extensible API & Web UI**: Offers a clean Python API for programmatic integration, a command-line interface, and an interactive web server for visual workflow management and result analysis.

## Installation

Requires Python 3.10+. Use Linux or WSL2 for backend development and bioinformatics tools.

### Install via pip (Recommended)

```bash
pip install uprobe
```

### Install from source (For development)

```bash
git clone https://github.com/UFISH-Team/U-Probe.git
cd U-Probe
pip install -e .
```

### Using conda environment

```bash
git clone https://github.com/UFISH-Team/U-Probe.git
cd U-Probe
conda env create -f environments.yaml
conda activate uprobe
pip install -e .
```

The conda environment also installs the command-line tools U-Probe calls (Bowtie2, Bowtie 1, BLAST+, MMseqs2, Jellyfish). With a pip install, install the ones you use from bioconda, e.g. `conda install -c bioconda bowtie2 bowtie jellyfish`. Bowtie 1 (bioconda package `bowtie`) runs the [exhaustive off-target search](#off-target-search), which short probes need: the smFISH, seqFISH+, MERFISH and HybISS templates use it, and so does the generated off-target attribute of a protocol without `attributes` when the genome enables it. `uprobe genome fetch` adds `bowtie` to a new genome's `align_index` when Bowtie 1 is installed.

## Use U-Probe with an AI agent

U-Probe ships an [Agent Skill](https://agentskills.io/specification) in [`uprobe/skills/u-probe/`](uprobe/skills/u-probe/). Its `SKILL.md` gives a coding agent the design workflow (clarify the request, choose the technology and template, fetch the genome, edit, validate, run, verify, deliver) and the rules behind it, and sends the agent to reference notes on each technology, on U-Probe's YAML and on probe-design concepts when it needs them; `scripts/check_probes.py` checks a result table before delivery. Install U-Probe and the aligners, then copy the skill to where your agent looks for skills:

```bash
pip install uprobe
conda install -c bioconda bowtie bowtie2 jellyfish
uprobe skill install --target codex            # ./.agents/skills/u-probe (run it at the repository root)
uprobe skill install --target claude --user    # ~/.claude/skills/u-probe
uprobe skill path                              # the copy inside the installed package
```

| Agent | Project (default, working directory) | `--user` | Source |
| --- | --- | --- | --- |
| Codex | `.agents/skills/u-probe` | `~/.agents/skills/u-probe` | [Codex docs: Agent Skills](https://developers.openai.com/codex/skills) |
| Claude Code | `.claude/skills/u-probe` | `~/.claude/skills/u-probe` | [Claude Code docs: Skills](https://code.claude.com/docs/en/skills) |

Without `--target` the skill is installed for both agents; `--target DIR` installs it as `DIR/u-probe` for any other agent that reads Agent Skills. After upgrading U-Probe, run the install again with `--force` so the skill matches the CLI. Then describe the experiment to the agent in plain language, e.g. "Design smFISH probes against human GAPDH mRNA for HeLa cells".

## Usage Guide

U-Probe provides three main entry points: **Command Line Interface (CLI)**, **Python API**, and the **HTTP/Web UI**.

For source installs, install the package and runtime dependencies before running examples:

```bash
pip install -e .
```

Most workflow commands require two YAML files:

1. **`genomes.yaml`**: genome paths such as FASTA and GTF.
2. **`protocol.yaml`**: probe design parameters.

### 1. Command Line Interface (CLI)

The CLI is perfect for running tasks quickly in the terminal or shell scripts. A typical session, also the path for scripts and AI agents, fetches a reference genome, starts from a published template, checks the protocol and runs it; each step is described below.

```bash
uprobe genome fetch --species mus_musculus --release latest --genomes genomes.yaml
uprobe templates show padlock_iss -o protocol.yaml    # then set genome, targets and encoding
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --json
```

#### 🤖 Using U-Probe with an AI Agent
The recommended way to use U-Probe with an AI agent such as Claude Code or Codex is the Agent Skill under [`uprobe/skills/u-probe/`](uprobe/skills/u-probe/); install it with `uprobe skill install` (see [Use U-Probe with an AI agent](#use-u-probe-with-an-ai-agent)). U-Probe no longer ships a terminal agent of its own: the `uprobe agent` command has been removed. The web server keeps its chat agent as an optional extra (see below).

#### 🌐 Start Web Server 
Run the FastAPI backend in WSL/Linux and the separate React frontend in Windows. First check `config.ini`: use local Linux paths under `[Paths]` and `frontend_url = http://localhost:5173` under `[Server]`.

```bash
cd /path/to/U-Probe
conda activate uprobe  # or your existing Python 3.10+ environment
uprobe server --env development --host 127.0.0.1 --port 8005 --workers 1
```

Development mode reloads Python changes automatically. Alternatively, run `APP_ENV=development PORT=8005 python -m uprobe.http.server` from the project root. Running without overrides uses `config.ini`, which may select production mode.

The web chat (the `/agent` endpoints) is the optional `web-agent` extra: install it with `pip install 'uprobe[web-agent]'`, or `pip install -e '.[web-agent]'` from source. Without it the server starts and those endpoints answer 503; deployments that serve the chat, such as app.u-probe.org, must install it.

In a separate Windows PowerShell terminal:

```powershell
cd /path/to/probe-web-ui
pnpm install  # first run or after dependency changes
pnpm dev --host 127.0.0.1 --port 5173
```

- UI: http://localhost:5173
- API docs: http://localhost:8005/docs
- Connection check: http://localhost:5173/api/ should return `{"message":"Hello, World!"}`.

The frontend proxies `/api` to backend port `8005`. No separate database service is needed; users and tasks are stored locally. Create a test account through `POST /auth/register` in the API docs; this endpoint needs no email verification. Full probe workflows also require FASTA/GTF data and the relevant bioinformatics tools.

#### 📥 Fetch a Reference Genome
Download an Ensembl genome FASTA and GTF, verify them against Ensembl's `CHECKSUMS`, and register them in `genomes.yaml`:

```bash
uprobe genome fetch --species homo_sapiens --release latest --out ./genomes --genomes genomes.yaml
uprobe genome fetch --species arabidopsis_thaliana --division plants --build-index -t 10
uprobe genome search mouse            # find Ensembl species names
uprobe genome list -g genomes.yaml    # entries, file and index status (--json for scripts)
```

- `--species` takes an Ensembl name (`homo_sapiens`), display name (`human`) or taxonomy id (`9606`). `--release latest` is resolved to a number, which is recorded in the entry.
- Sources: Ensembl (vertebrates plus yeast, worm and fly) and the Ensembl Genomes divisions `plants`, `fungi`, `metazoa` and `protists` (`--division`; they have their own release numbers). Ensembl Bacteria is not supported.
- The unmasked primary assembly FASTA is used when Ensembl provides one, otherwise the toplevel FASTA. Files are decompressed (U-Probe's tools need plain FASTA) into `<out>/<name>/`; the entry name defaults to `<species>_<assembly>_e<release>` (`_eg<release>` for Ensembl Genomes).
- `--soft-masked` fetches the soft-masked FASTA instead (Ensembl `dna_sm`: repeats in lowercase; the name gets the suffix `_sm`). Every template can use it, since the aligners and Jellyfish ignore case; the Oligopaint template needs it, because its `extracts.target_region.masked: drop` skips the genome windows with a soft-masked or N base, as OligoMiner mines a repeat-masked genome.
- Re-running skips files that are already present and verified, and resumes interrupted downloads. `--dry-run` prints the resolved URLs without downloading; `--mirror` replaces the Ensembl FTP base URL.
- A new entry gets `align_index: [bowtie2]`, plus `bowtie` when Bowtie 1 (`bowtie-build`) is installed, for the [exhaustive off-target search](#off-target-search) (the log says so when it is not), and `jellyfish: true` when Jellyfish is installed. `--align-index` and `--jellyfish/--no-jellyfish` set them; re-fetching keeps the settings of an existing entry.
- `--build-index` runs the same index build as `uprobe build-index`.

#### Protocol Templates

U-Probe ships complete protocol templates for published probe technologies: smFISH (Stellaris-style), MERFISH, seqFISH+, HCR v3.0, STARmap, padlock ISS (HybISS), MiP-seq, Oligopaint/OligoMiner and the double-hybridization RCA design. Each template declares the source protocol's [hybridization conditions](#hybridization-conditions); where a paper gives a Tm window together with the conditions it was computed under (MERFISH, Oligopaint, HybISS), the Tm attribute is computed under those conditions and the window is applied as a filter. A `meta` block cites the source papers, maps every parameter to its source and lists what you must replace (genome, targets, barcodes/readouts). Every template passes `uprobe validate`. Every RNA template designs on the reference transcript (`transcripts: mane`: MANE Select, the Ensembl canonical transcript without it) and takes the region its source designed on (`meta.parameter_sources` quotes the evidence): the whole mRNA for smFISH, MERFISH, HCR v3, HybISS and MiP-seq, whose published probes also lie in UTRs, and the CDS for seqFISH+ and STARmap, which used the coding region, and for the double-hybridization RCA design. Raj et al. 2008 tiled the coding region of each smFISH target; the smFISH template takes the whole mRNA, as the Stellaris designer suggests for expanding a set, so that genes sharing coding sequence with paralogs or pseudogenes still reach 48 probes (`source: CDS` for the coding region alone). Templates whose target windows are 30 nt or shorter (smFISH, seqFISH+, MERFISH, HybISS) count off-target genes with the exhaustive Bowtie 1 search, so their genome needs `bowtie` in `align_index`; the smFISH template keeps only probes without an off-target gene, exempts copies of the target (`exempt_copies`), recommends (commented out) an expression table, by transcript when you have one, with `min_fraction_of_target`, and offers `off_target_pool` as a commented-out option for genes with many paralogs; see [Off-target search](#off-target-search).

```bash
uprobe templates list
uprobe templates show hcr_v3 -o protocol.yaml
```

#### ✅ Validate a Protocol
Check a protocol (and optionally the genome config) before running it. Every problem is reported with its location in the YAML and a suggested fix; the exit status is 1 when there are errors.

```bash
uprobe validate -p protocol.yaml
uprobe validate -p protocol.yaml -g genomes.yaml --check-targets --json
```

Probe expressions are checked with the grammar of the safe evaluator that `uprobe run` uses, hybridization blocks with the same rules as the run, and the probes are then built on synthetic target regions, so expressions that would fail at run time are reported too. `-g` adds genome checks (files, aligner indexes required by the attributes, installed aligners), `--check-targets` looks the targets up in the GTF or FASTA, and `--json` prints the findings with a summary for scripts and agents.

#### 🌟 Complete Workflow (Recommended)
To run the entire pipeline from genome index construction to final probe generation in one go:

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --threads 10
```

Useful flags:

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --continue-invalid --raw --threads 10
```

- `--continue-invalid`: a gene target that is not a `gene_name` of the GTF stops the run before anything is designed (status `error`, exit 1), with an error that lists those targets. With `--continue-invalid`, they are skipped and reported instead: logged, listed in `summary.json` (`inputs.invalid_targets` and an `invalid_targets` warning), and the run ends `partial` (or `no_valid_targets` if no target is left). `uprobe validate --check-targets` finds them before a run. Tasks submitted through the web server skip invalid targets.
- `--raw`: also write every candidate with its attributes, before filtering (`<name>_<time>_raw.xlsx`).
- `--threads` (default 10): threads for building the indexes and for every threaded step of the attributes (Bowtie, Bowtie2, Jellyfish) whose own `threads` option is not set; an attribute's `threads` takes precedence. `uprobe post-process` takes `--threads` too; web-server tasks use the server's `task_threads`.

Every run writes `summary.json` to the output directory, also when it stops early or fails. It gives the outcome (`status`: `ok`, `partial`, `no_probes`, `no_valid_targets`, `no_target_sequences`, `no_probes_constructed` or `error`, with a `message`; `partial` means some targets have no probes, or fewer than `equal_space` requested), how many candidates each target kept after every post-processing stage and which stage removed the most, statistics of the summary attributes, the off-target genes found by `mapped_genes` attributes (`off_targets`), the transcripts each gene target's windows came from (`targets.<name>.transcripts`, see `transcripts` below), the output file names, warnings, and the inputs, including the Tm/fold conditions in effect and the genome's Ensembl `source` when it is recorded. The final probe table is written as `<name>_<time>.xlsx` and as plain tab-separated `<name>_<time>_probes.tsv`. Excel limits never fail a run or lose data: a cell longer than Excel's 32,767 characters (e.g. the off-target list of a repeat window) is cut in the `.xlsx` and ends with `...[truncated: N characters; full value in <file>.tsv]`, and a table with more rows than a worksheet holds (1,048,576 with the header) is written only as TSV; the raw table then gets a `<name>_<time>_raw.tsv` too, and `summary.json` lists the affected tables in an `excel_limits` warning (`output_files.raw_tsv`). Web-server result archives include the TSV files. With `--json`, `uprobe run` prints the summary to stdout and sends logs and tool output to stderr.

Exit status of `uprobe run` (`summary.json` and the `--json` output are written before it exits):

| Code | Meaning | `status` |
| --- | --- | --- |
| 0 | Probes were designed | `ok`, `partial` |
| 1 | The run failed | `error` |
| 2 | Invalid command-line usage (missing or bad options) | none; no summary is written |
| 3 | The run finished without any probe | `no_probes`, `no_probes_constructed`, `no_target_sequences`, `no_valid_targets` |

#### 🔧 Step-by-Step Execution
For advanced users who need intermediate results or custom workflows, you can execute each step individually:

```bash
# 1. Build genome index
uprobe build-index -p protocol.yaml -g genomes.yaml -t 10

# 2. Validate target genes against the GTF file
uprobe validate-targets -p protocol.yaml -g genomes.yaml --continue-invalid

# 3. Extract target region sequences
uprobe generate-targets -p protocol.yaml -g genomes.yaml -o ./results --continue-invalid

# 4. Construct initial probes from target sequences
uprobe construct-probes -p protocol.yaml -g genomes.yaml --targets ./results/target_sequences.csv -o ./results

# 5. Post-process probes (requires target + probe columns)
uprobe post-process -p protocol.yaml -g genomes.yaml --probes ./results/constructed_probes_combined.csv -o ./results

# 6. Generate HTML analysis report
uprobe generate-report -p protocol.yaml -g genomes.yaml --probes ./results/probes_*.csv -o ./results
```

#### 🧬 Generate Barcodes

```bash
uprobe generate-barcodes \
  --strategy max_orthogonality \
  --name barcodes \
  --num-barcodes 16 \
  --length 8 \
  --alphabet ACT \
  --output ./results/barcodes
```

This writes `barcodes.csv` (column `sequence`) and `barcodes.txt`.

### 2. Python API (Ideal for Backend Integration)

If you are developing a web backend or data analysis pipeline, we recommend directly using `UProbeAPI`. It returns Pandas DataFrames, making it easy to process further.

```python
import uprobe
from pathlib import Path

# Initialize API
api = uprobe.UProbeAPI(
    protocol_config=Path("protocol.yaml"),
    genomes_config=Path("genomes.yaml"),
    output_dir=Path("./results")
)

# --- Method 1: Complete Workflow ---
df_final = api.run_workflow(threads=10)

# --- Method 2: Step-by-Step Execution ---
api.build_genome_index(threads=10)
api.validate_targets()  # raises on a target not in the GTF; continue_on_invalid=True skips it
df_targets = api.generate_target_seqs()
df_probes = api.construct_probes(df_targets)

import pandas as pd
df_combined = pd.concat([df_targets, df_probes], axis=1)
df_final = api.post_process_probes(df_combined)

# Generate HTML report
reports = api.generate_report(df_final)
print(reports["html_reports"])

# Quick barcode generation
barcodes = api.quick_generate_barcodes(num_barcodes=16, length=8, alphabet="ACT")
```

## Configuration Details

U-Probe relies on two main configuration files to run. Here is a breakdown of how to structure them.

### 1. `genomes.yaml`
This file maps a genome name to its corresponding file paths. It tells U-Probe where to find the reference genome and annotation files, and which aligner indices to build.

```yaml
# Example genomes.yaml
GRCh38:
  fasta: "/path/to/GRCh38.fasta"
  gtf: "/path/to/GRCh38.gtf"
  align_index:
    - bowtie2
    - bowtie     # Bowtie 1 transcriptome index for mapped_genes with aligner: bowtie (exhaustive off-target search)
    - blast
  jellyfish: false  
```

`uprobe build-index` builds the indexes of the aligners in `align_index`; the transcriptome (`transcript.fa`, spliced from the GTF) is indexed for `bowtie2`, `blast` and `bowtie`, which has no genome index.

Entries written by `uprobe genome fetch` also record their origin under `source` (Ensembl division, release, assembly and URLs), which each run copies into `summary.json`; the verified checksums are kept in `genome_source.json` next to the files. The workflow ignores keys it does not use.

### 2. `protocol.yaml`
This is the core configuration file that defines all parameters for a specific probe design run. It is highly customizable.

```yaml
# Example protocol.yaml snippet
name: my_rna_probe_design
genome: GRCh38
targets:
  - CD4

extracts:
  target_region:
    source: mRNA  # genome / mRNA / exon / CDS / UTR
    length: 30
    overlap: 15

# Define how genes map to specific barcodes
encoding:
  CD4:
    BC1: ACGAGCCTTCCA
    BC2: CGGTAATGGACT

# Define the structure of your probes
probes:
  probe_1:
    template: "{part1}{part2}"
    parts:
      part1:
        expr: "rc(target_region[0:20])"
      part2:
        template: "CC{barcode1}TGCGTCTATTT{barcode2}TAGTGGAGCCT"
        parts:
          barcode1:
            expr: "encoding[target]['BC1']"
          barcode2:
            expr: "encoding[target]['BC2']"
  probe_2:
    template: "{part1}AGGCTCCACTA"
    parts:
      part1:
        expr: "rc(target_region[-10:])"

# Optional. If omitted or empty, the CLI auto-generates attributes,
# post_process, and summary based on source mode and probe structure.
attributes:
  target_gc:
    target: target_region
    type: gc_content
  target_tm:
    target: target_region
    type: annealing_temperature
  target_fold:
    target: target_region
    type: fold_score

post_process:
  filters:
    target_tm:
      condition: target_tm >= 50 & target_tm <= 90
  sorts:
    is_ascending:
      - target_gc
      - target_fold   # lower fold_score = less structure

summary:
  report_name: rna_report   # rna_report / dna_report
  attributes:
    - target_gc
    - target_tm
    - target_fold
```

#### Key Sections in `protocol.yaml`:
- **`name` & `genome`**: Basic metadata. The `genome` must match a key in your `genomes.yaml`.
- **`targets`**: A list of target gene names or IDs you want to design probes for.
- **`extracts`**: Parameters for extracting target sequences: the sliding window `length` and `step` (or `overlap`), and the `source` the windows are taken from:
  - `genome`: targets are regions such as `chr1:1000-2000`, tiled along the genomic sequence (the reference, + strand). Regions are 1-based and inclusive, as genome browsers and Ensembl show them: `chr1:1000-2000` is the 1,001 bases from position 1,000 to 2,000, and the first base of a chromosome is 1. `avoid_otp` `target_regions` use the same notation. A window's `sub_region` `s-e` is its position along the region, so it covers the genomic bases `start + s - 1` to `start + e - 1`. Earlier versions read region targets as 0-based (`chr1:1000-2000` was bases 1,001-2,000, while `avoid_otp` already read them 1-based), so a DNA design made from the same protocol now starts one base earlier; a start of 0 is an error.
  - `mRNA`: the whole mature RNA of each target gene, i.e. all its `exon` features: 5' UTR, CDS and 3' UTR of a protein-coding gene, the processed RNA of a non-coding gene. Use it when probes may bind anywhere on the mRNA, as smFISH probe sets do.
  - `exon`: the coding sequence of a protein-coding gene (its `CDS` features, without the UTRs), and all exons only of a gene without CDS. Despite its name it is not the whole mRNA of a protein-coding gene; `uprobe validate` warns and points to `mRNA`.
  - `CDS`: the same as `exon`, under the name that says what it extracts.
  - `UTR`: the UTR features (`UTR` in GENCODE GTFs, `five_prime_utr` and `three_prime_utr` in Ensembl GTFs).

  With `transcripts: all` (the default), each window lies within one exon (or CDS or UTR segment), so no window spans a splice junction. With `exon`, `CDS` and `UTR`, segments shorter than the window are left out, and a target left with none gets no windows (a log warning and `no_candidate_windows` in `summary.json`) while the other targets are designed as usual; with `mRNA`, short exons give no window but still count in the positions along the RNA (`target_start`, see [Probes per target](#probes-per-target)). When `transcripts` selects one transcript of a target, `mRNA`, `exon` and `CDS` tile its spliced sequence instead, so windows also span its exon junctions; see [Windows across exon junctions](#windows-across-exon-junctions).
  - `transcripts` (RNA sources): which transcripts of each gene the windows come from. `all` (the default) takes the exons of all its transcripts, the merged exon model, so windows can lie on exons that only rare isoforms, or retained-intron transcripts, contain. `canonical` takes the transcript the GTF tags `Ensembl_canonical`; `mane` the one tagged `MANE_Select` (human), and for a gene without one (other species, most non-coding genes) its `Ensembl_canonical` transcript, with a log note; a list of transcript ids, e.g. `transcripts: [ENST00000646664]`, takes those (with or without version), and every target needs at least one of them. Ensembl and GENCODE GTFs carry both tags. An unknown id, an id of a gene that is not a target, a target without a listed id, or a gene without the tag stops the run with an error; `uprobe validate --check-targets` reports them beforehand. With one transcript, the windows tile that spliced transcript (`mRNA`) or its spliced CDS (`exon`, `CDS`), across its exon junctions, `target_start` is the position on it, and `exempt_copies` judges copies against it. `summary.json` lists each target's transcripts: `selection`, `chosen_by` (`all`, `MANE_Select`, `Ensembl_canonical` or `listed`) and `ids`.
  - `masked` (`source: genome`): `keep` (the default) upper-cases the windows; `drop` skips every window with a soft-masked (lowercase) or N base, as OligoMiner mines a repeat-masked genome. It needs a soft-masked FASTA (`uprobe genome fetch --soft-masked`, or UCSC `hg38.fa`); the run log and `uprobe validate --check-targets` warn when a target has no lowercase base. The Oligopaint template uses `drop`.
- **`encoding`**: Mapping of specific genes to custom barcodes or identifiers (e.g., assigning `BC1` and `BC2` to `CD4`).
- **`probes`**: The core of the design, powered by a **Directed Acyclic Graph (DAG)** architecture. This allows for complex, modular probe construction where parts and probes can reference each other.
  - **`template`**: Construct sequences using placeholders (e.g., `{part1}{part2}`).
  - **`expr`**: Apply Python-like expressions. You can use built-in functions like `rc()` (reverse complement), slice sequences (`[0:20]`), or fetch from the `encoding` map (`encoding[target]['BC1']`). Expressions are evaluated safely (never with Python `eval`): only string/integer literals, the names `target_region`, `target`, `encoding`, `target_parts`, `probe_index` and probe/part names (e.g. `probe_1.part1`), subscripts and slices, `+` concatenation, `%` (remainder of two integers) and `rc()` are allowed; anything else (other function calls, attribute or method access, comprehensions, lambdas, ...) is rejected with an error naming the probe part.
  - **`probe_index`**: the number of the probe within its target, 0, 1, 2, ... along the target (5' to 3' on an RNA, along the region on DNA), so a part can vary from probe to probe, e.g. rotate readouts with an encoding list: `rc(encoding[target]['readouts'][(probe_index + 1) % 4])` (an encoding value may be a list of sequences). The final probes are numbered after post-processing, the parts that read `probe_index` are then rebuilt and their attributes recomputed, and the final table gets a `probe_index` column; filters and sorts see each candidate numbered in extraction order (`uprobe validate` warns about attributes of such parts). The `merfish` template uses it to put 3 of each gene's 4 readouts on every probe.
  - **DAG References**: Because of the DAG structure, subsequent probes or parts can dynamically reference the sequences of previously defined probes/parts in their expressions.
- **`attributes`**: Define the biochemical or physical properties you want to calculate for specific targets or probe parts. Available attribute types include:
  - `gc_content`: Calculates the GC ratio of the sequence.
  - `annealing_temperature`: Calculates the melting temperature (Tm) using Primer3, under the protocol's `hybridization` conditions. With `tm_method: multi_padlock_design` it uses the formula of the HybISS padlock design code instead (see [Hybridization conditions](#hybridization-conditions)).
  - `fold_score`: The negated minimum free energy, −MFE in kcal/mol (≥ 0), of the sequence folded alone as RNA by ViennaRNA, at `hybridization.temperature_c` when set. Higher means more stable secondary structure, so lower is better: filter with `<=` and sort it in `is_ascending`.
  - `self_match`: Counts 4-mer matches between the sequence and its reverse complement, a measure of self-dimerization; lower is better.
  - `mapped_sites`: Aligns the sequence to the genome (Bowtie2 `--very-sensitive-local -k 100`, both strands) and lists the sites, with their number in `<attribute>_num`; `<attribute>_num == 1` keeps sequences with one site, OligoMiner's Unique Mode. A window across an exon junction has no contiguous genomic site ([Windows across exon junctions](#windows-across-exon-junctions)).
  - `mapped_genes`: Counts the number of unique genes the sequence aligns to in the transcriptome, via Bowtie2 or, exhaustively within `max_mismatches` substitutions, via Bowtie 1. See [Off-target search](#off-target-search); with `exempt_copies`, genes that copy most of the target are not counted ([Copies of the target](#copies-of-the-target)).
  - `kmer_count`: Counts k-mer occurrences in the genome (via Jellyfish) to evaluate specificity: the highest genome count among the sequence's k-mers, so lower is better. Each `kmer_len` (default 35) has its own index, `jf_genome/<fasta stem>.k<kmer_len>.jf` next to the FASTA, built at its first use. `strand: both` (the default) counts canonical k-mers, i.e. on both strands; `strand: sense` counts the sequence's k-mers as written on the FASTA's forward strand, as OligoMiner's and PaintSHOP's dictionaries do (index `<fasta stem>.k<kmer_len>.stranded.jf`), and `antisense` their reverse complements. For windows across exon junctions see [Windows across exon junctions](#windows-across-exon-junctions).
- **`hybridization`** (optional): The experimental conditions (temperature, salt, probe concentration, formamide) that `annealing_temperature` and `fold_score` are computed under. See [Hybridization conditions](#hybridization-conditions).
- **`gene_expression`** (optional): An expression table of your sample, by gene or by transcript; `mapped_genes` attributes then ignore off-target genes that are not expressed, or optionally expressed far below the target. See [Expression-aware off-target counting](#expression-aware-off-target-counting).
- **`post_process`**: Define strict `filters` (e.g., Tm ranges) based on the calculated attributes, and `sorts` to rank the best probes (ascending or descending). Lower is better for `fold_score`, `self_match`, `kmer_count` and `mapped_genes`: list them under `is_ascending` (`uprobe validate` warns when they are sorted descending). When `post_process` is generated, every sort is ascending: GC, Tm, mapped genes, then the fold, self-match and k-mer scores. Its specificity filter is `target_mapped_genes <= 4` in RNA mode and, in DNA mode (`source: genome`), `target_kmer_count <= 100` on the generated `target_kmer_count` attribute: the highest genome count of the 35-mers of the genomic window (Jellyfish; the genome needs `jellyfish: true`).
- **`remove_overlap`**: Optional post-processing step for controlling spacing between probes. Put it under `post_process` when you need it, e.g. `post_process: { remove_overlap: { location_interval: 0 } }`. It visits each target's candidates in order of position and keeps one when it starts more than `location_interval` bases after the end of the last kept one. In RNA mode the positions are those along the target (`target_start`/`target_end`: on the gene's merged exon model, or on its one selected transcript; see [Probes per target](#probes-per-target)), so windows from alternative exons of different transcripts are compared too, and a site that several exon entries contain is kept once, with the transcripts of all of them in `transcript_names` (and their number in `n_trans`).
- **`off_target_pool`**: Optional post-processing step that lets at most `max_probes_per_gene` probes of a target's set hit the same off-target gene, e.g. `post_process: { off_target_pool: { attribute: target_mapped_genes, max_probes_per_gene: 2 } }`. See [Off-target pool limit](#off-target-pool-limit).
- **`equal_space`**: Optional post-processing step that keeps `number_desired` probes per target, spread evenly along the target, e.g. `post_process: { equal_space: { number_desired: 48 } }`, or per target: `equal_space: { GAPDH: { number_desired: 48 } }`. See [Probes per target](#probes-per-target).
- **`summary`**: Define the `report_name` (e.g., `rna_report` or `dna_report`) and the specific `attributes` you want to visualize and output in the final report.
- **`meta`** (optional): Free-form documentation that the workflow ignores; the built-in templates use it for citations and parameter sources.

The `post_process` steps always run in this order, whatever their order in the protocol: `filters`, `avoid_otp` (DNA mode), `remove_overlap`, `off_target_pool`, `equal_space`, `sorts`.

#### Windows across exon junctions

When `extracts.target_region.transcripts` selects exactly one transcript of a target (`canonical`, `mane`, or a list with one of its ids), U-Probe tiles the spliced sequence of that transcript, as one would tile its cDNA: `source: mRNA` joins all its exons, `exon` and `CDS` its CDS pieces (its exons if it has no CDS), 5' to 3', every one however short. Windows therefore also span its exon-exon junctions. They are valid windows of the spliced RNA and bind only the mature RNA, not the pre-mRNA. On Ensembl 113, the MANE Select mRNA of ACTB (ENST00000646664, 6 exons) gives 1,793 20-nt windows instead of the 1,698 that lie within one exon, exactly the 20-mers of its Ensembl cDNA; ACTB, EEF1A1, GAPDH, HSPA8, RPL13A and TUBB give 12,401 instead of 11,679 (19 more per junction). `transcripts: all` (the merged exon model), several selected transcripts of a target, and `source: UTR` keep one window per exon or segment.

- `target_start`/`target_end`: the window's position on the spliced transcript (or spliced CDS), which `remove_overlap`, `equal_space`, `off_target_pool` and `exempt_copies` use as for any window.
- `exon_name`: the exon (or CDS piece) the window lies on, `chrom_start_end_strand`, or, for a window across junctions, the exons it spans, 5' to 3', joined by `;`, e.g. `7_5530524_5530601_-;7_5529535_5529663_-`.
- `start`/`end` and `sub_region`: the window's position along the exons that `exon_name` names, joined; for a window within one exon, its position on that exon, as for every other window. ACTB's window at 70-89, on the exons above, covers the last 9 nt of the 78-nt first exon and the first 11 of the second. `exon_name` and `start`/`end` thus place every window on the genome.
- `transcript_names`: the transcript; `n_trans`: 1.
- Attributes: `mapped_genes` and `exempt_copies` search the spliced transcripts (`transcript.fa`), so they find a junction window on its own transcript like any other window; the sequence attributes (Tm, GC, folding, self-match) read the window itself. `mapped_sites` and `kmer_count` search the genome, where an intron interrupts a junction window: `mapped_sites` finds the window's own site in part (a local alignment) or not at all, and `kmer_count` finds its k-mers that cross the junction only in processed copies of the gene, so their values for junction windows do not compare with those of other windows. Use `mapped_genes` to judge the specificity of RNA probes.

#### Probes per target

`equal_space` works in DNA (`source: genome`) and RNA mode (`mRNA`, `exon`, `CDS`, `UTR`). It runs after the filters, `remove_overlap` and `off_target_pool`, orders each target's remaining candidates along the target and keeps every k-th one (evenly spaced ranks, first and last included), so the kept probes follow the candidates from one end of the target to the other. Use it together with `remove_overlap`; otherwise overlapping windows can be kept (`uprobe validate` warns).

- **Along the target**: in DNA mode, along the genomic region. In RNA mode, windows number `start`/`end` within their exon or UTR, so U-Probe also writes `target_start`/`target_end`: the window's position along the target. When `transcripts` selects one transcript of the target, that is the position on the spliced transcript (`mRNA`) or its spliced CDS (`exon`, `CDS`); see [Windows across exon junctions](#windows-across-exon-junctions). Otherwise it is the position on the gene's merged exon model, the union of the extracted exons (CDS or UTRs) of its transcripts with the introns removed, read 5' to 3': for a gene with one transcript and `source: mRNA` the position on the spliced RNA, while exons of different transcripts that overlap share positions, and with `exon`, `CDS` and `UTR`, segments shorter than the window length are not extracted and do not count. Direct sequences use their own coordinates.
- **Fewer candidates than requested**: the target keeps all of them. `summary.json` lists `requested_probes` and `final` for every target and a `fewer_probes_than_requested` warning with the counts, the run status is `partial`, and `uprobe run` logs the per-target counts as `final/requested`.

#### Hybridization conditions

By default, Tm is computed under PCR-like conditions (50 mM Na+, 0.6 mM dNTP, 50 nM oligo, no formamide; these match the fisheye reference pipeline). FISH hybridizations usually run in SSC buffer with formamide, so declare the conditions of your experiment with a top-level `hybridization` block:

```yaml
hybridization:
  temperature_c: 37     # hybridization temperature (°C); used by fold_score and auto-generated Tm filters
  na_mM: 390            # monovalent cations; 2x SSC (0.3 M NaCl + 0.03 M Na3-citrate) = 390 mM Na+
  mg_mM: 0              # divalent cations
  dntp_mM: 0            # dNTPs (chelate Mg2+; no effect when mg_mM is 0)
  probe_nM: 50          # concentration of each probe oligo
  formamide_pct: 30     # formamide, % (v/v)

attributes:
  probe_1_tm:
    target: probe_1.part1
    type: annealing_temperature          # uses the block above
  readout_tm:
    target: probe_2
    type: annealing_temperature
    hybridization:                       # per-attribute override, merged over the block
      formamide_pct: 10
```

| Field | Unit | Default when unset | Used by |
| --- | --- | --- | --- |
| `temperature_c` | °C | 37 (ViennaRNA default) | `fold_score`, auto-generated Tm filters |
| `na_mM` | mM | 50 | `annealing_temperature` |
| `mg_mM` | mM | 0 | `annealing_temperature` |
| `dntp_mM` | mM | 0.6 | `annealing_temperature` |
| `probe_nM` | nM | 50 | `annealing_temperature` |
| `formamide_pct` | % (v/v) | 0 | `annealing_temperature` |

- Every field is optional. Without the block (or with an empty one) results are identical to earlier versions. An `annealing_temperature` or `fold_score` attribute may carry its own `hybridization` mapping with the same fields, overriding the top-level values for that attribute only.
- Values must be numbers and `>= 0`; `na_mM` and `probe_nM` must be `> 0`, `formamide_pct` and `temperature_c` at most 100, and unknown field names are rejected.
- **Tm model**: Primer3 nearest-neighbor Tm (SantaLucia 1998, PNAS 95:1460-65) with SantaLucia salt correction; Mg2+ is converted to a Na+ equivalent as in von Ahsen et al. 2001 (Clin Chem 47:1956-61). Primer3 parameters describe DNA/DNA duplexes, so Tm against RNA targets is an approximation.
- **`tm_method: multi_padlock_design`** (on an `annealing_temperature` attribute): the Tm formula of Moldia/multi_padlock_design (`lib/screenseq.py`), the code that designed the HybISS padlock probes (Gyllborg et al. 2020, Nucleic Acids Res 48:e112, [doi:10.1093/nar/gkaa792](https://doi.org/10.1093/nar/gkaa792)): Sugimoto et al. 1996 DNA/DNA nearest-neighbor parameters (Nucleic Acids Res 24:4501-4505, [doi:10.1093/nar/24.22.4501](https://doi.org/10.1093/nar/24.22.4501)) with one initiation term, the probe concentration as strand concentration, + 16.6 log10([Na+] + 3.795 sqrt([Mg2+])) (mol/l) and - 0.72 °C per % formamide; `dntp_mM` is not used. It reproduces the repository's values (to 1e-12 °C) under its fixed conditions, 75 mM Na+, 10 mM Mg2+, 100 nM probe and 20% formamide, and lets the `padlock_iss` template apply the published 65-75 °C window: on the 1,195 published HybISS targets Primer3 gives 3.4 ± 1.5 °C less, and only 327 of them fall within 65-75 °C by Primer3.
- **Formamide**: U-Probe uses Primer3's correction, which lowers Tm by (2.88 − 0.453 × GC fraction) °C per mol/l formamide (Blake & Delcourt 1996, Nucleic Acids Res 24:2095-2103, [doi:10.1093/nar/24.11.2095](https://doi.org/10.1093/nar/24.11.2095)). `formamide_pct` is converted to mol/l with formamide's density (1.1334 g/ml at 20 °C) and molar mass (45.04 g/mol), so 1% (v/v) = 0.252 M and Tm drops 0.61-0.72 °C per percent depending on GC content, close to the linear 0.60 °C/% measured by Hutton 1977 (Nucleic Acids Res 4:3537-3555, [doi:10.1093/nar/4.10.3537](https://doi.org/10.1093/nar/4.10.3537)). It is preferred over a fixed linear factor because it accounts for GC content and is applied inside the same Primer3 call.
- **Folding**: `fold_score` folds at `temperature_c`; ViennaRNA does not model formamide, and salt is left at its default.
- A per-attribute override suits a published Tm window that was computed under other conditions than the experiment's. The MERFISH template, for example, hybridizes in 2x SSC with 30% formamide (top-level block) but computes `binding_tm` at 300 mM Na+, 5 nM probe and no formamide, the conditions of its 66-76 °C window.
- `uprobe validate` reports an invalid block or override with its path. The conditions in effect are logged for every Tm/fold attribute and recorded in `summary.json` (`inputs.hybridization_conditions`). When a protocol declares them, they are also shown in the HTML report header and stored in the result `.xlsx` as custom document properties (`hybridization:<attribute>`).
- **Auto-generated Tm filters** (when `post_process` is left out) follow the hybridization temperature: each Tm attribute keeps [`temperature_c`, `temperature_c` + 40] °C, the 40 °C width of the legacy window, since Tm under declared FISH conditions (e.g. 50% formamide) can fall tens of degrees below the legacy 50-90 °C. An attribute's own `temperature_c` wins; an override that changes salt, probe or formamide without its own `temperature_c` (such as the MERFISH design conditions above) keeps 50-90 °C. Conditions declared without `temperature_c` use the default 37 °C (37-77 °C), with a warning from `uprobe run` and `uprobe validate`. Without the block (or with an empty one) the filters stay at 50-90 °C. The chosen windows are logged.

#### Off-target search

A `mapped_genes` attribute aligns a sequence to the transcriptome and counts the distinct genes it hits, the target gene included, so `target_mapped_genes <= 1` keeps sequences that hit no other gene. Choose the search with `aligner`:

| `aligner` | Search | Needs |
| --- | --- | --- |
| `bowtie2` | Bowtie2 `--very-sensitive-local -k 100`: local alignments, which may cover only part of the sequence, found by a seed heuristic | `bowtie2` in `align_index` |
| `bowtie` | Bowtie 1 `-v <max_mismatches> -a`: every end-to-end site with at most `max_mismatches` substitutions (0-3, default 2); indels are not considered | `bowtie` in `align_index` and Bowtie 1 installed |

```yaml
attributes:
  target_mapped_genes:
    target: target_region
    type: mapped_genes
    aligner: bowtie
    max_mismatches: 2
    strand: sense
```

- `strand` (both aligners): `both` (default) counts sites on either transcript strand; `sense` only sites where the sequence itself occurs in a transcript, right for `target_region`, whose probes bind exactly there; `antisense` only sites of its reverse complement, right for a probe's binding sequence. A probe does not bind an RNA that contains its own sequence, so `both` also counts antisense transcripts the probe cannot bind.
- Bowtie2's local mode requires a score of at least 20 + 8 ln(length), more than a perfect match scores (2 per base) below 22 nt: shorter sequences never align and count 0 genes, so a filter on them keeps every candidate. `uprobe run` logs a warning and `uprobe validate` reports it; use `aligner: bowtie` for smFISH-length oligos. `uprobe validate` also warns when Bowtie2 aligns sequences of 22-30 nt, whose off-target sites it finds only in part (see below).
- On the Ensembl 113 human transcriptome (386,885 primary-assembly transcripts), for all 15,295 18-22-mers of the ACTB and GAPDH MANE Select transcripts: `aligner: bowtie` (`max_mismatches: 2`) reported exactly the sites a brute-force scan of every transcript window found, for a random sample of 1,000 of them (147,426 sites); `bowtie2` found 1.5% of the off-target genes with a site within 2 substitutions (none for 18-21-mers; for 30- and 40-mers it found 56% and 82%). The Bowtie 1 search took 4.6 s for the 15,295 sequences with 8 threads (2.5 s with `strand: sense`; Bowtie2 1.6 s); its index took about 5 minutes to build and 650 MB of disk.
- Per window, which is what a filter such as `target_mapped_genes <= 1` sees: for 1,500 random windows of each length from 400 random MANE Select mRNAs, searched in the 410,909 Ensembl 113 human cDNA and ncRNA transcripts (patches included, genes counted by name) with `strand: sense`, Bowtie2 flagged none of the 20-mers, 69% of the 25-mers, 88% of the 28-mers, 87% of the 30-mers and 97-99% of the 40-, 46- and 52-mers that have a site within 2 substitutions in another gene, listing on average 0%, 49%, 71%, 70% and 87-94% of those genes. It reports at most 100 alignments, so it misses most genes of windows in repeats; from 40 nt on it also reports partial matches (about 25 nt or more), indels and sites with 3 or more mismatches, which a search within 2 substitutions cannot see (3-5% of the windows had only such hits). The templates therefore use `aligner: bowtie` for windows of 30 nt or less (smFISH, seqFISH+, MERFISH, HybISS) and Bowtie2 for the 40-52-nt windows of MiP-seq, STARmap, the double-hybridization RCA design and HCR v3, all with `strand: sense`.
- Protocols without `attributes`: in RNA mode the generated `target_mapped_genes` uses `aligner: bowtie` (`max_mismatches: 2`, `strand: sense`) when the genome's `align_index` includes `bowtie` and Bowtie 1 is installed (the run builds the index if it is missing), and `bowtie2` otherwise; the log says which and why, with a warning when the windows are shorter than 22 nt. `uprobe validate -g` makes the same choice.
- Each `mapped_genes` attribute also writes `<attribute>_off_targets` to the result tables (`.xlsx` and `_probes.tsv`, next to the attribute): the genes other than the row's target that the sequence hits, each with its fewest mismatches, e.g. `ACTG1:0, POTEF:2` (names only with `bowtie2`). Genes are named after the GTF's `gene_name`; a name that several gene ids share is written `name(gene_id)`. Genes that an expression table rules out go to `<attribute>_unexpressed`, copies of the target to `<attribute>_copies`, and with a table of transcripts `<attribute>_hit_transcripts` gives the transcripts each off-target gene was valued by (see below).

#### Expression-aware off-target counting

A probe can only bind RNAs that are present. Supply an expression table of your sample, and `mapped_genes` attributes no longer count off-target genes that are not expressed, or are expressed far below the target:

```yaml
gene_expression:
  table: hek293_tpm.tsv   # your table: tab-separated (comma-separated when the name contains .csv), may be gzipped
  level: gene             # default; transcript for a table of transcripts (recommended, see below)
  gene_column: gene_id    # Ensembl gene ids (with or without version) or gene names; default: the first column
  value_column: TPM       # default: the only numeric column
  min_value: 1            # off-target genes with a lower value are not counted
  min_fraction_of_target: 0.05   # optional: nor those below 5% of the target's value
  absent: count           # default; ignore: off-target genes missing from the table are not counted
```

- Genes are looked up by Ensembl gene id, then by gene name; `gene_column` may list several columns, e.g. `[gene_id, gene_name]`. A gene listed more than once counts with its highest value. Leading `#` comment lines are skipped, and GTEx `.gct` files can be used directly (choose the tissue with `value_column`).
- An off-target gene below the threshold is not counted in `<attribute>` and is listed in `<attribute>_unexpressed` instead of `<attribute>_off_targets`. The threshold is `min_value`, or with `min_fraction_of_target` the higher of the two below. The target gene always counts, so `<attribute> <= 1` still means "no (expressed) off-target gene".
- `min_fraction_of_target` (optional, 0-1) makes the threshold relative to the target: an off-target gene of a target counts only if its value is at least max(`min_value`, `min_fraction_of_target` x the target's value). A spot needs many probes on one RNA, and an RNA at a few percent of the target's level adds at most a few percent of false spots, so low-level pseudogenes and paralogs of a highly expressed target no longer remove its probes, while a gene expressed at a similar level still does. The target's value is read from the same table (the highest of its gene ids). A target that the table does not list keeps `min_value` alone: the log warns, and `summary.json` lists it in `off_targets.gene_expression.targets_not_in_table` (warning `targets_not_in_expression_table`), next to the value and threshold of every listed target (`off_targets.gene_expression.targets`). Per attribute, `below_target_fraction` counts the off-target genes that only the fraction removed, for at least one target.
- Off-target genes that the table does not list count as expressed (`absent: count`, the conservative default). With `absent: ignore` they are not counted and are listed in `<attribute>_unexpressed`, for a table you trust to list every gene the sample expresses. Either way the log and `summary.json` (`off_targets`: `missing_from_table` per attribute, with a warning `off_target_genes_not_in_expression_table`) report how many there were, with examples, and which of them the final probes hit. A table that matches no gene of the annotation stops the run: its identifiers or species do not fit the genome.
- The table applies to every `mapped_genes` attribute. Relative paths are resolved from the directory where uprobe runs. `uprobe validate` checks the block, the file and its columns.

##### A table of transcripts

Use a table of transcripts when you have one, e.g. salmon's `quant.sf`, kallisto's `abundance.tsv` or RSEM's `isoforms.results`. False signal comes only from RNA molecules that carry a probe's site, so an off-target gene's value for a sequence is then the sum of the values of the gene's transcripts that the sequence hits, not the whole gene's. A gene expressed on isoforms without the site no longer counts.

```yaml
gene_expression:
  table: quant.sf
  level: transcript
  transcript_column: Name   # default: the first column
  value_column: TPM
  min_value: 1
  min_fraction_of_target: 0.05
```

- Transcripts are looked up by id with or without version (`ENSMUST00000073605.15` matches `ENSMUST00000073605`; a GENCODE FASTA header such as `ENST...|ENSG...|...` by its first field). The genome's GTF, from which the searched transcriptome is built, assigns them to genes, so `gene_column` does not apply.
- The target's value, for `min_fraction_of_target`, is the sum over all its transcripts in the table.
- A hit transcript missing from the table makes its gene count (`absent: count`) or adds nothing to the gene's value (`absent: ignore`). `missing_from_table` then counts the genes with such transcripts, and `missing_transcripts` the transcripts.
- `<attribute>_hit_transcripts` lists, for each off-target gene the table valued, the transcripts the sequence hits with their values (`NA` when not listed), e.g. `Mir100hg: ENSMUST00000212079=0.014, ENSMUST00000362798=0`. `<attribute>_off_targets` and `off_target_pool` still work by gene, and `summary.json` reports the level (`off_targets.gene_expression.level`).
- Both aligners report the transcripts they hit. Bowtie2 reports at most 100 alignments per sequence, so for genes with many isoforms some hit transcripts may be missing from the sum.
- A table whose identifiers fit the other level (gene ids with `level: transcript`, or transcript ids with `level: gene`) stops the run with a hint, and `uprobe validate` reports it.

Measured with template runs of the smFISH template as in [Copies of the target](#copies-of-the-target) (reference mRNA, 20-mers, Bowtie 1 within 2 substitutions, `exempt_copies`, `target_mapped_genes <= 1`, 48 probes per gene) with TPM >= 1 and `min_fraction_of_target: 0.05`. The tables are salmon quantifications of HeLa-S3 and NIH/3T3 on the Ensembl 113 transcriptomes, by gene and by transcript; a gene's TPM is the sum of its transcripts' TPM, so both tables come from one quantification. Final probes (genes with 48), candidate windows without a counted off-target gene, and off-target genes counted for at least one window:

| Table | Human ACTB, EEF1A1, GAPDH, HSPA8, RPL13A, TUBB | Human, 30 random genes | Mouse Actb, Eef1a1, Gapdh, Hspa8, Rpl13a, Rpl3, Tubb5 | Mouse, 30 random genes |
| --- | --- | --- | --- | --- |
| By gene | 288 (6); 9,314; 2,777 | 1,337 (26); 64,631; 12,894 | 313 (5); 9,015; 1,840 | 1,426 (28); 84,142; 11,245 |
| By transcript | 288 (6); 10,233; 1,078 | 1,350 (26); 74,107; 11,589 | 335 (6); 9,984; 724 | 1,427 (28); 91,248; 9,999 |

- Mouse Gapdh went from 26 to 48 probes. Its 20-mers hit two lncRNA genes mostly on isoforms that NIH/3T3 hardly expresses. Of the 689 windows that hit ENSMUSG00000141751 (90 TPM), 495 hit only isoforms that sum to 1.5 TPM or less; the 194 that also hit its 89-TPM isoform still count. Mir100hg (129 TPM over 169 isoforms) is hit by 272 windows, all on isoforms at 0.11 TPM or less.
- Elsewhere the table of transcripts counted 61% fewer off-target genes for the housekeeping genes and 10-11% fewer for the random genes, and 8-15% more windows passed, but few final probes changed, since most genes reach 48 either way. Among the random genes, TEX54 went from 13 to 23 probes.

#### Copies of the target

A processed pseudogene, a retrocopy, or a gene that overlaps the target (a readthrough gene, a duplicated locus) can hold most of the target's sequence. Every probe of the target then hits it too, so no probe set avoids it, and counting it as an off-target leaves such a target with few or no probes, although probe sets for it work in practice. With `exempt_copies`, a `mapped_genes` attribute does not count these genes:

```yaml
attributes:
  target_mapped_genes:
    target: target_region
    type: mapped_genes
    aligner: bowtie
    max_mismatches: 2
    strand: sense
    exempt_copies:      # or `exempt_copies: true` for these defaults
      kmer: 30          # k-mer length, 12-32
      min_shared: 0.5   # share of the target's k-mers that one transcript of the gene must hold
```

- An off-target gene is a copy of a target when one of its transcripts contains, exactly and in the same orientation, at least `min_shared` of the distinct `kmer`-mers of the target's sequence. A copy is not counted in `<attribute>` and is listed in `<attribute>_copies` (name and fewest mismatches, as in `_off_targets`) instead of `<attribute>_off_targets`, so `off_target_pool` does not limit it either. It is exempt whatever its expression; `gene_expression` applies to the other off-target genes. Paralogs share few exact 30-mers and stay counted.
- The target's sequence is what its candidate windows cover, placed by `target_start`/`target_end`: in RNA mode its merged exon model (the extracted segments of its transcripts, all of them or those `transcripts` selects, with the introns removed: the coordinates `remove_overlap` and `equal_space` use), so the share is that of the sequence the probes are drawn from, e.g. of the MANE Select mRNA with `source: mRNA` and `transcripts: mane`; for a direct sequence, the sequence itself. With `source: mRNA` that is the whole mRNA; with `exon` or `CDS`, U-Probe extracts the CDS of protein-coding genes, so their copies are judged on the coding sequence. Bases that no window covers are unknown, so keep `step` at most the window length (`uprobe validate` warns).
- The transcripts of the genes the alignment hit are read from `transcript.fa`, the FASTA the transcriptome index is built from (extracted again if it is missing), in one pass, with the k-mers packed into 64-bit integers. Checking every gene of the Ensembl 113 human transcriptome (410,909 transcripts, 597 Mb) against 4 targets took 9 s; in the runs below the check added 3-13 s to a run.
- It needs no expression table, since a copy is a matter of sequence, and works with either aligner, so it can be enabled on its own.
- The log and `summary.json` name every exempt gene: per attribute, `copy_genes`, `copies` (each target's copies with their shares) and `final_copy_genes` (those the final probes hit). A target whose windows give no k-mer keeps all its off-target genes and is reported (`targets_without_sequence`, warning `copies_not_checked`).
- Near-identical gene families cannot be told apart by any probe: SPDYE8-15 hold 84-89% of the 30-mers of SPDYE17, so with `exempt_copies` the SPDYE17 probes count them as copies and detect the whole family, as `<attribute>_copies` shows.

Shares of the target's 30-mers on Ensembl 113, for the canonical mRNA and for the merged CDS that `source: exon` extracts; at `min_shared: 0.5` the paralog and the human GAPDH pseudogenes stay counted:

| Target | Off-target gene | Canonical mRNA | Merged CDS |
| --- | --- | --- | --- |
| ACTB | ACTG1 (paralog) | 4.8% | 7.4% |
| GAPDH | GAPDHP1 (the closest of its ~70 pseudogenes) | 23.5% | 29.9% |
| EEF1A1 | EEF1A1P5 (processed pseudogene) | 35.6% | 81.2% |
| RPL13A | RPL13AP5 (processed pseudogene) | 44.9% | 64.9% |
| mouse Gapdh | Gm13394 (retrocopy, 668 TPM in NIH/3T3) | 75.6% | |

Measured with template runs of the smFISH template on Ensembl 113 (GRCh38 and GRCm39 primary assemblies with their GTFs): windows from the whole reference mRNA, across its exon junctions (`source: mRNA`, `transcripts: mane`: MANE Select in human, Ensembl canonical in mouse), 20-mers, Bowtie 1 within 2 substitutions, `target_mapped_genes <= 1`, `remove_overlap` 2, 48 probes per gene, with HeLa-S3 and NIH/3T3 tables (TPM, `min_value: 1`). Final probes (genes with 48); each row adds its option to the row above, except `exempt_copies` alone, which is the template's default:

| Setting | Human ACTB, EEF1A1, GAPDH, HSPA8, RPL13A, TUBB | Human, 30 random genes | Mouse Actb, Eef1a1, Gapdh, Hspa8, Rpl13a, Rpl3, Tubb5 | Mouse, 30 random genes |
| --- | --- | --- | --- | --- |
| No expression table | 119 (1) | 1,198 (24) | 120 (1) | 1,288 (23) |
| `exempt_copies` alone (template default) | 119 (1) | 1,231 (24) | 125 (1) | 1,300 (23) |
| Table, TPM >= 1 | 180 (2) | 1,255 (25) | 125 (1) | 1,394 (27) |
| + `min_fraction_of_target: 0.05` | 193 (2) | 1,285 (25) | 203 (2) | 1,394 (27) |
| + `exempt_copies` (recommended) | 193 (2) | 1,333 (26) | 275 (4) | 1,394 (27) |
| + `absent: ignore` | 264 (5) | 1,337 (26) | 335 (6) | 1,426 (28) |

- The relative threshold helps highly expressed genes whose pseudogenes are expressed at a few percent of their level: human GAPDH 7 -> 14 probes and ACTB 38 -> 44, mouse Eef1a1 17 -> 48 and Rpl13a 17 -> 46. It changed little for random genes.
- Copies matter where a copy is expressed: mouse Rpl3 (Rpl3-ps1 holds 85% of its 30-mers) went from 7 to 48 probes and Hspa8 from 18 to 48 with both options; among the random genes the copies were lncRNA genes overlapping NDUFV2 and EIF3L on the same strand, the SPDYE family, and mouse Vmn2r92 for Vmn2r91 (73%). Human GAPDH has no copy above 30%; it reached 48 probes only when, on top of the fraction, its pseudogenes missing from the HeLa table were ignored. Mouse Gapdh has 35 retrocopies holding 50-92% of its 30-mers, all exempt, but others holding less are missing from the NIH/3T3 table and kept it at 6 probes until they were ignored. `absent: ignore` trusts the table to list every expressed gene.

Where the windows come from matters as much. Up to U-Probe 0.2 the template used `source: exon`, the CDS of all coding isoforms; `mRNA` with `transcripts: all` takes the exons of every isoform, retained introns included. Probes (genes with 48) and the share of final probes whose site lies on the reference mRNA, without an expression table (`exempt_copies` on) and with the recommended settings (TPM >= 1, `min_fraction_of_target: 0.05`, `exempt_copies`):

| Windows | Expression table | Human ACTB, EEF1A1, GAPDH, HSPA8, RPL13A, TUBB | Human, 30 random genes | Mouse Actb, Eef1a1, Gapdh, Hspa8, Rpl13a, Rpl3, Tubb5 | Mouse, 30 random genes |
| --- | --- | --- | --- | --- | --- |
| `exon` (CDS), all transcripts (template up to 0.2) | none | 11 (0), 45% | 1,021 (13), 95% | 38 (0), 68% | 1,082 (17), 98% |
| `exon` (CDS), all transcripts (template up to 0.2) | recommended | 125 (2), 94% | 1,119 (16), 95% | 222 (3), 95% | 1,209 (19), 99% |
| `mRNA`, all transcripts | none | 288 (6), 16% | 1,334 (27), 68% | 332 (5), 27% | 1,337 (26), 76% |
| `mRNA`, all transcripts | recommended | 288 (6), 31% | 1,376 (28), 72% | 336 (7), 38% | 1,400 (28), 77% |
| `mRNA`, `transcripts: mane` (template) | none | 119 (1), 100% | 1,231 (24), 100% | 125 (1), 100% | 1,300 (23), 100% |
| `mRNA`, `transcripts: mane` (template) | recommended | 193 (2), 100% | 1,333 (26), 100% | 275 (4), 100% | 1,394 (27), 100% |

The reference mRNA gave the six human housekeeping genes 119 probes against 11 from the CDS of all isoforms without a table, and 193 against 125 with it; paralogs and pseudogenes share their coding sequence more than their UTRs. All transcripts reach 48 probes for nearly every housekeeping gene, but mostly on exons of other isoforms and retained introns: only 16-38% of these probes lie on the reference mRNA.

#### Off-target pool limit

In smFISH a false spot needs many probes of the set on the same off-target RNA, so a few probes that hit one off-target gene do no harm, while many do. `off_target_pool` limits how many probes of one target's set may hit the same off-target gene:

```yaml
post_process:
  filters:
    target_mapped_genes:
      condition: target_mapped_genes <= 2   # at most one off-target gene per probe
  remove_overlap:
    location_interval: 2
  off_target_pool:
    attribute: target_mapped_genes   # a mapped_genes attribute, or a list of them
    max_probes_per_gene: 2
  equal_space:
    number_desired: 48
  sorts:
    is_ascending: [target_mapped_genes]
```

- For each target, probes are visited from best to worst rank, the order `sorts` defines (the current order when there are no `sorts`), and a probe is dropped when one of its off-target genes (`<attribute>_off_targets`, all listed attributes together) is already hit by `max_probes_per_gene` kept probes. With a `gene_expression` table only expressed off-target genes are listed, so only they are limited; copies of the target (`exempt_copies`) are not listed either. `max_probes_per_gene: 0` keeps only probes without off-target genes.
- It runs after `remove_overlap`, so the limit applies to probes that can end up in the same set, and before `equal_space`, which then chooses the final N among the probes that respect it (any subset still does). Per-probe filters such as `target_mapped_genes <= 3` still apply first.
- `summary.json` reports the stage in `stages` and, in `off_targets.pool`, how many probes it removed and the genes that reached the limit most often.
- A probe in a repeat can hit hundreds of genes, each of them by no other probe of the set, so the limit alone keeps it; combine the limit with a per-probe filter such as `target_mapped_genes <= 2` (at most one off-target gene per probe).
- `equal_space` picks the final probes by position, not by their off-target genes, so when a target has more candidates than it needs, the limit mostly swaps clean probes for probes with an off-target gene instead of adding probes (measured below). It would pay off with an `equal_space` that prefers probes without off-target genes, which U-Probe does not have yet.
- For padlock/RCA probes a single probe makes signal, so allow no expressed off-target per probe with a filter instead: `target_mapped_genes <= 1` (with `gene_expression`).

Measured on Ensembl 113 human genes with template runs of the smFISH template (the MANE Select mRNAs, `exempt_copies`, off-targets within 2 substitutions by Bowtie 1, `remove_overlap` 2, 48 probes per gene; the genes of [Copies of the target](#copies-of-the-target)): `target_mapped_genes <= 1` alone, against `target_mapped_genes <= 2` with `max_probes_per_gene: 2`.

| Genes | Expression table | `<= 1`: probes (genes with 48) | `<= 2` + pool: probes (genes with 48) | Final probes with an off-target gene, `<= 2` + pool |
| --- | --- | --- | --- | --- |
| 30 random protein-coding genes | none | 1,231 (24) | 1,249 (24) | 47% |
| 30 random protein-coding genes | HeLa, TPM >= 1 | 1,303 (26) | 1,312 (26) | 34% |
| ACTB, EEF1A1, GAPDH, HSPA8, RPL13A, TUBB | none | 119 (1) | 132 (1) | 45% |
| ACTB, EEF1A1, GAPDH, HSPA8, RPL13A, TUBB | HeLa, TPM >= 1 | 180 (2) | 195 (2) | 36% |

The limit added 9-18 probes per gene set while a third to almost half of the final probes came to carry an off-target gene (each such gene hit by at most 2 probes of a set). An expression table gained more: 61 and 72 probes with `<= 1`. The smFISH template therefore keeps `target_mapped_genes <= 1`, recommends `gene_expression` and leaves `off_target_pool` as a commented-out option for genes with many paralogs or pseudogenes.

For more detailed examples and advanced configurations, please refer to the [`tests/data/*.yaml`](https://github.com/UFISH-Team/U-Probe/tree/main/tests/data "Click to visit here") directory and the built-in [protocol templates](#protocol-templates) (`uprobe templates list`).

## Community & Support

- 📖 **Documentation**: [uprobe.readthedocs.io](https://uprobe.readthedocs.io/)
- 💬 **GitHub Discussions**: [Ask questions and share ideas](https://github.com/UFISH-Team/U-Probe/discussions)
- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/UFISH-Team/U-Probe/issues)
- 🚀 **Contributing**: See our [contributing guide](https://uprobe.readthedocs.io/en/latest/contributing.html)

## Citation

If you use U-Probe in your research, please cite our paper (or software):

```bibtex
@software{uprobe2026,
  title={U-Probe: Universal Probe Design Tool},
  author={Zhang, Qian and Xu, Weize and Cai, Huaiyuan},
  year={2025},
  url={https://github.com/UFISH-Team/U-Probe},
  version={1.0.0}
}
```
*(Note: Update the citation format with actual journal details once published.)*

## License

U-Probe is released under the [MIT License](LICENSE). See the `LICENSE` file for details.

## Acknowledgments

We thank the bioinformatics community for valuable feedback during development, and the authors of the following tools that U-Probe integrates:

- [PantheonOS](https://github.com/aristoteleo/PantheonOS) - The multi-agent framework powering the web server's chat agent
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - Fast and memory-efficient sequence alignment
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - Sequence similarity search  
- [MMseqs2](https://github.com/soedinglab/MMseqs2) - Ultra-fast and sensitive sequence search and clustering
- [Jellyfish](https://github.com/gmarcais/Jellyfish) - Fast k-mer counting
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) - RNA secondary structure prediction
- [Primer3](https://primer3.org/) - Primer and probe design algorithms
- [FastAPI](https://fastapi.tiangolo.com/) & [React](https://react.dev/) - Powering our interactive Web UI




