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
- **Agentic AI Workflow**: Integrates hierarchical LLM-based AI agents (via PantheonOS) for conversational design. Users can describe experimental goals in plain language or provide scRNA-seq data, and the agents will autonomously construct configurations, select parameters, and execute the design pipeline.
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

The conda environment also installs the command-line tools U-Probe calls (Bowtie2, Bowtie 1, BLAST+, MMseqs2, Jellyfish). With a pip install, install the ones you use, e.g. `conda install -c bioconda bowtie2 bowtie jellyfish`; Bowtie 1 (`bowtie`) is only needed for the [exhaustive off-target search](#off-target-search).

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

#### 🤖 AI Smart Assistant 
U-Probe includes an interactive AI Agent powered by Pantheon. You can design probes through natural language conversation without writing YAML files manually.

```bash
uprobe agent
```

CLI Agent sessions set `UPROBE_OUTPUT_DIR` automatically. New generated files should be under:

```text
<workspace>/outputs/agent/<user>/<session>/agent_runs/<run>/
```

Use `uprobe agent --force` after upgrading if you need to refresh the installed Pantheon team template.

#### 🌐 Start Web Server 
Run the FastAPI backend in WSL/Linux and the separate React frontend in Windows. First check `config.ini`: use local Linux paths under `[Paths]` and `frontend_url = http://localhost:5173` under `[Server]`.

```bash
cd /path/to/U-Probe
conda activate uprobe  # or your existing Python 3.10+ environment
uprobe server --env development --host 127.0.0.1 --port 8005 --workers 1
```

Development mode reloads Python changes automatically. Alternatively, run `APP_ENV=development PORT=8005 python -m uprobe.http.server` from the project root. Running without overrides uses `config.ini`, which may select production mode.

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
- Re-running skips files that are already present and verified, and resumes interrupted downloads. `--dry-run` prints the resolved URLs without downloading; `--mirror` replaces the Ensembl FTP base URL.
- `--build-index` runs the same index build as `uprobe build-index`.

#### Protocol Templates

U-Probe ships complete protocol templates for published probe technologies: smFISH (Stellaris-style), MERFISH, seqFISH+, HCR v3.0, STARmap, padlock ISS (HybISS), MiP-seq, Oligopaint/OligoMiner and the double-hybridization RCA design. Each template declares the source protocol's [hybridization conditions](#hybridization-conditions); where a paper gives a Tm window together with the conditions it was computed under (MERFISH, Oligopaint, HybISS), the Tm attribute is computed under those conditions and the window is applied as a filter. A `meta` block cites the source papers, maps every parameter to its source and lists what you must replace (genome, targets, barcodes/readouts). Every template passes `uprobe validate`.

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

Probe expressions are checked with the grammar of the safe evaluator that `uprobe run` uses, hybridization blocks with the same rules as the run, and the probes are then built on synthetic target regions, so expressions that would fail at run time are reported too. `-g` adds genome checks (files, aligner indexes required by the attributes), `--check-targets` looks the targets up in the GTF or FASTA, and `--json` prints the findings with a summary for scripts and agents.

#### 🌟 Complete Workflow (Recommended)
To run the entire pipeline from genome index construction to final probe generation in one go:

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --threads 10
```

Useful flags:

```bash
uprobe run -p protocol.yaml -g genomes.yaml -o ./results --continue-invalid --raw --threads 10
```

Every run writes `summary.json` to the output directory, also when it stops early or fails. It gives the outcome (`status`: `ok`, `partial`, `no_probes`, `no_valid_targets`, `no_target_sequences`, `no_probes_constructed` or `error`, with a `message`; `partial` means some targets have no probes, or fewer than `equal_space` requested), how many candidates each target kept after every post-processing stage and which stage removed the most, statistics of the summary attributes, the off-target genes found by `mapped_genes` attributes (`off_targets`), the output file names, warnings, and the inputs, including the Tm/fold conditions in effect and the genome's Ensembl `source` when it is recorded. The final probe table is written as `<name>_<time>.xlsx` and as plain tab-separated `<name>_<time>_probes.tsv`. With `--json`, `uprobe run` prints the summary to stdout and sends logs and tool output to stderr.

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
api.validate_targets()
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
    - bowtie     # optional: Bowtie 1 transcriptome index for mapped_genes with aligner: bowtie
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
    source: exon  # genome / exon / CDS / UTR
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
    is_descending:
      - target_fold

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
- **`extracts`**: Parameters for extracting target sequences. The `source` can be `genome`, `exon`, `CDS`, or `UTR`. You can also define the sliding window `length` and `overlap`.
- **`encoding`**: Mapping of specific genes to custom barcodes or identifiers (e.g., assigning `BC1` and `BC2` to `CD4`).
- **`probes`**: The core of the design, powered by a **Directed Acyclic Graph (DAG)** architecture. This allows for complex, modular probe construction where parts and probes can reference each other.
  - **`template`**: Construct sequences using placeholders (e.g., `{part1}{part2}`).
  - **`expr`**: Apply Python-like expressions. You can use built-in functions like `rc()` (reverse complement), slice sequences (`[0:20]`), or fetch from the `encoding` map (`encoding[target]['BC1']`). Expressions are evaluated safely (never with Python `eval`): only string/integer literals, the names `target_region`, `target`, `encoding`, `target_parts` and probe/part names (e.g. `probe_1.part1`), subscripts and slices, `+` concatenation and `rc()` are allowed; anything else (other function calls, attribute or method access, comprehensions, lambdas, ...) is rejected with an error naming the probe part.
  - **DAG References**: Because of the DAG structure, subsequent probes or parts can dynamically reference the sequences of previously defined probes/parts in their expressions.
- **`attributes`**: Define the biochemical or physical properties you want to calculate for specific targets or probe parts. Available attribute types include:
  - `gc_content`: Calculates the GC ratio of the sequence.
  - `annealing_temperature`: Calculates the melting temperature (Tm) using Primer3, under the protocol's `hybridization` conditions.
  - `fold_score`: Calculates the RNA folding minimum free energy (MFE) using ViennaRNA (lower is more stable), at `hybridization.temperature_c` when set.
  - `self_match`: Calculates the potential for self-dimerization.
  - `mapped_sites`: Aligns the sequence to the genome (via Bowtie2) and counts off-target mapped sites.
  - `mapped_genes`: Counts the number of unique genes the sequence aligns to in the transcriptome, via Bowtie2 or, exhaustively within `max_mismatches` substitutions, via Bowtie 1. See [Off-target search](#off-target-search).
  - `kmer_count`: Counts k-mer occurrences in the genome (via Jellyfish) to evaluate specificity.
- **`hybridization`** (optional): The experimental conditions (temperature, salt, probe concentration, formamide) that `annealing_temperature` and `fold_score` are computed under. See [Hybridization conditions](#hybridization-conditions).
- **`gene_expression`** (optional): An expression table of your sample; `mapped_genes` attributes then ignore off-target genes that are not expressed. See [Expression-aware off-target counting](#expression-aware-off-target-counting).
- **`post_process`**: Define strict `filters` (e.g., Tm ranges) based on the calculated attributes, and `sorts` to rank the best probes (ascending or descending).
- **`remove_overlap`**: Optional post-processing step for controlling spacing between probes. Put it under `post_process` when you need it, e.g. `post_process: { remove_overlap: { location_interval: 0 } }`. It visits each target's candidates in order of position and keeps one when it starts more than `location_interval` bases after the end of the last kept one. In RNA mode the positions are those on the gene's merged exon model (`target_start`/`target_end`, see [Probes per target](#probes-per-target)), so windows from alternative exons of different transcripts are compared too, and a site that several exon entries contain is kept once, with the transcripts of all of them in `transcript_names` (and their number in `n_trans`).
- **`off_target_pool`**: Optional post-processing step that lets at most `max_probes_per_gene` probes of a target's set hit the same off-target gene, e.g. `post_process: { off_target_pool: { attribute: target_mapped_genes, max_probes_per_gene: 2 } }`. See [Off-target pool limit](#off-target-pool-limit).
- **`equal_space`**: Optional post-processing step that keeps `number_desired` probes per target, spread evenly along the target, e.g. `post_process: { equal_space: { number_desired: 48 } }`, or per target: `equal_space: { GAPDH: { number_desired: 48 } }`. See [Probes per target](#probes-per-target).
- **`summary`**: Define the `report_name` (e.g., `rna_report` or `dna_report`) and the specific `attributes` you want to visualize and output in the final report.
- **`meta`** (optional): Free-form documentation that the workflow ignores; the built-in templates use it for citations and parameter sources.

The `post_process` steps always run in this order, whatever their order in the protocol: `filters`, `avoid_otp` (DNA mode), `remove_overlap`, `off_target_pool`, `equal_space`, `sorts`.

#### Probes per target

`equal_space` works in DNA (`source: genome`) and RNA mode (`exon`, `CDS`, `UTR`). It runs after the filters, `remove_overlap` and `off_target_pool`, orders each target's remaining candidates along the target and keeps every k-th one (evenly spaced ranks, first and last included), so the kept probes follow the candidates from one end of the target to the other. Use it together with `remove_overlap`; otherwise overlapping windows can be kept (`uprobe validate` warns).

- **Along the target**: in DNA mode, along the genomic region. In RNA mode, exon, CDS and UTR windows number `start`/`end` within their exon or UTR, so U-Probe also writes `target_start`/`target_end`: the window's position on the gene's merged exon model, the union of the extracted exons (CDS or UTRs) of all its transcripts with the introns removed, read 5' to 3'. For a gene with one transcript this is the position on the spliced RNA; exons of different transcripts that overlap share positions. Segments shorter than the window length are not extracted and do not count. Direct sequences use their own coordinates.
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
- **Formamide**: U-Probe uses Primer3's correction, which lowers Tm by (2.88 − 0.453 × GC fraction) °C per mol/l formamide (Blake & Delcourt 1996, Nucleic Acids Res 24:2095-2103, [doi:10.1093/nar/24.11.2095](https://doi.org/10.1093/nar/24.11.2095)). `formamide_pct` is converted to mol/l with formamide's density (1.1334 g/ml at 20 °C) and molar mass (45.04 g/mol), so 1% (v/v) = 0.252 M and Tm drops 0.61-0.72 °C per percent depending on GC content, close to the linear 0.60 °C/% measured by Hutton 1977 (Nucleic Acids Res 4:3537-3555, [doi:10.1093/nar/4.10.3537](https://doi.org/10.1093/nar/4.10.3537)). It is preferred over a fixed linear factor because it accounts for GC content and is applied inside the same Primer3 call.
- **Folding**: `fold_score` folds at `temperature_c`; ViennaRNA does not model formamide, and salt is left at its default.
- A per-attribute override suits a published Tm window that was computed under other conditions than the experiment's. The MERFISH template, for example, hybridizes in 2x SSC with 30% formamide (top-level block) but computes `binding_tm` at 300 mM Na+, 5 nM probe and no formamide, the conditions of its 66-76 °C window.
- `uprobe validate` reports an invalid block or override with its path. The conditions in effect are logged for every Tm/fold attribute and recorded in `summary.json` (`inputs.hybridization_conditions`). When a protocol declares them, they are also shown in the HTML report header and stored in the result `.xlsx` as custom document properties (`hybridization:<attribute>`).
- **Auto-generated Tm filters** (when `post_process` is left out) follow the hybridization temperature: each Tm attribute keeps [`temperature_c`, `temperature_c` + 40] °C, the 40 °C width of the legacy window, since Tm under declared FISH conditions (e.g. 50% formamide) can fall tens of degrees below the legacy 50-90 °C. An attribute's own `temperature_c` wins; an override that changes salt, probe or formamide without its own `temperature_c` (such as the MERFISH design conditions above) keeps 50-90 °C. Conditions declared without `temperature_c` use the default 37 °C (37-77 °C), with a warning from `uprobe run` and `uprobe validate`. Without the block (or with an empty one) the filters stay at 50-90 °C. The chosen windows are logged.

#### Off-target search

A `mapped_genes` attribute aligns a sequence to the transcriptome and counts the distinct genes it hits, the target gene included, so `target_mapped_genes <= 1` keeps sequences that hit no other gene. Choose the search with `aligner`:

| `aligner` | Search | Needs |
| --- | --- | --- |
| `bowtie2` (templates and generated attributes) | Bowtie2 `--very-sensitive-local -k 100`: local alignments, which may cover only part of the sequence, found by a seed heuristic | `bowtie2` in `align_index` |
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
- Bowtie2's local mode requires a score of at least 20 + 8 ln(length), more than a perfect match scores (2 per base) below 22 nt: shorter sequences never align and count 0 genes, so a filter on them keeps every candidate. `uprobe run` logs a warning and `uprobe validate` reports it; use `aligner: bowtie` for smFISH-length oligos.
- On the Ensembl 113 human transcriptome (386,885 primary-assembly transcripts), for all 15,295 18-22-mers of the ACTB and GAPDH MANE Select transcripts: `aligner: bowtie` (`max_mismatches: 2`) reported exactly the sites a brute-force scan of every transcript window found, for a random sample of 1,000 of them (147,426 sites); `bowtie2` found 1.5% of the off-target genes with a site within 2 substitutions (none for 18-21-mers; for 30- and 40-mers it found 56% and 82%). The Bowtie 1 search took 4.6 s for the 15,295 sequences with 8 threads (2.5 s with `strand: sense`; Bowtie2 1.6 s); its index took about 5 minutes to build and 650 MB of disk.
- Each `mapped_genes` attribute also writes `<attribute>_off_targets` to the result tables (`.xlsx` and `_probes.tsv`, next to the attribute): the genes other than the row's target that the sequence hits, each with its fewest mismatches, e.g. `ACTG1:0, POTEF:2` (names only with `bowtie2`). Genes are named after the GTF's `gene_name`; a name that several gene ids share is written `name(gene_id)`.

#### Expression-aware off-target counting

A probe can only bind RNAs that are present. Supply an expression table of your sample, and `mapped_genes` attributes no longer count off-target genes that are not expressed:

```yaml
gene_expression:
  table: hek293_tpm.tsv   # your table: tab-separated (comma-separated when the name contains .csv), may be gzipped
  gene_column: gene_id    # Ensembl gene ids (with or without version) or gene names; default: the first column
  value_column: TPM       # default: the only numeric column
  min_value: 1            # off-target genes with a lower value are not counted
```

- Genes are looked up by Ensembl gene id, then by gene name; `gene_column` may list several columns, e.g. `[gene_id, gene_name]`. A gene listed more than once counts with its highest value. Leading `#` comment lines are skipped, and GTEx `.gct` files can be used directly (choose the tissue with `value_column`).
- An off-target gene below `min_value` is not counted in `<attribute>` and is listed in `<attribute>_unexpressed` instead of `<attribute>_off_targets`. The target gene always counts, so `<attribute> <= 1` still means "no (expressed) off-target gene".
- Off-target genes that the table does not list count as expressed. The log and `summary.json` (`off_targets`, with a warning `off_target_genes_not_in_expression_table`) report how many there were, with examples, and which of them the final probes hit. A table that matches no gene of the annotation stops the run: its identifiers or species do not fit the genome.
- The table applies to every `mapped_genes` attribute. Relative paths are resolved from the directory where uprobe runs. `uprobe validate` checks the block, the file and its columns.

#### Off-target pool limit

In smFISH a false spot needs many probes of the set on the same off-target RNA, so a few probes that hit one off-target gene do no harm, while many do. `off_target_pool` limits how many probes of one target's set may hit the same off-target gene:

```yaml
post_process:
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

- For each target, probes are visited from best to worst rank, the order `sorts` defines (the current order when there are no `sorts`), and a probe is dropped when one of its off-target genes (`<attribute>_off_targets`, all listed attributes together) is already hit by `max_probes_per_gene` kept probes. With a `gene_expression` table only expressed off-target genes are listed, so only they are limited. `max_probes_per_gene: 0` keeps only probes without off-target genes.
- It runs after `remove_overlap`, so the limit applies to probes that can end up in the same set, and before `equal_space`, which then chooses the final N among the probes that respect it (any subset still does). Per-probe filters such as `target_mapped_genes <= 3` still apply first.
- `summary.json` reports the stage in `stages` and, in `off_targets.pool`, how many probes it removed and the genes that reached the limit most often.
- For padlock/RCA probes a single probe makes signal, so allow no expressed off-target per probe with a filter instead: `target_mapped_genes <= 1` (with `gene_expression`).

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

- [PantheonOS](https://github.com/aristoteleo/PantheonOS) - The multi-agent framework powering our AI agent design
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - Fast and memory-efficient sequence alignment
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - Sequence similarity search  
- [MMseqs2](https://github.com/soedinglab/MMseqs2) - Ultra-fast and sensitive sequence search and clustering
- [Jellyfish](https://github.com/gmarcais/Jellyfish) - Fast k-mer counting
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) - RNA secondary structure prediction
- [Primer3](https://primer3.org/) - Primer and probe design algorithms
- [FastAPI](https://fastapi.tiangolo.com/) & [React](https://react.dev/) - Powering our interactive Web UI




