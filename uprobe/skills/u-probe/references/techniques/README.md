# Technique references: index

Read this first when a user names an imaging technology or describes an experiment and you must pick
a probe architecture and a U-Probe template. Each file covers one technology: architecture, published
parameters with DOIs, target region, specificity, the U-Probe template and its departures, pitfalls,
a checklist.

## Decision table

| Experimental goal | Technology | File | U-Probe template |
| --- | --- | --- | --- |
| Count / localize 1-3 RNAs per cell at single-molecule resolution with commercial dye-labelled oligos (cells, yeast, worms, thin sections) | smFISH (Raj 2008 / Stellaris) | `smfish.md` | `smfish_stellaris` |
| Same readout, cheaper and dye-flexible: unlabelled primaries + one labelled FLAP secondary | smiFISH (Tsanov 2016) | `smifish.md` | none; protocol YAML in `smifish.md` (adapted from `smfish_stellaris`) |
| Enzyme-free amplification for embryos, thick or autofluorescent tissue, low-abundance RNA; 1-5 targets (10 with spectral imaging); quantitative or single-molecule | HCR v3.0 split-initiator pairs | `hcr-v3.md` | `hcr_v3` |
| Hundreds to ~10,000 genes per cell, error-robust binary barcodes, cell typing in tissue | MERFISH (encoding + readout probes, MHD4 codebook) | `merfish.md` | `merfish` |
| ~1,000-10,000 genes at sub-cellular density by pseudocolour dilution (cells, thin tissue) | seqFISH+ (primary probes with 4 readout sites) | `seqfish-plus.md` | `seqfish_plus` |
| Targeted in situ sequencing after reverse transcription; RCA, one probe per site gives a spot | Padlock ISS / HybISS | `padlock-iss-hybiss.md` | `padlock_iss` |
| RT-free in situ sequencing in intact thick tissue (hydrogel), padlock + primer pair (SNAIL) | STARmap | `starmap.md` | `starmap` |
| RT-free padlock + RCA initiator primer with pairwise barcodes, sequencing-by-ligation readout | MiP-seq | `mip-seq.md` | `mipseq` |
| RNA-templated circle probe + amplification primer (U-Probe RCA configuration) | Double-hybridization RCA | `double-hyb-rca.md` | `double_hyb_rca` |
| DNA FISH: chromosome paints, loci, chromatin tracing | Oligopaint / OligoMiner | `oligopaint-dna-fish.md` | `oligopaint` |

## Choosing between the RNA methods

| Question | Then | Why |
| --- | --- | --- |
| 1-3 genes, absolute counts, simplest assay? | smFISH or smiFISH | 24-48 unamplified oligos per RNA; no barcodes. |
| RNA shorter than ~600 nt, or sample too bright/thick for unamplified spots? | HCR v3 (or an RCA method) | Amplified signal; 5-20 pairs suffice. |
| More than ~10 genes in the same cells? | MERFISH or seqFISH+ | Combinatorial barcodes, readout probes, codebooks. |
| Genes with very high expression in a multiplex panel? | image them sequentially or drop them | Optical crowding breaks combinatorial decoding (`merfish.md` section 5). |
| Does one bound probe (or pair) make a spot? | HCR, padlock/RCA: yes | Each probe or pair must be specific; filter per probe. |
| | smFISH, smiFISH, MERFISH, seqFISH+: no | A false spot needs many probes on one off-target RNA; a few promiscuous probes cost background. |

## Probe-level numbers at a glance (sources in each file)

| Technology | Binding length | Probes per RNA | Readout / initiator part | Primary hybridization |
| --- | --- | --- | --- | --- |
| smFISH | 20 nt (17-22) | 25-48 | none (3' dye) | 10 % formamide, 2x SSC, 37 C |
| smiFISH | 26-32 nt (20 in some protocols) | 24 (12-48) | 28-nt FLAP site | 15 % formamide, 1x SSC, 37 C |
| HCR v3 | 2 x 25 nt per pair, 2-nt gap | 20 pairs (>= 25-30 for single molecules) | 18-nt initiator half + 2-nt spacer | 30 % formamide, 5x SSC, 37 C |
| MERFISH | 30 nt | 92 (48-92 in later work) | 3 x 20-nt readouts (2 x 30 in Chen 2015) + 2 x 20-nt primers | 30 % formamide, 2x SSC, 37 C |
| seqFISH+ | 28 nt | 24 | 4 x 15-nt readout sites | 40 % formamide, 2x SSC, 37 C |

## Rules for every RNA design

- Probe binding sequence = reverse complement of the RNA (`rc(target_region)`); use `strand: sense`
  in `mapped_genes` on `target_region`. Readout/initiator/FLAP sites = rc of the labelled oligo.
- Design on the reference transcript (`transcripts: mane`) unless the user names isoforms.
- Give every Tm with its conditions (Na+, formamide, probe concentration); a Tm window is meaningless
  without them. U-Probe computes Tm under the protocol's `hybridization` block.
- Windows of 30 nt or less need `aligner: bowtie` (Bowtie 1) for off-targets; add `gene_expression`
  so that unexpressed genes do not remove probes.
- Run `uprobe validate -p protocol.yaml -g genomes.yaml --check-targets` before `uprobe run`; read
  `summary.json` (`partial` = fewer probes than requested for some target).
- Each file lists where its template departs from the papers; tell the user when a departure matters.
