# Probe Architecture

## Overview

Scope: V01, original 2025 method, including dual-RNA and sgRNA implementations. The live-cell probe is an RNA–protein complex; DNA oligos are manufacturing templates. Fixed-cell DNA FISH and smiFISH validation reagents are separate assay components. Reported: S01, Results and STAR Methods; S02, Table S1; S03, Table S2.

## Canonical Architecture

All diagrams below are **schematic**, not exact nucleotide sequences. RNA order and termini are reported in S02, Table S1(6), PDF p.3 (printed p.2), V01. The mature RNA does not contain the upstream DNA T7 promoter or PCR adapter.

```text
Dual-RNA implementation (two separate strands):
crRNA:    5′ — [target-complementary spacer] — [crRNA common region] — 3′
tracrRNA: 5′ — [dye, in the universal dye-tracrRNA implementation] — [tracrRNA] — 3′
                  crRNA common region pairs with tracrRNA; both associate with dCas9

sgRNA implementation (one strand):
5′ — [target-complementary spacer] — [sgRNA common region/scaffold] — 3′
      terminal dye position depends on the labeling implementation

Manufacturing DNA, Twist pool (coding/sense arrangement inferred from T7 transcription):
5′ — [forward adapter] — [T7 promoter] — [spacer DNA] — [crRNA/sgRNA common DNA] — 3′
IDT crRNA forward oligo omits the separate adapter and uses a shorter common-DNA overlap.
```

RNA 5′→3′ orientation is explicit; DNA sense-strand orientation is inferred from Table S1(1–4) and the transcription workflow, not a strand label printed on every DNA template. crRNA and tracrRNA pair as separate molecules; no covalent crRNA–tracrRNA junction is implied. The spacer is complementary to its bound DNA strand; the published table does not identify genomic strand for each entry. Do not reverse-complement Table S2 blindly. PAM belongs to genomic DNA. S01, design methods; S02, Table S1; V01.

### Published sequence records

Exact fixed-region strings below are reported in S02, Table S1(1–6), V01. DNA primers/adapters retain the printed T alphabet; mature RNA retains U. No reverse complements have been substituted. Lengths are **derived** by counting letters in the displayed exact strings. Fixed means constant within these published constructs, not all future implementations.

| Record | Exact published string, joined across printed line wraps | Molecule / orientation / representation |
|---|---|---|
| crRNA common region | `GUUUAAGAGCUAUGCUGUUUUG` | RNA; 5′→3′; exact fixed region; 22 nt derived |
| tracrRNA | `GGAACCAUUCAAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU` | RNA; explicitly 5′→3′; exact sequence; dye is specified separately |
| sgRNA common region | `GUUUAAGAGCUAUGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU` | RNA; 5′→3′; exact fixed region; no target spacer included |
| T7 promoter | `TAATACGACTCACTATA` | DNA; exact fixed region; orientation inferred as coding/sense arrangement |
| crRNA reverse primer | `CAAAACAGCATAGCTCTTAAAC` | DNA; exact published primer notation; synthesis direction inferred from primer use |
| sgRNA reverse primer | `AAAAAAAGCACCGACTCGGTGCCAC` | DNA; exact published primer notation; synthesis direction inferred from primer use |
| Adapter 1 | `CTTGGTCAGACGAGTGC` | DNA; exact library adapter; orientation inferred from forward-primer use |
| Adapter 2 | `TGGTACGGGAACAGCAC` | DNA; exact library adapter; orientation inferred from forward-primer use |
| Adapter 3 | `CTCACCGCTCTTGTAGC` | DNA; exact library adapter; orientation inferred from forward-primer use |
| Adapter 4 | `AAGACGGCTGAGCCATC` | DNA; exact library adapter; orientation inferred from forward-primer use |

The Table S1 strings containing N are **degenerate/schematic templates**, not exact probes. Its prose specifies 12–19 variable N positions after an initial G in relevant manufacturing constructs and 13–20 nt spacers in RNA pool schematics; the drawn number of N letters must not override those ranges. Repetitive-region examples are separately listed and must not be forced into the non-repetitive pool length rule. S02, Table S1(1–7).

Table S2 supplies exact spacer strings in DNA alphabet, not complete chemically labeled RNAs. For example, `GTTCTGTTCCCCTGTGAACC` is the published spacer at D2 of `Chr3_ACAP2-nonrep_P1-P4 18-20bp`; it is an exact **sequence example**, variable across targets, with genomic strand not supplied in that sheet. No assembled full probe is asserted from it. S03, named sheet, row 2; V01.

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region.crRNA | spacer | DNA recognition | 5′ part of crRNA | RNA 5′→3′; complementary to bound DNA | 13–20 nt for pool schematics | Variable; exact spacers in S03, DNA notation | variable | yes | NGG-associated genomic sites; see design table | reported | S02 Table S1(6); S01 design methods; V01 |
| other_functional_sequence.crRNA | crRNA common region | tracrRNA pairing / guide structure | 3′ of spacer | RNA 5′→3′ | 22 nt, counted | Exact fixed region above | fixed | yes | Published construct includes A–U flip | reported | S02 Table S1(6); S01 RNA preparation; V01 |
| other_functional_sequence.tracrRNA | dCas9 tracrRNAs | Partner RNA in guide complex | Separate strand | Explicit RNA 5′→3′ | See exact sequence; no source length stated | Exact sequence above | fixed | conditional | Required for dual-RNA implementation | reported | S02 Table S1(6); V01 |
| other_functional_sequence.sgRNA | crRNA/sgRNA common region | Single-guide scaffold | 3′ of spacer in sgRNA | Explicit RNA 5′→3′ | See exact sequence; no source length stated | Exact fixed region above | fixed | conditional | Alternative to dual-RNA arrangement | reported | S02 Table S1(3,6); V01 |
| fluorophore.tracrRNA | 5′-dye-tracrRNAs | Fluorescence | 5′ end of separate tracrRNA | Terminal conjugate, not nucleotide | Not applicable | Cy3 example; no nucleotide identity | sample_specific | conditional | Commercially synthesized dye-tracrRNA | reported | S01 RNA preparation and Results; V01 |
| chemical_modification.terminal | 3′-azido-modified RNA | Dye attachment | RNA 3′ extension | 3′ chemistry | Extension length not specified | Azido adenosine additions and DBCO-dye reaction | sample_specific | conditional | PAP / click labeling alternative | reported | S01 Dye labeling of RNAs; V01 |
| adapter.forward | Adapter sequence | Library-selective PCR | Upstream of T7 promoter, DNA only | Forward/sense inferred | 17 nt each, counted | Exact adapters 1–4 above | sample_specific | conditional | Twist templates; assignment not resolved for every library | reported | S02 Table S1(1,3,5); V01 |
| promoter.T7 | T7 promoter | IVT initiation | Upstream of spacer DNA | Sense arrangement inferred | 17 nt, counted | Exact DNA string above | fixed | yes | Manufacturing template only | reported | S02 Table S1(1–4); V01 |
| primer.reverse.crRNA | Reverse primer | Template PCR | Separate DNA oligo; anneals to common template region | Primer direction inferred from function | 22 nt, counted | Exact DNA string above | fixed | conditional | crRNA template manufacture | reported | S02 Table S1(1–2); V01 |
| primer.reverse.sgRNA | Reverse primer | Template PCR | Separate DNA oligo | Primer direction inferred from function | 25 nt, counted | Exact DNA string above | fixed | conditional | sgRNA template manufacture | reported | S02 Table S1(3–4); V01 |

Evidence column describes component presence/identity; counted lengths and inferred orientations are explicitly qualified in their cells. dCas9 is a required protein partner, not a nucleotide component. S01, Results design subsection; V01.

## Target-Binding Design

All reported rules below refer to S01, STAR Methods, “Design tiling crRNAs/sgRNAs targeting non-repetitive regions”, V01, unless a different locator is given. They describe that design pipeline, not every published repetitive guide.

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | Genomic DNA | V01 | reported | S01 design methods |
| Target region | hg38 regions on Chr3, Chr8, Chr13, Chr14, ChrX; MYC/FOS enhancer selection uses ABC predictions | V01 | reported | S01 design methods; S02 Table S1(7) |
| Probe length | Spacer 13–20 bp in method design; compared 13–17 and 18–20 nt groups; website advertises 18–20 bp | V01 / online interface scoped separately | reported | S01 Figure 2 text; S04 homepage |
| GC constraints | 35%–80% | V01 non-repetitive guides | reported | S01 design methods |
| Tm constraints | Not reported in examined design section | V01 | unavailable | S01 design methods |
| Specificity filtering | Cas-OFFinder against hg38: ≤1 mismatch for <16 bp, ≤2 for ≥16 bp; remove guides with >100 off-target loci; remove guides involved in regions with ≥2 off-target loci separated by <2000 bp | V01 paper pipeline | reported | S01 design methods |
| Repeat masking | No explicit repeat-mask program or threshold specified | V01 | unavailable | S01 design methods |
| Exon selection | Not applicable as a transcript-exon selection requirement; genomic loci selected for imaging | V01 | reported | S01 design methods |
| Isoform handling | No transcript-isoform selection rule specified for genomic guide design | V01 | unavailable | S01 design methods |
| Sequence complexity | Reject TTTT, 10×G, CTCGAG, GCTNAGC, CCANNNNNNTGG; initial G for IVT efficiency | V01 | reported | S01 design methods |
| Spacing | Remove guides too close to previous guide (<10 bp); exact operational meaning should be checked in code | V01 paper pipeline | reported | S01 design methods; S05 remove_overlap_3rd.py |
| Overlap | See spacing; no additional universal overlap threshold stated | V01 | unavailable | S01 design methods |
| Number of probes per target | ACAP2 experiment compares 12, 24, 48, 96, 192; these are experimental sets, not a universal requirement | V01 | reported | S01 Figure 2 and Results |
| Chromatin context | Guide abundance and accessibility dominate measured performance; targeting near ATAC-seq/H3K27ac peaks is used | V01 | reported | S01 Figure 2 Results |

Public-code interpretation (**inferred**, S05 pinned commit): `generate_possible_guides.py`, function `generate_sgRNA_seq`, scans candidate lengths from longest to shortest to find a G-start spacer; it reverse-complements minus-strand candidate DNA before output. `01_run_generate_possible_guides.sh` supplies hg38.fa and GC bounds 0.35/0.8. `remove_overlap_3rd.py` compares interval gaps against 10. These observations do not prove Table S2 was produced by this exact commit.

Preserve paper/code differences: `04_filter_off_targets_ifnoHF.sh` applies `$6<2` to filenames containing `18` and excludes `chr3_KI.*alt`; the paper’s blanket ≥16-bp mismatch statement does not describe those special cases. Do not silently replace one specification with the other. S01 design methods (reported); S05 named script (inferred); V01/current public code.

## Barcode System

No barcode/codeword layer is specified in the examined V01 workflow. Identity follows targeting-pool/channel assignment. Barcode length, barcode-library size, codebook, Hamming distance, error correction and sequencing-round assignment are not applicable to this direct imaging architecture; no claim of their universal absence from future implementations is intended. The PCR adapter is a manufacturing selection handle, not an imaging barcode. S01, Results/assembly; S02, Table S1; V01 (reported architecture, inferred functional distinction).

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | No separate sequential-readout handle specified | V01 | unavailable | S01 design/assembly; S02 Table S1 |
| Readout probe | Labeled guide RNA complex supplies signal directly | V01 | reported | S01 Figure 1 and Results |
| Fluorophore | Universal 5′-Cy3-tracrRNA example; Cy3/Cy5 terminal guide-labeling variants | V01 examples | reported | S01 Results and Dye labeling; S02 Table S1 |
| Imaging round | Continuous/time-lapse imaging; no combinatorial round codebook | V01 | reported | S01 Figure 5 and microscopy methods |
| Sequencing primer | Not applicable; imaging does not sequence guides | V01 | reported | S01 workflow |
| Decoding role | Fluorescence localization and tracking | V01 | reported | S01 image analysis |

Fixed-cell validation smiFISH uses FLAP-Cy5 and target-specific primary probes (S01 smiFISH section; S03 RNA FISH sheets). Those handles must not be inserted into the live-cell guide architecture.

## Primer / Linker / Adapter Components

DNA manufacture uses the published forward adapters, T7 promoter, and common reverse primers listed above. Table S1(2) prints the IDT forward oligo with a shorter common-region overlap than the complete RNA tail in Table S1(6); they are different molecular stages. Do not equate that forward oligo with mature crRNA or append primers to it as an asserted synthesized full probe. No separate universal linker is established. Reported / unavailable: S02, Table S1(1–6); V01.

## Amplification-Related Components

PCR primers and T7 promoter are manufacturing components. The imaging workflow does not use padlock ligation, rolling-circle amplification, or a sequencing amplification handle. S01, DNA preparation, RNA synthesis and assembly; S02, Table S1; V01 (reported workflow).

## Chemical Modifications

Universal dye-tracrRNAs are synthesized by IDT. Alternative IVT RNA labeling includes 3′ azido-nucleotide addition with yeast PAP followed by DBCO-dye click chemistry, and 5′ EDC-mediated labeling. Dye-UTP incorporation was tested but gave inefficient labeling in the reported comparison. Exact 3′ extension length/dye stoichiometry per RNA is not established. Reported / unavailable: S01, Results and “Dye labeling of RNAs”; V01.

## Probe-Design Workflow

Paper sequence: select hg38 regions → enumerate NGG-associated, G-start guides in specified length range → apply GC/motif filters → Cas-OFFinder screen → remove frequent/clustered off-target candidates → enforce spacing → choose pools and manufacture RNA. Preserve chromatin-accessibility context when selecting targets. S01, design methods and Figure 2 Results; V01 (reported).

Public software accepts a BED region input and genome FASTA, produces candidate guides, requires a Cas-OFFinder stage, then filters to a final spacer BED. It does not by itself document complete dye-labeled oligo manufacture. S05, README and code (reported documentation / inferred code behavior).

## Official Probe-Design Software

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Oligo-LiveFISH online design tool | https://oligo-livefish.org/ | Endorsed by S01 and linked to S05 on official homepage | No numbered release inspected | Backend language unresolved | Cell line, gene, selected non-repetitive region | Spacer pool according to instructions | Advertises 18–20-bp tiling spacers | Online service; backend dependency versions unresolved | Homepage/instructions reachable 2026-10-15; successful design job not tested | Not established | reported | S01 availability; S04 homepage/instructions; V01 |
| Oligo-LiveFISH crRNA/sgRNA design | https://github.com/LinBiologyLab/Oligo-LiveFISH_crRNA_sgRNA_design | Official website links directly to repository | Commit 210b34c2c54fd0f7e282a67e8eca144d96f44ad7; numbered release not established | Python and shell | BED regions, hg38.fa, guide length bounds, Cas-OFFinder results | Filtered spacer coordinates/sequences; not full synthesized labeled probes | See design workflow and paper/code differences | Biopython; BEDTools; Cas-OFFinder; shell utilities; scripts reference module environment/CUDA | Public source fetched 2026-10-15; not executed end to end | No maintenance commitment established | reported | S04 homepage GitHub link; S05 README, code/*.py, code/*.sh; current public code |

Reproduction cautions from code inspection (**inferred**, S05): `run_all.sh` comments out stages 1–4; paths/environment modules require configuration; `04_filter_off_targets_ifnoHF.sh` references `$PIPELINE_DIR`. The README presents a simplified final-output example, whereas the writer embeds chromosome/start/end/strand/spacer into BED column 4. This is not evidence of a turnkey reproducible release. Cas-OFFinder is a third-party dependency, not an Oligo-LiveFISH-authored program.

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | Published Table S2, publisher mmc2.xlsx | V01 | Eleven genomic spacer sheets; FOS/MYC RNA FISH and qRT-PCR sheets separately | Actual XLSX downloaded and all sheets inspected, 2026-10-15 | S03, sheet inventory below |
| D02 | Published Table S1, publisher mmc1.pdf | V01 | RNA scaffold sequences, degenerate templates, four adapters, primers, repetitive spacers and hg38 regions | Actual PDF downloaded; pp.2–4 inspected, 2026-10-15 | S02 Table S1(1–7) |
| D03 | S05 README/example_data at pinned commit | Current public design-code snapshot associated with V01 | Example input BED and README output rows; no complete experimental probe library established | Source inspected, 2026-10-15 | S05 README; example_data/chr3_repeatregion_lf20K.bed |

Table S2 sequence-row counts (**derived** by counting nonempty A/C/G/T spacer cells in column D, excluding header/formatting-only rows):

| Sheet | Spacer records | Observed string-length range (nt; derived) |
|---|---:|---|
| Chr3_ACAP2-nonrep_P1-P4 18-20bp | 96 | 18–20 |
| Chr3_ACAP2-nonrep_P5-P8 18-20bp | 96 | 18–20 |
| Chr3_ACAP2-nonrep_P5-P8 13-17bp | 96 | 13–17 |
| Chr3_ACAP2-nonrep_TwistBio | 136 | 13–20 |
| Chr13_TUBGCP3 | 96 | 13–20 |
| ChrX_PLCXD1 | 96 | 14–20 |
| MUC4 gene | 96 | 18–20 |
| FOS_promoter | 237 | 13–20 |
| FOS_enhancer | 343 | 13–20 |
| MYC_promoter | 340 | 13–20 |
| MYC_enhancer | 375 | 13–20 |

S03 columns A–E contain chromosome, start, end, spacer, length formula. No explicit strand field is present in these sheets. Two RNA FISH sheets have 48 named primary sequences each; qRT-PCR has four named primers. They are assay-validation resources, not live-cell guide libraries. S03, named sheets (reported cells; counts derived).

## Version Differences

Within V01, the usual universal 5′-dye-tracrRNA implementation differs from dye-labeled crRNA and terminally labeled sgRNA alternatives. SL-Oligo-LiveFISH changes localization/imaging, not necessarily the guide sequence architecture. LiveFISH is predecessor context; Oligopaint is not adopted as a version. S01, Introduction/Results; S02 Table S1. No later developer method version was established in the scoped search.

## Unresolved Information

- Table S2 has no explicit genomic-strand column or coordinate-origin declaration. The code emits 0-based BED and 1-based TXT alternatives; do not transfer a convention automatically to the spreadsheet.
- Eight rows in `Chr3_ACAP2-nonrep_P5-P8 13-17bp` have `End − Start != len(Spacer Sequence)`: Excel rows 6, 14, 21, 27, 45, 52, 57, 65. This is a derived internal-consistency observation, not a corrected coordinate assignment. S03 columns B–D.
- The examined workbook has no separately named P1–P4 13–17-bp sheet; do not invent it to complete an experimental matrix.
- Per-library adapter choice, full synthesized DNA assemblies, strand-to-guide mapping for every row, and exact chemical-label stoichiometry are not established.
- The public code commit is not demonstrated to reproduce every published table exactly; paper/code filter differences remain explicit.

These gaps are Research-stage findings, not independent validation. V01; S02, S03, S05.
