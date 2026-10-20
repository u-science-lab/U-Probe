# Probe Architecture

## Overview

Two molecular roles must remain distinct: encoding probes bind target RNA; readout probes bind complementary handles. Binary codewords describe signal patterns, not a contiguous nucleotide segment. Ex situ manufacturing primers are not in situ sequencing primers (S01, Fig. 1 and Methods; S03 §3–4; V01; reported).

## Canonical Architecture

All diagrams below are **schematic sequences**, not nucleotide strings. `RC(x)` denotes reverse complement, never the unchanged target strand. DNA hybridizes antiparallel to RNA. The left/right ordering is based on S03 §3.4–3.5 and Table 2 (reported); extension to the V01 extra priming region uses S01 Probe design and Probe synthesis (reported placement, synthesis interpretation noted below).

### V01 — encoding core, template, and finished probe

```text
Encoding core, ssDNA:
5′ — [RC(readout probe i), 30 nt] — [RC(target RNA region), 30 nt] — [RC(readout probe j), 30 nt] — 3′

Array template, 1001-gene design:
5′ — [Index Primer 1 sequence] — [encoding core] — [RC(Index Primer 2)] — 3′

Array template, 140-gene design:
5′ — [Index Primer 1 sequence] — [common RT priming region, 20 nt] — [encoding core] — [RC(Index Primer 2)] — 3′

Readout oligo:
5′ — [readout probe, 30 nt] — [3′ Cy5] — 3′
```

The T7 promoter is introduced by the reverse PCR primer; it is not a universal in situ encoding-probe segment (S01, Probe synthesis; reported). RT initiates with the common primer in the 140-gene libraries and Index Primer 1 in the 1001-gene library. Thus a whole array template should not automatically be called the final probe. Exact finished-probe termini for every member remain unresolved without Table S5 and Fig. S3 sequence-level inspection. The 90-nt core length is **derived** as 30+30+30, not the full probe length. Only two of each gene’s four on-bit readouts occur in an individual encoding core (S01, Probe design; reported).

### V02 — three-readout design

```text
5′ — [priming region] — [one or two 20-nt readout sequences] — [30-nt target-binding region]
   — [remaining two or one 20-nt readout sequences] — [priming region] — 3′
Readout oligo: 5′ — [dye–disulfide linkage] — [20-nt readout probe] — 3′
```

Adenosine spacers occur at junctions; the exact placement and finished-probe termini are library-dependent. Diagram omits these junction details rather than asserting zero-length spacers (S04, SI Encoding Probe Design; reported).

### V03/V04 — anchoring molecule, distinct from encoding probes

V03 anchor: 5′ — [acrydite] — [15-nt alternating dT / thymidine-LNA region] — 3′. It pairs with poly(A) RNA and anchors it to the gel (S05, SI Staining with Encoding Probes; reported). V04 instead reports a 20-nt alternating dT/LNA region plus a 20-nt reverse complement of a readout sequence and a 5′ acrydite. Relative order of the two nucleotide regions is unresolved in the inspected prose; do not reuse the V03 length (S10, Encoding probe staining; S11 corrected concentration; reported components).

### V05 — bDNA molecular layers

```text
Primary amplifier:   5′ — [(20-nt secondary-binding site + 2-nt spacer) × 5] — [20-nt encoding-handle-binding region] — 3′
Secondary amplifier: 5′ — [(20-nt final-readout-binding site + 2-nt spacer) × 5] — [20-nt primary-repeat-binding region] — 3′
Final readout DNA:   5′ — [20-nt final readout probe] — 3′
+ fluorescent dye/disulfide linkage; exact attachment end not independently resolved here
```

Amplifier 3′ binding-region placement is explicitly reported (S08, Amplifier staining). The five repeated blocks and 130-nt total length are independently **derived** from S12 Table S1 by removing formatting whitespace: 5×(20+2)+20 = 130 nt. Primary and secondary binding is antiparallel. S12 headers “Original readout sequences” are retained verbatim: all 16 entries match the V02 readout **probe** sequences and the primary 3′ tail, not the encoding-probe handle (derived equality; inferred strand interpretation supported by S04 Table S1 and S08 Fig. 1). Do not place that column unchanged onto encoding probes. The 4×4 and 9×9 tests in S08 are separate configurations; S12 is the 16-bit 5×5 library.

### V06

Encoding probes contain a 30-nt targeting sequence and three readout sequences selected from 69 identities. The inspected main text does not resolve all finished-oligo terminal/junction details; do not transfer V02 primers or spacer sequences (S09, Results first subsection; reported components).

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region | targeting sequence / targeting region | RNA recognition | central core | antisense; 5′→3′ DNA | 30 nt | Variable region; VCAN example below | variable | yes | V01 OligoArray rules below | reported | S01 Probe design; S03 §3.4; V01 |
| readout_sequence_left | readout sequence | presents one on-bit | 5′ flank of target | RC of readout probe | 30 nt V01; 20 nt V02 | Library sequences; not abstract codeword | barcode_library | yes | Use assigned on-bit; avoid same sequence twice in V01 core | reported | S01 Probe design; S03 §3.4; S04 SI; V01/V02 |
| readout_sequence_right | readout sequence | presents another on-bit | 3′ flank of target | RC of readout probe | 30 nt V01; 20 nt V02 | Library sequences | barcode_library | yes | V02 has third readout on one side | reported | S01 Probe design; S04 SI Encoding Probe Design; V01/V02 |
| primer_index_1 | Index Primer 1 | selective template amplification / V01 MHD2 RT | 5′ end of array template | sequence as primer | 20 nt | Library-specific; example below | sample_specific | conditional | Manufacturing; not an RNA-targeting arm | reported | S01 Probe design and synthesis; V01 |
| primer_index_2_site | Index Primer 2 | reverse-primer recognition | 3′ template end | RC of primer’s index segment | 20 nt | Library-specific | sample_specific | conditional | T7 promoter is added by PCR reverse primer | reported | S01 Probe design; S03 §3.5; V01 |
| primer_common_rt | common priming region | reverse transcription during probe production | between index 1 and first readout | primer-oriented on design strand | 20 nt | CGGGTTTAGCGCCGGAAATG, RT primer exact sequence | sample_specific | conditional | 140-gene libraries; absent from 1001-gene design | reported | S01 Probe design and Probe synthesis Step 3; V01 |
| fluorophore_readout | Cy5 | fluorescent detection | readout oligo 3′ | not a nucleotide | Not applicable | /3Cy5Sp/ source notation | version_specific | yes | V01 photobleaching workflow | reported | S01 readout probe table; V01 |
| cleavage_site_readout | disulfide bond | reductive removal of dye | readout 5′ dye linkage | not a nucleotide | Not applicable | Exact linker chemical structure unresolved | version_specific | yes | TCEP cleavage | reported | S04 Table S1 note; V02 |
| spacer_junction | adenosine nucleotide spacers | avoid G-rich junctions | readout/target junctions | 5′→3′ | Library-dependent | Not transcribed from unavailable Dataset S1 | sample_specific | conditional | Avoid terminal G triplets / G quadruplets | reported | S04 SI Encoding Probe Design; V02 |
| anchor_sequence_polya | anchor probe | retain poly(A) RNA in gel | separate oligo | antiparallel to poly(A) | 15 nt V03 | Alternating dT and dT+; exact start phase not asserted | version_specific | conditional | 5′ acrydite, gel-clearing workflow | reported | S05 SI Staining with Encoding Probes; V03 |
| other_functional_sequence_primary | primary amplifier | bind encoding handle and recruit secondaries | binding segment at 3′; repeats 5′ | 5′→3′ | 130 nt in 5×5 library | Exact strings retained in S12 workbook | barcode_library | conditional | Three-letter primary sequences; five sites | reported | S08 Fig. 1, Amplifier staining; S12 Sheet1; V05 |
| other_functional_sequence_secondary | secondary amplifier | recruit final readout | binding segment at 3′; repeats 5′ | 5′→3′ | 130 nt in 5×5 library | Exact strings retained in S12 workbook | barcode_library | conditional | Three-letter secondary sequences; five sites | reported | S08 Fig. 1, Amplifier staining; S12 Sheet1; V05 |

### Exact published example, not a universal sequence

The VCAN template in S03 Table 2 is a DNA **sequence example** printed as the following five blocks in 5′→3′ design order (orientation explicitly explained in §3.4–3.5). These are the source strings, with no complement transformation:

| Role | Published DNA block |
|---|---|
| 5′ priming region | `CGCGGGCTATATGCGAACCG` |
| First readout handle | `TTAGTCGTCGCAACGCCAGGCAACTCATGC` |
| RNA-targeting region | `TAAAGAAATTAGATAGGCTGGAAATGCTTA` |
| Second readout handle | `AAATTGCGTGACGGACCTGGGCCATTGGCC` |
| 3′ priming region | `GCGTTGTATGCCCTCCACGC` |

These lengths sum to 130 nt (**derived**, 20+30+30+30+20). This is a **template**, not proof of all final probe termini. It has no common internal RT primer site and must not be generalized to the 140-gene design (S03 Table 2; S01 Probe design).

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | Cellular RNA; human IMR90 library | V01 | reported | S01 Probe design |
| Target region | Most abundant isoform; transcript models divided into 1-kb regions for computation | V01 | reported | S01 Probe design |
| Probe length | 30-nt targeting sequence; not full encoding-probe length | V01/V02 | reported | S01 Probe design; S04 SI |
| GC constraints | No numeric target-GC cutoff established in inspected V01 text; V02 43–63% | V01/V02 | reported | S01 Probe design; S04 SI Encoding Probe Design |
| Tm constraints | V01 target duplex >70 °C; exclude cross-hybridization targets >72 °C and predicted internal structure >76 °C. V02 66–76 °C | V01/V02 | reported | S01 Probe design; S04 SI |
| Specificity filtering | V01 reject mapping to >1 cellular RNA species; keep multiple sites on same RNA. Assembled-probe homology filters: >14 nt rRNA/tRNA, >17 nt genes with FPKM >10,000 | V01 | reported | S01 Probe design |
| Repeat masking | No named repeat-mask program established in inspected design sections; this does not mean masking is absent in every implementation | V01/V02 | unavailable | S01 Probe design; S04 SI Encoding Probe Design |
| Exon selection | Explicit exon-selection rule not established in inspected sections | V01/V02 | unavailable | S01 Probe design; S04 SI |
| Isoform handling | V01 most abundant isoform. V02 abundance-weighted 17-nt penalty tables; isoform and gene specificity indices each 75–100% | V01/V02 | reported | S01 Probe design; S04 SI |
| Sequence complexity | V01 exclude contiguous repeats of ≥6 identical bases. V02 readout constraints differ from target constraints | V01 | reported | S01 Probe design |
| Spacing | V02 select nonoverlapping valid regions starting at 5′ end of each isoform; additional gap not specified | V02 | reported | S04 SI Encoding Probe Design |
| Overlap | V06 nonoverlap for 9,050 genes; up to 20-nt overlap for 1,000 shorter/less targetable genes | V06 | reported | S09 Results first subsection |
| Number of probes per target | V01 ~192 (SD 2) in 140-gene experiments; ~94 (SD 6) in 1001-gene. V02 92 putative probes/gene. V06 up to 48 for nonoverlap design | Library-specific | reported | S01 Probe design; S04 SI; S09 Results |

**Source conflict retained:** S03 §3.1 Step 4 describes a “lower bound” of 72 °C for off-target hybridization and says no internal secondary structure with Tm “lower than” 76 °C; S01 instead excludes off-targets above 72 °C and predicted internal structures above 76 °C. S03 also says no run “longer than six”, whereas S01 excludes six or more. These are not interchangeable rules. V01 extraction above follows S01, with S03 wording retained here; no silent correction or runnable parameter equivalence is asserted.

V02 specificity calculations use 17-nt isoform/transcriptome tables weighted by measured abundance; the separate ncRNA table uses 15-nt sequences. Nearest-neighbor Tm calculations assume 300 mM monovalent salt and 5 nM probe concentration. A release-84 GRCh38 ncRNA URL is given explicitly; the cDNA annotation release is not thereby established (S04 SI Encoding Probe Design; reported).

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | 16 bits MHD4 or 14 bits MHD2; V06 69 bits | V01/V06 | reported | S01 Results; S09 Results |
| Number of barcodes | 140 MHD4; 1,001 MHD2; 12,903 in V06 selected code | V01/V06 | reported | S01 Results; S09 Results |
| Codebook | S01 Tables S1/S3 identified but not parsed; author-hosted examples linked by S07 | V01 | reported | S01 Probe design; S07 Additional Data and Resources |
| Hamming distance | 4 for MHD4; 2 for MHD2; weight 4 in both | V01 | reported | S01 Results |
| Error correction | MHD4 corrects one bit and detects any two-bit errors; MHD2 detects single errors but does not uniquely correct them | V01 | reported | S01 Results and Construction of measured words |
| Round assignment | One bit per round V01; two per round V02; three per round V06 | V01/V02/V06 | reported | S01 imaging Methods; S04 MERFISH Imaging; S09 Results |
| Channel assignment | V01 Cy5; V02 bit-specific Cy5/Alexa750 table below; V06 three-color assignment not fully transcribed | Version-specific | reported | S01 readout table; S04 Table S1; S09 Results |
| Fixed vs experiment-specific barcode sets | Codeword-to-gene assignment is experiment-specific; V01 includes two shuffled MHD4 codebooks. V01 controls: 5 blank + 5 no-target for MHD4; 11 blank + 5 no-target for MHD2 | V01 | reported | S01 Probe design |

## Readout System

The following are **exact published readout-probe strings**, not encoding handles. V01 source notation `/3Cy5Sp/` is preserved; it denotes 3′ Cy5. The source’s oligo/3′ modification context and S03 §3.2/3.4 establish 5′→3′ readout orientation. Bit order is the S01 table, not a reconstructed codebook. All rows: reported, V01, S01 Materials and Methods → Probe design → readout probe table.

| Bit | Readout probe (DNA, 5′→3′, source notation) |
|---|---|
| 1 | CGCAACGCTTGGGACGGTTCCAATCGGATC/3Cy5Sp/ |
| 2 | CGAATGCTCTGGCCTCGAACGAACGATAGC/3Cy5Sp/ |
| 3 | ACAAATCCGACCAGATCGGACGATCATGGG/3Cy5Sp/ |
| 4 | CAAGTATGCAGCGCGATTGACCGTCTCGTT/3Cy5Sp/ |
| 5 | GCGGGAAGCACGTGGATTAGGGCATCGACC/3Cy5Sp/ |
| 6 | AAGTCGTACGCCGATGCGCAGCAATTCACT/3Cy5Sp/ |
| 7 | CGAAACATCGGCCACGGTCCCGTTGAACTT/3Cy5Sp/ |
| 8 | ACGAATCCACCGTCCAGCGCGTCAAACAGA/3Cy5Sp/ |
| 9 | CGCGAAATCCCCGTAACGAGCGTCCCTTGC/3Cy5Sp/ |
| 10 | GCATGAGTTGCCTGGCGTTGCGACGACTAA/3Cy5Sp/ |
| 11 | CCGTCGTCTCCGGTCCACCGTTGCGCTTAC/3Cy5Sp/ |
| 12 | GGCCAATGGCCCAGGTCCGTCACGCAATTT/3Cy5Sp/ |
| 13 | TTGATCGAATCGGAGCGTAGCGGAATCTGC/3Cy5Sp/ |
| 14 | CGCGCGGATCCGCTTGTCGGGAACGGATAC/3Cy5Sp/ |
| 15 | GCCTCGATTACGACGGATGTAATTCGGCCG/3Cy5Sp/ |
| 16 | GCCCGTATTCCCGCTTGCGAGTAGGGCAAT/3Cy5Sp/ |

V01 MHD4 uses bits 1–16; MHD2 uses bits 1–14 (S01 Probe design; reported). The handles on encoding probes are reverse complements; none are silently substituted above. V01 readout is fluorescence detection, not primer extension or sequencing.

V02 uses different 20-nt DNA probes. The dye is attached at the 5′ end through a disulfide bond. All rows: reported, V02, S04 SI Table S1. These are exact source strings in probe orientation, not three-letter encoding-handle sequences.

| Bit | Readout probe name | Sequence (5′→3′) | Dye |
|---|---|---|---|
| 1 | RS0015 | ATCCTCCTTCAATACATCCC | Cy5 |
| 2 | RS0083 | ACACTACCACCATTTCCTAT | Alexa750 |
| 3 | RS0095 | ACTCCACTACTACTCACTCT | Alexa750 |
| 4 | RS0109 | ACCCTCTAACTTCCATCACA | Cy5 |
| 5 | RS0175 | ACCACAACCCATTCCTTTCA | Cy5 |
| 6 | RS0237 | TTTCTACCACTAATCAACCC | Alexa750 |
| 7 | RS0247 | ACCCTTTACAAACACACCCT | Cy5 |
| 8 | RS0255 | TCCTATTCTCAACCTAACCT | Alexa750 |
| 9 | RS0307 | TATCCTTCAATCCCTCCACA | Alexa750 |
| 10 | RS0332 | ACATTACACCTCATTCTCCC | Cy5 |
| 11 | RS0343 | TTTACTCCCTACACCTCCAA | Cy5 |
| 12 | RS0384 | TTCTCCCTCTATCAACTCTA | Alexa750 |
| 13 | RS0406 | ACCCTTACTACTACATCATC | Cy5 |
| 14 | RS0451 | TCCTAACAACCAACTACTCC | Alexa750 |
| 15 | RS0468 | TCTATCATTACCCTCCTCCT | Alexa750 |
| 16 | RS0548 | TATTCACCTTACAAACCCTC | Cy5 |

V02 design uses A/T/G encoding handles (per-base generation probabilities 25%/25%/50%), retained GC 40–50%, no internal G stretch >3 nt, and no cross-homology or human-transcriptome homology >11 contiguous nt (S04 SI Encoding Probe Design; reported). Their complementary readout probes therefore use A/T/C; this alphabet relationship is **derived**, not an author-defined barcode nucleotide alphabet.

V05 final readouts bind secondary amplifiers. S12 provides a different final readout sequence for each bit; reusing the V02 sequence as the V05 final readout would be incorrect. Fluorescent dye/linkage metadata are not encoded in the S12 nucleotide strings (S08 Fig. 1; S12 Table S1).

## Primer / Linker / Adapter Components

See component records. S01 lists experiment-specific PCR primers and the common RT primer. The reverse PCR primer adds the T7 promoter; RT primers manufacture the probe, not read out RNA in situ. V02 adds adenosine spacers. No universal linker/adapter string across versions is established (S01 Probe design/synthesis; S04 SI Encoding Probe Design; reported).

## Amplification-Related Components

Ex situ V01 amplification: indexed PCR → T7 RNA transcription → DNA reverse transcription → alkaline RNA removal (S01 Probe synthesis; reported). V05 in situ signal amplification: encoding handle → primary amplifier → secondary amplifier → fluorescent readout (S08 Fig. 1; reported). Theoretical 5×5 gain is 25 binding sites; measured brightness gain must not be equated with this theoretical value. S08 reports 10.5-fold signal gain for the tested 5×5 configuration (Results, Properties of bDNA amplification; reported).

## Chemical Modifications

V01 readout 3′ Cy5; V02 readout 5′ dye–disulfide; V03 anchor 5′ acrydite with alternating DNA/LNA T residues; V04 anchor has distinct length/components. These are different molecules and versions, not a list of modifications present on every encoding probe (S01 Probe design; S04 Table S1; S05 SI; S10 Encoding probe staining; reported).

## Probe-Design Workflow

V01: quantify isoform abundance with Cufflinks v2.1 and ENCODE RNA data using GENCODE v18; choose most abundant isoform; design 30-nt targets using OligoArray2.1; screen specificity; assign codewords/readouts; concatenate templates and re-screen abundant RNA homology; synthesize indexed pools (S01 Probe design; reported). Exact genome FASTA checksum and coordinate conventions remain unresolved.

V02: load hg38 transcriptome and abundances; build 17-nt and 15-nt penalty tables and thermodynamic annotations; filter target windows; choose nonoverlapping sites; randomly assign three of four on-bit readouts and one of two arrangements; add spacers/primers; screen assembled oligos (S04 SI Encoding Probe Design; reported). S06 example code explicitly reverse-complements target and readout input sequences during assembly (lines 332–375; **inferred** algorithm behavior from inspected code). The example is not proof of identical historical outputs or random seeds.

## Official Probe-Design Software

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MERFISH_analysis | [ZhuangLab/MERFISH_analysis](https://github.com/ZhuangLab/MERFISH_analysis) | Author-endorsed by S07; readme links original and V02 papers | eba1735a4c350dbff31c06be1ab2df9d153e5c3d | MATLAB; legacy Python analysis dependencies | Transcript FASTA, isoform FPKM, ncRNA FASTA, readouts FASTA, codebook CSV | TargetRegions objects, probe/primer FASTA, codebooks, analysis outputs | Precompute off-target tables; filter target regions; assemble antisense/readout components; design primers | MATLAB; matlab-storm, storm-analysis, export-fig; README lists Python 2.7 for legacy storm-analysis. Example uses Bioinformatics functions and parallel pool | README and selected code fetched 2026-10-15 | README records updates through August 2021; no abandonment inference | reported | S06 README and example_scripts/library_design_example.m; later pipeline associated with V02 |
| OligoArray | Referenced external package | Third-party software used by original authors; not MERFISH-specific official authorship | 2.1 in S01; 2.0 terminology in S03 | Java | Transcript sequences / BLAST database | Candidate targeting oligos | Thermodynamic target, off-target and secondary-structure filtering | OligoArrayAux; legacy BLAST as described by S03 | Historical use reported; executable not installed or tested | Not assessed | reported | S01 Probe design; S03 §3.1; V01 |

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | NIHMS726733 supplementary Tables S1/S3/S5 | V01 | Code assignments and array template sequences identified from paper; full sequence files not parsed | Linked; retrieval failed 2026-10-15 | S01 Probe design; S02 supplementary inventory |
| D02 | pnas.1612826113.sd01.xlsx | V02 | Encoding-probe template library reported; readout sequences directly inspected in article Table S1 | Dataset link identified; binary retrieval failed 2026-10-15 | S04 SI Encoding Probe Design; S13 |
| D03 | 41598_2019_43943_MOESM2_ESM.xlsx | V05 | 16 bit-mapped primary/secondary amplifier pairs, original and final readout columns | Workbook parsed 2026-10-15 | S12 Sheet1 rows 2–18 |
| D04 | pnas.1912459116.sd01.csv / sd02.csv | V06 | Codebook / probe resources identified by main-text Dataset S1/S2 citations; files not parsed | Article inventory inspected 2026-10-15 | S09 Results first subsection, Supplementary Material |
| D05 | MERFISH_Examples.zip / MERFISH_Examples2.zip | V01/V02 software examples | Codebooks, raw/analyzed data and fast-design example inputs described on author page | Page and code inspected; archives not accessible 2026-10-15 | S07 Additional Data and Resources; S06 README/example script |

Sequence collections are linked with source-level provenance; they are not represented as fully reconstructed final synthesized probes. Grades and reference limitations are in benchmarks.md.

## Version Differences

The architecture, probe numbers, handle lengths, fluorescence attachment position, signal-removal chemistry and target rules differ by version. The original MHD4/MHD2 alternatives belong to V01; later clearing/expansion/bDNA features are not original features. See Canonical Architecture and description.md. V06 uses overlapping probes for a subset, not as a universal MERFISH rule (S01/S04/S05/S08/S09/S10; locators above).

## Unresolved Information

Original Table S5 and Fig. S3 binary contents; complete V02 Dataset S1; all exact final encoding-probe termini; full transcript/cDNA reference snapshots and target coordinates; V04 anchor subregion order; V06 complete sequence-level assembly; historical design random seeds and reproducible dependencies. S12’s “Original readout sequences” label is ambiguous about molecular strand: functional interpretation is explicitly inferred above. No unresolved sequence is invented.
