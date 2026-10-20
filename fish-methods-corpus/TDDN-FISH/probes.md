# Probe Architecture

## Overview

TDDN-FISH employs a multi-component probe system consisting of:
1. Bifunctional primary probes (target-binding + readout sequence)
2. Tetrahedral DNA Dendritic Nanostructure (TDDN) amplification system with three monomer types (T0, T1, T2)
3. Fluorophore-labeled detection oligos

The TDDN system provides hierarchical signal amplification through layer-by-layer self-assembly of tetrahedral DNA monomers. Each monomer consists of four single-stranded DNA oligonucleotides that anneal to form a tetrahedron, with sticky ends for programmable assembly (reported: S01, Figure 1a and Figure 2a).

## Canonical Architecture

### Primary Probe Architecture (5' → 3')

The primary probe is a linear, single-stranded DNA oligo containing two functional regions:

```
5' — [Target-specific sequence (35 nt)] — [Readout sequence (30 nt)] — 3'
```

**Note**: The exact 5'→3' orientation of the two regions is not explicitly stated in the examined sources. The reported orientation shows target-binding on the 5' side and readout sequence on the 3' side based on the schematic in Figure 2c, but this requires verification against the actual sequences in Supplementary Data 1.

### TDDN Component Architecture

**T0 Monomer (Core):**
- Assembled from 4 ssDNA strands: T0a, T0b, T0c, T0d
- 17 base pair side length
- 4 sticky ends:
  - T0a sticky end: binds to primary probe readout sequence
  - T0b, T0c, T0d sticky ends: identical; bind to T1 sticky ends

**T1 Monomer (First Dendritic Layer):**
- Assembled from 4 ssDNA strands: T1a, T1b, T1c, T1d
- 17 base pair side length
- 4 sticky ends:
  - T1a sticky end: binds to T0 sticky ends
  - T1b, T1c, T1d sticky ends: identical; bind to T2 sticky ends

**T2 Monomer (Second Dendritic Layer):**
- Assembled from 4 ssDNA strands: T2a, T2b, T2c, T2d
- 17 base pair side length
- 4 sticky ends:
  - T2a sticky end: binds to T1 sticky ends
  - T2b, T2c, T2d sticky ends: identical; bind to fluorophore-labeled detection oligos

### Complete Probe Assembly (Schematic)

```
RNA target
    ↓ (hybridization)
Primary probe — [Target-specific region] — [Readout sequence]
                                            ↓ (sticky end binding)
                                       T0 monomer
                                            ↓ (sticky end binding)
                                       T1 monomer
                                            ↓ (sticky end binding)
                                       T2 monomer
                                            ↓ (sticky end binding)
                                 Fluorophore-labeled detection oligos
```

### Disulfide Cleavage Site

The primary probe contains a disulfide bond (–S–S–) that is cleavable by TCEP, enabling probe removal between imaging rounds. The exact position of the disulfide bond within the primary probe is not explicitly reported in the examined sources (reported: S01, Figure 1c and STAR Methods).

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Primary probe | Primary probe; primer probe | Target detection | N/A | N/A | 65 nt total (35 nt target + 30 nt readout) | Exact sequences not provided | variable | yes | Target region: 35 nt; Readout: 30 nt (based on MERFISH) | reported | S01, STAR Methods |
| Target-specific region | Target sequence | RNA hybridization | 5' region of primary probe | 5'→3' | 35 nt | Exact sequences not examined | variable | yes | Designed with OligoArray 2.0 | reported | S01, STAR Methods |
| Readout sequence | Reading sequence; readout sequence | TDDN binding site | 3' region of primary probe | 5'→3' | 30 nt | Exact sequences not examined | fixed (across method) | yes | Designed based on MERFISH | reported | S01, STAR Methods |
| Disulfide bond | Disulfide bond; –S–S– | TCEP-cleavable cleavage site | Within primary probe backbone | N/A | N/A | N/A | fixed | yes | Enables probe removal | reported | S01, STAR Methods |
| T0 monomer | Tetrahedral DNA (T0); Tn monomer | Core amplification scaffold | N/A | N/A | 17 bp side length; hydrodynamic diameter ~5.8 nm | Exact sequences not examined | fixed | yes | Assembled from 4 complementary ssDNA strands | reported | S01, STAR Methods and Figure 2a |
| T1 monomer | Tetrahedral DNA (T1); Tn monomer | First dendritic layer | N/A | N/A | 17 bp side length; hydrodynamic diameter ~5.8 nm | Exact sequences not examined | fixed | yes | Assembled from 4 complementary ssDNA strands | reported | S01, STAR Methods and Figure 2a |
| T2 monomer | Tetrahedral DNA (T2); Tn monomer | Second dendritic layer | N/A | N/A | 17 bp side length; hydrodynamic diameter ~5.8 nm | Exact sequences not examined | fixed | yes | Assembled from 4 complementary ssDNA strands | reported | S01, STAR Methods and Figure 2a |
| T0 sticky ends | Sticky end | Programmable binding sites | T0 ends | N/A | Not reported | Exact sequences not examined | fixed | yes | 1 sticky end for primary probe; 3 for T1 binding | reported | S01, STAR Methods |
| T1 sticky ends | Sticky end | Programmable binding sites | T1 ends | N/A | Not reported | Exact sequences not examined | fixed | yes | 1 sticky end for T0; 3 for T2 binding | reported | S01, STAR Methods |
| T2 sticky ends | Sticky end | Fluorophore binding sites | T2 ends | N/A | Not reported | Exact sequences not examined | fixed | yes | 1 sticky end for T1; 3 for detection oligo binding | reported | S01, STAR Methods |
| Fluorophore-labeled detection oligo | Dye-modified single-stranded nucleic acids; fluorophore-labeled probe | Fluorescence readout | T2 sticky ends | N/A | Not reported | Exact sequences not examined | fixed | yes | Binds to T2 sticky ends | reported | S01, STAR Methods |

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | RNA (mRNA, lncRNA, miRNA, viral RNA) | V01 | reported | S01, Abstract |
| Target region | 35 nt target-specific sequence per probe | V01 | reported | S01, STAR Methods |
| Probe length | 65 nt total (35 nt target + 30 nt readout) | V01 | reported | S01, STAR Methods |
| GC constraints | Not explicitly reported | V01 | unavailable | Not examined |
| Tm constraints | Tm > 70°C for properly hybridized probe; minimum Tm 72°C for off-target hybridization | V01 | reported | S01, STAR Methods |
| Specificity filtering | Minimum Tm 72°C for hybridization to potential off-target sequences | V01 | reported | S01, STAR Methods |
| Repeat masking | Not explicitly reported | V01 | unavailable | Not examined |
| Exon selection | Not explicitly reported | V01 | unavailable | Not examined |
| Isoform handling | Not explicitly reported | V01 | unavailable | Not examined |
| Sequence complexity | No contiguous runs of the same nucleotide longer than 6 | V01 | reported | S01, STAR Methods |
| Spacing | Not explicitly reported | V01 | unavailable | Not examined |
| Overlap | Not explicitly reported | V01 | unavailable | Not examined |
| Number of probes per target | 3 primary probes for ACTB mRNA detection; 1 probe for miR-21 (72 nt) | V01 | reported | S01, STAR Methods and Figure 2 |
| Internal secondary structure | No internal secondary structures with Tm < 76°C | V01 | reported | S01, STAR Methods |
| Probe design software | OligoArray 2.0 for target sequences; MERFISH-derived readout sequences | V01 | reported | S01, STAR Methods |

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | 2 digits (one per imaging round in demonstrated experiments) | V01 | reported | S01, Methods |
| Number of barcodes | F^N combinations (F fluorescence channels × N imaging rounds) | V01 | reported | S01, Introduction |
| Codebook | Not explicitly reported; barcode determined by fluorescence pattern across rounds | V01 | unavailable | Not examined |
| Hamming distance | Not explicitly reported | V01 | unavailable | Not examined |
| Error correction | Not incorporated in current version | V01 | reported | S01, Discussion |
| Round assignment | Sequential; each round assigned unique fluorophore combinations | V01 | reported | S01, STAR Methods |
| Channel assignment | Cy5 = "1", Cy3 = "2", FAM = "3" | V01 | reported | S01, STAR Methods |
| Fixed vs experiment-specific barcode sets | Experiment-specific; different targets assigned different fluorophore combinations across rounds | V01 | reported | S01, STAR Methods |
| Barcode example | 11 (Cy5 round 1 + Cy5 round 2), 12 (Cy5 round 1 + Cy3 round 2), etc. | V01 | reported | S01, STAR Methods |

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | 30 nt sequence derived from MERFISH | V01 | reported | S01, STAR Methods |
| Readout probe | Fluorophore-labeled detection oligos | V01 | reported | S01, STAR Methods |
| Fluorophore | Cy5, Cy3, FAM (three-channel system demonstrated) | V01 | reported | S01, STAR Methods |
| Imaging round | Sequential rounds; 2 rounds demonstrated for 9-gene detection; up to 15 rounds tested | V01 | reported | S01, Figure 3e and STAR Methods |
| Sequencing primer | Not applicable | V01 | not_applicable | N/A |
| Decoding role | Color code assignment based on highest normalized fluorescence intensity per spot per round | V01 | reported | S01, STAR Methods |

## Primer / Linker / Adapter Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T0a | T0a | T0 strand with primary probe binding site | T0 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky end binds to readout sequence of primary probe | reported | S01, STAR Methods |
| T0b-T0d | T0b, T0c, T0d | T0 strands with T1 binding sites | T0 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky ends bind to T1a | reported | S01, STAR Methods |
| T1a | T1a | T1 strand with T0 binding site | T1 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky end binds to T0 sticky ends | reported | S01, STAR Methods |
| T1b-T1d | T1b, T1c, T1d | T1 strands with T2 binding sites | T1 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky ends bind to T2a | reported | S01, STAR Methods |
| T2a | T2a | T2 strand with T1 binding site | T2 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky end binds to T1 sticky ends | reported | S01, STAR Methods |
| T2b-T2d | T2b, T2c, T2d | T2 strands with fluorophore binding sites | T2 monomer | N/A | Not reported | Exact sequence not examined | fixed | yes | Sticky ends bind to fluorophore-labeled detection oligos | reported | S01, STAR Methods |

## Amplification-Related Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TDDN monomer | Tetrahedral DNA monomer; Tn | Self-assembling amplification unit | N/A | N/A | 17 bp side length; ~5.8 nm hydrodynamic diameter | Exact sequences not examined | fixed | yes | 4 ssDNA strands anneal to form tetrahedron | reported | S01, Figure 2a and STAR Methods |
| TDDN assembly | Layer-by-layer assembly | Hierarchical amplification | N/A | N/A | Not reported | N/A | N/A | yes | 3:1 molar ratio (Tn:Tn-1) for each layer | reported | S01, STAR Methods |

## Chemical Modifications

| Component | Original term | Function | Position | Modification type | Evidence status | Source / locator |
|---|---|---|---|---|---|---|
| Disulfide bond | Disulfide; –S–S– | TCEP-cleavable cleavage site | Within primary probe | Chemical crosslink | reported | S01, STAR Methods |
| Fluorophore | Cy5, Cy3, FAM | Fluorescent label | Detection oligos | Covalent dye conjugation | reported | S01, STAR Methods |
| DAPI | DAPI | Nuclear stain | N/A | Fluorescent dye | reported | S01, STAR Methods |

## Probe-Design Workflow

The probe-design workflow involves:

1. **Target selection**: Identify RNA targets of interest
2. **Reference preparation**: Obtain reference genome/transcriptome sequences
3. **Target region selection**: Use OligoArray 2.0 to select 35 nt target regions with:
   - Tm > 70°C for proper hybridization
   - Minimum Tm 72°C for off-target binding
   - No internal secondary structures with Tm < 76°C
   - No contiguous runs of same nucleotide > 6
4. **Readout sequence assignment**: Attach 30 nt readout sequence (derived from MERFISH)
5. **TDDN component synthesis**: Synthesize 4 ssDNA strands per TDDN monomer (T0a-d, T1a-d, T2a-d)
6. **Fluorophore-labeled detection oligo design**: Design oligos complementary to T2 sticky ends
7. **Barcode assignment**: Assign fluorophore combinations across imaging rounds

Software used: OligoArray 2.0 (target design), MERFISH readout sequences (readout sequences) (reported: S01, STAR Methods).

## Official Probe-Design Software

Not yet examined. No dedicated TDDN-FISH probe design software identified. The paper reports using OligoArray 2.0 for target sequence design and MERFISH-derived readout sequences.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | Supplementary Data 1 (Nature Communications article) | V01 | Primary probe sequences, TDDN component sequences, fluorophore-labeled detection oligo sequences | Available via article supplementary materials (2025-10-20) | S01, Supplementary Data 1 |

Note: The sequences in Supplementary Data 1 were not examined in this research stage.

## Version Differences

Not applicable. Only V01 has been identified.

## Unresolved Information

- **Exact primary probe sequences**: Available in Supplementary Data 1 but not examined
- **Exact TDDN component sequences**: Available in Supplementary Data 1 but not examined
- **Exact fluorophore-labeled detection oligo sequences**: Available in Supplementary Data 1 but not examined
- **Primary probe 5'→3' orientation**: Not explicitly verified against actual sequences
- **Disulfide bond position**: Reported but exact location within primary probe not specified
- **T0, T1, T2 sticky end sequences**: Not explicitly reported
- **Official probe-design software**: No dedicated software identified
- **TDDN-FISH codebook documentation**: Barcode assignments appear experiment-specific, no centralized codebook provided
