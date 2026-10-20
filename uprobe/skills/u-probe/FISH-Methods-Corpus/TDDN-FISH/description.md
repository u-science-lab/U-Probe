# Method Overview

## Purpose

TDDN-FISH is designed for rapid, enzyme-free, multiplexed spatial RNA detection and imaging. It enables single-cell and subcellular resolution mapping of RNA distributions in cultured cells and tissue sections, including short RNA targets such as miRNAs (reported: S01, Abstract).

## Core Principle

The method employs hierarchical self-assembling tetrahedral DNA dendritic nanostructures (TDDNs) as signal amplification units. TDDNs are assembled in a layer-by-layer fashion from three tetrahedral DNA monomers (T0, T1, T2), each containing four single-stranded overhangs (sticky ends) for programmable hierarchical conjugation. A bifunctional primary probe hybridizes to the target RNA and carries a readout sequence that serves as the attachment site for TDDN assembly. Fluorophore-labeled detection oligos bind to the outermost TDDN layer, producing amplified fluorescence signals (reported: S01, Figure 1a and Results).

## Overall Workflow

1. **Sample preparation**: Fixation with 4% formaldehyde, proteinase K treatment, permeabilization
2. **Primary probe hybridization**: Overnight hybridization of bifunctional primary probes (target-specific + readout sequence)
3. **TDDN assembly**: Sequential hybridization of T0, T1, and T2 tetrahedral monomers to the primary probe readout sequence
4. **Fluorophore binding**: Fluorophore-labeled detection oligos hybridize to T2 sticky ends
5. **Imaging**: Confocal fluorescence microscopy
6. **Probe removal**: TCEP (tris(2-carboxyethyl)phosphine) cleavage of disulfide bonds removes TDDN probes
7. **Iterative rounds**: Repeat steps 2-6 for additional targets using different fluorophore assignments

Each round of TDDN-FISH requires approximately 1 hour post-hybridization (reported: S01, Figure 1d and e).

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `nanostructure`.
- `readout_strategy`: `sequential_fish`.

Structural evidence: the primary recognition probe is linear (35-nt target + 30-nt readout); it is recognized by a hierarchical tetrahedral DNA dendritic nanostructure (T0-T1-T2) that acts as a signal amplifier, not as the recognition element (see `probes.md` Canonical Architecture and Probe Components). The architecture is therefore `nanostructure` (signal amplifier attached to a primary probe) rather than `linear` or `padlock`: there is no padlock circularization and the nanostructure does not provide target specificity. The F^N encoding is decoded by sequential imaging rounds with different fluorophore assignments per round, so readout is `sequential_fish` (not a single-round fluorescent barcode). Source: S01, Fig. 1a and Fig. 2 (TDDN assembly), Introduction, and "Multi-round TDDN-FISH imaging".

## Target Recognition Strategy

Target recognition is achieved through a bifunctional primary probe that contains:
- A 35-nucleotide target-specific sequence
- A 30-nucleotide readout sequence

The target-specific portion hybridizes directly to the RNA molecule of interest. The readout sequence provides a universal binding site for the TDDN amplification system (reported: S01, STAR Methods "Primer probe design and synthesis").

## Encoding Strategy

TDDN-FISH uses an F^N cyclic encoding-decoding framework:
- **F**: Number of fluorescence channels (up to 3 demonstrated: Cy5, Cy3, FAM)
- **N**: Number of imaging rounds

Each target RNA is assigned a unique barcode based on its fluorescence pattern across multiple imaging rounds. The barcode is read by observing which fluorophore channels are active in each round for each detected spot (reported: S01, Introduction and Methods).

Additionally, for tissue cell typing applications, the method uses a co-detection approach where multiple marker genes sharing the same barcode identify a cell type, inspired by FISH-nCHIPs (reported: S01, Discussion).

## Amplification Strategy

Signal amplification is achieved through hierarchical TDDN assembly:
- **T0 monomer**: Core structure with 4 sticky ends (1 for primary probe, 3 for T1 binding)
- **T1 monomer**: First dendritic layer with 4 sticky ends (1 for T0, 3 for T2 binding)
- **T2 monomer**: Second layer with 4 sticky ends (1 for T1, 3 for fluorophore-labeled detection oligos)

The modular assembly yields exponential signal amplification through the dendritic branching structure. Stepwise fluorescence enhancement is observed with each TDDN layer addition (reported: S01, Figure 2d).

TDDN monomers are synthesized by thermal annealing (95°C denaturation, 4°C cooling) of four complementary oligonucleotides in TM buffer. Assembly proceeds at room temperature with defined molar ratios (3:1 Tn:Tn-1 for each layer addition) (reported: S01, STAR Methods).

## Imaging / Readout Strategy

The readout system uses fluorophore-labeled single-stranded DNA detection oligos that bind to the sticky ends of T2 monomers. Three fluorophores are demonstrated:
- Cy5 (red)
- Cy3 (yellow/orange)
- FAM (green)

Detection uses confocal laser scanning microscopy with:
- 640 nm laser (Cy5 excitation)
- 561 nm laser (Cy3 excitation)
- 488 nm laser (FAM excitation)
- 405 nm laser (DAPI for nuclear staining)

Emission is collected through bandpass filters (697/58 nm, 615/40 nm, 525/50 nm, 447/60 nm) (reported: S01, STAR Methods "Multi-round TDDN-FISH imaging").

## Decoding Strategy

Decoding involves:
1. **Image registration**: Align images across rounds using MATLAB scripts to correct positional shifts
2. **Spot localization**: Identify individual RNA spots using TrackMate plugin (ImageJ)
3. **Color assignment**: Assign each spot a code based on which fluorophore channel has the highest normalized intensity in each round
4. **Barcode mapping**: Assign unique two-digit codes (e.g., 11, 12, 13, 21, 22, 23, etc.) based on color combinations across rounds

For tissue samples, Cellpose is used for cell segmentation, and fluorescence intensity normalization is applied within each cell mask region (reported: S01, STAR Methods "Imaging analysis").

## Major Versions or Variants

V01: TDDN-FISH as described in the 2025 Nature Communications publication. No variants or subsequent versions identified in the examined sources.

## Strengths

- **Rapid detection**: ~1 hour per imaging round, approximately 8-fold faster than HCR-FISH
- **High sensitivity**: Stronger fluorescence signals than smFISH, with only 3 primary probes for ACTB detection
- **Short RNA detection**: Successfully images miR-21 (72 nucleotides) with a single primary probe
- **Enzyme-free**: Eliminates batch variability from DNA ligase or polymerase activity fluctuations
- **Modular design**: TDDN components self-assemble through programmed sticky-end interactions
- **Scalable multiplexing**: F^N encoding enables detection of F^N RNA types with F channels and N rounds
- **Subcellular resolution**: Demonstrated in both cultured cells and tissue sections up to 40 μm thickness
- **Demonstrated target range**: mRNAs (e.g., ACTB, POLR2A, HSP70, mTOR, B2M), lncRNAs (MALAT1, NEAT1), miRNAs (miR-21), viral RNAs (influenza vRNAs)

Evidence: S01, Abstract, Results, and Discussion.

## Limitations

- **Signal overlap**: Ultra-high-density RNA populations (e.g., highly abundant mRNAs comparable to ACTB) may suffer from signal overlap
- **Tissue thickness**: Diffusion limitations affect probe penetration in specimens thicker than 100 μm
- **Resolution trade-off**: Lower magnification imaging enables faster whole-slide scanning but at reduced spatial precision
- **No error correction in encoding**: The method does not currently incorporate Hamming-distance-encoded probe sets with error correction
- **Enzyme-dependent reference**: The method optionally integrates with scRNA-seq data for cell typing, which requires enzyme-based library preparation for the reference data

Evidence: S01, Discussion.

## Unresolved or Version-Specific Points

- Exact primary probe sequences: Available in Supplementary Data 1 but exact sequences not examined in this research
- Full TDDN component sequences: Available in Supplementary Data 1 but exact sequences not examined
- Whether subsequent publications describe method improvements or new versions
- Whether the method has been applied by independent research groups
- Probe design software: No official TDDN-FISH-specific design tool identified; target sequences designed using OligoArray 2.0
