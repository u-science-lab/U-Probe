# Experimental Protocol

## Protocol Scope

This protocol covers TDDN-FISH for spatial RNA detection in cultured cells and mouse tissue sections. The method includes: sample fixation, primary probe hybridization, TDDN assembly and fluorophore binding, iterative multi-round imaging with TCEP probe removal, and image analysis. The protocol is based on S01 STAR Methods sections "TDDN-FISH procedure on cultured cells and tissue sample" and "Multi-round TDDN-FISH imaging".

## Sample Preparation

### Cell Culture

HeLa, A549, and MDCK cell lines were obtained from the Shanghai Cell Bank, Chinese Academy of Sciences. Cells were cultured in DMEM supplemented with 10% FBS and 1% penicillin-streptomycin at 37°C under 5% CO₂. For FISH assays, cells were plated onto glass-bottom dishes and grown to 70-80% confluence prior to fixation (reported: S01, STAR Methods "Cell culture").

### Tissue Sample Preparation

Six-week-old female Balb/c mice were used. Brain tissue was rapidly collected after euthanasia, frozen in optimal cutting temperature compound (Tissue-Tek O.C.T. Compound), and stored at -80°C. Fresh frozen brain samples were sectioned into 10 μm-thick slices using a cryostat and mounted onto functionalized coverslips. Sections were air-dried for 5 min at room temperature (reported: S01, STAR Methods "Mouse tissue sample preparation").

### Coverslip Functionalization

40 mm coverslips were cleaned by immersion in potassium dichromate and concentrated sulfuric acid solution for 2 h, followed by thorough rinsing with water. Coverslips were soaked in 95% ethanol for 12 h and air-dried. A working solution of 50× APES (Anti-Slice Escaping Agentia) was freshly prepared (1:50 dilution in acetone). Dried coverslips were immersed for 20-30 s, rinsed three times with distilled water, and stored dry (reported: S01, STAR Methods "Coverslip functionalization").

## Fixation

Cells and tissue samples were washed twice with PBS and fixed with 4% formaldehyde in PBS at room temperature for 10 min. Formaldehyde was aspirated and samples washed three times with PBS (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

## Permeabilization

Samples were incubated with 5 μg/mL Proteinase K in diethyl pyrocarbonate (DEPC)-treated PBS buffer for 20 min at room temperature. Samples were then transferred to permeabilization buffer (2× SSC, 0.5% Triton X-100) at room temperature for 15 min, followed by three washes with 2× SSC buffer. Samples were transferred to 70% ethanol at -20°C overnight (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

## Reverse Transcription

Not applicable. TDDN-FISH does not involve reverse transcription; it directly detects RNA targets through DNA probe hybridization.

## Probe Hybridization

### Primary Probe Hybridization

Samples were incubated in hybridization buffer A with the following composition:
- 10% dextran sulfate
- 2× SSC
- 30% formamide
- 20 mM RVC (RNA venier complex)
- 1 mg/mL yeast tRNA
- 0.1% RNase-free BSA
- Primary probes at 50 nM for each RNA target

Incubation conditions: 8-12 h at 37°C or 42°C (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

### TDDN Probe Hybridization

After primary probe hybridization and washing, samples were incubated in hybridization buffer B:
- 10% dextran sulfate
- 2× SSC
- 5% formamide
- 20 mM RVC
- 1 mg/mL yeast tRNA
- 0.1% RNase-free BSA
- TDDN-FISH probes

Incubation: 45 min to 1 h at room temperature (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

## Washing

### Post-primary probe washing

Hybridization wash buffer A (2× SSC, 30% formamide, 20 mM RVC) at 37°C for 10 min, three times, followed by three washes with 2× SSC for 5 min at room temperature. Samples stored in 2× SSC at 4°C if not proceeding immediately (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

### Post-TDDN probe washing

Hybridization wash buffer B (2× SSC, 5% formamide, 20 mM RVC) at room temperature, three times (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

## Ligation

Not applicable. TDDN-FISH does not involve ligation reactions.

## Amplification

Signal amplification is achieved through TDDN hierarchical assembly rather than enzymatic amplification:
- T0 monomers bind to primary probe readout sequences
- T1 monomers hybridize to T0 sticky ends
- T2 monomers hybridize to T1 sticky ends
- Fluorophore-labeled detection oligos bind to T2 sticky ends

TDDN monomers were mixed step by step at a 3:1 molar ratio (Tn:Tn-1) at room temperature for 1 h. Fluorophore-labeled T2 was pre-incubated with detection oligos at 1:3 molar ratio before adding to TDDN assembly (reported: S01, STAR Methods "Assembly of tetrahedral DNA dendritic nanostructures (TDDNs)").

## Polymerization

Not applicable. TDDN-FISH does not involve polymerization reactions.

## Sequencing / Readout

Not applicable. TDDN-FISH uses fluorescence imaging for readout, not sequencing.

## Imaging

### Equipment

Nikon Ti2 inverted microscope with:
- Confocal laser scanning
- Yokogawa CSU-X1 spinning disk unit
- Andor DU-897X EMCCD camera

Objectives:
- SR HP APO TIRF 100× oil immersion (NA 1.49) for high-resolution cell imaging
- Apo TIRF 60× oil immersion (NA 1.49) for single-particle analysis in tissues
- Plan Apo VC 20× air objective (NA 0.75) for large-area single-cell scans

Four-channel sequential excitation:
- 640 nm laser (Cy5)
- 561 nm laser (Cy3)
- 488 nm laser (FAM)
- 405 nm laser (DAPI)

Emission bandpass filters:
- 697/58 nm (Cy5)
- 615/40 nm (Cy3)
- 525/50 nm (FAM)
- 447/60 nm (DAPI)

Focal stability maintained via Nikon's Perfect Focus System (reported: S01, STAR Methods "Multi-round TDDN-FISH imaging").

### Nuclear Staining

DAPI (1 mg/mL) staining for 5 min (reported: S01, STAR Methods "TDDN-FISH procedure on cultured cells and tissue sample").

## Decoding

Image analysis pipeline:
1. Image registration using custom MATLAB scripts (correct positional shifts across rounds)
2. Spot localization using TrackMate plugin in ImageJ (for cell samples)
3. Cell segmentation using Cellpose algorithm (for tissue samples)
4. Fluorescence intensity normalization
5. Color code assignment based on highest normalized intensity channel
6. Two-digit barcode assignment per spot (e.g., 11 = Cy5 round 1 + Cy5 round 2)

Cy5 = "1", Cy3 = "2", FAM = "3" (reported: S01, STAR Methods "Imaging analysis").

## Critical Parameters

### Hybridization Temperature

Optimized between 37-42°C for primary probe hybridization (reported: S01, STAR Methods and Supplementary Figure 1d-g).

### Formamide Concentration

Optimized between 10-30% in hybridization buffers (reported: S01, STAR Methods).

### TCEP Cleavage

0.25 μM TCEP dissolved in 2× SSC, pH adjusted to 7.0 using NaOH, incubated for 15 min for probe removal between rounds (reported: S01, STAR Methods "Multi-round TDDN-FISH imaging").

### Buffer Composition

- TM buffer for TDDN assembly: 50 mM Tris-HCl, 8 mM MgSO4, pH 7.5
- TDDN annealing: 95°C denaturation for 5 min, then cooled to 4°C for 30 min
- Final TDDN concentration: 1 μM (reported: S01, STAR Methods).

### Flow Chamber Imaging

For multi-round experiments, samples were assembled in Bioptechs FCS2 flow chamber with temperature control. Fluidics controlled via peristaltic pump at 500 μL/min (10 rpm) (reported: S01, STAR Methods "Multi-round TDDN-FISH imaging").

## Version-Specific Differences

Not applicable. Only V01 has been identified.

## Unreported Parameters

The following parameters were not reported in the examined sources:
- Exact primary probe sequences (contained in Supplementary Data 1, not examined)
- Exact TDDN component sequences (contained in Supplementary Data 1, not examined)
- Exact fluorophore-labeled detection oligo sequences
- Detailed imaging acquisition parameters (exposure times, laser powers)
- Cell-type-specific marker gene selection criteria beyond scRNA-seq analysis
- Specific genome/transcriptome references used for probe design

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fixation | Formaldehyde | 4% | Room temperature | 10 min | PBS | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| Proteinase K treatment | Proteinase K | 5 μg/mL | Room temperature | 20 min | DEPC-treated PBS | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| Permeabilization | Triton X-100 | 0.5% | Room temperature | 15 min | 2× SSC | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| Primary probe hybridization | Primary probes | 50 nM per RNA | 37-42°C | 8-12 h | Hybridization buffer A (10% dextran sulfate, 2× SSC, 30% formamide, 20 mM RVC, 1 mg/mL yeast tRNA, 0.1% RNase-free BSA) | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| TDDN probe hybridization | TDDN probes | Not reported | Room temperature | 45 min - 1 h | Hybridization buffer B (10% dextran sulfate, 2× SSC, 5% formamide, 20 mM RVC, 1 mg/mL yeast tRNA, 0.1% RNase-free BSA) | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| TDDN assembly | T0, T1, T2 monomers | 1 μM each strand | Room temperature | 1 h per layer | TM buffer (50 mM Tris-HCl, 8 mM MgSO4) | 7.5 | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| Post-primary wash | - | - | 37°C | 10 min | Hybridization wash buffer A (2× SSC, 30% formamide, 20 mM RVC) | Not reported | 3 | - | - | Cells and tissue | reported | S01, STAR Methods |
| Post-TDDN wash | - | - | Room temperature | Not reported | Hybridization wash buffer B (2× SSC, 5% formamide, 20 mM RVC) | Not reported | 3 | - | - | Cells and tissue | reported | S01, STAR Methods |
| TCEP cleavage | TCEP | 0.25 μM | Not reported | 15 min | 2× SSC | 7.0 | - | - | - | Multi-round | reported | S01, STAR Methods |
| Nuclear staining | DAPI | 1 mg/mL | Not reported | 5 min | Not reported | Not reported | - | - | - | Cells and tissue | reported | S01, STAR Methods |
| TDDN annealing | - | 1 μM per strand | 95°C denaturation, 4°C cooling | 5 min denaturation, 30 min cooling | TM buffer | 7.5 | - | - | - | TDDN synthesis | reported | S01, STAR Methods |
