# Experimental Protocol

## Protocol Scope

Research-stage extraction, not an executable protocol. Unresolved: the publisher full methods and supplementary files could not be retrieved. The inspected author repositories do not supply a complete wet-laboratory protocol. This is an access/extraction gap, not evidence that the paper omits the information (S01; S04 README opening notice; V01). V01 denotes the original 2026 publication; code observations refer specifically to the inspected author implementation snapshot, not necessarily every experiment in that paper.

## Sample Preparation

Mouse brain is an evidenced application (S07 dataset description). Thickness, storage and preparation conditions remain unresolved here.

## Fixation

Unresolved; see Protocol Scope.

## Permeabilization

Unresolved; see Protocol Scope.

## Reverse Transcription

Applicability is unresolved. The software supports mRNA and cDNA design modes; this does not establish that reverse transcription is part of the published assay (S04 05_Probe_design cell 2 and cell 10).

## Probe Hybridization

Target and padlock design settings are documented in probes.md. In silico temperatures and salt/formamide inputs must not be substituted for verified wet-laboratory conditions (S04 05_Probe_design cell 10; 03_Padlock_join cell 3).

## Washing

Unresolved; see Protocol Scope.

## Ligation

Padlocks are explicitly designed by the author code. Ligase, reaction partners, actual ligation orientation, concentrations, temperature and duration remain unresolved (S04 03_Padlock_join; S01 full methods inaccessible).

## Amplification

RCA is reported by S04 README; enzyme, buffer, duration and temperature remain unresolved.

## Polymerization

The publisher Introduction describes post-MSI hydrogel embedding and clearing. Composition and polymerization conditions were not retrieved (S01 Introduction; V01 integrated workflow).

## Sequencing / Readout

Fluorescence readout is supported by RO design and channel-coded image decoding. Do not infer a sequencing chemistry from the presence of padlocks (S04 01_RO_design; 06_Decode). Readout exchange and stripping conditions remain unresolved.

## Imaging

The decoder notebook requests projected images, recommends avoiding prior stitching, and requires channel names to match the codebook. Its pixel scalers and channel offsets are calibration examples; they are not universal microscope parameters (reported; S04 06_Decode cells 2–4).

## Decoding

Follow author notebook 06 for preprocessing, registration, Spotiflow spot detection, Cellpose segmentation, decoding and matrix generation. Parameters are sample-dependent; the code was inspected but not executed (S04 06_Decode cells 4–10).

## Critical Parameters

No verified wet-laboratory parameter records were extracted. Blank table means unresolved values, not zero and not absence of these steps.

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Version-Specific Differences

Standalone OpenFISH and the integrated MALDI workflow are separate experimental scopes. The latter requires post-MSI processing; complete conditions and preprint differences remain unresolved (S01 Introduction; S02 inaccessible full text).

## Unreported Parameters

Unextracted because of source access: reagent concentrations, times, temperatures, buffers/pH, section thickness, ligase/polymerase, hydrogel chemistry, fluorophore conjugation and readout removal. These must not be described as absent from the publication. The author workflow README says a step-by-step wet protocol will be available; only computational instructions were accessible (S04 README opening notice).
