# Experimental Protocol

## Protocol Scope

All reported conditions below are from S01, STAR Methods, in the named sample/variant subsection. RT means room temperature as stated by the source; A549 fixation defines RT as 22°C. Missing values are not filled from common practice. V01 covers A549 endogenous RNA and mouse tissue; V02 covers Perturb-RAEFISH. This is an evidence summary, not a merged executable protocol.

## Sample Preparation

V01 A549-Cas9 cells are seeded on coverslips 3–4 days before fixation. Mouse C57BL/6 tissue is OCT embedded, rapidly frozen and stored at −80°C; liver and E13.5 placenta sections are 10 µm, lymph node sections 8 µm. The cryostat is −20°C. V02 uses doxycycline-inducible A549-Cas9i screen cells; Cas9 is induced with 1 µg/mL doxycycline for one week before seeding. Reported: S01, Cell lines, Mouse tissue section and CRISPR screen library construction.

## Fixation

V01 A549: 4 mL 4% PFA in DPBS, 22°C, 20 min. Mouse: 4 mL 4% PFA in PBS, RT, 30 min. Both are fixed after RCA with 3 mL 4% PFA in PBS for 30 min at RT. V02 initial fixation: 4% PFA in DPBS, RT, 20 min. Reported: S01, corresponding RAEFISH protocol subsections.

## Permeabilization

V01 A549: 3 mL 0.5% Triton X-100 in DPBS for 10 min; mouse: 3 mL 0.5% Triton X-100 in DPBS for 20 min at RT. V02: 0.5% Triton X-100 in DPBS for 10 min. Reported: S01, corresponding RAEFISH protocol subsections.

## Reverse Transcription

Specimen RNA reverse transcription is not part of the reported RAEFISH detection workflow. Probe manufacturing separately uses limited-cycle PCR, in vitro transcription, reverse transcription, alkaline hydrolysis and purification. Maxima H Minus RT: padlock/priming libraries 55°C for 15 min; encoding/MERFISH libraries 50°C for 60 min. A phosphate-labeled RT primer produces ligatable padlock probes. Reported: S01, Probe synthesis.

## Probe Hybridization

V01 A549: 50 µL, 2×SSC, 50% formamide, 20 U murine RNase inhibitor (MRI), each padlock/priming oligo 3 nM, 37°C for 16 h. Mouse: 25 µL, 2×SSC, 50% formamide, 10% dextran, 0.1 mg/mL yeast tRNA, 10 U MRI, each oligo 3 nM, 37°C for 16 h. A549 encoding: 50 µL, 2×SSC, 40% formamide, 10% dextran sulfate, 20 U MRI, 0.1 mg/mL yeast tRNA, 3 nM per encoding oligo, 37°C for 16 h. Mouse encoding: 25 µL, same stated constituents but 400 U MRI, 37°C for 16 h. Preserve this reported MRI difference. V02 primary: 50 µL, 2×SSC, 30% formamide, 200 U MRI, 5 nM per padlock and 100 nM shared priming oligo, 37°C for 16 h. V02 co-hybridization: 3 nM per RAEFISH encoding probe and 0.5 nM per MERFISH probe in 50 µL 2×SSC, 40% formamide, 10% dextran sulfate, 37°C for 16 h. Reported: S01, three RAEFISH protocol subsections.

## Washing

V01 A549 primary wash: 2 mL 50% formamide in 2×SSC, 0.1% Tween 20 and 120 U MRI, three times at 37°C for 15 min each; then DPBS twice, 2 min each. Mouse primary wash: three 2-mL washes with 50% formamide, 400 U MRI and 0.1% Tween 20 at 37°C, 15 min each; do not silently supply an unstated wash-buffer constituent. After ligation: 2 mL 20% formamide in 2×SSC at 37°C for 15 min, then two washes of 2 mL 20% formamide/2×SSC plus 800 nM toehold oligo at 37°C, 15 min each; DPBS three times, 5 min each. Encoding wash: twice in 2 mL 40% formamide/2×SSC/0.1% Tween at 37°C for 15 min, then 3 mL 40% formamide/2×SSC at RT for 2 min. V02 primary wash is 45% formamide in 2×SSC, 0.1% Tween 20, 120 U MRI, 3 mL three times for 15 min; its temperature is not specified in that sentence. Reported: S01, corresponding protocols.

## Ligation

V01 A549: 100 µL with 50 U T4 ligase, 1 µM splint, 1 mM DTT, 1 mM ATP, 0.2 mg/mL BSA, 40 U MRI, 1×T4 buffer, RT 2 h. Mouse: 100 µL with 4000 U Hi-T4 ligase, same reported splint/DTT/ATP/BSA/MRI amounts, 1×Hi-T4 buffer, RT 16 h. V02: 50 µL with 2000 U Hi-T4, 1 µM splint, 1 mM DTT, 1 mM ATP, 0.2 mg/mL BSA, 20 U MRI, 1×Hi-T4 buffer, RT 2 h. Reported: S01, corresponding protocols.

## Amplification

V01 A549: 50 µL RCA mixture, 50 U Phi29 polymerase, 0.5 mM dNTP solution mix as stated, 1 mM DTT, 0.2 mg/mL BSA, 20 U MRI, 1×Phi29 buffer, 30°C 2 h. Mouse: 100 µL, 100 U Phi29, 0.5 mM dNTP, 1 mM DTT, 0.2 mg/mL BSA, 40 U MRI, 1×Phi29 buffer, 30°C 3 h. V02 refers to A549 RCA. Reported: S01, RAEFISH protocols. Do not infer whether the reported dNTP concentration denotes each nucleotide or total.

## Polymerization

Hydrogel embedding/polymerization is not included in the reported RAEFISH sample protocol. DNA polymerization occurs during RCA; see Amplification. S01, RAEFISH protocols and Image system.

## Sequencing / Readout

V01: 47 dual-color rounds. Each round: 10 nM each Alexa Fluor 750 and ATTO 647N common oligo plus 6.7 nM each of two adapters in 20% v/v ethylene carbonate/2×SSC, RT 20 min; wash without oligos for 7 min. Remove/block signals using 65% formamide/2×SSC with 1 mM blocking oligos for 24 min, preserving the source-reported unit. V02: 14 bits/7 rounds, directly labeled Alexa Fluor 750/Alexa Fluor 647 oligos at 3 nM each. Paired RNA MERFISH is 26 bits/13 rounds and imaged before V02. Reported: S01, Sequential secondary probe hybridization and imaging.

## Imaging

Imaging buffer: 50 mM Tris-HCl pH 8.0, 10% w/v glucose, 2 mM Trolox, 0.5 mg/mL glucose oxidase, 40 µg/mL catalase in 2×SSC. Far-red/near-infrared channels: 647/638 nm and 750/749 nm; fiducial beads: 488/477 nm. A549: 15 z positions at 0.6 µm, depth 8.4 µm; tissue: 21 positions at 0.6 µm, depth 12 µm. Nikon 60× oil objective, NA 1.40; 107.9-nm pixels. DAPI after V01 readout at 1:2000. Reported: S01, Sequential secondary probe hybridization and imaging and Image system.

## Decoding

Fiducial drift correction, morphological background removal (disk radius 4 or 5 pixels), normalized pixel-code matching to 25,281 possible codes; distance threshold 0.5176 in A549 and 0.65 in tissue. Foci are filtered using area, code distance and vector magnitude with empty-code estimates. S01 inconsistently calls these “94 rounds” and a “96-dimensional” vector despite specifying 94 bits and 47 physical rounds elsewhere. The inspected mouse archive decoder sets RoundsOfHybs=94 and uses that vector size (inferred: S03, RAEFISH/runM2_2_HPC.m); the manuscript wording remains discrepant and is not silently rewritten. V02 uses a separate segmentation/decoding workflow. Reported: S01, Image analysis; S03, code inventory.

## Critical Parameters

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Primary recognition | Padlock and priming probes | 3 nM each | 37°C | 16 h | 2×SSC, 50% formamide | Not reported | 1 | Not applicable | Not applicable | V01 A549 | reported | S01, A549 protocol |
| Ligation | T4 ligase; splint | 50 U in 100 µL; 1 µM | RT | 2 h | 1×T4 ligase buffer; additives above | Not reported | 1 | Not applicable | Not applicable | V01 A549 | reported | S01, A549 protocol |
| RCA | Phi29 | 50 U in 50 µL | 30°C | 2 h | 1×Phi29 buffer; additives above | Not reported | 1 | Not applicable | Not applicable | V01 A549 | reported | S01, A549 protocol |
| RCA | Phi29 | 100 U in 100 µL | 30°C | 3 h | 1×Phi29 buffer; additives above | Not reported | 1 | Not applicable | Not applicable | V01 mouse | reported | S01, mouse protocol |
| Readout hybridization | Common dyes; adapters | 10 nM each dye; 6.7 nM each adapter | RT | 20 min | 20% ethylene carbonate in 2×SSC | Not reported | 47 rounds | 2 channels/round | Not applicable | V01 | reported | S01, Sequential secondary probe hybridization and imaging |
| Readout hybridization | Directly labeled oligos | 3 nM each | RT as shared readout procedure | 20 min as shared procedure | Secondary hybridization buffer | Not reported | 7 rounds | 2 channels/round | Not applicable | V02 | reported | S01, Sequential secondary probe hybridization and imaging |

## Version-Specific Differences

V02 changes primary formamide, oligo concentrations, ligase conditions and readout scheme. Mouse V01 differs from A549 V01 in fixation, permeabilization, hybridization volume, ligation and RCA. These are separate sample protocols within V01, not evidence for a chronological method release. S01, STAR Methods.

## Unreported Parameters

GC/Tm design cutoffs, full probe-manufacturing PCR cycle/reaction details and wash pH are not provided in the examined S01 procedure; synthesis cites an earlier workflow. No defaults were inserted. The 94/96 decoding dimensionality and splint length inconsistency are explicit open questions. See probes.md.
