# Experimental Protocol

## Protocol Scope

Numerical execution conditions below come from the inspected November 2024 bioRxiv v2 [S05], the V01 precursor. They are not asserted to be the final 2025 Methods. Final-publication oligo/pooling records are independently available in S02/S03 for V01/V02. No unified protocol is synthesized across these scopes.

## Sample Preparation

S05, Mice / Human: C57BL/6J mice older than 8 weeks; embryos E12.5, E13.5 and E14.5; human HCC. Thin samples were O.C.T.-embedded, snap frozen, stored at −80 °C and sectioned at 10 µm. Thick brain was perfusion-fixed, sectioned at 100 µm and stored in 70% ethanol at 4 °C. Human tissue was initially placed in RPMI 1640 before O.C.T. processing. [reported; V01 precursor]

## Fixation

Separate harvest, section and post-RCA fixation conditions appear in Critical Parameters. Preserve the source term “Glutaric dialdehyde”. [reported; S05, Mice / Sample treatment for thin tissue; V01 precursor]

## Permeabilization

Pepsin/HCl treatment applies to thin sections. Thick sections use ethanol rehydration; do not transfer the thin pepsin treatment to them. [reported; S05, thin/thick tissue sections; V01 precursor]

## Reverse Transcription

Not applicable to the described direct-RNA padlock workflow: no cDNA conversion step is used in S05, Sample treatment for thin/thick tissue. [reported workflow scope; V01 precursor]

## Probe Hybridization

Oligo-dT blocking precedes padlock hybridization. Thin samples use staged temperatures; thick samples use prolonged incubation. See parameter records. [reported; S05, thin/thick tissue; V01 precursor]

## Washing

PBST is 0.25% Tween in 1X PBS. The defined washing buffer contains 10% formamide and 2X SSC. Thin post-hybridization washing is three times for 10 min, followed by two PBST rinses; thick is three times for 20 min, followed by three PBST washes. See S05 for intervening washes around ligation, RCA and staining. [reported; S05, thin/thick tissue; V01 precursor]

## Ligation

SplintR ligates the padlock nick after RNA-target hybridization. Ligation duration differs between thin and thick tissue. [reported; S05, thin/thick tissue; V01 precursor]

## Amplification

Phi29 RCA uses external primers and aminoallyl-dUTP. Post-amplification BS(PEG)9 fixation is followed by formamide washing; thick tissue additionally receives Cubic L treatment. [reported; S05, thin/thick tissue; V01 precursor]

## Polymerization

Not applicable: S05 describes enzymatic DNA polymerization during RCA, already recorded under Amplification, but no hydrogel polymerization step. [workflow scope; V01 precursor]

## Sequencing / Readout

No sequencing cycles are used for identity readout. Competitive fluorescent/unlabeled imaging probes stain RCA products in one round. S02 supplies final V01/V02 pair recipes. Its theoretical and actual mixtures are different; use the actual recipe appropriate to the instrument and configuration. [reported; S01, Fig. 1c; S02, imaging sheets; V01/V02]

## Imaging

S05, Imaging: thin tissue used a 40X NA 0.60 objective, 9 z-planes at approximately 1-µm intervals, 4% tile overlap, 0.1625 × 0.1625 µm/pixel. Cy3/Cy5 and AF488/TxRed/DAPI channels were serially acquired. Thick tissue used Zeiss LSM880/980 with 0.207 × 0.207 × 0.700 µm/voxel sampling. [reported; V01 precursor; these are instrument examples, not universal requirements]

## Decoding

S05, Image processing / Gene decoding: focal stacking, illumination correction, inter-channel registration and stitching precede spot calling. Thin tissue uses top-hat/local maxima; thick tissue uses Airlocalize 3D fitting. Normalize channels and classify in color space. The optional Ch4-only barcode is handled separately because the other-channel denominator is zero. [reported; V01 precursor]

Coordinate conventions differ: S05 gives x=(Ch1−Ch2)/sum, y=2Ch3/sum−1; S01 Extended Data Fig. 2 gives X=2Ch3′−1, Y=Ch2′−Ch1′. These are different axis definitions and must not be mixed. [reported; V01]

## Critical Parameters

All records preserve source wording/units; blank experimental values were not filled. In particular, the preprint DAPI concentration is retained as printed and has not been confirmed against the final Methods.

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Harvest fixation | PFA | 4% (w/v) | 4 °C | Overnight | Not reported | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Mice |
| Section fixation | PFA; Glutaric dialdehyde | 4%; 0.1% | Room temperature | 15 min | Not reported | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Permeabilization | Pepsin | 0.01% | 37 °C | 2 min; 10 min for mouse embryo | 0.1 M HCl | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Dehydration | Ethanol | 80%; pure | Not reported | 10 min; 2 min respectively | Not applicable | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Blocking | Oligo-dT; KCl; BSA; Yeast tRNA; Ribolock | 100 nM; 50 mM; 20 µg/mL; 20 µg/mL; 1 U/µL | Room temperature | 10 min | Ampligase buffer; 20% formamide | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Padlock hybridization | Each padlock; KCl; BSA; Yeast tRNA; Ribolock | 200 nM; 50 mM; 20 µg/mL; 20 µg/mL; 1 U/µL | 55 °C then 45 °C | 10 min then 1 h 50 min | Ampligase buffer; 20% formamide | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Post-hybridization wash | Formamide; SSC | 10%; 2X | Not reported | 10 min per wash | Washing buffer | Not reported | 3 washes | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Ligation | SplintR ligase; BSA; Ribolock | 2.5 U/µL; 20 µg/mL; 1 U/µL | 37 °C | 2 h | SplintR buffer | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| RCA | Phi29; dNTP; aminoallyl-dUTP; glycerol; BSA; RCA primer | 0.25 U/µL; 250 µM; 50 µM; 10%; 20 µg/mL; 600 nM | 30 °C | Overnight | Phi29 polymerase buffer | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Post-RCA fixation | BS(PEG)9 | 10 µg/µL | Room temperature | 15 min | PBST | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Post-RCA formamide wash | Formamide | 65% | Room temperature | 5 min per wash | Not otherwise specified | Not reported | 3 washes | Not applicable | Not applicable | V01 precursor / 10-µm tissue | reported | S05, Materials and Methods: Sample treatment for thin tissue |
| Rehydration | PBST | 0.25% Tween in 1X PBS | Not reported | 20 min each | PBST | Not reported | 3 incubations | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Blocking | Oligo-dT mixture | Same composition as thin-tissue block | Room temperature | 30 min | As thin-tissue block | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Padlock hybridization | Padlock mixture | As thin-tissue mixture | 45 °C | 24 h; gentle shaking | As thin-tissue mixture | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Post-hybridization wash | Washing buffer | 10% formamide; 2X SSC | Not reported | 20 min each | Washing buffer | Not reported | 3 washes; then 3 PBST washes | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Ligation | Ligation mixture | As thin-tissue mixture | 37 °C | Overnight; gentle shaking | SplintR buffer | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| RCA | RCA mixture | As thin-tissue mixture | 30 °C | 12–24 h | Phi29 polymerase buffer | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Post-RCA fixation | BS(PEG)9 | As thin-tissue mixture | Not separately specified | 30 min | PBST | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Lipid removal | Cubic L (TCI T3740) | Not reported | Room temperature | 1 h | Not reported | Not reported | Not reported | Not applicable | Not applicable | V01 precursor / 100-µm mouse brain | reported | S05, Materials and Methods: Sample treatment for thick tissue |
| Thin imaging-probe stain | 17 imaging-probe pairs | 120 nM per pair | 50 °C then 37 °C | 9 min then 21 min | 20% formamide; 2X SSC | Not reported | Not applicable | 1 | Not applicable | V01 precursor / thin tissue | reported | S05, Imaging probe staining |
| Thick imaging-probe stain | Imaging probes | As thin-tissue mixture | 37 °C | 2 h | As thin-tissue mixture | Not reported | Not applicable | 1 | Not applicable | V01 precursor / 100-µm tissue | reported | S05, Imaging probe staining |
| Thin nuclear stain | DAPI | 5 µg/µL (source notation) | Room temperature | 5 min | Not reported | Not reported | Not applicable | 1 | Not applicable | V01 precursor / thin tissue | reported | S05, Imaging probe staining |
| Thick nuclear stain | DAPI | Not separately specified | Room temperature | 30 min | Not reported | Not reported | Not applicable | 1 | Not applicable | V01 precursor / 100-µm tissue | reported | S05, Imaging probe staining |
| Refractive-index matching | Iohexol | 0.9 M | Room temperature | 1 h | 2X SSC | Not reported | Not applicable | 1 | Not applicable | V01 precursor / 100-µm tissue | reported | S05, Imaging probe staining |

S02, RCA-primer sheets: each primer is prepared at 10 µM and pooled as an RCA-primer stock [reported; V01/V02]. This does not establish a final reaction concentration for the final paper.

## Version-Specific Differences

Final V01 and V02 have different imaging/primer libraries and empirical pooling recipes [S02]. The 2024 v2 protocol is not silently promoted to V02. S03 supplies sample-specific padlock dilution ratios; it is not a general recipe. [reported; S02/S03]

## Unreported Parameters

Final 2025 full Methods were subscription-restricted in the inspected page. Final-paper reaction conditions, complete 64-panel execution conditions, reaction volumes, pH, source-confirmed resolution of the preprint DAPI units, and publication-time design settings remain unresolved. The separate FFPE demonstration in S04 Fig. 39 does not establish a complete FFPE protocol here.
