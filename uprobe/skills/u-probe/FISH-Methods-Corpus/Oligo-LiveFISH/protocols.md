# Experimental Protocol

## Protocol Scope

V01, original 2025 Oligo-LiveFISH. The live-cell workflow is separated from sequential fixed-cell DNA FISH specificity validation and smiFISH expression assays. Reported: S01, STAR Methods. Values below are selected reported execution conditions, not substituted defaults or an independently validated operating procedure.

## Sample Preparation

U2OS/HeLa are cultured in DMEM with GlutaMax and 10% FBS. Human T cells and hESC-derived neurons require their own preparation workflows. Live-cell fRNP delivery uses 2–4 × 10^5 suspended cells and the Neon 10-µL system. Reported: S01, Cell lines and Assembly and live cell delivery; V01.

## Fixation

Not applicable to the core live-cell tracking procedure. In the separate sequential DNA FISH experiment, 4–6 h after electroporation cells are washed with PBS, fixed with 4% formaldehyde for 10 min at room temperature, washed twice, and imaged for guide signals before subsequent DNA FISH. This assay must not be presented as live imaging. S01, Sequential (Oligo-)LiveFISH-FISH; V01 (reported).

## Permeabilization

Electroporation delivers core live-cell fRNPs; no fixed-cell permeabilization is prescribed for that step. The separate sequential DNA FISH assay uses 70% ethanol at −20°C overnight after recording guide signals. Reported: S01, Assembly and Sequential FISH; V01.

## Reverse Transcription

Not applicable to genomic guide labeling. T7 IVT makes RNA from DNA; it is not reverse transcription. qRT-PCR is a separate expression-validation assay. S01, RNA synthesis and qRT-PCR sections; V01.

## Probe Hybridization

Dual-RNA annealing: 20 mM HEPES, pH 7.5, 150 mM KCl; 2:1–5:1 molar ratio with unlabeled component in excess, or equimolar if both labeled. Incubate 95°C 5 min, 70°C 10 min, gradually cool to room temperature, add 1 mM MgCl2, incubate 40°C 5 min, cool gradually. sgRNA folding omits the 95°C step and starts at 70°C 10 min. fRNP assembly uses room-temperature incubation for 10 min before delivery. S01, RNA preparation and Assembly; V01 (reported).

## Washing

IVT products are DNase-I-treated for 20 min and purified with RNA Clean & Concentrator-25. Labeled RNA is purified after modification/click reactions. These are reagent purification steps. Core cell-medium exchange and the T-cell single wash are distinct from stringent fixed-cell FISH washes. S01, RNA synthesis, Dye labeling, Assembly; V01 (reported).

## Ligation

Not applicable to the described guide-pool labeling workflow; no probe ligation step. S01, preparation through assembly; V01.

## Amplification

DNA oligo pools from Twist or IDT are PCR-amplified using Q5 Hot Start High-Fidelity 2X Master Mix and purified from agarose gels. HiScribe T7 Quick High Yield RNA Synthesis Kit produces guide RNA. Full PCR cycling and IVT incubation settings are not supplied in the examined preparation text; do not import kit defaults. S01, DNA oligo pools and RNA synthesis; S02 Table S1; V01 (reported/unavailable).

## Polymerization

Yeast poly(A) polymerase is used for optional 3′ labeling, not tissue embedding: 1–10 µM RNA, 2 mM 2′-Azido-2′-dATP, 1200 U PAP, 1× PAP reaction buffer, 37°C overnight. After purification, react with 1 mM DBCO-dye and 1 U/µL RNase inhibitor at 37°C for 1 h, then purify. S01, Dye labeling; V01 (reported).

## Sequencing / Readout

Direct fluorescence from guide complexes; no sequencing cycles. Assembly uses 0.5–10 pmol labeled crRNA:tracrRNA or sgRNA, equal molar dCas9-EGFP or 1.2–1.5× molar dCas9. U2OS and HeLa: 1400 V, 15 ms, 4 pulses; activated T cells: 1400 V, 10 ms, 3 pulses; hESC iNs: 1300 V, 20 ms, 3 pulses. S01, Assembly and live cell delivery; V01 (reported).

## Imaging

Exchange U2OS/HeLa growth medium to BrightCell MEMO Photostable Media with 10% FBS. Hoechst 33342 is reported at 0.1 µg/mL, 6–8 h after transfection or immediately before imaging. Microscope A1 uses 37°C, humidified 5% CO2; A2 uses 37°C stage incubation, approximately 3–10 W/cm² illumination and 18- or 50-ms exposures in specified experiments. DHPSF microscope C provides two-color 3D tracking; do not assign A2 illumination to C. S01, Assembly and Microscopy imaging; S02 Table S3; V01 (reported).

## Decoding

Localization and trajectory analysis use Fiji, PYME and MATLAB. The paper’s punctate-signal SNR metric and trajectory filters matter for its performance comparisons. No barcode decoding algorithm is required. S01, Image processing and quantitative analysis; V01 (reported).

## Critical Parameters

Reported execution records; “NR” means not reported in the cited step, “NA” means inapplicable to that operation. All rows are V01. Do not combine alternative rows into one synthetic recipe.

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RNA annealing | crRNA + tracrRNA | Unlabeled:labeled 2:1–5:1; equimolar if both labeled | 95°C then 70°C | 5 min then 10 min | 20 mM HEPES, 150 mM KCl | 7.5 | NA | NA | NA | V01 dual RNA | reported | S01 RNA preparation |
| Post-cooling folding | RNA | 1 mM MgCl2 added | 40°C | 5 min then gradual cooling | Folding buffer | 7.5 | NA | NA | NA | V01 dual RNA / sgRNA | reported | S01 RNA preparation |
| 3′ modification | RNA, azido-dATP, PAP | RNA 1–10 µM; nucleotide 2 mM; PAP 1200 U | 37°C | Overnight | 1× PAP buffer | NR | NA | NA | NA | V01 optional labeling | reported | S01 Dye labeling |
| Click labeling | DBCO-dye, RNase inhibitor | 1 mM dye; 1 U/µL inhibitor | 37°C | 1 h | NR | NR | NA | NA | NA | V01 optional labeling | reported | S01 Dye labeling |
| fRNP assembly | Labeled guide + dCas9 | Guide 0.5–10 pmol; dCas9-EGFP 1× or dCas9 1.2–1.5× molar | Room temperature | 10 min | NR | NR | NA | NA | NA | V01 | reported | S01 Assembly |
| Electroporation | fRNPs, 2–4 × 10^5 cells | Delivered in Neon 10-µL system | NR | Pulse 15 ms | Kit conditions, exact composition NR | NR | 4 pulses, 1400 V | NA | NA | V01 U2OS/HeLa | reported | S01 Assembly |
| Electroporation | fRNPs | Same system; dose optimized | NR | Pulse 10 ms | Kit conditions | NR | 3 pulses, 1400 V | NA | NA | V01 activated T cells | reported | S01 Assembly |
| Electroporation | fRNPs | Same system; dose optimized | NR | Pulse 20 ms | Kit conditions | NR | 3 pulses, 1300 V | NA | NA | V01 hESC iNs | reported | S01 Assembly |
| FOS stimulation | ZnSO4·7H2O | 400 µM | Culture conditions | Added 3.5 h before imaging | Cell medium | NR | NA | NA | NA | V01 FOS experiment only | reported | S01 Assembly |
| Sequential DNA FISH fixation | Formaldehyde | 4% | Room temperature | 10 min | Solution composition otherwise NR | NR | NA | NA | NA | V01 fixed validation only | reported | S01 Sequential FISH |

For the alternative 5′ EDC labeling procedure, S01 specifies 120 mM GMP solution in IVT, Cy3-amine/imidazole chemistry or an amination/AF647 NHS-ester route, each followed by purification. These alternatives are not the universal dye-tracrRNA synthesis procedure. Reaction-volume and multistep precipitation details remain in that source; no combined 3′/5′ protocol is inferred.

## Version-Specific Differences

All conditions above are V01. Dual RNA versus sgRNA changes folding and labeling; cell type changes electroporation settings; SL imaging changes optical acquisition. Sequential DNA FISH and smiFISH are validation branches, not versions of live-cell guide labeling. S01, relevant STAR Methods subsections; S02 Tables S1/S3.

## Unreported Parameters

The cited core preparation text does not specify complete PCR cycling, IVT incubation duration/temperature, every reaction volume, exact fRNP assembly buffer, cooling rates, or single-RNA dye stoichiometry. Do not invent these from reagent names. The concentration range is not a universal optimum. S01, preparation/labeling/assembly and Limitations; V01 (unavailable within examined scope).
