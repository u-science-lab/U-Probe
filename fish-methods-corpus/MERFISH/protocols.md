# Experimental Protocol

## Protocol Scope

Evidence extraction for V01 original IMR90 measurements, with separate V02 high-throughput, V03 clearing, V04 expansion, V05 bDNA and V06 high-plex scopes. This is not a synthetic protocol assembled across publications. S03 is an additional detailed developer protocol; its values do not silently replace S01. Table records below are reported unless explicitly noted. Full equipment and reagent purchasing inventories remain in the original methods.

## Sample Preparation

V01 uses human primary fibroblasts IMR90 (ATCC), plated at 350,000 cells per 22-mm #1.5 coverslip and cultured in Eagle’s Minimum Essential Medium at 37 °C, 5% CO2 for 48–96 h (S01, Sample preparation and labeling; reported). V02 uses U-2 OS cells for the high-throughput library; IMR-90 is also used for tests (S04 SI Cell Culture). V03 includes adult mouse hypothalamus cryosections, 10 µm thick; it is not an IMR90-only protocol (S05 Fig. 4 and SI Tissue Preparation; reported).

## Fixation

V01: 4% paraformaldehyde in 1× PBS, 20 min at room temperature, then sodium borohydride reduction and washes. After encoding hybridization and fiducial addition, postfix with 4% v/v paraformaldehyde in 2× SSC for 30 min (S01 Sample preparation and labeling; reported). S03 §5.2 instead specifies 15 min initial fixation; retain this as its protocol-specific value, not a correction to S01.

## Permeabilization

V01: 0.5% v/v Triton in PBS for 2 min at room temperature, followed by three ice-cold PBS washes (S01 Sample preparation; reported). V02: 10 min in 0.5% v/v Triton X-100 (S04 SI Cell Culture; reported). V03 tissue adds SDS/ethanol treatments; do not apply the cultured-cell timing to tissue (S05 SI Tissue Preparation; reported).

## Reverse Transcription

Not applicable to target RNA conversion in the unamplified in situ workflow: target RNA is directly hybridized by DNA probes. Reverse transcription occurs during ex situ encoding-probe manufacture; S01 Step 3 and the table specify its conditions. V01 RNA template is subsequently removed by alkaline hydrolysis. This is distinct from an in situ cDNA intermediate (S01 Probe synthesis; reported).

## Probe Hybridization

V01 encoding stain: 10 µL pool, 100 µM for 140-gene libraries or 200 µM for 1001-gene library; 18–36 h at 37 °C. Readout stain: 10 nM, 15 min per round. Full buffer compositions are in Critical Parameters (S01 Sample preparation / MERFISH imaging; reported). V02 uses ~200 µM encoding pool, batch titration, 36–48 h, followed by two-color 1 nM readout probes (S04 SI; reported).

## Washing

V01 post-encoding: three 10-min washes at 47 °C in encoding wash buffer. Per-round readout wash: 3 min in 20% v/v formamide/2× SSC/2 mM VRC. V02 post-encoding: two 30-min washes at 47 °C; readout uses an ethylene-carbonate buffer (S01 imaging/preparation; S04 SI; reported). V05 amplifier wash sequences are separate from encoding and readout washes (S08 Amplifier staining).

## Ligation

Not applicable to the described V01 direct RNA hybridization/fluorescence workflow; no target-dependent ligation step is used (S01 Fig. 1 and Methods). bDNA binding is hybridization assembly, not padlock ligation (S08 Fig. 1).

## Amplification

V01 encoding library manufacture uses PCR → IVT → RT → hydrolysis ex situ. PCR cycle count is not a universal fixed number; S03 §4.2 recommends following amplification kinetics. V05 adds primary and secondary amplifier staining in situ; it uses 5×5 amplifiers for multiplexed measurements. S08’s 4×4 and 9×9 tests are separate brightness experiments (S01 Probe synthesis; S03 §4; S08 Results and Methods; reported).

## Polymerization

Not applicable to original V01. V03 embeds RNA-anchor-stained samples in polyacrylamide and clears them by proteinase K/detergent digestion; detailed conditions are in the table (S05 SI Sample Embedding and Clearing; reported). V04 uses an expandable gel and re-embedding, a different chemistry from V03’s nonexpanding gel. S10 Methods describes 7.7% wt/wt sodium acrylate, 4% v/v 19:1 acrylamide/bis, 2 M NaCl and 60 mM Tris-HCl pH 8 in its monomer solution, followed by expansion in 0.5× SSC plus 0.2% v/v proteinase K with buffer exchange every 30 min until expansion stops, typically ~2 h (reported).

## Sequencing / Readout

No nucleotide sequencing cycles apply to these fluorescence barcode measurements. V01 reads one bit per fluorescent hybridization, then photobleaches; 16 rounds MHD4 or 14 rounds MHD2. V02 reads two bits per round over eight rounds and removes dyes by TCEP cleavage (S01 imaging Methods; S04 SI MERFISH Imaging; reported). V06 reads 69 bits with three-color imaging over 23 rounds (S09 Results; reported).

## Imaging

V01: automated wide-field oblique-incidence imaging, 100× 1.45-NA objective, 167-nm sample-plane pixels and 40×40-µm regions. Approximately 75–100 regions per experiment are imaged in the described workflow. Fiducial beads register rounds. Its single focal plane does not imply complete volume capture (S01 MERFISH imaging; reported). V03 tissue acquires seven ~1.5-µm optical sections per FOV for 10-µm sections (S05 Fig. 4; reported).

## Decoding

V01: multi-Gaussian spot fits; bead-based rotation/translation; connect spots within one pixel (167 nm); construct on/off words. Accept exact or one-bit corrected MHD4 words, but only exact MHD2 words (S01 Construction of measured words; reported). V06 uses voxel-based intensity decoding and aggregation instead; do not universalize the original spot-link threshold (S09 Materials and Methods; reported).

## Critical Parameters

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fixation | Paraformaldehyde | 4% | Room temperature | 20 min | 1× PBS | Not specified | Not applicable | Not applicable | Not applicable | V01 IMR90 | reported | S01 Sample preparation and labeling |
| Background reduction | Sodium borohydride | 0.1% w/v | Not specified | 5 min | Water | Not specified | Not applicable | Not applicable | Not applicable | V01 IMR90 | reported | S01 Sample preparation and labeling |
| Permeabilization | Triton | 0.5% v/v | Room temperature | 2 min | 1× PBS | Not specified | Not applicable | Not applicable | Not applicable | V01 IMR90 | reported | S01 Sample preparation and labeling |
| Encoding hybridization | Encoding probe pool | 100 µM (140-gene); 200 µM (1001-gene) | 37 °C | 18–36 h | 10 µL; 2× SSC, 30% v/v formamide, 2 mM VRC, 1 mg/mL yeast tRNA, 10% w/v dextran sulfate | Not specified | One encoding incubation | Not applicable | Not applicable | V01 IMR90 | reported | S01 Sample preparation and labeling |
| Encoding washes | Encoding wash buffer | 2× SSC; 30% v/v formamide; 2 mM VRC | 47 °C | 10 min each | Encoding wash buffer | Not specified | 3 washes | Not applicable | Not applicable | V01 IMR90 | reported | S01 Sample preparation and labeling |
| Postfixation | Paraformaldehyde | 4% v/v | Room temperature | 30 min | 2× SSC | Not specified | One | Not applicable | Not applicable | V01, after bead staining | reported | S01 Sample preparation and labeling |
| Readout hybridization | Cy5 readout probe | 10 nM | Heated objective 37 °C; solution temperature not separately given | 15 min | 1 mL; 2× SSC, 10% v/v formamide, 10% w/v dextran sulfate, 2 mM VRC | Not specified | Per round | 16 MHD4 / 14 MHD2 | Not applicable | V01 | reported | S01 MERFISH imaging |
| Readout wash | Readout wash buffer | 2× SSC; 20% v/v formamide; 2 mM VRC | Not separately specified | 3 min | 2 mL wash buffer | Not specified | Per round | 16 / 14 | Not applicable | V01 | reported | S01 MERFISH imaging |
| Imaging buffer | Tris-HCl; glucose; Trolox; glucose oxidase; catalase | 50 mM; 10% w/v; 2 mM; 0.5 mg/mL; 40 µg/mL | Objective 37 °C | Not specified | 2 mL; 2× SSC | 8 (Tris-HCl) | Per round | 16 / 14 | Not applicable | V01 | reported | S01 MERFISH imaging |
| Photobleaching | 641-nm excitation | 200 mW at microscope backport | Not separately specified | 3 s per region | 2× SSC; 2 mM VRC | Not specified | Per round | 16 / 14 | Not applicable | V01 | reported | S01 MERFISH imaging |
| Probe manufacture: IVT | T7 polymerase; template; NTPs | 2 µL enzyme; ~1 µg template; 10 mM each NTP | 37 °C | 4 h | 20 µL; 1× reaction buffer; 1× RNase inhibitor | Not specified | One reaction | Not applicable | Not applicable | V01 ex situ | reported | S01 Probe synthesis Step 2 |
| Probe manufacture: RT | Maxima H−; RT primer; dNTPs | 300 units; 2 nmol; 1.6 mM each dNTP | 50 °C; then 85 °C | 45 min; then 5 min inactivation | 50 µL; 1× RT buffer; 60 units RNasin | Not specified | One reaction | Not applicable | Not applicable | V01 ex situ | reported | S01 Probe synthesis Step 3 |
| Probe manufacture: hydrolysis | EDTA and NaOH mixture | Added 20 µL of 0.25 M EDTA / 0.5 N NaOH | 95 °C | 10 min | Added to RT reaction; not asserted as final concentration | Not specified | One | Not applicable | Not applicable | V01 ex situ | reported | S01 Probe synthesis Step 4 |
| Permeabilization | Triton X-100 | 0.5% v/v | Room temperature | 10 min | 1× PBS | Not specified | One | Not applicable | Not applicable | V02 U-2 OS | reported | S04 SI Cell Culture and Encoding Probe Staining |
| Encoding hybridization | Encoding probe pool | ~200 µM; titrated per batch | 37 °C | 36–48 h | 30 µL; 2× SSC, 30% formamide, 2 mM VRC, 0.1% w/v yeast tRNA, 1% v/v murine RNase inhibitor, 10% w/v dextran sulfate | Not specified | One | Not applicable | Not applicable | V02 U-2 OS | reported | S04 SI Cell Culture and Encoding Probe Staining |
| Encoding washes | Encoding wash buffer | As above without supplements | 47 °C | 30 min each | Encoding wash buffer | Not specified | 2 washes | Not applicable | Not applicable | V02 U-2 OS | reported | S04 SI Cell Culture and Encoding Probe Staining |
| Readout hybridization | Cy5 and Alexa750 probes | 1 nM each | Not separately specified in procedural paragraph | 5 min exchange + 6 min continuous flow | 2 mL + 2 mL; 2× SSC, 1% v/v ethylene carbonate, 10% w/v dextran sulfate, 2 mM VRC | Not specified | Per round | 8 two-color rounds | Not applicable | V02 | reported | S04 SI MERFISH Imaging |
| Readout wash | Ethylene carbonate | 10% v/v | Not specified | 9 min | 2 mL; 2× SSC; 2 mM VRC | Not specified | Per round | 8 | Not applicable | V02 | reported | S04 SI MERFISH Imaging |
| Reductive cleavage | TCEP hydrochloride | 50 mM | Not specified | 4 min initial flow; then 15 min at 0.1 mL/min | 2× SSC; initial volume unresolved source wording | Not specified | Per round | 8 | Not applicable | V02 | reported | S04 SI MERFISH Imaging |
| Gel anchoring | Alternating dT/LNA anchor | 1 µM where appropriate | 37 °C encoding incubation | 36 h encoding incubation | Encoding buffer as S05; encoding probes 100 µM unless specified | Not specified | One | Not applicable | Not applicable | V03 cultured cells | reported | S05 SI Staining with Encoding Probes |
| Gel polymerization | 19:1 acrylamide/bis; APS; TEMED | 4% v/v; 0.03% w/v; 0.15% v/v | Room temperature | 1.5 h | 60 mM Tris-HCl; 0.3 M NaCl; fiducial beads | 8 | One | Not applicable | Not applicable | V03 cultured cells | reported | S05 SI Sample Embedding and Clearing |
| Clearing digestion | Proteinase K | 1% v/v stock (source catalog P8107S) | 37 °C | >12 h | 0.8 M guanidine-HCl; 50 mM Tris-HCl; 1 mM EDTA; 0.5% v/v Triton X-100 | 8 | One | Not applicable | Not applicable | V03 cultured cells | reported | S05 SI Sample Embedding and Clearing |
| Corrected encoding stain | Encoding probe pool; anchor | ~50 µM; 3.3 µM | See S10 staining protocol | See S10 staining protocol | 30 µL encoding hybridization solution | Not specified | One | Not applicable | Not applicable | V04; corrected from ~300 µM | reported | S11 Author Correction, Encoding probe staining |
| Primary amplifier staining | Each primary amplifier | 5 nM | 37 °C | 30 min for MERFISH | 50 µL; 2× SSC, 10% v/v formamide, 0.1% w/v yeast tRNA, 1% v/v murine RNase inhibitor, 10% w/v dextran sulfate | Not specified | One | Before readout | Not applicable | V05 gel-cleared samples | reported | S08 Amplifier staining |
| Secondary amplifier staining | Each secondary amplifier | 5 nM | 37 °C | 30 min for MERFISH | 50 µL; same amplifier hybridization buffer | Not specified | One | Before readout | Not applicable | V05 gel-cleared samples | reported | S08 Amplifier staining |

## Version-Specific Differences

V01 uses photobleaching; V02 uses disulfide cleavage and two-color imaging. V03 adds matrix anchoring/clearing; V04 expands samples; V05 adds bDNA; V06 uses more barcode bits and overlapping probes for a subset. V04’s published encoding concentration ~300 µM is corrected explicitly to ~50 µM by S11. Both values are retained here as original and corrected, not competing unlabelled conditions.

S03 §3.1 threshold wording conflicts with S01 and must not be silently converted to executable design constraints (see probes.md). S04 SI MERFISH Imaging says “Two microliters” for an initial cleavage flush followed by a flow rate of 0.1 mL/min; this initial volume is retained as unresolved source wording, not changed to milliliters by assumption.

## Unreported Parameters

Unknown values in the table mean the cited procedural passage does not specify a separate value. No pH, exposure, fixed PCR cycle number, reagent stock activity, or thermal condition is supplied from common practice. Complete V04/V06 supplemental protocol extraction is unfinished; V04 supplementary PDF was downloaded but not fully examined, and V06 SI binary access was not completed. See sources.md and research_manifest.md for exact inspection scope.
