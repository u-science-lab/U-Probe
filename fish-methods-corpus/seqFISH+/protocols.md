# Experimental Protocol

## Protocol Scope

V01: Eng et al. (2019), 10,000-gene RNA seqFISH+ in NIH/3T3 cells and mouse brain. S01 Methods is the principal experimental record. S05 identifies the author supplementary protocol; its procedure was inspected only through a mirrored PDF. The current official V.1 page returned no readable procedure during Fix. Mirror-only claims remain inferred pending official-body authentication. Conditions below follow S01 unless S05 is explicitly named. These are source-grounded records, not a reconstructed replacement for the complete original protocols.

## Sample Preparation

NIH/3T3 cells were grown on functionalized coverslips to approximately 80–90% confluence. Brain experiments used male C57BL/6J mice, P23 for cortex and P40 for olfactory bulb; distinguish the female 10-week-old mice used for separate smFISH validation. Brain was cryoprotected with 30% sucrose in PBS at 4°C until sinking, embedded/frozen in OCT and stored at −80°C. The paper used 5-µm cryosections. [S01, Methods, Cell culture experiment, Mice and Tissue slices experiment; reported; V01]

Coverslip preparation: cells, high-setting plasma cleaning for 5 min; 1% bind-silane in pH 3.5, 10% acidic ethanol for 30 min at room temperature; three ethanol rinses; >90°C for 30 min. Tissue coverslips instead received 1 M HCl for 1 h, then water, 1 M NaOH for 1 h, three water rinses and 1% bind-silane for 1 h at room temperature. Remaining steps followed the cell branch. The source reports poly-D-lysine at 100 µg/µl for >1 h; preserve that printed unit pending experimental clarification rather than silently replacing it. [S01, Methods, Coverslip functionalization; reported; V01]

## Fixation

Cells: freshly prepared 4% formaldehyde in 1× PBS for 10 min at room temperature. Tissue: perfuse 8 min at 4°C with 10 U/ml heparin and 0.5% (w/v) NaNO₂ in 0.1 M PBS, then 8 min with 4% PFA in 0.1 M PBS at 4°C; post-fix extracted brain in 4% PFA for 2 h at room temperature. [S01, Methods, Cell culture experiment and Tissue slices experiment; reported; V01]

## Permeabilization

Cells: 70% ethanol for 1 h at room temperature. Tissue slices: 70% ethanol at 4°C for >1 h, followed by 8% SDS in 1× PBS for 30 min at room temperature. This pre-hybridization tissue clearing is separate from subsequent proteinase K digestion. [S01, Methods, cell/tissue experiments; reported; V01]

## Reverse Transcription

Not applicable as an in situ RNA-recognition step. RT is used during probe manufacture after IVT, with a uracil-containing forward primer. USER enzyme at 1:30 dilution removes the primer over approximately 24 h at 37°C. RNA templates are removed with 1 M NaOH at 65°C for 15 min, followed by 1 M acetic acid neutralization. [S01, Methods, Primary-probe construction; reported; V01]

## Probe Hybridization

Primary buffer: 40% formamide, 2× SSC and 10% (w/v) dextran sulfate. S01 mixes primaries with 2 µM LNA polyT15 and 2 µM LNA polyT30. Cells use approximately 1 nM per primary, 24 probes/gene, for approximately 36–48 h at 37°C in a humid chamber. Tissue primary hybridization is approximately 60 h at 37°C. S01 does not independently restate a tissue per-probe concentration in the Tissue slices experiment paragraph; the S05 mirror gives a general 1 nM/probe recommendation (inferred pending official-body authentication). [S01, Methods, Primary-probe construction and cell/tissue experiments; S05 mirror, Primary probes hybridization; reported (S01); inferred (mirror-only recommendation); V01]

Readout buffer: 10% ethylene carbonate, 10% dextran sulfate, 4× SSC and 1:100 SUPERase In RNase inhibitor; Alexa Fluor 647 and Cy3B readouts at 50 nM each, Alexa Fluor 488 readout at 100 nM. Use 100 µl per serial hybridization through a flow cell requiring approximately 25 µl; incubation at room temperature is 17 min for cells and 20 min for tissue. [S01, Methods, seqFISH+ imaging; reported; V01]

## Washing

After primary hybridization, cells receive 40% formamide in 2× SSC for 30 min at 37°C and three 1-ml 2× SSC rinses. Tissue receives the reported 40% formamide wash for 30 min at 37°C, three 2× SSC rinses and one PBS rinse. The S05 mirror explicitly includes 0.1% Triton X-100 (inferred pending official-body authentication) in the primary wash; do not erase this protocol-versus-paper formulation difference. [S01, Methods, cell/tissue experiments; S05 mirror, Primary probes wash; reported (S01); inferred (mirror-only wash addition); V01]

After each readout hybridization, flow approximately 300 µl of 10% formamide/0.1% Triton X-100 in 2× SSC; rinse with approximately 200 µl 4× SSC containing 1:1,000 SUPERase In. Stain with 10 µg/ml DAPI in 4× SSC plus the same inhibitor dilution for approximately 15 s. After imaging, flow 55% formamide/0.1% Triton X-100 in 2× SSC for 1 min and incubate 1 min, then rinse with 4× SSC. [S01, Methods, seqFISH+ imaging; reported; V01]

## Ligation

Not applicable to the reported RNA seqFISH+ recognition/readout workflow. No ligation or padlock step is specified. Do not transfer DNA seqFISH+ ligation procedures to this method. [S01, complete cell/tissue/imaging workflow; V01]

## Amplification

No in situ enzymatic amplification step. Manufacture: Twist oligoarray pools → limited-cycle PCR → QIAquick purification → IVT → RT → USER cleavage → alkaline RNA removal → ethanol precipitation/phenol–chloroform extraction/Zeba desalting → formulation. PCR cycle count, primer strings and IVT/RT reaction recipes are not restated in S01. Readouts were coupled from 5′-amine DNA: 5 nmol DNA with 25 µg dye in 0.5 M sodium bicarbonate/10% dimethyl fumarate for at least 6 h at 37°C, then purified by ethanol precipitation, HPLC and columns; a 500 nM working stock was prepared. Preserve the solvent name as printed. [S01, Methods, Primary-probe construction and Readout-probe synthesis; reported; V01]

## Polymerization

Cells: 1:1,000 TetraSpeck beads in 2× SSC for 5–10 min at room temperature; rinse; preincubate in degassed 4% acrylamide/2× SSC for 5 min. Apply 10 µl 4% hydrogel containing 0.2% APS and 0.2% TEMED in 2× SSC, sandwich and polymerize 1 h at room temperature under nitrogen. The paper writes crosslinking as 1:19; the S05 mirror identifies acrylamide:bis as 19:1 (inferred pending official-body authentication). Retain the source notation rather than treating them as opposing formulations. [S01, Methods, Cell culture experiment; S05 mirror, Hydrogel preparation; reported (S01); inferred (mirror-only ratio notation); V01]

Clear cells for >1 h at 37°C with 1:100 proteinase K, 50 mM Tris-HCl pH 8, 1 mM EDTA, 0.5% Triton X-100 and 500 mM NaCl. Rinse in 2× SSC, treat with 1:10 Label-IT for 30 min at 37°C, post-fix with 4% PFA in PBS for 15 min at room temperature and quench in 1 M Tris-HCl pH 8 for 10 min at room temperature. [S01, Methods, Cell culture experiment; reported; V01]

Tissue: treat with 0.1 mg/ml Acryoloyl-X succinimidyl ester in PBS for 30 min at room temperature; preincubate in 4% acrylamide/2× SSC for 30 min. Apply 20 µl hydrogel containing 0.05% APS/0.05% TEMED; incubate under nitrogen at 4°C for 30 min, then 37°C for 2.5 h. Digest for >3 h at 37°C with the cell digestion buffer plus 1% SDS. Rinse, treat with 0.1 mg/ml label-X for 45 min at 37°C and re-embed; gelation may be shortened to 2 h. Store prepared cell/tissue samples in 4× SSC plus 2 U/µl SUPERase In at 4°C for no longer than 6 h, or image immediately. [S01, Methods, Tissue slices experiment; reported; V01]

## Sequencing / Readout

Not nucleotide sequencing. Four base-20 barcoding rounds correspond to 80 serial hybridizations, each followed by washing, imaging and readout stripping. Three readout channels are acquired per serial hybridization. The 81st hybridization repeats the first for retention QC. [S01, Methods, encoding strategy and seqFISH+ imaging; Extended Data Fig. 2; reported; V01]

## Imaging

Use anti-bleaching buffer: 10% (w/v) glucose, 1:100 catalase, 0.5 mg/ml glucose oxidase, 0.02 U/µl SUPERase In, 50 mM Tris-HCl pH 8 in 4× SSC. Reported platform: Leica DMi8, Yokogawa CSU-W1, Andor Zyla 4.2 Plus, 63× oil/1.40 NA objective, motorized stage and custom µManager integration. Methods specify two z-slices per FOV, 0.35-µm steps, at 647, 561, 488 and 405 nm. The main text describes the brain measurement as a single optical plane of 0.75-µm thickness; retain this acquisition/analysis wording difference without inventing a projection rule. [S01, main text adjacent to Fig. 3 and Methods, seqFISH+ imaging; reported; V01]

Final segmentation staining: 10 µg/ml DAPI, 50 nM LNA T20–Alexa Fluor 647 and 1:100 Nissl stain in PBS for 30 min at room temperature. [S01, Methods, seqFISH+ imaging; reported; V01]

## Decoding

DAPI plus TetraSpeck beads register cell images; tissue uses DAPI alone. Deconvolve with a bead-derived PSF; manual ImageJ ROIs define cells. Detect local maxima and refine with radialcenter.m. Search symmetric nearest neighbours within 1 pixel for cells and 1.4 pixels for tissue, resolve competing patterns by residual distance and Hamming-distance-one matching, reject remaining ambiguity and require agreement in at least three of four seed rounds. Decode channels separately and use unused barcodes for off-target evaluation. [S01, Methods, Image registration, Image processing and Barcode calling; reported; V01]

## Critical Parameters

All rows are reported V01 values from S01; manufacture, cell/tissue preparation and comparison-assay records are explicitly scoped. S05 mirror-only guidance remains qualified in prose and is not promoted to an authenticated S01 parameter. “Not stated” means the cited paragraph does not give the field; “Not applicable” identifies a different type of operation. Buffer definitions and additional handling are above.

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Fixation | Formaldehyde | 4% | Room temperature | 10 min | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Permeabilization | Ethanol | 70% | Room temperature | 1 h | Not stated | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Permeabilization | Ethanol | 70% | 4°C | >1 h | Not stated | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Pre-clear | SDS | 8% | Room temperature | 30 min | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Primary hybridization | Primary probes | Approximately 1 nM/probe | 37°C | Approximately 36–48 h | 40% formamide, 2× SSC, 10% dextran sulfate | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Primary-probe construction; Cell culture experiment |
| Primary hybridization | Primary probes | Not restated in tissue paragraph | 37°C | Approximately 60 h | Primary buffer described above | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Primary wash | Formamide | 40% | 37°C | 30 min | 2× SSC explicitly given for cells | Not stated | Not applicable | Not applicable | Not applicable | V01 cells/tissue | reported | S01, cell/tissue experiments |
| Cell gelation | Acrylamide; APS; TEMED | 4%; 0.2%; 0.2% | Room temperature | 1 h | 2× SSC; 10 µl | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Tissue gelation | Acrylamide; APS; TEMED | 4%; 0.05%; 0.05% | 4°C then 37°C | 30 min then 2.5 h | 2× SSC; 20 µl | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Digestion | Proteinase K | 1:100 | 37°C | >1 h cells; >3 h tissue | Digestion buffer above; tissue adds 1% SDS | 8 (Tris) | Not applicable | Not applicable | Not applicable | V01 separate sample durations | reported | S01, cell/tissue experiments |
| Readout hybridization | Alexa Fluor 647/Cy3B/Alexa Fluor 488 readouts | 50/50/100 nM | Room temperature | 17 min cells; 20 min tissue | 10% EC, 10% dextran sulfate, 4× SSC, inhibitor 1:100 | Not stated | Not applicable | 1–80; repeat 1 in 81 | Not applicable | V01 | reported | S01, seqFISH+ imaging |
| Readout wash | Formamide; Triton X-100 | 10%; 0.1% | Not separately stated | Not separately stated | 2× SSC; approximately 300 µl | Not stated | Not applicable | After each readout | Not applicable | V01 | reported | S01, seqFISH+ imaging |
| Stripping | Formamide; Triton X-100 | 55%; 0.1% | Not separately stated | 1 min flow +1 min incubation | 2× SSC | Not stated | Not applicable | Between readouts | Not applicable | V01 | reported | S01, seqFISH+ imaging |
| Primer cleavage in manufacture | USER | 1:30 | 37°C | Approximately 24 h | Not restated | Not stated | Not applicable | Not applicable | Not applicable | V01 probe preparation | reported | S01, Primary-probe construction |
| RNA template hydrolysis | NaOH | 1 M | 65°C | 15 min | Neutralize with 1 M acetic acid | Not stated | Not applicable | Not applicable | Not applicable | V01 probe preparation | reported | S01, Primary-probe construction |

| Coverslip cleaning | Plasma | High setting | Not stated | 5 min | Not applicable | Not applicable | Not applicable | Not applicable | Not applicable | V01 cell coverslips | reported | S01, Coverslip functionalization |
| Coverslip cleaning, acid | HCl | 1 M | Room temperature | 1 h | Not stated | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue coverslips | reported | S01, Coverslip functionalization |
| Coverslip cleaning, base | NaOH | 1 M | Room temperature | 1 h | Not stated | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue coverslips | reported | S01, Coverslip functionalization |
| Silanization | Bind-silane | 1% | Room temperature | 30 min | 10% acidic ethanol | 3.5 | Not applicable | Not applicable | Not applicable | V01 cell coverslips | reported | S01, Coverslip functionalization |
| Silanization | Bind-silane | 1% | Room temperature | 1 h | Source refers to cell formulation | 3.5 (referenced cell formulation) | Not applicable | Not applicable | Not applicable | V01 tissue coverslips | reported | S01, Coverslip functionalization |
| Coverslip heat curing | Not applicable | Not applicable | >90°C | 30 min | Not applicable | Not applicable | Not applicable | Not applicable | Not applicable | V01 coverslips | reported | S01, Coverslip functionalization |
| Adhesion treatment | Poly-D-lysine | 100 µg/µl (printed unit) | Room temperature | >1 h | Water | Not stated | Not applicable | Not applicable | Not applicable | V01 coverslips | reported | S01, Coverslip functionalization |
| Coverslip storage | Functionalized coverslips | Not applicable | 4°C | No longer than 2 weeks | Dry storage | Not applicable | Not applicable | Not applicable | Not applicable | V01 coverslips | reported | S01, Coverslip functionalization |
| Cell culture endpoint | NIH/3T3 cells | Approximately 80–90% confluence | Not restated | Not stated | Culture medium not restated | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Perfusion before fixation | Heparin; NaNO₂ | 10 U/ml; 0.5% (w/v) | 4°C | 8 min | 0.1 M PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Perfusion fixation | PFA | 4% | 4°C | 8 min | 0.1 M PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Brain post-fixation | PFA | 4% | Room temperature | 2 h | Not restated | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Cryoprotection | Sucrose | 30% | 4°C | Until brain sinks | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Frozen storage / sectioning | OCT-embedded brain; cryosections | 5-µm sections | −80°C storage | Not stated | OCT embedding medium | Not applicable | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Primary anchoring reagent formulation | LNA polyT15; LNA polyT30 | 2 µM each | Not separately stated | Not separately stated | 40% formamide; 2× SSC; 10% dextran sulfate | Not stated | Not applicable | Not applicable | Not applicable | V01 primary-probe formulation | reported | S01, Primary-probe construction |
| Fiducial bead incubation | TetraSpeck beads | 1:1,000 | Room temperature | 5–10 min | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Pre-gel incubation | Degassed acrylamide | 4% | Room temperature | 5 min | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Post-digestion stabilization | Label-IT | 1:10 | 37°C | 30 min | Not stated | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Post-stabilization fixation | PFA | 4% | Room temperature | 15 min | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Quenching | Tris-HCl | 1 M | Room temperature | 10 min | Tris-HCl | 8 | Not applicable | Not applicable | Not applicable | V01 cells | reported | S01, Cell culture experiment |
| Pre-gel functionalization | Acryoloyl-X succinimidyl ester | 0.1 mg/ml | Room temperature | 30 min | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Pre-gel incubation | Acrylamide | 4% | Room temperature | 30 min | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Post-digestion stabilization | label-X | 0.1 mg/ml | 37°C | 45 min | Stock preparation referenced, not restated | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Re-embedding | Acrylamide; APS; TEMED | 4%; 0.05%; 0.05% (previous gel step referenced) | Previous gel step referenced: 4°C then 37°C | Previous step referenced; gelation may be shortened to 2 h | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | V01 tissue | reported | S01, Tissue slices experiment |
| Prepared sample storage | SUPERase In | 2 U/µl | 4°C | No longer than 6 h, or image immediately | 4× SSC | Not stated | Not applicable | Not applicable | Not applicable | V01 cells/tissue | reported | S01, Cell culture experiment; Tissue slices experiment |
| Readout dye conjugation | 5′-amine DNA; Alexa Fluor 647/Cy3B/Alexa Fluor 488 dyes | 5 nmol DNA; 25 µg dye | 37°C | At least 6 h | 0.5 M sodium bicarbonate; 10% dimethyl fumarate (published wording) | Not stated | Not applicable | Not applicable | Not applicable | V01 readout manufacture | reported | S01, Readout-probe synthesis |
| Readout working stock | Dye-conjugated readouts | 500 nM | −20°C storage | Not stated | Resuspended in water | Not stated | Not applicable | Not applicable | Not applicable | V01 readout manufacture | reported | S01, Readout-probe synthesis |
| Primary stock storage | Primary probes | Not stated | −20°C | Until use | Primary hybridization buffer | Not stated | Not applicable | Not applicable | Not applicable | V01 probe manufacture | reported | S01, Primary-probe construction |
| Readout wash rinse | SUPERase In | 1:1,000 | Not stated | Not stated | 4× SSC; approximately 200 µl | Not stated | Not applicable | After each readout wash | Not applicable | V01 imaging | reported | S01, seqFISH+ imaging |
| DAPI during serial imaging | DAPI; SUPERase In | 10 µg/ml; 1:1,000 | Not separately stated | Approximately 15 s | 4× SSC | Not stated | Not applicable | Each serial hybridization | Not applicable | V01 imaging | reported | S01, seqFISH+ imaging |
| Anti-bleaching imaging buffer | Glucose; catalase; glucose oxidase; SUPERase In; Tris-HCl | 10% (w/v); 1:100; 0.5 mg/ml; 0.02 U/µl; 50 mM | Not separately stated | Not stated | 4× SSC | 8 (Tris-HCl) | Not applicable | Each serial hybridization | Not applicable | V01 imaging | reported | S01, seqFISH+ imaging |
| Optical acquisition | 647/561/488/405-nm channels | Two z-slices; 0.35-µm steps | Not applicable | Not stated | Imaging buffer above | 8 (Tris-HCl) | Not applicable | 1–81; 81 repeats 1 | Not applicable | V01 imaging | reported | S01, seqFISH+ imaging |
| Final segmentation staining | DAPI; LNA T20–Alexa Fluor 647; Nissl | 10 µg/ml; 50 nM; 1:100 | Room temperature | 30 min | PBS | Not stated | Not applicable | After serial readout | Not applicable | V01 final imaging | reported | S01, seqFISH+ imaging |
| smFISH comparison hybridization | Primary probes | 10 nM/probe | 37°C | Overnight | 40% hybridization buffer | Not stated | Not applicable | Not applicable | Not applicable | V01 comparison assay, cells | reported | S01, smFISH |
| smFISH comparison hybridization | Primary probes | 10 nM/probe | 37°C | >16 h | 40% hybridization buffer | Not stated | Not applicable | Not applicable | Not applicable | V01 comparison assay, tissue | reported | S01, smFISH |

## Version-Specific Differences

V01 sample branches differ in fixation, clearing, initiator concentrations and stabilization as above. The S05 mirror permits 5–20-µm cryosections while S01's experiment used 5 µm. S05 gives typical digestion durations of 1 h/3 h for cells/tissue, whereas S01 specifies >1 h/>3 h. S05 gives a general 15–20-min readout incubation; S01 specifies 17/20 min. Preserve source scope rather than merge these into a single exact protocol. [S01, cell/tissue/imaging methods; S05 mirror, Tissue section preparation, Digestion/Clearing, Readouts hybridization; reported (S01); inferred (mirror-only guidance, official-body equivalence unresolved)]

The S01 smFISH comparison uses two readout sites and 10 nM primary per probe, overnight at 37°C for cells or >16 h for tissue; these are comparison-assay parameters, not replacements for the seqFISH+ panel. [S01, Methods, smFISH; reported; V01]

## Unreported Parameters

Exact probe/primer sequences, PCR cycle count and source-specific manufacture reaction recipes are unresolved in the examined original paper/protocol/repository resources. Exact reference release and design-software execution are also unresolved. Readout wash durations and temperature are not independently stated in the S01 imaging paragraph. Printed unusual units/solvent names above are preserved and require clarification before operational reuse. These gaps are not filled with routine laboratory defaults.
