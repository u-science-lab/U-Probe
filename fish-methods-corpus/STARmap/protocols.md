# Experimental Protocol

## Protocol Scope

V01 is original STARmap (S01, Wang 2018); V02 is STARmap PLUS (S11, Zeng 2023). S03 is a later developer protocol. Its probe-design section explicitly includes STARmap and STARmap PLUS, but its numerical procedure is recorded here as **S03 / V02 protocol conditions**, not as evidence of the conditions used in 2018. RIBOmap/TEMPOmap-specific splint and metabolic-labeling steps are excluded. All S03 conditions below are reported unless a calculation is explicitly labeled derived. This is a sourced summary; consult the complete source procedure before execution.

## Sample Preparation

V01: S01 reports imaging 8-µm-thick volumes for the 160-gene experiment and 150-µm-thick tissue slices for volumetric STARmap (Results, spatial cell typing and scaling to 3D). Imaging volume does not by itself establish cryostat section thickness or a fixation recipe.

S03 / V02: Step 9B describes perfused mouse tissue, OCT freezing and **20-µm cryosections**, mounted one per well. This later thickness must not replace the V01 experimental scope.

## Fixation

S03 / V02, Step 9B: perfuse with 20 mL PBS then 20 mL 4% PFA at approximately 1 mL per 5 s. Optional post-fixation is 4% PFA for 3–4 h at 4°C, followed by 30% sucrose–PBS overnight at 4°C; the source warns that post-fixation may reduce detection sensitivity. Fix mounted tissue slices with 600 µL 4% PFA in 1× PBS per well, room temperature, 15 min. Cultured adherent cells instead receive 400 µL 1.6% PFA in 1× PBS, room temperature, 15 min (Step 9A). These sample-specific options are not combined.

## Permeabilization

S03 / V02, Step 10A: add 600 µL prechilled methanol and incubate at −20°C or −80°C for 1 h; equilibrate at room temperature for 5 min, remove methanol, and quench with 600 µL buffer A for 5 min at room temperature. Alternative Step 10B preserves endogenous fluorescence: 600 µL quenching buffer B, room temperature, 10 min. Buffer B contains glycine, yeast tRNA, RVC and Triton X-100 in PBSTR; its recipe uses 1 mL 10% Triton X-100 per 10 mL final volume (derived: 1% Triton). Briefly rinse with PBSTR (Step 11). See S03 Reagent setup for the complete buffer recipes.

## Reverse Transcription

Not applicable to the SNAIL recognition/amplification workflow: probes bind RNA directly and bypass reverse transcription. The primer is the ligation template and the RCA primer; the circularized **padlock is the RCA template** (reported: S01 Results, SNAIL; S03 Overview of the protocol).

## Probe Hybridization

S03 / V02, Steps 15–18: heat primer/padlock mixtures at 90°C for 5 min and immediately place on ice. Use 400 µL hybridization mixture per well: 2× SSC, 10% formamide, 20 mM RVC, 1% Tween-20, 0.1 mg/mL tRNA, 0.2 U/µL SUPERase·In, and 100 nM **per probe**. For more than 100 genes use 1–5 nM per probe, optimized for tissue type; Tween-20 may be replaced by the same concentration of Triton X-100 for tissue slices. Incubate at 40°C for 12–36 h with shaking, parafilm seal and a water reservoir. RIBOmap splint reagent in the shared source table is not included in STARmap.

## Washing

S03 / V02, Steps 19–20: wash with 600 µL PBSTR twice, then 600 µL high salt buffer once, 37°C for 20 min each; briefly rinse with 600 µL PBSTR at room temperature. PBST is PBS with 0.1% Tween-20; PBSTR adds RNase inhibitor. High salt buffer is prepared from **2 mL 20× SSC plus 8 mL PBSTR** (Reagent setup; derived SSC concentration: 4×, not 2×). Post-ligation washes: PBSTR twice for 5 min (Step 23). Post-RCA washes: PBST twice for 5 min (Step 27).

## Ligation

S03 / V02, Steps 21–22: 400 µL mixture with 1× T4 DNA ligase buffer, 0.5 mg/mL BSA, 0.2 U/µL SUPERase·In and 0.1 U/µL T4 DNA ligase. For more than 100 genes the source specifies 2.5-fold higher ligase concentration (derived: 0.25 U/µL). Incubate at room temperature for 2 h under parafilm. Thaw ligase buffer at room temperature (source critical step). The primer splints the padlock ends for circularization; RNA is the recognition target, not the ligation splint.

## Amplification

S03 / V02, Steps 24–26: 400 µL mixture containing 1× Phi29 buffer, 250 µM dNTP mixture (source notation retained), 20 µM 5-(3-aminoallyl)-dUTP, 0.5 mg/mL BSA, 0.2 U/µL SUPERase·In and 0.2 U/µL Phi29 DNA polymerase. For more than 100 genes use 2.5-fold higher polymerase concentration (derived: 0.5 U/µL). Preincubate at 4°C for 30 min, then incubate at 30°C for 2 h with shaking in a humidified oven under parafilm. These source conditions do not support an 8–14 h RCA duration.

## Polymerization

V01: S01 Results (HTC) reports acrylic acid NHS ester modification and acrylamide copolymerization; numerical V01 coupling conditions are not reported in the examined main text. S02 supplementary methods have not been inspected here.

S03 / V02, Steps 28–33: optional protein comapping places primary antibody binding after RCA and before MA-NHS treatment. Block with 5 mg/mL BSA in PBS-Triton at room temperature for 30 min; primary-antibody incubation depends on the antibody. Use 400 µL modification mixture with 25 mM MA-NHS and 100 mM sodium bicarbonate (stock pH approximately 8.4), room temperature, 1 h.

Steps 34–40: monomer buffer is 4% acrylamide, 0.2% bis, 2× SSC and 0.2% TEMED (source concentrations). Incubate with monomer buffer at room temperature for 25 min. Mix 100 µL monomer buffer with 2 µL 10% ammonium persulfate (source: approximately 0.2% final). Polymerize under a coated coverslip in a nitrogen-purged desiccator at room temperature for 1 h. The source does not document a second 2% acrylamide / 0.05% bis re-embedding step.

Steps 42–43: where protein detection is not required, digest with 400 µL 0.2 mg/mL proteinase K in 2× SSC and 1% SDS at 37°C for 1 h; wash with 600 µL PBST three times, 5 min each. **Skip both steps for STARmap PLUS protein comapping** or other experiments requiring protein signal preservation.

## Sequencing / Readout

S03 / V02, Steps 44–48: dephosphorylate before SEDAL with 0.25 U/µL Antarctic phosphatase, 0.5 mg/mL BSA and 1× phosphatase buffer, 400 µL, 37°C for 1 h. Wash with PBST three times for 5 min; strip twice at room temperature for 10 min and wash with PBST three times for 5 min.

Step 50: use 400 µL sequencing mixture containing 1× T4 ligase buffer, 0.5 mg/mL BSA, 10 µM reading probe, **5 µM total fluorescent decoding-probe mixture** and 0.1875 U/µL T4 DNA ligase; incubate at room temperature for at least 3 h. The decoding mixture is prepared by equal-volume mixing of 16 stocks (Reagent setup). Step 51: three 600 µL washes in washing/imaging buffer, room temperature, 10 min each. Step 60: strip twice for 10 min, then PBST three times for 10 min. Stripping buffer recipe: 6 mL formamide, 100 µL 10% Triton X-100, 3.9 mL water (Reagent setup). Room temperature is explicit for the initial stripping step (47); Step 60 does not separately restate temperature.

The S03 overview describes sequential **overlapping** dinucleotides: reading position shifts by one base each round, with a final quality-control round for a 1,000-gene experiment. Table 4 gives cycle-specific reading probes; Table 5 gives the dinucleotide-to-dye assignments (see probes.md). These are not three independent two-base chunks.

## Imaging

V01: fluorescent Nissl staining after six sequencing rounds was used to segment cell bodies (S01 Results, 160-gene experiment). This does not establish a universal V02 stain.

S03 / V02: DAPI at 100 ng/mL in 400 µL washing/imaging buffer for 3 h (Step 52; some organs may require overnight staining). DAPI is acquired during the imaging workflow, not only after all sequencing rounds. Steps 55–59 describe Leica TCS SP8 four-color acquisition plus DAPI; typical objectives are 40× NA 1.3 oil or 63× NA 1.4 oil at high amplicon densities; recommended x–y pixels <0.15 µm and z steps <0.35 µm. After all SEDAL rounds, secondary-antibody incubation for STARmap PLUS is at room temperature for at least 12 h, with 1:100 dilution suggested but antibody dependent (Steps 62–63).

## Decoding

V01: S01 Results reports two-base encoding, six rounds for the 160-gene experiment and approximately 1.8% versus 29.4% error in the crowded four-gene comparison. The reported filtered fractions are 20% (four-gene) and 40% (1,020-gene). These are experiment-specific measures, not guarantees that every erroneous read is rejected.

S03 / V02: convert the round-by-round colors into a barcode using the two-base scheme, then match the barcode–gene codebook (Overview; Fig. 4; Steps 64 onward). Sequential rounds interrogate overlapping dinucleotides. Exact experiment codebooks require their own verified resources; the generic Table 5 dye mapping is not a complete gene codebook.

## Critical Parameters

Reported conditions retain source scope. Unknown values are not zeros.

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Tissue perfusion | PFA | 4%; 20 mL after 20 mL PBS | Not stated | Flow approximately 1 mL/5 s | PBS | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 9B |
| Optional post-fixation | PFA | 4% | 4°C | 3–4 h | Not separately stated | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 9B(ii) |
| Mounted tissue fixation | PFA | 4%; 600 µL | Room temperature | 15 min | 1× PBS | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 9B(vi) |
| Regular permeabilization | Methanol | 600 µL | −20°C or −80°C | 1 h | Not applicable | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 10A |
| Fluorescence-preserving option | Quenching buffer B | 600 µL; recipe in Reagent setup | Room temperature | 10 min | PBSTR-based | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 10B; Reagent setup |
| Probe denaturation | Primer and padlock mixtures | Stock mixtures | 90°C then ice | 5 min | Stock buffer | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 15 |
| Hybridization | Primer + padlock | 100 nM per probe; >100 genes: 1–5 nM per probe | 40°C | 12–36 h | 2× SSC; 10% formamide; 20 mM RVC; 1% Tween-20; 0.1 mg/mL tRNA; 0.2 U/µL SUPERase·In | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 16–18 |
| Post-hybridization wash | PBSTR then high salt buffer | 600 µL per wash | 37°C | 2 PBSTR + 1 high salt wash; 20 min each | High salt: 2 mL 20× SSC + 8 mL PBSTR | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 19–20; Reagent setup |
| Padlock ligation | T4 DNA ligase | 0.1 U/µL; >100 genes: 2.5-fold higher | Room temperature | 2 h | 1× ligase buffer; 0.5 mg/mL BSA; 0.2 U/µL SUPERase·In; 400 µL | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 21–22 |
| RCA | Phi29; dNTP mixture; aminoallyl-dUTP | 0.2 U/µL Phi29 (>100 genes: 2.5-fold higher); 250 µM dNTP mixture; 20 µM aminoallyl-dUTP | 4°C then 30°C | 30 min then 2 h | 1× Phi29 buffer; 0.5 mg/mL BSA; 0.2 U/µL SUPERase·In; 400 µL | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 24–26 |
| Amine modification | MA-NHS | 25 mM; 400 µL | Room temperature | 1 h | 100 mM sodium bicarbonate | Stock approximately 8.4 | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 31–32 |
| Monomer incubation | Acrylamide; bis; TEMED | 4%; 0.2%; 0.2% | Room temperature | 25 min | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 34–35 |
| Polymerization | Monomer buffer + APS | 100 µL buffer + 2 µL 10% APS; source approximately 0.2% APS | Room temperature; nitrogen-purged desiccator | 1 h | Monomer buffer | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 37–40 |
| Optional tissue clearing | Proteinase K; SDS | 0.2 mg/mL; 1%; 400 µL | 37°C | 1 h | 2× SSC | Not stated | Not applicable | Not applicable | Not applicable | S03; omit for V02 protein comapping | reported | S03 Step 42 |
| Post-clearing wash | PBST | 600 µL | Not separately stated | 3 × 5 min | PBST | Not stated | Not applicable | Not applicable | Not applicable | S03; omit with Step 42 for V02 protein comapping | reported | S03 Step 43 |
| Dephosphorylation | Antarctic phosphatase | 0.25 U/µL; 400 µL | 37°C | 1 h | 1× phosphatase buffer; 0.5 mg/mL BSA | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Steps 44–45 |
| SEDAL ligation | Reading probe; decoding mixture; T4 ligase | 10 µM; 5 µM total mixture; 0.1875 U/µL | Room temperature | At least 3 h | 1× ligase buffer; 0.5 mg/mL BSA; 400 µL | Not stated | Each SEDAL round | Each SEDAL round | Each SEDAL round | S03 / V02 | reported | S03 Step 50 |
| DAPI stain | DAPI | 100 ng/mL; 400 µL | Not separately stated | 3 h; overnight may be needed | Washing/imaging buffer | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 52 |
| Initial stripping | Formamide-based buffer | 6 mL formamide + 100 µL 10% Triton + 3.9 mL water | Room temperature | 2 × 10 min | Stripping buffer | Not stated | Not applicable | Not applicable | Not applicable | S03 / V02 | reported | S03 Step 47; Reagent setup |
| Between-round stripping | Formamide-based buffer | Same recipe | Not separately restated | 2 × 10 min; then PBST 3 × 10 min | Stripping buffer; then PBST | Not stated | Each SEDAL round | Each SEDAL round | Each SEDAL round | S03 / V02 | reported | S03 Step 60 |
| V01 coupling conditions | Acrylic acid NHS ester | Not reported in examined main text | Not reported in examined main text | Not reported in examined main text | Not reported in examined main text | Not stated | Not applicable | Not applicable | Not applicable | V01 | unavailable | S01 Results, HTC; S02 not inspected |
| V01 segmentation | Fluorescent Nissl stain | Not extracted | Not extracted | After six sequencing rounds | Not extracted | Not stated | Not applicable | Not applicable | Not applicable | V01; 160-gene experiment | reported | S01 Results, 160-gene experiment |

| Application HTC coupling | Acrylic acid NHS ester | 20 mM | Not separately stated for this treatment | 1 h | PBST | Not stated | Not applicable | Not applicable | Not applicable | S11 / V02 (2023 application) | reported | S11 Methods, STARmap PLUS procedure |
| Application RNA sequencing | SEDAL | See S11 procedure | Room temperature for ligation | At least 3 h per ligation | S11 mixture | Not stated | 8 for main panel; 4 for validation | 8; 4 | 8; 4 | S11 / V02; 2,766-gene main / 64-gene validation | reported | S11 Methods, sequencing and imaging |
| Application segmentation | Propidium iodide | Manufacturer instructions in source | Not specified in inspected passage | After sequencing/protein staining | See S11 procedure | Not stated | Not applicable | Post-sequencing | Not applicable | S11 / V02 | reported | S11 Methods, final staining |

## Version-Specific Differences

S01 supports V01 acrylic-NHS HTC, RNA recognition, RCA and SEDAL. S03 explicitly includes both STARmap and STARmap PLUS in its probe-design section, while its numerical procedure is a later protocol implementation. V02 adds protein comapping: primary antibody before MA-NHS, secondary antibody after sequencing, and omission of proteinase K when retaining protein signals (S03 Steps 28–43, 62–63). Neither S03 tissue fixation nor its 2 h RCA is asserted to reconstruct the 2018 experiments. RIBOmap and TEMPOmap remain sibling methods.


### S11 application procedure versus later S03 protocol

Reported: S11 Methods, “Tissue collection and sample preparation for STARmap PLUS” and “STARmap PLUS” procedure. The 2023 application uses fresh-frozen brain cut into 20-µm sections; mounted slices receive 4% PFA in PBS at room temperature for 15 min, then cold-methanol permeabilization. After RCA/primary antibody treatment, it uses **20 mM acrylic acid NHS ester in PBST for 1 h**, not the later S03 MA-NHS recipe. Its gel polymerization is APS/TEMED at room temperature under nitrogen for 1 h. S11 uses **PI** for segmentation, **eight sequencing rounds for 2,766 genes**, and a separate **four-round 64-gene** validation experiment. Later S03 uses MA-NHS and DAPI; these changes are protocol-specific within V02, not features retroactively assigned to every STARmap PLUS experiment. This preserves differing reported reagents instead of silently selecting one.

V01 thick tissue is also a separate readout configuration: S01 Fig. 6 and scaling Results report 28 genes read linearly, four genes per cycle across seven cycles. It is not seven cycles of the five-base barcode decoding scheme.
## Unreported Parameters

V01 acrylic-NHS concentration, temperature and duration are not reported in the examined S01 main text and are not supplied by S03, which uses MA-NHS. Original supplementary procedures remain unexamined, so this is not a claim of absence from all official sources. Original panel reference assemblies, annotation releases and per-gene designs remain unresolved pending direct S02 inspection. Antibody-specific conditions remain sample dependent in S03. Do not replace these gaps with common-practice defaults.
