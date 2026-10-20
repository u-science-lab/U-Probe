# Experimental Protocol

> Original experiment conditions (S02, V01-V03) and later generalized targeted recommendations (S11-S14, V04) are separate scopes.

## Protocol Scope

Original methods are in S02: V01 untargeted library preparation/SOLiD (supplementary pp. 5-8), V02 targeted brain experiments/four-round SBS (pp. 18-23), and V03 targeted metastatic breast cancer/SOLiD (pp. 29-31). Later generalized targeted recommendations are V04: S11 Tissue Preparation v2, S12 Probe Generation v2, S13 Sequencing Library Preparation v1, and S14 In Situ Sequencing (Illumina Chemistry) v1, routed by S03. These are meaningful chemistry/encoding scopes, not author-numbered releases.

Evidence status is reported for the following explicitly cited reagent/parameter steps unless identified as derived or inferred. Unless explicitly labeled otherwise, the detailed targeted steps and parameter rows below describe **V04 recommendations**, not an exact reconstruction of every original experiment. S11 and S12 abstracts explicitly call their protocols generalizations. The examined protocols.io collection has a dedicated Illumina module; original detailed SOLiD experimental protocols are in S02, not absent.

## Sample Preparation

- **Mouse brain (Option A)**: Transcardial perfusion with ice-cold 1× PBS (5–10 mL) followed by ice-cold 4% PFA in 1× PBS (10–20 mL); post-fixation in 4% PFA / 1× PBS at 4 °C overnight; vibratome sectioning at 50 µm thickness into 70% ethanol for permeabilization (≥1 h at RT, or store at 4 °C for up to months); rehydrate in 1× PBS 15 min ×2 at RT before use. [S11 Steps 5–7]
- **Mouse brain (Option B)**: PFA perfusion as in Option A, followed by 30% sucrose cryoprotection at 4 °C until sunk, freezing, cryosectioning and 70% ethanol permeabilization; the post-sectioning formalin step is skipped for this PFA-perfused branch. [S11 Steps 5.5, 8–10]
- **Fresh frozen tissue (Option C; generalized recommendation)**: Snap-freeze in OCT / dry ice / isopentane, store at −80 °C, cryosection onto SuperFrost Plus slides, fix with ice-cold 10% buffered formalin for 15 min at RT for 15-20 µm sections (longer for thicker sections; in fume hood), rinse in 1× PBS 3× 5 min at RT, permeabilize in 70% ethanol (1 h at RT or 4 °C overnight). [S11 Steps 8–10]

## Fixation

- **Mouse brain (transcardial PFA)**: 4% PFA in 1× PBS at 4 °C overnight (typically 20–30 mL per brain). [S11 Step 5.3]
- **Fresh frozen tissue (post-sectioning)**: ice-cold 10% buffered formalin for 15 min at RT for 15-20 µm sections (longer for thicker sections; in fume hood); 4–5 mL per slide in 4-well plates. [S11 Step 9.2]
- **Quenching of free aldehydes** (optional, after fixation): 100 mM glycine in 1× PBS at 4 °C overnight; then wash in 1× PBS 30 min at 4 °C. [S11 Step 5.4]

## Permeabilization

- 70% ethanol, ≥1 h at RT or 4 °C overnight; storage in 70% ethanol at 4 °C for up to months. [S11 Steps 6, 9]

## Reverse Transcription

**V01 only, reported (S02 supplementary pp. 5-6):** following expansion/re-embedding/passivation, DNase I treatment removes endogenous DNA. RT uses modified random-octamer primers with sample-specific constant regions (see probes.md), SuperScript IV, RNase inhibitor, DTT, dNTP and inosine. The culture and 50-µm hippocampal slice include 40 µM aminoallyl-dUTP in RT, followed by formaldehyde fixation; the 15-µm slice did not have cDNA anchoring. RNase/RNase H and Endonuclease V processing precede cDNA circularization with CircLigase II (3 U/µL, 1 M betaine, 2.5 mM MnCl₂, 1× buffer; 3 h at 60 °C). RCA primers differ by sample (S02 p. 6). Original wording `250 µM dNTP` is not expanded into an unsupported per-nucleotide interpretation here.

**V02-V04 targeted:** reverse transcription is not applicable because padlocks bind and are ligated directly on RNA [S02, supplementary pp. 20-22, 30; S13 Step 9; reported].

## Probe Hybridization

- **Pre-hybridization**: wash gels with wash buffer (2× SSC, 20% formamide in water) for 30 min at RT. [S13 Step 8.1]
- **Hybridization mix**: 2× SSC, 20% formamide (v/v), pooled padlock probes at variable concentration (suggested 100 nM per probe); total volume 200 µL. [S13 Step 8.2]
- **Incubation**: overnight at 37 °C. [S13 Step 8.3]
- **Washes**: wash buffer (2× SSC / 20% formamide) for 30 min ×3 at 37 °C; then 1× PBS for 30 min at 37 °C. [S13 Steps 8.4–8.5]
- **Pooled probe concentration**: pooled library up to 200 µM total probe; 100 nM per probe as a starting point; for equal numbers of probes per gene, per-probe stock concentration = total stock concentration / (number of genes × probes per gene), a derived dilution relation. This stock concentration is distinct from the working hybridization concentration. [S13 Step 4 + Step 8.2]

## Washing

- Wash buffer: 2× SSC, 20% formamide (v/v) in water. [S13 Step 1.3]
- Wash-10 buffer: 2× SSC, 10% formamide (v/v) in water (used for amplicon detection probe hybridization). [S13 Step 1.4]
- Stripping solution: 80% formamide (v/v), 0.01% Triton X-100 (v/v) in water (used at 80 °C for detection probe stripping). [S13 Step 1.5]
- All wash volumes are 200 µL for PCR tubes; 500–1000 µL for microcentrifuge tubes. [S13 Guidelines]

## Ligation

- **Buffer pre-equilibration**: wash gels with 1× SplintR Ligase Buffer for 30 min at RT. [S13 Step 9.1]
- **Pre-ligation**: SplintR Ligase (1× buffer + 1.25 U/µL enzyme in 200 µL) for 6 h at 4 °C (extended pre-ligation allows the enzyme to diffuse into the sample without ligation activity; the protocol notes that 3 h to overnight can also work). [S13 Step 9.2]
- **Ligation**: fresh ligation mix (1× SplintR buffer + 1.25 U/µL SplintR ligase) for 6 h at 37 °C (can be extended overnight; 6 h is recommended to minimize RNA degradation). [S13 Step 9.3]
- **Ligation termination**: wash with 2× SSC for 30 min at RT. [S13 Step 9.4]
- SplintR Ligase: NEB M0375L, 25 U/µL stock, 1.25 U/µL final. [S13 Materials + Step 9.2]

## Amplification

**Rolling-circle amplification (RCA) by Phi29:**

- **Buffer pre-equilibration**: 1× Phi29 Buffer wash for 15 min at RT. [S13 Step 11.1]
- **Pre-amplification**: Phi29 DNA polymerase (1 U/µL) + 1× Phi29 Buffer in 200 µL, 6 h at 4 °C (extended pre-amplification enables uniform enzyme diffusion; can be shortened to 3 h or extended to overnight). [S13 Step 11.2]
- **RCA mix**: 1× Phi29 Buffer, 250 µM each dNTP, 40 µM 5-(3-aminoallyl)-dUTP (aadUTP), 1 U/µL Phi29 in 200 µL total. [S13 Step 11.3]
- **Incubation**: overnight at 30 °C. [S13 Step 11.3]
- **Termination**: 1× PBS for 30 min at RT. [S13 Step 11.4]
- **BS(PEG)9 cross-linking**: 5 mM BS(PEG)9 in 1× PBS, 2 h at RT; quench with 1 M Tris pH 8 for 15 min; wash in 1× PBS for 30 min. [S13 Step 12]

Phi29 DNA polymerase: Enzymatics P7020-HC-L, 100 U/µL stock, 1 U/µL final. [S13 Materials + Step 11.2/11.3]

## Polymerization

### Primary expansion gel, V04

S11 Step 3.1 (PDF pp. 9-10) distinguishes precursor stocks (38% w/v sodium acrylate, 50% w/v acrylamide, 2% w/v bis) from **StockX**, which combines them with NaCl, PBS and water. The source table reports:

| Ingredient | Stock | Volume | Source-listed StockX concentration |
|---|---|---|---|
| Sodium acrylate | 38 g/100 mL | 2.25 mL | 8.6 g/100 mL |
| Acrylamide | 50 g/100 mL | 0.5 mL | 2.5 g/100 mL |
| Bis | 2 g/100 mL | 0.75 mL | 0.15 g/100 mL |
| NaCl | 5 M | 4 mL | 11.7 g/100 mL |
| PBS | 10× | 1 mL | 1× |
| Water | Not applicable | 0.9 mL | Not applicable |
| Total | | 9.4 mL | |

**Reported source arithmetic ambiguity:** the listed volumes sum to 9.4 mL, while the listed final concentrations use a different apparent denominator. Both are preserved; this record does not silently substitute a corrected StockX recipe. MES pH 6.5 is a passivation-buffer stock, not the expansion-gel buffer.

For primary gelation, combine StockX:0.5% 4HT:10% TEMED:10% APS at 47:1:1:1, adding 4HT, TEMED, then APS; the 200-µL example is 188/4/4/4 µL (final 4HT 0.01%, TEMED/APS each 0.2%). Pre-incubate 30 min at 4 °C and gel 2 h at 37 °C [S11 Steps 14.1-14.5; reported].

### First re-embedding and passivation, V04

First re-embedding solution per 5 mL: 375 µL 40% acrylamide/bis 19:1, 4525 µL water, 25 µL 1 M Tris base, 37.5 µL each 10% TEMED and APS. Source-listed final concentrations are approximately 2.85% acrylamide/0.15% bis, 5 mM Tris, 0.075% each TEMED/APS. Pre-incubate 30 min at RT, nitrogen purge 10 min, then gel 1.5 h at 37 °C; the source reports a locked 3.3× linear expansion factor [S11 Steps 17.1-17.4; reported].

Passivation activation: equal volumes of 300 mM NHS in 4 M ethanolamine-HCl and 300 mM EDC in 200 mM MES pH 6.5; final 150 mM NHS, 150 mM EDC, 2 M ethanolamine, 100 mM MES, 2 h at RT. Coupling: equal volumes 4 M ethanolamine-HCl and 125 mM borate pH 8.5; final 2 M ethanolamine and 62.5 mM borate, 40 min at RT. Wash 1× PBS 20 min ×3 at RT [S11 Steps 18.1-18.4, PDF pp. 32-34; reported]. Buffer pH values are those reported by the source, not independent measurements of the final mixtures.

### Second re-embedding, V04

The distinct glass-anchoring gel uses approximately 3.8% acrylamide/0.2% bis with beads, as specified under Sequencing / Readout [S14 Step 7; reported]. Do not substitute this for the first re-embedding formulation.

## Sequencing / Readout

### V04 generalized Illumina SBS

The detailed steps below are the later generalized Illumina module [S14]. Original V02 brain SBS acquired four rounds and used a 60 °C cleavage incubation (S02 supplementary pp. 22-23); the later module uses 50 °C. Preserve these separate source/experiment conditions.

**Sample preparation for sequencing:**

- **Glass-bottom plate preparation**: Bind-Silane (10 µL) + 100% ethanol (4000 µL) + glacial acetic acid (100 µL) + water (890 µL) reaction mix, 45 min at RT, then wash 2× in 100% ethanol and dry. [S14 Step 6]
- **Second re-embedding gelling solution** (per 1 mL): Acrylamide/Bis 19:1 40% (w/v), 100 µL; water, 875 µL; 1 M Tris base, 5 µL (5 mM final); 10% TEMED, 5 µL (0.05%); 10% APS, 5 µL (0.05%); 0.2 µm TetraSpeck beads, 10 µL. [S14 Step 7]
- **Nitrogen purge**: 10 min; gelation at 37 °C for 1.5 h. [S14 Step 7.2]
- **Post-gelation wash**: 1× PBS for 20 min ×3 at RT. Plate can be stored at 4 °C for up to two weeks. [S14 Step 7.3]

**Dideoxy capping of cellular 3′ DNA ends** (to block background from cellular DNA during SBS):

- **Pre-incorporation**: TdT buffer (1×) + 250 µM CoCl₂ + 50 µM each ddNTP in 200 µL, 20 min at RT. [S14 Step 9.1]
- **TdT reaction**: TdT buffer (1×) + 250 µM CoCl₂ + 50 µM each ddNTP + 0.4 U/µL TdT in 200 µL, 1.5 h at 37 °C. [S14 Step 9.2]
- TdT enzyme: NEB M0315S, 20 U/µL stock, 0.4 U/µL final. [S14 Materials + Step 9.2]

**Sequencing primer hybridization:**

- **Primer**: TCTCGGGAACGCTGAAGACGGC (22 nt), HPLC-purified, ≥100 nmol synthesis scale. [S14 Materials / Step 10]
- **Hybridization mix**: 100 µM primer stock diluted 1:40 in 4× SSC, giving 2.5 µM primer. [S14 Step 10.1]
- **Incubation**: 1 h at 37 °C. [S14 Step 10.1]
- **Washes**: 4× SSC for 10 min ×3 at 37 °C, then PR2 for 10 min ×3 at RT. [S14 Steps 10.2–10.3]

**Sequencing cycles** (repeated for each base of the barcode; 7 rounds for a 7-base barcode):

- **Incorporation pre-wash + reaction**: ISS Incorporation Mix (0.5× Illumina IMT + 2.5 mM MgCl₂ + 1× Taq Reaction Buffer), 15 min wash ×2 then 10 min at 50 °C. [S14 Steps 12.1–12.2]
- **Post-incorporation washes**: 2% Zwittergent (w/v) in PR2 for 15 min ×2 at 50 °C, then PR2 for 15 min at 50 °C. [S14 Step 12.3]
- **DAPI staining**: 1 mg/L DAPI in PR2 for 15 min at RT. [S14 Step 12.4]
- **Imaging buffer wash**: USM (Illumina scan mix) for 15 min ×2 at RT; image in USM. [S14 materials and Step 13.2]
- **Imaging**: see Imaging section below. [S14 Step 13.2]
- **Cleavage**: CMS (Illumina) for 10 min wash ×2 at RT, then 20 min at 50 °C; post-cleavage wash in PR2 (10 min ×2 at 50 °C; 10 min ×2 at RT). [S14 Step 14]
- **Iteration**: if more bases remain to be sequenced, return to incorporation step. [S14 Step 15]

**Sequencing primer / RCA primer note:** the 22-nt sequencing primer is identical to the default padlock backbone and binds the complementary RCA product. The separate 18-nt RCA primer binds backbone bases 1-18 on the circle; reverse complementation of those bases gives its letters (derived). See probes.md for sequence notation [S13 Step 5; S14 materials/Step 10; S04 assembly lines 83-91; V04].

### V01 original untargeted and V03 original cancer SOLiD

S02 supplementary pp. 7-8 describes automated untargeted ligation sequencing: five offset primers, stripping, primer hybridization, T4 DNA ligation, washes, imaging and cleavage. The V03 targeted protocol is separately specified on pp. 30-31: 2.5 µM sequencing primer in 5× SASC for 20 min; 6 U/µL T4 DNA ligase and 1:40 SOLiD oligos for 1 h; extended four-wash/3 h total manual washes, imaging, Quick CIP treatment and two-step cleavage. All sample reactions are reported at RT, with specified reservoirs held at 80 °C or 4 °C. V01 and V03 wash/cycle schedules differ and must not be combined. The original source contains execution details beyond S12's barcode-construction explanation [S02 named sections; reported].


## Imaging

- **Microscope**: spinning-disk confocal, with 405, 488, 561, 640 and 685 nm laser lines. [S14 materials and Step 13.2]
- **Recommended objective**: Nikon 40× water-immersion long-working-distance, NA 1.15, CFI ApoLWD Lambda S 40XC WI (MRD77410). [S14 Step 13.2; S13 Step 14.6]
- **Imaging parameters** (per the Alon 2021 protocols.io protocol):

| Channel | Excitation (nm) | Laser power (mW) | Emission range (nm) | Exposure (ms) | Template base | Added base |
|---|---|---|---|---|---|---|
| A | 685 | 30 | 705–845 | 400 | G | C |
| B | 640 | 160 | 663–737 | 200 | T | A |
| C | 561 | 150 | 575–590 | 200 | A | T |
| D | 488 | 150 | 500–550 | 400 | C | G |
| DAPI | 405 | 100 | 440–460 | 300 | n/a | n/a |

(All lasers at 100% output power; tile overlap 2%; z-spacing 0.4 µm; image longer-to-shorter wavelength.) [S14 Step 13.2]

## Decoding

V02 original brain: retain perfect matches across four logical bases [S02, supplementary p. 23; reported]. V03 cancer: allow one mismatch in seven-symbol barcodes [S02 p. 31; reported]. V04: use the selected full-length or embedded codebook and experiment's gene assignment; the supplied library is not an original-brain default [S12 Steps 3-7; reported]. ExSeqProcessing documents six analysis stages [S10, Stage Overview; reported].

V01 untargeted: match short in situ reads to the ex situ Illumina library from extracted amplicons. Uniqueness is 92% of matches overall and 97% of non-rRNA matches, not an overall read-matching yield [S01, matching paragraph; reported]. Ex situ annotation used GRCm38/mm10 and UCSC RefSeq downloaded 2018-08-14; this is an analysis reference, not an established targeted-design database [S02 supplementary p. 11; reported].

## Critical Parameters

| Step / condition | Reagent | Concentration | Temperature | Duration | Buffer | pH | Cycle number | Imaging round | Sequencing cycle | Version / sample scope | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Tissue fixation (mouse brain) | 4% PFA in 1× PBS | 4% (w/v) | 4 °C | Overnight | 1× PBS | 7.4 | n/a | n/a | n/a | V04 (mouse brain) | reported | S11 Step 5.3 |
| Tissue fixation (fresh frozen) | 10% buffered formalin | 10% | RT | 15 min for 15-20 µm sections; longer for thicker sections | n/a | n/a | n/a | n/a | n/a | V04 fresh-frozen tissue recommendation | reported | S11 Step 9.2 |
| Permeabilization | 70% ethanol | 70% (v/v) | 4 °C / RT | ≥1 h (or up to months at 4 °C) | water | n/a | n/a | n/a | n/a | V04 | reported | S11 Steps 6, 9 |
| LabelX anchoring | LabelX | 0.05–0.1 mg/mL (1:20 to 1:10 dilution of 1 mg/mL stock) | 37 °C | Overnight (≥6 h) | 1× Buffer A or 20 mM MOPS | 7.7 | n/a | n/a | n/a | V04 (mouse brain & human tumor tested at 0.1 mg/mL; 0.05 mg/mL works in thin tumor sections) | reported | S11 Step 12 |
| Expansion gel StockX | Mixed monomers, NaCl and PBS; see Polymerization | Source table lists 8.6% sodium acrylate / 2.5% acrylamide / 0.15% bis; arithmetic ambiguity retained | 4 °C pre-incubation; 37 °C gelation | 30 min; 2 h | PBS / NaCl | Not reported as final measured pH | Not applicable | Not applicable | Not applicable | V04 | reported | S11 Steps 3.1, 14 |
| Gelation accelerator + initiator | TEMED; APS | 0.2% each (10% stock at 1:50) | added immediately before gelation | — | — | — | n/a | n/a | n/a | V04 | reported | S11 Step 14 |
| Mild Proteinase K | Proteinase K in 1× digestion buffer + 500 mM NaCl | 8 U/mL | RT, 37 °C or 60 °C (alternatives) | Overnight at RT, 3 h at 37 °C, or 1 h at 60 °C | 1× digestion buffer | n/a | n/a | n/a | n/a | V04 (mouse brain) | reported | S11 Steps 15.1, 15.3 |
| Robust Proteinase K | Proteinase K in 1× digestion buffer + 800 mM guanidine HCl | 8 U/mL | 60 °C then 37 °C | 1 h at 60 °C, then overnight at 37 °C; replace mix at elapsed 24 h, continue overnight; 48 h total | 1× digestion buffer | n/a | n/a | n/a | n/a | V04 robust-digestion recommendation | reported | S11 Steps 15.1, 15.4 |
| Expansion | water | — | RT | 3× 20 min washes | water | Not reported | n/a | n/a | n/a | V04 | reported | S11 Step 16 |
| Re-embedding | 19:1 acrylamide:bis-acrylamide (40%) | ~2.85% acrylamide / 0.15% bis | 37 °C | 1.5 h | 5 mM Tris base | Not reported | n/a | n/a | n/a | V04 (locks the linear expansion factor at ~3.3×) | reported | S11 Step 17 |
| Passivation activation | NHS; EDC; ethanolamine-HCl | 150 mM NHS; 150 mM EDC; 2 M ethanolamine | RT | 2 h | 100 mM MES | 6.5 (source buffer) | Not applicable | Not applicable | Not applicable | V04 | reported | S11 Steps 18.1-18.2 |
| Passivation coupling | Ethanolamine-HCl | 2 M | RT | 40 min | 62.5 mM sodium borate | 8.5 (source buffer) | Not applicable | Not applicable | Not applicable | V04 | reported | S11 Step 18.3 |
| Padlock hybridization | pooled padlock probes | 100 nM per probe working concentration (suggested; distinct from 200 µM total stock) | 37 °C | Overnight | 2× SSC / 20% formamide | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 8 |
| SplintR ligation pre-diffusion | SplintR ligase | 1.25 U/µL | 4 °C | 6 h (or overnight) | 1× SplintR buffer | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 9.2 |
| SplintR ligation | SplintR ligase | 1.25 U/µL | 37 °C | 6 h (or overnight) | 1× SplintR buffer | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 9.3 |
| RCA primer hybridization | universal RCA primer (TCT TCA GCG TTC CCG A*G*A, *PS*) | 500 nM | 37 °C | 2 h | 2× SSC / 20% formamide | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 10 |
| RCA pre-amplification | Phi29 DNA polymerase | 1 U/µL | 4 °C | 6 h (or overnight) | 1× Phi29 buffer | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 11.2 |
| RCA | Phi29 DNA polymerase + dNTPs + aadUTP | 1 U/µL Phi29; 250 µM each dNTP; 40 µM aadUTP | 30 °C | Overnight | 1× Phi29 buffer | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 11.3 |
| BS(PEG)9 cross-linking | BS(PEG)9 | 5 mM | RT | 2 h | 1× PBS | 7.4 | n/a | n/a | n/a | V04 | reported | S13 Step 12 |
| Detection probe hybridization | /5Cy3/TCTCGGGAACGCTGAAGA | 100 nM | 37 °C | 1 h | 2× SSC / 10% formamide (Wash-10) | n/a | n/a | n/a | n/a | V04 (optional QC) | reported | S13 Step 14 |
| Detection probe stripping | 80% formamide / 0.01% Triton X-100 | 80% (v/v) formamide | 80 °C (pre-heated; sample at RT) | 3× 10 min | water | n/a | n/a | n/a | n/a | V04 | reported | S13 Step 14.7 |
| Second re-embedding | 19:1 acrylamide:bis 40% | ~3.8% acrylamide / 0.2% bis | 37 °C | 1.5 h (after 10 min N₂ purge) | 5 mM Tris | Not reported | n/a | n/a | n/a | V04 | reported | S14 Step 7 |
| TdT ddNTP capping | TdT + ddNTPs | 0.4 U/µL TdT; 50 µM each ddNTP | 37 °C | 1.5 h | 1× TdT buffer + 250 µM CoCl₂ | n/a | n/a | n/a | n/a | V04 (blocks cellular 3′ ends) | reported | S14 Step 9 |
| Sequencing primer hybridization | TCTCGGGAACGCTGAAGACGGC | 2.5 µM (1:40 of 100 µM stock in 4× SSC) | 37 °C | 1 h | 4× SSC | n/a | n/a | n/a | n/a | V04 | reported | S14 Step 10 |
| SBS incorporation | 0.5× IMT + 2.5 mM MgCl₂ | 0.5× IMT | 50 °C | 10 min | 1× Taq Reaction Buffer | n/a | n/a | n/a | n/a | V04 | reported | S14 Step 12.2 |
| SBS imaging | spinning-disk confocal; 4 colors (685/640/561/488 nm) | — | RT | tile + z-stack | USM (Illumina) | n/a | n/a | 1 of 7 | 1 of 7 (one base per cycle) | V04 (template base G→T→A→C; added base C→A→T→G; 7 cycles for a 7-base barcode) | reported | S14 Step 13.2 |
| SBS cleavage | CMS (Illumina) | As supplied | 50 °C | 20 min | CMS; PR2 is the subsequent wash | n/a | n/a | n/a | n/a | V04 (cleaves dye and unblocks reversible terminator) | reported | S14 Step 14.3 |
| DAPI staining (in sequencing) | DAPI | 1 mg/L | RT | 15 min | PR2 | n/a | n/a | n/a | n/a | V04 | reported | S14 Step 12.4 |
| Antibody signal recovery | Streptavidin-Alexa 488 (or alternative) | 10 µg/mL | 4 °C | Overnight | PR2 | n/a | n/a | n/a | n/a | V04 (optional, post-sequencing) | reported | S14 Step 16 |
| Untargeted ExSeq in situ sequencing | SOLiD ligation; 17 rounds shown | n/a | n/a | n/a | n/a | n/a | 17 (rounds) | 17 (rounds) | 17 (cycles) | V01 (Fig. 1Biii example) | reported | S01 Fig. 1Biii |

## Version-Specific Differences

- V01: modified random-octamer RT, sample-dependent cDNA fixation, enzymatic processing, cDNA circularization/RCA and SOLiD sequencing [S02 supplementary pp. 5-8; reported].
- V02: original brain four-round checksum encoding; cortex maximum 8 probes/transcript, Tm threshold initially 55 °C lowered to minimum 52 °C when needed; hippocampus maximum 16, threshold 50 °C. Original search starts the next successful window 5 nt past the previous end. The original Tm inequality wording is ambiguous and preserved in probes.md [S02 pp. 19-20, 23; reported].
- V03: cancer length-7 distance-3 codes, 9-16 probes/gene, dedicated original SOLiD protocol; do not replace with the V04 embedded library or generalized SBS recipes [S02 pp. 29-31; reported].
- V04: generalized targeted execution and configurable design defaults, rather than an exact retrospective record of every panel [S11/S12 abstracts; reported].

V01 expression correlation r=0.89 uses full ex situ ExSeq data versus an adjacent RNAseq slice; ten acquired fields give r=0.47 and detect 3,039 genes [S01 Fig. 2Giv-vii; reported].

## Unreported Parameters

The original S9/S10 workbook cells are not yet examined; gene/accession/barcode and probe-sequence resources are explicitly public for brain and cancer panels [S02 supplementary pp. 20, 29, 90]. Historical targeted off-target database releases remain unresolved in the examined design sources. The V01 analysis reference is separately documented above.

The StockX volume/concentration inconsistency in S11 Step 3.1 and original Tm wording in S02 p. 19 remain source ambiguities, not invented corrected values. Exact original reagent lot numbers and complete historical execution environments have not been reconstructed. Current software maintenance remains unresolved. Detailed passivation and original SOLiD protocols are available and must not be listed as unreported.
