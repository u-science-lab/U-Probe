# Experimental Protocol

## Protocol Scope

MiP-seq protocol for simultaneous detection of DNA, RNA, proteins, and biomolecules in cultured cells and tissue sections. This protocol covers RNA detection as the primary application, with notes on adaptations for DNA and protein detection. Species demonstrated include mouse, rat, zebrafish, and human cell lines. Reference genome/annotation versions were not explicitly stated in the examined sources.

## Sample Preparation

### Cell Culture

Cell types used:
- HEK-293T cells
- HBMEC (human brain microvascular endothelial cells)
- Bend3 (mouse brain endothelial cells)
- C6 (rat glioma cells)
- PK-15 cells (porcine kidney)

Cell seeding density: 9×10³ cells/well in 96-well plates (PerkinElmer, 6055500).

Culture media:
- DMEM (Gibco, 11965118) with 10% FBS (Gibco, 10093) for HEK-293T
- DMEM with 10% FBS, 2 mM L-glutamine, 1 mM sodium pyruvate, 1× MEM nonessential amino acids, 1× MEM amino acid solution, 1× MEM vitamins, 1× penicillin-streptomycin for HBMEC
- Neurobasal with 1.25% FBS, 2% B-27 supplement, 0.5 mM L-glutamine, 0.5% penicillin-streptomycin for primary neurons

Cell culture conditions: 37°C, 5% CO₂.

### Tissue Preparation

Mouse/rat brains:
1. Transcardial perfusion with 4% PFA in PBS
2. Post-fixation in 4% PFA for 12 hours at 4°C
3. Cryoprotection in 30% sucrose in PBS
4. Embedding in OCT at −80°C
5. Sectioning at 8 µm or 10 µm thickness with Leica cryostat
6. Storage at −80°C

Zebrafish brains:
1. Fixation with 4% PFA in PBS for 12 hours at 4°C
2. Cryoprotection in 30% sucrose in PBS
3. Dehydration series: 25%, 50%, 75%, 100% methanol in PBS (5 min each)
4. Storage in 100% methanol at −20°C

### Pretreatment Before MiP-seq

For cultured cells:
1. Wash twice with PBS
2. Fix with 4% PFA for 10 min at room temperature
3. Wash three times with PBSTR (PBS, 0.1% Tween-20, 0.1 U/µL RNase inhibitor) for 5 min each
4. Methanol pre-cooled at −20°C treatment for 15 min at −80°C
5. Seal with Secure-Seal hybridization chambers

For brain slices:
1. Fix with 4% PFA for 10 min (pre-cooled)
2. Wash three times with DEPC-PBSTR for 5 min each
3. Methanol treatment at −80°C for 15 min
4. Pepsin treatment: 2 mg/mL pepsin in 0.1 M HCl at 37°C for 90 s
5. Wash three times with DEPC-PBSTR

For zebrafish:
1. Rehydration series: 75%, 50%, 25% methanol in PBS (5 min each)
2. Wash three times with DEPC-PBSTR
3. Proteinase K (20 mg/mL) treatment at room temperature for 30 min
4. Wash three times with DEPC-PBSTR

## Fixation

4% PFA (paraformaldehyde) in PBS.

## Permeabilization

Methanol at −20°C or −80°C for 15 min (cultured cells and brain slices).

Pepsin (2 mg/mL in 0.1 M HCl) at 37°C for 90 s (brain slices).

Proteinase K (20 mg/mL) at room temperature for 30 min (zebrafish).

## Reverse Transcription

Not applicable. MiP-seq detects RNA directly without reverse transcription.

## Probe Hybridization

### Padlock Probe Phosphorylation

Padlock probes are phosphorylated at 5' end with T4 Polynucleotide Kinase (Vazyme, N102-01) before hybridization.

### Hybridization Conditions

For RNA detection:
- Padlock probe concentration: 30 nM/probe
- RCA initiator concentration: 30 nM/probe (equimolar to padlock probes)
- Hybridization buffer: 2× SSC, 10% formamide, 20 mM RVC
- Incubation: overnight at 37°C

For DNA detection:
- Padlock probe concentration: 100 nM
- Initiator primer concentration: 100 nM
- Hybridization solution: 10% formamide in 2× SSC
- Incubation: overnight at 37°C

## Washing

After hybridization:
- Cultured cells: Wash twice with DEPC-PBSTR for 20 min each, then once with 4× SSC in DEPC-PBSTR
- Brain slices: Same as cultured cells
- DNA detection: Wash once with 25% formamide in 2× SSC, then three times with PBST at room temperature

## Ligation

### RNA Detection (SplintR Ligase)

- Enzyme: SplintR ligase (NEB, M0375L)
- Concentration: 1 U/µL
- Buffer: 1× SplintR buffer
- RNase inhibitor: 0.2 U/µL
- Incubation: 2 hours at 37°C
- Followed by washing twice with DEPC-PBSTR for 20 min each

### DNA Detection (T4 DNA Ligase)

- Enzyme: T4 DNA ligase (Thermo Fisher, EL0011)
- Concentration: 0.2 U/µL
- Buffer: 1× ligation buffer with 0.2 µg/µL BSA
- Incubation: 3 hours at 16°C
- Followed by three washes with PBST

## Amplification

### Rolling-Circle Amplification (RCA)

- Enzyme: Phi29 polymerase (Vazyme, N106-01)
- Concentration: 1 U/µL
- Buffer: 1× RCA buffer
- dNTPs: 0.25 µM
- BSA: 0.2 µg/µL
- Glyceryl: 5%
- Incubation temperature: 30°C
- Duration: 2 hours (RNA) or 6 hours (DNA)
- Followed by washing twice with DEPC-PBSTR or PBST for 20 min each

## Polymerization

Not applicable.

## Sequencing / Readout

### Sequencing Ligation Reaction

- Enzyme: T4 DNA ligase (Thermo Fisher Scientific, EL0012)
- Concentration: 0.2 U/µL
- Buffer: 1× T4 DNA ligase buffer with 1× BSA
- Anchor primer concentration: 1 µM
- Query primer concentration: 1 µM
- Incubation: 2 hours at 25°C

### Stripping Between Rounds

- Buffer: 60% formamide in 2× SSC
- Duration: 10 min at room temperature
-重复: Two times

### Sequential Dilution Strategy (Optional)

- Different sequencing primer sets for different gene subsets
- Reference genes (e.g., Gad2, Gfap) sequenced in every round
- Signals reconstructed by alignment of reference gene signals

## Imaging

### Microscope

Leica TCS SP8 confocal microscope with 60× water-immersion objective (NA 1.2).

### Channels

DAPI, Alexa Fluor 488, Alexa Fluor 546, Alexa Fluor 594, and Alexa Fluor 647.

### Imaging Buffer

- Glucose: 10% (w/v)
- Catalase: 40 µg/mL (Sigma C30)
- Glucose oxidase: 0.5 mg/mL (Sigma G2133)
- RNase inhibitor: 0.02 U/µL
- Tris-HCl pH 8: 50 mM
- Buffer: 2× SSC

## Decoding

1. **Image registration**: Affine transformation matching subsequent rounds to first-round image.
2. **Spot detection**: White top-hat filtering and h_maxima detection.
3. **Intensity extraction**: 5×5×3 cube extraction around each spot; max intensity identification.
4. **Channel combination**: I_ch1,ch2 = sqrt(I_ch1 × I_ch2) for dual-channel signals.
5. **Filtering**: Signal-to-noise ratio thresholding; border spot removal.
6. **Barcode matching**: Combined codes matched to gene codebook.
7. **Quality control**: Signal validation per spot.

## Critical Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Padlock probe arm length | 12-16 bp | Complementary to target |
| Initiator primer arm length | 12-16 bp | Complementary to target |
| Padlock probe concentration | 30 nM/probe | For RNA detection |
| Initiator concentration | 30 nM/probe | Equimolar to padlock |
| SplintR ligase concentration | 1 U/µL | For RNA ligation |
| RCA temperature | 30°C | For Phi29 amplification |
| RCA duration | 2 hours (RNA) / 6 hours (DNA) | Longer for DNA |
| Sequencing ligation temperature | 25°C | For T4 DNA ligase |
| Sequencing ligation duration | 2 hours | |
| Formamide stripping | 60% in 2× SSC | Between rounds |
| Image Z-resolution | Interpolated to match reference | For registration |

## Version-Specific Differences

No significant version differences were identified for V01. The preprint (2023) and published paper (2024) describe the same workflow.

## Unreported Parameters

- Specific buffer pH values (not explicitly stated)
- Exact sodium chloride concentrations beyond 2× SSC
- Detailed washing volumes
- Specific Cat/Dog/Lamp/Pax/etc. gene selection criteria beyond "exon close to 5' end"
- Complete barcode library sequences
- Exact sequence of anchor primers and query probes
- Minimum number of probes per gene for robust detection
- Specificity filtering criteria
- Off-target alignment parameters
