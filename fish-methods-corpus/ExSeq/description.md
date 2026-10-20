# Method Overview

## Purpose

ExSeq combines expansion microscopy with in situ RNA sequencing to map transcripts in intact specimens and relate their positions to subcellular morphology [S01, abstract and Figs. 1-6; reported, V01-V03].

## Core Principle

RNA is anchored to a swellable hydrogel with LabelX. Digestion, expansion, re-embedding and carboxylate passivation make the expanded specimen compatible with enzymatic library preparation. Untargeted ExSeq generates cDNA circles after random-primed reverse transcription; targeted ExSeq circularizes DNA padlocks directly on RNA with SplintR, bypassing reverse transcription. Phi29 RCA produces amplicons for in situ sequencing [S01, Materials and methods, Figs. 1 and 4; S02, supplementary pp. 3-8, 18-23; reported, V01-V03].

## Overall Workflow

1. Anchor RNA, gel the specimen, digest proteins, expand, re-embed and passivate the hydrogel.
2. Prepare either a random-primed cDNA library (V01) or an RNA-splinted padlock library (V02-V04), then circularize and RCA-amplify it.
3. Stabilize amplicons, immobilize the specimen for repeated chemistry, and perform sequencing and fluorescence imaging.
4. Register images and decode spatially localized reads. Untargeted experiments additionally use ex situ Illumina sequencing of extracted amplicons to interpret short in situ reads.

Reported workflow: S01 Materials and methods; S02 supplementary pp. 5-11, 20-23, 29-31. Generalized targeted execution is documented separately by S11, S13 and S14 (V04); those recommendations must not substitute for original experiment conditions. Polymerization recipes, including the distinct first and second re-embedding gels, are in protocols.md.

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `in_situ_sequencing`.

Structural evidence: V02-V04 are RNA-targeted padlocks with terminal homology arms that pair antiparallel with adjacent portions of a contiguous RNA region; SplintR closes a nick to circularize the padlock (see `probes.md` Canonical Architecture; Probe Components include `padlock_arm` and a ligation site); amplification is Phi29 RCA. V01 uses random-octamer RT to generate cDNA circles (still a circular RCA template, not a linear primary probe). Readout is SOLiD sequencing-by-ligation (V01/V03) or Illumina SBS (V02), both per-base in situ sequencing reactions on the RCA product (see `probes.md` Readout System and the Imaging / Readout Strategy section above). Source: S02 supplementary pp. 5-8, 18-23, 29-31; S12, S13, S14 (V04). Untargeted V01's lack of a predefined gene barcode (in situ/ex situ matching) does not change the architecture or readout labels.

## Target Recognition Strategy

Targeted padlocks carry two terminal homology arms that pair antiparallel with adjacent portions of a contiguous RNA region. The default 32-nt recognition region is split into two 16-nt arms. SplintR closes a nick between the padlock 3'-OH and 5'-phosphate; no extra RNA nucleotide is inserted between the paired halves [S02, supplementary pp. 18-20; S12, Fig. 1A and Step 9.2; reported; S04, homology selection and assembly lines 463-482, inferred strand interpretation; V02-V04]. Untargeted recognition uses modified random-octamer RT primers, with sample-specific constant regions [S02, supplementary pp. 5-6; reported, V01].

## Encoding Strategy

- V02 original visual cortex/hippocampus: four logical symbols in {0,1,2,3}, with R4=(R1+R2+R3) mod 4. Physical SOLiD and Illumina barcode regions are each seven nucleotides, but the brain experiments acquired four SBS rounds. The last three physical positions do not imply seven decoded rounds [S02, supplementary pp. 19-20, 23; reported].
- V03 original cancer: length-7, minimum Hamming-distance-3 logical barcodes generated with DNABarcodes; the analysis permits one substitution mismatch. The two physical barcode regions follow the experiment-specific implementation in Fig. S12 [S02, supplementary pp. 29-31; reported].
- V04 later generalized design: the supplied 326-entry library has length 7 and minimum distance 3, with a 58-entry subset whose first four symbols have distance 2. Barcode conversion depends on the selected readout specification; this library is not attributed to every original experiment [S12, Steps 3-4, 12; reported; S04 CSV counts/distance, derived].

V01 has no predefined gene barcode: biological cDNA sequence supplies identity through in situ/ex situ matching [S01, Fig. 1C; reported].

## Amplification Strategy

Phi29 performs RCA on cDNA circles (V01) or circularized padlocks (V02-V04). In the generalized targeted protocol the RCA primer is the 18-nt DNA oligo `TCT TCA GCG TTC CCG A*G*A`, with phosphorothioate bonds denoted by `*`. It binds the padlock backbone; the 22-nt SBS sequencing primer instead binds its complement on the RCA product. The RCA mix contains 250 µM **each** dNTP plus 40 µM aminoallyl-dUTP, followed by 5 mM BS(PEG)9 cross-linking. These solution concentrations do not establish an incorporation fraction [S13, Steps 5, 10-12; S14, Step 10; reported, V04]. Original untargeted oligos and conditions differ by sample [S02, supplementary pp. 5-6; reported, V01].

## Imaging / Readout Strategy

Original untargeted ExSeq uses SOLiD sequencing by ligation [S02, supplementary pp. 7-8; reported, V01]. Original brain panels use four Illumina SBS rounds [S02, pp. 22-23; reported, V02]; cancer uses SOLiD ligation [S02, pp. 30-31; reported, V03]. FISSEQ means fluorescent in situ sequencing; it does not specify synthesis chemistry.

The generalized protocols.io collection includes a dedicated Illumina module: incorporation, four-color imaging, cleavage, and repetition for each barcode base. Its full length-7 implementation uses seven rounds [S14, Steps 11-15; reported, V04]. S12 Step 12.2 explains SOLiD barcode implementation. Detailed original SOLiD experimental conditions are independently available in S02, rather than only in that barcode-design explanation.

## Decoding Strategy

Original brain analysis retains perfect matches across all four logical bases [S02, supplementary p. 23; reported, V02]. Cancer analysis allows one mismatch in seven-symbol barcodes [S02, p. 31; reported, V03]. V04 decoding must use the selected barcode specification and experiment's gene-to-code assignment [S12, Steps 3-7, 12; reported].

For V01, short in situ reads are matched to the ex situ library: 92% of **matches** overall and 97% of matches aligned against non-rRNA were unique. These are not percentages of all acquired in situ reads. RNAseq was performed on an adjacent hippocampal slice: r=0.89 compares full ex situ ExSeq expression with RNAseq (Fig. 2Giv), whereas ten imaged fields yield r=0.47 (Fig. 2Gv). The ten fields detected 3,039 genes [S01, matching paragraph, Biological validations, Fig. 2Giv-vii; reported].

## Major Versions or Variants

V01-V04 are internal corpus identifiers, not official release names. V01 distinguishes random RT/cDNA architecture; V02 distinguishes original targeted brain checksum encoding and four-round SBS; V03 distinguishes the original cancer length-7 distance-3 system and SOLiD protocol; V04 distinguishes the generalized protocol/code recommendations, including the embedded four-symbol library. V01-V03 were all published in Alon et al. 2021; same publication year does not imply identical architecture or encoding [S02, supplementary pp. 5-8, 18-23, 29-31; S12, abstract and Steps 3-4, 12; reported distinctions, inferred version grouping].

## Strengths

- Expansion supports spatial localization of reads relative to subcellular morphology [S01, Figs. 1-3; reported, V01].
- Targeted ExSeq yielded approximately 62% of HCRv3.0-ExFISH counts in the reported same-cell HeLa comparison (slope 0.62; Pearson r=0.991, 60 cells). This is a relative comparison, not equality of sensitivity [S01, Fig. 4B; reported].
- Demonstrations include untargeted specimens and targeted mouse brain/cancer panels [S01, Figs. 1-6; reported, V01-V03].
- The developers provide probe-design code and analysis software [S12, abstract/Step 1; S01, Data and materials availability; reported, V04 and original analysis respectively].

## Limitations

- Untargeted expression comparisons depend on the imaged volume and on ex situ matching; ten-field and full ex situ correlations are distinct [S01, Fig. 2G; reported, V01].
- Multiround chemistry requires repeated imaging and registration. Original SOLiD experiments and generalized SBS have different cycle schedules [S02, supplementary pp. 7-8, 30-31; S14, Steps 11-15; reported].
- The generalized design pipeline uses legacy BLAST/MATLAB interfaces; the protocol warns about `blastlocal` compatibility [S12, Step 1.1; reported, V04].
- 74 nt and 16 probes per target are generalized defaults, not universal restrictions. Original visual cortex uses at most 8 probes/transcript; hippocampus at most 16 [S12, Steps 9.2, 13.2; S02, supplementary p. 20; reported].
- Original brain and cancer probe resources are publicly directed to Tables S9/S10, but their workbook cells remain uninspected in this corpus; exact selection reproducibility is not established [S02, supplementary pp. 20, 29, 90; reported resource scope].

## Unresolved or Version-Specific Points

The supplementary workbook's individual cells, historical off-target RefSeq database releases, and exact original selection reproducibility remain pending. The untargeted GRCm38/mm10 analysis reference must not be assumed to be a targeted-design reference. Source conflicts in original Tm wording and StockX arithmetic are retained in probes.md and protocols.md. Software availability does not establish active maintenance. No universal RCA copy count is assigned [S02, supplementary pp. 11, 19; S11, Step 3.1; S04, inspected revision; examined scopes].
