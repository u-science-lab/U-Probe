# Method Overview

Unless otherwise stated, this file concerns V01, original RNA seqFISH+ (2019). Bracketed evidence applies to the associated paragraph.

## Purpose

Measure spatial expression and subcellular positions of mRNAs from a large targeted gene panel while preserving cell and tissue context. The original study profiled 10,000 genes in NIH/3T3 cells and in cortex, subventricular zone and olfactory bulb of mouse brain. This is targeted hybridization, despite its broad transcriptome coverage. [S01, abstract, Figs. 1–3; reported]

## Core Principle

Optical crowding is reduced by assigning genes to sparse sequential readouts. Twenty temporal pseudocolours are acquired in each of three fluorescence channels. Transcripts can be localized in the sparse images, then combined computationally. A pseudocolour denotes the hybridization identity within a barcoding round, not another dye or a nucleotide. [S01, Fig. 1a and Methods, seqFISH+ encoding strategy; reported]

## Overall Workflow

Fix and permeabilize samples; hybridize primary probes; embed in hydrogel and clear; stabilize retained nucleic acids; perform repeated readout hybridization, washing, imaging and stripping; register images, locate dots and decode the temporal patterns. Cell and tissue stabilization differ. [S01, Methods, Cell culture experiment, Tissue slices experiment and seqFISH+ imaging; Extended Data Fig. 1; reported]

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `linear`.
- `readout_strategy`: `sequential_fish`.

Structural evidence: each gene is targeted by 24 linear ssDNA primary probes that contain a 28-nt RNA-binding region and four 15-nt readout-binding overhangs; there is no in situ circularization or RCA template (see `probes.md` Canonical Architecture; Probe Components; Probe Components do not include a `padlock_arm` or ligation site). Identity is decoded across 80 sequential rounds of fluorescent-oligo hybridization, imaging, and formamide stripping (see `probes.md` Readout System and the Imaging / Readout Strategy section above). Source: S01, Fig. 1a and Methods, Primary-probe construction, Readout-probe design, seqFISH+ imaging and encoding strategy. The DNA seqFISH+ extension (Takei et al., 2021; S09, S10) is a separately classified method and is not included here.

## Target Recognition Strategy

Each gene is targeted by 24 primary probes. Target-binding sequences are 28 nt and selected from coding exons first, adding UTR exons when necessary. The final manufactured probes are ssDNA; RNA targets are not reverse-transcribed in the sample. [S01, Methods, Primary-probe design/construction and encoding strategy; reported]

## Encoding Strategy

A gene has a four-unit base-20 codeword, read entirely within one fluorescence channel. Three informative units provide 20³ = 8,000 combinations per channel, with a fourth unit for error correction. The authors report a capacity of 24,000 error-correctable barcodes across three channels; the actual experiment uses 10,000 genes. Twenty serial hybridizations per barcoding round give 80 serial hybridizations overall. The three channels are acquired in parallel across that schedule. [S01, Methods, seqFISH+ encoding strategy; reported]

The public codebook contains 3,334 genes in the 647-nm channel and 3,333 each in the 561- and 488-nm channels. These sheet-specific counts are obtained by counting gene rows, not by assuming channel order from an unlabelled list. [S02, Codebook-647nm, Codebook-561nm, Codebook-488nm, rows 2 onward; derived]

## Amplification Strategy

There is no in situ enzymatic target or signal amplification in this workflow. Multiple primary probes contribute signal per transcript. PCR, in vitro transcription and reverse transcription are manufacturing steps performed before sample hybridization. Do not describe them as tissue RNA amplification. [S01, Methods, Primary-probe construction and sample/imaging sections; reported workflow]

## Imaging / Readout Strategy

Fluorescent 15-nt DNA readouts bind complementary overhangs on primary probes, then are removed with formamide between acquisitions. Alexa Fluor 647, Cy3B and Alexa Fluor 488 provide the three readout channels. A spinning-disk confocal microscope acquires the sparse images. Hybridization 81 repeats hybridization 1 to assess signal retention; it is not a fifth barcoding round. [S01, Methods, Readout-probe synthesis and seqFISH+ imaging; Extended Data Fig. 2; reported]

## Decoding Strategy

Register rounds, deconvolve, find local maxima and refine dot positions with radial-centre localization. Match spatial neighbours across barcoding rounds and compare to the codebook; ambiguous calls are filtered using spatial residuals and Hamming-distance-one matching, then rejected if still ambiguous. Retain calls supported when at least three of four rounds are used as decoding seeds. This seed-consistency criterion is distinct from a claim that any arbitrary one-symbol error is always correctable. [S01, Methods, Image registration, Image processing and Barcode calling; reported]

## Major Versions or Variants

V01 is the 2019 RNA method. Cell and tissue preparations are protocol branches, not separately named method versions. The 60-gene smFISH comparison uses two readout-binding sites per primary probe and is a validation assay, not the four-site canonical seqFISH+ panel. [S01, Methods, smFISH; reported]

DNA seqFISH+ is a later developer extension reported by Takei et al. (2021). It images genomic loci and uses distinct probe and encoding strategies. Its publication and official tissue repository were inspected for lineage/resource separation only; no DNA probe rules, assembly, sequences or benchmark grade are imported into V01. A full DNA extension record remains outside this original RNA draft. [S09, abstract/Fig. 1; S10, README and probes directory listing; reported]

## Strengths

The demonstrated panel combines large targeted gene coverage with localized RNA measurements in cells and tissue. Coding within a single channel avoids requiring cross-channel colocalization for an individual barcode. The pseudocolour scheme reduces simultaneous dot density and permits sub-diffraction localization using a confocal platform. [S01, main text and Fig. 1; reported]

## Limitations

The assay requires many serial hybridizations and image registration. The original 60-gene smFISH comparison estimated 49% detection efficiency in the tested cell experiment; this is not a universal sensitivity. Some signal is lost across 80 rounds, as shown by the hybridization-1/81 comparison. These limitations are experiment-specific. [S01, Fig. 2d, Extended Data Fig. 2 and Methods; reported]

For design benchmarking, the inspected original resources expose a codebook and processed measurements but do not resolve a complete sequence-level probe set or exact reference release. [S01, Primary-probe design and data availability; S02–S04; S06 inventory; unavailable within examined scope]

## Unresolved or Version-Specific Points

The full 5′→3′ order of primary-probe components, spacer identities, manufacturing primer strings, final probe strings, readout nucleotide library, UCSC assembly/annotation release and exact design implementation remain unresolved. See probes.md and benchmarks.md. The S05 mirror gives 5–20-µm sectioning guidance, whereas S01 used 5-µm sections. Official-body equivalence of the mirror remains unresolved. [S01, Methods, reported; S05 mirror, Tissue section preparation, inferred; V01]
