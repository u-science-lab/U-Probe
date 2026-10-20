# Method Overview

## Purpose

Targeted spatial RNA profiling with conventional fluorescence microscopy. Mouse and human tissue applications are reported. [reported; S01, Abstract; V01/V02]

## Core Principle

Each padlock carries four physical barcode segments. Imaging probes have sequences identical to these segments and bind their complementary copies in RCA products. Labeled/unlabeled competitors encode intensity levels. The four-symbol color-intensity codeword is distinct from the physical DNA sequence. [reported; S01, Fig. 1a–c; V01/V02]

## Overall Workflow

Padlock binding → ligation → RCA → pooled imaging-probe staining → one round of multichannel imaging → intensity-space gene calling → cell assignment. [reported; S01, Fig. 1b–f; V01/V02]

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `fluorescent_barcode`.

Structural evidence: each padlock carries four physical barcode segments flanked by terminal RNA-binding arms; ligation on RNA circularizes the padlock (see `probes.md` Canonical Architecture; Probe Components include `padlock_arm_5prime`, `padlock_arm_3prime`, `chemical_modification_padlock` referencing 5'-phosphate, and an `amplification_primer` row). Phi29 RCA produces the in situ amplicon. Readout is a single imaging round of competitive labeled/unlabeled imaging-probe binding to RCA products, with identity encoded by intensity ratios and color combinations across channels (Fig. 1c, Supplementary Fig. 2); there is no per-base sequencing chemistry or sequential stripping cycle, so readout is `fluorescent_barcode` rather than `in_situ_sequencing` or `sequential_fish`. Source: S01, Fig. 1a–f and Extended Data Fig. 7; S02 imaging sheets.

## Target Recognition Strategy

Terminal padlock arms recognize an RNA binding site. The inspected 2024 preprint describes SplintR ligation directly on RNA and 40-nt targets split into two 20-nt segments; this is preprint-scoped evidence. [reported; S05, Probe and barcode design / Sample treatment for thin tissue; V01 precursor]

## Encoding Strategy

Radius-vector selection separates intensity ratios despite differences in RCA-product brightness. V01 uses 15 three-channel combinations, doubled with a fourth channel, with an optional fourth-channel-only codeword. V02 expands intensity grading. [reported; S04, Supplementary Fig. 2; S01, Extended Data Fig. 7; V01/V02]

## Amplification Strategy

Circularized probes undergo RCA. The preprint specifies Phi29, external RCA primers and aminoallyl-dUTP, followed by BS(PEG)9 fixation. [reported; S05, Sample treatment for thin tissue; V01 precursor; S02, RCA-primer sheets, V01/V02]

## Imaging / Readout Strategy

V01 is represented by 17 imaging-probe sequence identities, supplied as labeled/unlabeled competitors; the final V02 workbook includes extra sequence families for high-Cy3 barcodes. One imaging round contains multiple spectral acquisitions and may contain z-stacks; it is not a single exposure. [reported/derived; S02, imaging sheets; S05, Imaging; V01/V02]

## Decoding Strategy

Spot intensities are normalized and assigned in color space. The inspected preprint describes Gaussian fitting and manually delimited boundaries; the current PRISM-Code supports a post-stitching pipeline and should not be assumed identical to publication-time code. [reported; S05, Gene decoding; S06, README; V01 precursor/current software]

## Major Versions or Variants

V01 and V02 are configurations in the 2025 paper, not separate later publications. V02 changes intensity levels and associated barcode/readout libraries. The 2024 preprint methods are separately scoped wherever final-method confirmation is unavailable. [reported; S01, Extended Data Fig. 7; S02, 64-panel sheets; S05]

## Strengths

A single staining/readout cycle avoids repeated decoding-round registration. The public workbook supplies full oligo records and experimental pooling recipes. [reported; S01, Fig. 1c–f; S02, imaging/padlock sheets; V01/V02]

## Limitations

Probe-pool calibration depends on optical setup; the workbook provides empirical mixtures and notes a FRET-like effect in the high-green V02 layer. RNA-site selection cannot be reproduced strictly from the examined resources because original transcript/reference releases and original design-run settings remain unresolved. [reported; S02, 64 imaging sheet row 71; S05, Probe and barcode design; inferred benchmark consequence]

## Unresolved or Version-Specific Points

The preprint says two binding targets per mouse gene; final mouse worksheets list one padlock per gene. The source scopes are preserved, not reconciled by assumption. See probes.md for backbone-length, HBX-arm, channel-order and nomenclature qualifications. The subscription-restricted final Methods were not available for direct comparison.
