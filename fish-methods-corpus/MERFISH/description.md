# Method Overview

## Purpose

Measure the identities, copy numbers, and spatial distributions of many selected RNA species in single cells. This is targeted imaging, requiring a designed probe library (S01, Abstract, Fig. 1, Probe design; V01; reported).

## Core Principle

A gene is assigned an abstract binary codeword. Multiple encoding probes bound along each RNA collectively supply the readout sequences for its on-bits. A different fluorescent readout probe interrogates each bit; spatially registered images yield the measured word (S01, Fig. 1; V01; reported).

## Overall Workflow

Select targets and codebook; design and manufacture encoding probes; fix and permeabilize the sample; hybridize encoding probes; iteratively hybridize readout probes, image, and extinguish fluorescence; register signals and decode. Clearing, expansion, and bDNA are separate later modifications (S01, Materials and Methods; S04/S05/S08/S10; reported, version-specific).

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `linear`.
- `readout_strategy`: `sequential_fish`.

Structural evidence: encoding probes are linear ssDNA with 30-nt target-binding regions; Probe Components do not include `padlock_arm` or a ligation site, and there is no in situ circularization/RCA in V01 (see `probes.md` Canonical Architecture and Probe Components). Readout is iterative hybridization of fluorescent readout probes, imaging, and signal removal by photobleaching (V01) or chemical cleavage (V02); decoding reads binary codewords across many rounds (see `probes.md` Readout System and the Imaging / Readout Strategy section above). Source: S01, Fig. 1, Probe design and imaging Methods; S04, Table S1 (V02); S08, Fig. 1 (V05 bDNA amplifier is a later signal-amplification addition, not an architecture change). V03-V06 retain linear probes and sequential readout; classification applies at the method level and individual version overrides are documented under "Major Versions or Variants" if they arise.

## Target Recognition Strategy

The original encoding probes contain 30-nt regions complementary to target RNA. Many different sites on one RNA provide sufficient brightness. Target and isoform selection are experiment-specific (S01, Probe design; V01; reported).

## Encoding Strategy

V01 uses a 16-bit modified Hamming-distance-4 (MHD4) code of weight 4 with 140 codewords, or a 14-bit MHD2 code of weight 4 with 1,001 codewords. The experiments called “140-gene” and “1001-gene” include controls: 130 and 985 biological RNA species, respectively. A molecule does not carry one contiguous 16- or 14-nucleotide barcode; its identity is the multi-round signal pattern (S01, Results, Probe design; reported; distinction is functional interpretation).

## Amplification Strategy

The aggregate amplification metadata is null because amplification differs by version. V01 does not amplify target RNA in situ. PCR, in vitro transcription, and reverse transcription manufacture encoding-probe libraries ex situ. V05 adds a primary/secondary bDNA assembly to amplify readout signal; it does not convert the workflow into target RNA amplification or sequencing (S01, Probe synthesis; S08, Fig. 1 and Methods; reported).

## Imaging / Readout Strategy

V01 uses 3′ Cy5 readout probes and photobleaching between single-color rounds. V02 uses 5′ disulfide-linked Cy5/Alexa750 probes, reductive cleavage, and two-color readout. Consequently, 16 bits require 16 rounds in the original MHD4 experiment but eight rounds in V02 (S01, imaging Methods; S04, Table S1 and MERFISH Imaging; reported).

## Decoding Strategy

V01 registers rounds with fiducial beads, links nearby spots, and constructs binary words. MHD4 permits assignment of an exact or one-bit-different word; MHD2 accepts exact matches and cannot correct a one-bit error uniquely. Later workflows use normalized pixel/voxel intensity vectors and connected components (S01, Construction of measured words; S10, MERFISH data analysis; S09, Materials and Methods; reported).

## Major Versions or Variants

| ID | Meaningful change | Evidence |
|---|---|---|
| V01 | Original two-readout encoding probes; MHD4/MHD2 alternatives | S01, Fig. 1 and Methods; reported |
| V02 | Three 20-nt readouts per encoding probe, fast three-letter sequences, chemical cleavage and multicolor imaging | S04, SI Encoding Probe Design and Table S1; reported |
| V03 | Poly(A)-binding anchors, polyacrylamide imprinting and protein/lipid removal; cultured cells and mouse hypothalamus | S05, SI Methods and Fig. 4; reported |
| V04 | Expandable gel, digestion, expansion and re-embedding to reduce molecular crowding | S10, Results and Methods; S11 correction; reported |
| V05 | Primary and secondary bDNA amplifiers add repeated detection sites | S08, Fig. 1; reported |
| V06 | 69-bit HD4 weight-4 code, 10,050 genes, expansion and three-color readout over 23 rounds | S09, Results, first subsection; reported |

## Strengths

Code redundancy can identify and correct errors within the chosen code’s limits, and imaging retains spatial context. V01 comparisons with conventional smFISH and RNA-seq support quantitative measurements for the tested libraries (S01, Figs. 2–3; reported). These results are scoped to those experiments, not universal performance guarantees.

## Limitations

Optical crowding, nonspecific binding/background, incomplete probe binding, and errors across rounds limit detection. V01 MHD2 exchanges error correction for greater coding capacity. V05 empirically rejects amplifier pairs with high background or low gain, showing that sequence screening alone does not guarantee performance (S01, Discussion; S08, Fig. 4 and Discussion; reported).

## Unresolved or Version-Specific Points

Detailed design outputs in the original Table S5 and V02 Dataset S1 were identified but not parsed because access returned challenge/error pages. Do not treat the synthesis template as an independently verified final ssDNA molecule. S03 contains thermodynamic threshold wording that conflicts with S01; see probes.md. All later-version values remain scoped. No independent validation is claimed.
