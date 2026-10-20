# Method Overview

## Purpose

RAEFISH profiles large predefined transcript panels with spatially resolved fluorescence imaging. “Whole-genome” describes panel scale, not unbiased discovery of arbitrary RNA sequences. V01 targets 23,312 human or 21,955 mouse genes; V02 reads engineered gRNAs (reported: S01, Summary and STAR Methods, Probe design).

## Core Principle

A reverse-padlock has invariant ends and an internal RNA-binding segment. A second RNA-binding oligo primes RCA when bound adjacent to the padlock. A splint assists ligation away from the RNA-interacting segment and is removed using a toehold oligo. Encoding oligos bind the amplicon and present readout-binding sites (reported: S01, Figure 1A and Probe design).

## Overall Workflow

Prepare probes by pooled synthesis and amplification; fix/permeabilize specimens; hybridize padlock and priming probes; ligate; remove splint; perform RCA; hybridize encoding probes; sequentially hybridize/read/strip fluorescent oligos; register images, decode molecules and assign them to cells (reported: S01, STAR Methods).

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `sequential_fish`.

Structural evidence: a reverse-padlock DNA is circularized by a splint on adjacent RNA targets, and a linear priming oligo supplies the RCA primer (see `probes.md` Canonical Architecture; Probe Components include `padlock_arm`, `ligation_sequence_splint`, an `amplification_sequence` row, and a `chemical_modification_padlock` row referencing 5'-phosphate). Phi29 RCA produces the in situ amplicon. Readout is sequential fluorescent-oligo hybridization across many rounds (47 dual-color rounds for V01; 7 for V02) with stripping between cycles (see `probes.md` Readout System and the Imaging / Readout Strategy section above). Source: S01, Fig. 1A, Probe design, Probe synthesis and Sequential secondary probe hybridization and imaging.

## Target Recognition Strategy

V01 uses two adjacent 30-nt RNA-binding regions with no overlap and a 0- or 1-nt gap. Human targets receive one probe set; most mouse targets receive two. V02 uses a 40-nt reverse-padlock recognition region and a shared 45-nt scaffold-binding priming region (reported: S01, Probe design).

## Encoding Strategy

V01: 94-bit binary codes, weight 4, minimum Hamming distance 4. Four physical 20-nt binding sites on each encoding probe implement the four positive bits. V02: 14-bit, weight-4, distance-2 codes shared by three gRNAs against the same perturbation gene. Molecular readout sequences and abstract codewords are distinct (reported: S01, Probe design).

## Amplification Strategy

Phi29 RCA amplifies a circularized reverse-padlock DNA template in situ. Separate PCR/transcription/reverse-transcription steps manufacture probe pools before the assay; they do not imply reverse transcription of specimen RNA (reported: S01, Probe synthesis and RAEFISH protocols).

## Imaging / Readout Strategy

V01 uses 47 dual-color imaging rounds for 94 bits, with adapter oligos and two common dye-labeled oligos. V02 uses 7 dual-color rounds with directly dye-labeled readout probes. Matched MERFISH expression profiling in the perturbation experiment is a separate 26-bit/13-round assay (reported: S01, Sequential secondary probe hybridization and imaging).

## Decoding Strategy

Endogenous RAEFISH uses fiducial registration, normalized pixel-to-code matching, connected-pixel grouping and blank-code-calibrated foci selection. V02 uses a separate gRNA decoding pipeline (reported: S01, Image analysis; S03 archive inventory).

## Major Versions or Variants

V01 and V02 are separated here because recognition lengths, coding distance, imaging cycles and readout reagents differ. Human and mouse V01 libraries differ in target coverage and number of probe sets. S06 is the preprint of S01; its supplementary-table numbering must not be substituted for final Table S1.

## Strengths

Demonstrated compatibility with transcriptome-scale predefined panels, subcellular RNA localization and short engineered gRNA targets. Invariant ends permit pool amplification; RCA supports fluorescent detection in thin tissue sections (reported: S01, Results, Figure 1 and Figure 7). These observations do not establish universal efficiency for all targets.

## Limitations

Target-dependent sensitivity and long sequential imaging remain limitations. S01, Discussion, relates under-detection of some targets to probe sequence features but presents secondary-structure/cross-hybridization explanations as hypotheses. The method is targeted; results in 8–10-µm sections do not establish performance in arbitrary thick intact organs (reported: S01, Discussion, Mouse tissue section and Image system).

## Unresolved or Version-Specific Points

GC/Tm acceptance thresholds, alignment filters, isoform-selection algorithm and exact design-reference accessions remain unresolved in the examined official material. Probe tables provide sequences and assignments; S03 MAT resources also contain transcript sequences and numeric ENSG/ENST identifiers. Their stable versioned accession interpretation remains unresolved. Duplicate symbols must not be collapsed. See probes.md and benchmarks.md.
