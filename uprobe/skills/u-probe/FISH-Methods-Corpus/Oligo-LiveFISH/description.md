# Method Overview

## Purpose

V01 tracks endogenous genomic DNA loci in live cells, including non-repetitive regions. Demonstrations include human U2OS and HeLa cells, primary human T cells, and hESC-derived induced neurons. Reported: S01, Summary and Results, diverse-cell-type experiments.

## Core Principle

A pool of fRNPs places multiple fluorescent guide complexes across a locus. Each dual-RNA complex contains nuclease-dead dCas9, target-complementary crRNA, and a universal tracrRNA; fluorescent sgRNA is an alternative. Reported: S01, Results, “Design of Oligo-LiveFISH…”; S02, Table S1(6); V01.

## Overall Workflow

Select genomic regions and guides → synthesize DNA templates → PCR → in vitro transcription (IVT) → purify/fold and label RNA → assemble dCas9 fRNPs → electroporate → image and track. A dye-labeled universal tracrRNA avoids separately labeling every crRNA species. Reported: S01, Results and STAR Methods; V01.

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `crispr`.
- `readout_strategy`: `fluorescent_barcode`, `live_imaging`.

Structural evidence: the recognition reagent is dCas9 programmed by a crRNA:tracrRNA duplex or sgRNA; there is no contiguous oligo with target-binding arms that ligates into a circular template (see `probes.md` Canonical Architecture; Probe Components list guide-RNA species and not `padlock_arm`/ligation site). The assay is documented for live-cell imaging (S01 microscopy methods; multiple human cell types; SL-Oligo-LiveFISH), so `live_imaging` applies. Readout is direct multicolor fluorescence imaging of fluorescently labeled guides (Cy3-tracrRNA, 5′-dye sgRNA, etc.) without sequential hybridization/cleavage cycles; identity is conveyed by channel and locus combination (`fluorescent_barcode`). Source: S01, Results and STAR Methods; S02, Table S1.

## Target Recognition Strategy

The spacer directs dCas9 to genomic DNA; computational selection uses NGG PAM sites and screens clustered off-target binding. PAM is a target-DNA feature, not an appended RNA barcode. Reported: S01, STAR Methods, “Design tiling crRNAs/sgRNAs targeting non-repetitive regions”; V01.

## Encoding Strategy

Locus identity is assigned by the delivered guide pool and fluorescence channel. The documented workflow is direct multicolor imaging rather than a sequential combinatorial codebook. No physical barcode component or Hamming-distance design is specified for V01. Reported workflow / unavailable codebook: S01, Figures 1G, 6 and STAR Methods; S02, Table S1.

## Amplification Strategy

PCR amplifies reagent-manufacturing DNA templates; T7 IVT produces RNA. Multiple bound guides accumulate fluorescence at a locus. An in situ enzymatic amplification reaction is not part of the described live-cell labeling workflow. Reported: S01, preparation/assembly sections; V01.

## Imaging / Readout Strategy

The paper uses widefield, spinning-disk confocal, and double-helix point-spread-function microscopy. Its SL-Oligo-LiveFISH implementation localizes sparse labeled loci; reported performance depends on microscope, target, and exposure. The headline 20-nm/50-ms result is not a universal specification. Reported: S01, Figure 5 and microscopy methods; S02, Table S3; V01.

## Decoding Strategy

Image localization and trajectory analysis replace sequencing-code decoding. Quantification includes signal-to-noise ratio, localization fitting, mean-square displacement, and velocity correlations. Reported: S01, Quantification and statistical analysis; V01.

## Major Versions or Variants

V01 includes crRNA:tracrRNA pools, fluorescent sgRNA pools, alternative terminal labeling chemistries, and 2D/3D SL imaging. LiveFISH is explicitly described as earlier work, with a different usual labeling arrangement. No independent Oligopaint-to-Oligo-LiveFISH version lineage is established by examined evidence. S01, Introduction and Results; S02, Table S1.

## Strengths

The reported labeling procedure avoids engineering a stable CRISPR imaging line for each locus and supports non-repetitive regions and several human cell types. These are demonstrated capabilities within the study, not guarantees for every target. Reported: S01, Results and Discussion; V01.

## Limitations

Labeling depends on guide abundance and chromatin accessibility; excess delivered fRNPs can raise background. The paper recommends concentration optimization and filters nonspecific localizations. The DHPSF implementation supports up to two simultaneous colors in Table S3. Reported: S01, Figure 2 and Limitations; S02, Table S3; V01.

## Unresolved or Version-Specific Points

Table S2 lacks strand columns and explicit coordinate-origin declarations. The short-spacer ACAP2 sheet contains coordinate/sequence-length discrepancies. Full library-specific synthesized oligo assemblies and exact adapter assignments are not established for every spacer. See probes.md and benchmarks.md. Derived / unavailable: S03, specified sheets; V01.
