# Method Overview

## Purpose

MiP-seq (Multi-omics in situ Pairwise Sequencing) enables simultaneous detection of DNA, RNA, proteins, and biomolecules at subcellular spatial resolution. The method combines dual-barcode padlock probe technology with rolling-circle amplification (RCA) and pairwise sequencing-by-ligation to achieve high-throughput multiplexed detection with fewer imaging rounds than conventional in situ sequencing methods.

## Core Principle

MiP-seq employs targeted padlock probes that hybridize to specific nucleic acid sequences or antibody-conjugated oligonucleotides. Upon target binding, the padlock probe is circularized by a target-dependent ligation step, creating an RCA template. Rolling-circle amplification generates rolling-circle products (RCPs) that contain multiple copies of the barcode sequence. Readout occurs via pairwise sequencing-by-ligation, where fluorescent query probes interrogate paired bases in the dual barcodes, enabling efficient decoding with 10N code capacity from N rounds of sequencing.

## Overall Workflow

1. **Sample preparation**: Tissue sections or cells are fixed, permeabilized, and pretreated.
2. **Probe hybridization**: Padlock probes and initiator primers are hybridized to target sequences (DNA, RNA, or antibody-conjugated oligonucleotides).
3. **Ligation**: SplintR ligase (for RNA) or T4 DNA ligase (for DNA) circularizes the padlock probe, creating an RCA template.
4. **Rolling-circle amplification (RCA)**: Phi29 polymerase amplifies the circularized probe, generating RCPs.
5. **Sequencing-by-ligation**: Multiple rounds of anchor primer and fluorescent query probe ligation decode the dual barcodes.
6. **Imaging**: Confocal microscopy captures fluorescence signals after each sequencing round.
7. **Decoding**: Image analysis identifies barcode signals and assigns gene identities.
8. **Cell segmentation and analysis**: Cellpose segments cells; decoded transcripts are assigned to cells.

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `in_situ_sequencing`.

Structural evidence: each padlock probe has two terminal target-binding arms (12-16 bp each); SplintR (RNA) or T4 (DNA) ligase circularizes the padlock to form the RCA template (see `probes.md` Canonical Architecture; Probe Components include `Target-binding arm 1`/`Target-binding arm 2`, `Barcode B1`/`B2`, and anchor/query probes). Phi29 RCA produces rolling-circle products containing tandem barcode copies. Readout is pairwise sequencing-by-ligation: anchor primers bind adjacent to the dual barcodes and fluorescent query probes interrogate paired bases per round, with N rounds reading 10^N codes (see `probes.md` Readout System and the Imaging / Readout Strategy section above). Source: S01, Fig. 1B and Materials and Methods (V01). The multi-omics capability (antibody-conjugated oligonucleotides, SNP/m6A detection) does not change the padlock + sequencing-by-ligation classification.

## Target Recognition Strategy

- **RNA detection**: Direct padlock probes hybridize to target RNA sequences; SplintR ligase circularizes the probe upon adjacent binding.
- **DNA detection**: Padlock probes hybridize to denatured genomic DNA; T4 DNA ligase performs circularization.
- **Protein/biomolecule detection**: Primary antibodies are conjugated with DBCO-modified oligonucleotides; a corresponding padlock probe detects the oligo, enabling protein detection via the same workflow.
- **Single-nucleotide specificity**: Ligase-dependent circularization requires precise base-pairing, enabling discrimination of SNPs and RNA modifications (e.g., m6A).

## Encoding Strategy

MiP-seq uses a **dual-barcode combinatorial encoding system**:
- Each padlock probe contains two barcode regions (B1 and B2).
- Pairwise sequencing interrogates paired bases (one from each barcode) in each round.
- Ten distinct dual-base combinations can be decoded per imaging channel pair.
- **10N codes** can be interrogated with N rounds of sequencing (versus 4N for conventional single-barcode methods).

## Amplification Strategy

**Rolling-circle amplification (RCA)** is employed:
- Phi29 polymerase extends from the circularized padlock probe.
- The amplified RCP contains multiple tandem copies of the barcode sequence.
- RCPs appear as fluorescent dots that are localized and counted.
- The initiator primer (IP) provides additional specificity by requiring three-way junction formation between IP, padlock probe, and target RNA.

## Imaging / Readout Strategy

- **Sequencing chemistry**: Sequencing-by-ligation using T4 DNA ligase.
- **Readout probes**: Four fluorescently labeled query probes interrogate paired bases.
- **Fluorophores**: Alexa Fluor 488, Alexa Fluor 546, Alexa Fluor 594, and Alexa Fluor 647.
- **Imaging system**: Leica TCS SP8 confocal microscope with 60× water-immersion objective (NA 1.2).
- **Signal stripping**: 60% formamide in 2× SSC removes query probes between rounds.
- **Channel acquisition**: DAPI, Alexa Fluor 488, Alexa Fluor 546, Alexa Fluor 594, and Alexa Fluor 647.

## Decoding Strategy

1. Raw images from multiple sequencing rounds are registered to the first-round image using affine transformations.
2. Candidate spots are identified by intensity thresholding and local maxima detection.
3. For each spot, fluorescence intensities from four channels across all rounds are extracted.
4. Dual-channel signals are combined (sqrt of product) to determine paired-base identity.
5. Barcodes are matched to the gene codebook to assign gene identity.
6. Cell segmentation (via Cellpose) assigns decoded transcripts to individual cells.

## Major Versions or Variants

Only one version (V01) has been documented. The method was first introduced as a preprint (bioRxiv, January 2023) and published in peer-reviewed form (Nature Biomedical Engineering, May 2024). The workflow has remained consistent between preprint and publication.

The **sequential dilution strategy** is described as a variant within V01 to resolve optically crowded signals: genes are allocated to different barcode sets with unique dilution primers, and reference genes are sequenced in every round for alignment.

## Strengths

- **High decoding capacity**: 10N codes achievable with N rounds (vs. 4N for single-barcode methods).
- **Reduced sequencing rounds**: Approximately 50% fewer rounds required compared to conventional in situ sequencing.
- **Single-nucleotide specificity**: Enables SNP detection, allele-specific expression, and RNA modification discrimination.
- **Multi-omics capability**: Detects DNA, RNA, proteins, and neurotransmitters in the same workflow.
- **Subcellular resolution**: Cytoplasmic and nuclear transcript localization is achievable.
- **Integration with functional imaging**: Compatible with calcium imaging and Raman imaging for multi-dimensional tissue atlases.
- **Spatial transcriptomics at scale**: Demonstrated with 100+ genes in single experiments and 217 genes in infection studies.

## Limitations

- **Targeted approach**: Requires prior knowledge of target sequences; not suitable for transcriptome-wide discovery.
- **Probe design dependency**: Performance depends on padlock probe design quality and target accessibility.
- **Fixed detection chemistry**: Sequencing-by-ligation chemistry limits flexibility in readout probe design.
- **Optical crowding**: High-expression genes may produce crowded signals; addressed by sequential dilution strategy but with increased complexity.
- **Protein detection indirect**: Requires antibody-oligo conjugation; not applicable to all targets without appropriate antibodies.
- **Specialized reagents**: Requires SplintR ligase, specific fluorophores, and custom software.

## Unresolved or Version-Specific Points

- Whether subsequent publications describe updated probe architectures or workflow modifications.
- Complete sequence information for the fixed barcode library and readout probes.
- Detailed primer design rules beyond the brief description (12-16 bp arms, exon-proximal selection).
- Whether the method has been adapted by other laboratories or commercialized.
