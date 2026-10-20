# Probe Architecture

## Overview

Not yet examined.

## Canonical Architecture

Not yet examined.

<!-- Neutral schematic only, never an exact sequence:
5' — [Component A] — [Component B] — [Component C] — 3'
Record each molecular species, topology, strand, pairing, orientation, version, evidence and source locator. Unknown order/components stay explicit. -->

<!-- Classification support: this section is the primary structural evidence for `metadata.yaml.classification.probe_architecture`. A circular topology plus a ligation site plus RCA amplification supports `padlock`. A contiguous 5'→3' oligo with no mandatory ligation supports `linear`. A guide-RNA species (sgRNA / crRNA / tracrRNA) with a Cas-family protein supports `crispr`. A hierarchical or dendritic nanostructure (e.g., tetrahedral DNA dendrimer) attached to a primary recognition element supports `nanostructure`. Cite the version and source locator. -->

## Probe Components

Not yet examined.

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|

<!-- Source / locator: a single canonical column name across every table in this file.
Locator must include page, section, figure, table, file, line, or accession. Add version.
Sequence class: fixed | variable | barcode_library | sample_specific | version_specific | unresolved.
Required: yes | no | conditional | unclear.
Evidence status: reported | derived | inferred | unavailable.
Source / locator includes locator and version. Preserve original author terms.
For every sequence record, identify representation type, molecule, direction, source notation and transformations. A schematic is never an exact sequence. -->

<!-- Classification support: presence of `padlock_arm`, `ligation_site`, and an `amplification_sequence` row referencing rolling-circle amplification supports `padlock`. Presence of guide-RNA species (sgRNA / crRNA / tracrRNA) and absence of canonical target-binding arms supports `crispr`. Presence of predefined nanostructure monomer components (e.g., tetrahedral DNA dendrimer layers) attached to a primary recognition probe supports `nanostructure`. None of those features with contiguous 5'→3' recognition supports `linear`. -->

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | Not yet examined. | | | |
| Target region | Not yet examined. | | | |
| Probe length | Not yet examined. | | | |
| GC constraints | Not yet examined. | | | |
| Tm constraints | Not yet examined. | | | |
| Specificity filtering | Not yet examined. | | | |
| Repeat masking | Not yet examined. | | | |
| Exon selection | Not yet examined. | | | |
| Isoform handling | Not yet examined. | | | |
| Sequence complexity | Not yet examined. | | | |
| Spacing | Not yet examined. | | | |
| Overlap | Not yet examined. | | | |
| Number of probes per target | Not yet examined. | | | |

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | Not yet examined. | | | |
| Number of barcodes | Not yet examined. | | | |
| Codebook | Not yet examined. | | | |
| Hamming distance | Not yet examined. | | | |
| Error correction | Not yet examined. | | | |
| Round assignment | Not yet examined. | | | |
| Channel assignment | Not yet examined. | | | |
| Fixed vs experiment-specific barcode sets | Not yet examined. | | | |

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | Not yet examined. | | | |
| Readout probe | Not yet examined. | | | |
| Fluorophore | Not yet examined. | | | |
| Imaging round | Not yet examined. | | | |
| Sequencing primer | Not yet examined. | | | |
| Decoding role | Not yet examined. | | | |

<!-- Classification support: this section is the primary structural evidence for `metadata.yaml.classification.readout_strategy`. Multiple imaging rounds with stripping/cleavage between cycles support `sequential_fish`. A sequencing-by-ligation or SBS reaction with a sequencing primer / query probe supports `in_situ_sequencing`. A single-round multichannel intensity/color codeword supports `fluorescent_barcode`. A documented live-cell compatible workflow in `protocols.md` supports `live_imaging`. Cite the version and source locator. -->

## Primer / Linker / Adapter Components

Not yet examined.

## Amplification-Related Components

Not yet examined.

## Chemical Modifications

Not yet examined.

## Probe-Design Workflow

Not yet examined.

## Official Probe-Design Software

Not yet examined.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|

<!-- Official status requires evidence of endorsement. Date availability and maintenance observations. -->

## Public Probe Sets

Not yet examined.

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|

## Version Differences

Not yet examined.

## Unresolved Information

Not yet examined.

