# STARmap

## Official Name

STARmap

## Full Name

spatially-resolved transcript amplicon readout mapping (S01 abstract)

## Aliases

- "Spatially-resolved Transcript Amplicon Readout Mapping" (capitalization variant used in some sources)
- "spatially-resolved transcript amplicon readout mapping" (lowercase variant used in the Science abstract)
- "3D STARmap" (the thick-tissue volumetric variant; see S01 main text)
- "SNAIL-STARmap" (when emphasizing the SNAIL probe-amplification component; see S01)
- "STARmap PLUS" (the V02 protein-co-detection extension; see S11 — Zeng et al. 2023 *Nature Neuroscience*)

## Method Category

In situ sequencing (ISS) using a two-probe SNAIL (specific amplification of nucleic acids via intramolecular ligation) padlock system, T4-DNA-ligase circularization, rolling-circle amplification (RCA) with Phi29 DNA polymerase, and SEDAL (sequencing with error-reduction by dynamic annealing and ligation) sequencing-by-ligation readout, all performed in tissue-hydrogel matrices (hydrogel-tissue chemistry, HTC). See `description.md` for details.

## Primary Application

3D intact-tissue sequencing of single-cell transcriptional states at subcellular resolution in mouse brain (V1 visual cortex, mPFC, hippocampus) and, in STARmap PLUS, integrated protein + RNA profiling in mouse models of neurodegeneration. STARmap has been demonstrated for 160 to 1,020 genes simultaneously at the cellular scale, and for tens of thousands of cells in cubic-millimeter 3D volumes.

## Original Publication

Wang X., Allen W.E., Wright M.A., Sylwestrak E.L., Samusik N., Vesuna S., Evans K., Liu C., Ramakrishnan C., Liu J., Nolan G.P., Bava F.-A., Deisseroth K. (2018). Three-dimensional intact-tissue sequencing of single-cell transcriptional states. *Science* 361(6400):eaat5691. DOI: 10.1126/science.aat5691. PMID 29930089. PMCID PMC6339868.

## Brief Summary

STARmap combines paired RNA recognition, padlock RCA, covalent hydrogel retention and SEDAL barcode sequencing. V01 demonstrated panels up to 1,020 genes and volumetric tissue mapping (reported: S01 Results; Fig. 1). STARmap PLUS adds protein comapping and improved probe design (S03 Overview; S11). Detailed later-protocol conditions are version scoped in protocols.md; uninspected benchmark candidates remain unassessed. See the production files below for selective retrieval.

## Known Versions

| Version | Year | Reference | Status | Brief note |
|---|---|---|---|---|
| V01 | 2018 | Wang et al. *Science* 361:eaat5691 | original | Foundational STARmap: acrylic acid NHS ester HTC; 5-nt barcode on padlock probe; 160- to 1020-gene panels; 8-µm imaged volumes and 150-µm tissue slices in mouse cortex. |
| V02 | 2023 | Zeng et al. *Nat Neurosci* 26(3):430–446 (STARmap PLUS) | variant | STARmap PLUS: RNA/protein comapping; 2,766-gene AD panel (S11); later S03 uses MA-NHS. |

The joint STARmap PLUS / RIBOmap / TEMPOmap protocol (S03) covers V01 and V02 design rules; RIBOmap and TEMPOmap are sibling methods on the same SNAIL-STARmap platform but are not STARmap versions and are not recorded here. See `metadata.yaml` for the version records.

## Corpus Files

- [index.md](index.md) — this file
- [metadata.yaml](metadata.yaml) — machine-readable routing and structured fields
- [description.md](description.md) — conceptual method overview
- [protocols.md](protocols.md) — experimental protocol
- [probes.md](probes.md) — probe architecture and design
- [benchmarks.md](benchmarks.md) — benchmark candidate assessment
- [sources.md](sources.md) — source provenance
- [research_manifest.md](research_manifest.md) — audit file (excluded from ordinary retrieval)

<!-- Corpus Files: link index.md, metadata.yaml, description.md, protocols.md, probes.md, benchmarks.md, sources.md. Keep routing concise; no detailed protocol/design. -->
