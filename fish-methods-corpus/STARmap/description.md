# Method Overview

## Purpose

STARmap (Spatially-resolved Transcript Amplicon Readout Mapping) is a targeted in situ sequencing method that retrieves high-content gene-expression information while preserving the three-dimensional positional anatomy of intact tissue at single-cell resolution. The published method was designed to map 160 to 1,020 genes simultaneously in mouse brain sections and in thick-tissue volumes, with the goal of integrating gene expression and spatial organization in complex biological tissues (S01 abstract).

## Core Principle

STARmap combines four interlocking principles:

1. **SNAIL (specific amplification of nucleic acids via intramolecular ligation)**: a two-probe padlock system in which a primer probe and a padlock probe must both hybridize to the same RNA molecule before T4 DNA ligase can circularize the padlock probe. The primer probe itself is not circularized; it serves as the ligation template and the primer for rolling-circle amplification (S01 abstract and main text: "a pair of primer and padlock probes (fig. S1A)"; "only when both probes hybridize to the same RNA molecule can the padlock probe be circularized and rolling-circle-amplified").
2. **Rolling-circle amplification (RCA)**: Phi29 DNA polymerase extends the primer probe around the circularized padlock, generating a DNA nanoball amplicon at the site of the original RNA. Aminoallyl-dUTP is incorporated into the amplicon to install primary amine groups for downstream chemistry (S01 main text).
3. **Hydrogel-tissue chemistry (HTC)**: the amine-functionalized amplicons are coupled to an acrylamide hydrogel via N-hydroxysuccinimide (NHS) ester chemistry (acrylic acid NHS ester in S01 and S11; methacrylic acid NHS ester / MA-NHS in the later S03 protocol), covalently anchoring the amplicon to the tissue-hydrogel matrix. In V01, protein digestion and lipid removal improve transparency; the later protein-comapping workflow skips proteinase K to preserve protein signals (S01 main text; S03 protocol).
4. **SEDAL (sequencing with error-reduction by dynamic annealing and ligation)**: a sequencing-by-ligation readout that reads a 5-nucleotide combinatorial barcode on the padlock probe in 6 imaging cycles (4 fluorescence channels; the separate 28-gene thick-tissue configuration uses seven linear readout rounds) using two-base (dinucleotide) encoding. Dynamic annealing/ligation supports specificity, while two-base decoding rejects inconsistent reads, giving a reported per-read error rate of ~1.8% vs. ~29.4% for cPAL methods (S01 main text, Fig. 1E).

Together these four components allow STARmap to (a) target a defined set of genes with high specificity (because both probes must hit the same RNA), (b) generate a localized, optically resolvable amplicon per RNA target, (c) clear the tissue for high-resolution imaging, and (d) read out a gene-unique barcode in a small number of imaging cycles.

## Overall Workflow

V01 conceptual workflow (reported: S01 Results; Fig. 1): RNA recognition by primer/padlock pairs, padlock circularization, Phi29 RCA, amine functionalization and acrylic-NHS hydrogel embedding, tissue clearing, cyclic SEDAL imaging and barcode decoding. Fluorescent Nissl staining after six sequencing rounds supports cell-body segmentation in the 160-gene experiment. Numerical conditions from the later S03 procedure are kept separately in protocols.md; they do not reconstruct the 2018 experiment. No second embedding recipe is attributed to S03.

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `in_situ_sequencing`.

Structural evidence: the recognition step uses a SNAIL primer/padlock pair that binds adjacent sites on the RNA target; the padlock's 5' and 3' ends pair with 6 nt each of the primer, and T4 DNA ligase joins them on that DNA splint, so the primer, not the RNA, templates the ligation (S02, Fig. S1A legend; see `probes.md` Canonical Architecture; Probe Components include the padlock 5′/3′ ends, the ligation substrates, and the primer-pairing region that splints the ligation and primes RCA); amplification is Phi29 rolling-circle amplification. Readout is SEDAL sequencing-by-ligation on the RCA amplicon (see `probes.md` Readout System), which reads a 5-nt barcode two bases at a time across 6 cycles. Source: S01 main text and Fig. 1; S02 Fig. S1A (V01). The same classification applies to V02 (STARmap PLUS); the protein-comapping addition does not change padlock architecture or SEDAL readout.

## Target Recognition Strategy

STARmap uses **two separate DNA probes that must both hybridize to the same RNA molecule** for signal generation. This dual-probe requirement provides target specificity at the recognition step. The original Wang 2018 paper reports a single gene-targeting segment on each probe. The S01 main text does not enumerate per-component lengths; its supplementary materials (S02) do: each primer or padlock has 19–25 nt that hybridize to the RNA with a designed Tm of 60 °C, and primer and padlock pair over only 6 nt on each arm, with a Tm below room temperature so that primer–padlock hybridization is negligible during DNA–RNA hybridization at 40 °C but allows T4 DNA ligation in the following step (S02, Fig. S1A legend, p. 10). The two RNA-binding halves are separated by a 0–2 nt gap in V01 (S02, STARmap probe design, p. 3; the later S03 protocol recommends a 1–2 nt interspace). The 6-nt arms of the padlock anneal to the primer probe and bring the padlock's 5' and 3' ends together for T4 ligation.

**Method version: V01.** S03 explicitly describes components for both STARmap and STARmap PLUS, but exact original-panel sequences and ambiguous barcode placement remain unresolved. See probes.md for schematic architecture and evidence limits.

## Encoding Strategy

STARmap uses a **5-nucleotide combinatorial nucleotide barcode** (A, T, C, G) embedded in the padlock probe, with a theoretical library size of 1,024 (4^5). The barcode is read by SEDAL two-base (dinucleotide) encoding: each SEDAL imaging cycle reads 2 bases simultaneously using 16 fluorescent decoding probes (one per dinucleotide combination) across 4 fluorescence channels. Inconsistent color sequences can be rejected by two-base decoding (S01 main text: "five-base barcode (library size of 1024)"; S03: 16 decoding probes, 4 channels, 60% formamide stripping). The reported per-read error rate is ~1.8%.

**Method version: V01.** S01 reports the encoding scheme and experimental error comparison. S03 cycle parameters are a later implementation; longer barcode options are described in probes.md.

## Amplification Strategy

Phi29 extends the primer probe around the circularized padlock template to produce tandem padlock-complement copies. Amine-modified nucleotides enable covalent hydrogel retention (reported: S01 Results, HTC). S03 / V02 Step 24 specifies 20 µM aminoallyl-dUTP and 250 µM dNTP mixture; Steps 25–26 specify 30 min at 4°C followed by 2 h at 30°C. The primer is not the RCA template. See protocols.md for the complete later-protocol conditions.

## Imaging / Readout Strategy

SEDAL ligates a fluorescent decoding probe to a phosphorylated reading probe on the DNA amplicon, then images and strips the ligation product. S01 Fig. 1E explicitly reports the reading-probe 5′ phosphate. S03 Tables 4–5 specify cycle-dependent reading oligos and 16 degenerate decoding oligos with four dyes. The reading position shifts by one base each round, interrogating overlapping dinucleotides (reported: S03 Overview). V01 uses Nissl segmentation after sequencing; S11 uses PI, while the later S03 workflow uses DAPI and, for STARmap PLUS, secondary-antibody imaging after sequencing (S01 160-gene Results; S03 Steps 52, 62–63).

## Decoding Strategy

Successive color observations are decoded through the two-base scheme and then matched to a barcode–gene codebook (reported: S01 SEDAL Results; S03 Overview and Fig. 4). The overlapping dinucleotides provide consistency checks; they are not three disjoint two-base chunks. Error rejection removes inconsistent reads but does not guarantee removal of every error. Table 5 gives decoding-oligo dinucleotide-to-dye groups; experiment-specific gene assignments remain unverified pending direct supplementary-table inspection.

## Major Versions or Variants

| Version | Reference | Distinctive features | Status |
|---|---|---|---|
| V01 | S01 (Wang 2018 *Science*) | Acrylic acid NHS ester HTC; 5-nt barcode on padlock probe; 160- to 1020-gene panels; 8-µm imaged volumes / 150-µm thick slices in mouse cortex; SEDAL sequencing-by-ligation | original |
| V02 | S11 (Zeng 2023 *Nat Neurosci*) | STARmap PLUS: protein comapping and improved probe design; 2,766-gene AD panel. S11 uses acrylic NHS; later S03 uses MA-NHS | variant |

V02 is retained as a version because S03 Overview explicitly describes STARmap PLUS as STARmap with protein detection and improved probe design. Its chemistry and workflow differences justify this routing choice. The S03 protocol document covers both V01 and V02 probe architecture and additionally covers RIBOmap and TEMPOmap, which are sibling methods on the same SNAIL-STARmap platform (different probe designs) and are not STARmap versions. See `metadata.yaml` for the version records.

## Strengths

- **High specificity by dual-probe requirement**: the primer + padlock pair must both hybridize to the same RNA before the padlock can be circularized (S01 main text).
- **Low reported per-read error rate**: ~1.8% for SEDAL vs. ~29.4% for cPAL methods (S01 main text), in the crowded four-gene comparison (S01 SEDAL Results), with two-base decoding used for error rejection.
- **Targeted gene panel scalability**: 160- to 1020-gene panels demonstrated (S01 abstract).
- **3D volumetric capability**: thick-tissue STARmap demonstrated at 150 µm section thickness with >30,000 cells per cubic-millimeter volume (S01 abstract, Fig. 3).
- **Covalent amplicon retention via HTC**: the amine-functionalized amplicons are covalently coupled to the hydrogel, preventing amplicon loss or migration during the multiple SEDAL cycles.
- **In-tissue cell segmentation via Nissl counterstain**: cell bodies are defined after SEDAL sequencing, enabling per-cell read assignment.

These are conservatively reported strengths from the S01 main text, not promotional interpretations.

## Limitations

- **Targeted, not untargeted**: STARmap requires a pre-designed gene panel (160–1020 genes in S01); it is not a transcriptome-wide method like untargeted FISSEQ. The 5-nt barcode library size is 1,024 (4^5), capping the multiplexed gene count (S01 main text: "five-base barcode (library size of 1024)").
- **No reverse transcription step**: STARmap operates on RNA directly, not on cDNA. This is a design choice, not a limitation per se, and does not by itself quantify sensitivity.
- **Multiple SEDAL cycles impose imaging time**: Six rounds are reported for the five-base 160-gene configuration. The separate thick-tissue experiment reads 28 genes in seven linear rounds, four genes per round (S01 Fig. 6); those are not seven rounds to decode the same five-base barcode. The 2018 paper reports 20% reads filtered in 4-gene experiments and 40% filtered in 1020-gene experiments (S01 main text).
- **Probe design software maintenance**: the later S03 design pipeline uses Picky 2.2, Dedupe/BBMap, NUPACK and BLAST; S05 is a related Bava-lab SNAIL utility whose equivalence to that pipeline and maintenance status are unverified (see `probes.md`).
- **Data availability for V01**: no GEO/SRA accession is reported in the S01 paper; the only public raw-image data deposit traceable is the ClusterMap RawImage directory (S08). A full benchmark-grade assessment of V01 is constrained by the uninspected supplementary resources and unverified reference context (see `benchmarks.md`).

These are conservatively reported limitations from the S01 main text and the supplementary materials; they are not promotional.

## Unresolved or Version-Specific Points

The S02 workbooks (Tables S1–S2) remain uninspected; per-gene sequences, barcode–gene assignments and original reference releases are therefore not confirmed here. The S02 PDF was read for the V01 probe-design rules (STARmap probe design, p. 3; Fig. S1A legend, p. 10; see probes.md); later S03 design rules are not evidence of the exact parameters used for the 2018 panels. Software maintenance and complete dataset inventories remain unverified. V01 uses the article identifier eaat5691, not a missing print page range. See probes.md and benchmarks.md for scoped gaps.
