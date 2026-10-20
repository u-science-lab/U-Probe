# Research Manifest

Audit only; exclude from ordinary retrieval. Documents actual source inspection for the TDDN-FISH Research stage.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date | Inspected revision |
|---|---|---|---|---|---|---|---|---|
| S01 | Primary publication | Nature Communications article and PDF | Abstract, Introduction, Results (Establishment of TDDN-FISH method, Verification of universality and efficiency, Mapping subcellular atlas of influenza A virus RNA, Mapping cellular neighborhoods in mouse brain, Visualization of organization of cell types in mouse primary visual cortex), Discussion, STAR Methods (Animal samples, Tetrahedral DNA sequence design, Synthesis of tetrahedral DNAs, Assembly of TDDNs, Primer probe design and synthesis, Cell culture, TDDN-FISH procedure on cultured cells and tissue sample, Multi-round TDDN-FISH imaging, Imaging analysis), Figures 1-6 | Method overview, probe architecture (T0/T1/T2 monomer design with sticky ends, primary probe 35 nt target + 30 nt readout), target-binding design rules (OligoArray 2.0, Tm constraints, complexity rules), encoding strategy (F^N cyclic encoding), amplification strategy (TDDN hierarchical assembly), imaging strategy (confocal microscopy, 4 lasers, TCEP cleavage), decoding strategy (TrackMate, Cellpose, MATLAB scripts), experimental protocol parameters, demonstrated applications (cultured cells, tissue sections, influenza vRNA, mouse brain), strengths and limitations | Supplementary Data 1 (XLSX with all probe sequences) not directly examined; exact sequences not extracted | Exact primary probe sequences; exact TDDN component sequences; exact fluorophore-labeled detection oligo sequences; genome build used; transcriptome version; cell-type-specific marker gene selection criteria; imaging acquisition parameters (exposure times, laser powers) | 2026-10-15 | 2025-10-20 |
| S02 | Official repository (image processing code) | Zenodo | Repository description, abstract, authors, topics | Image processing pipeline components (MATLAB scripts, ImageJ/TrackMate, Cellpose), publication date | Code itself not downloaded or examined in detail; no version control history examined | Specific code structure, dependencies, installation requirements, code documentation completeness | 2026-10-15 | 2025-09-30 |

## Unresolved Research Questions

1. **Exact primary probe sequences**: Supplementary Data 1 contains all 65 nt primary probe sequences (35 nt target + 30 nt readout) for all targets used in the study (ACTB, POLR2A, HSP70, mTOR, B2M, MALAT1, NEAT1, miR-21, Satb2, Gad1, Slc17a6, Lmnb1, Slc17a7, PPIA, PA vRNA, PB1 vRNA, and others for 60-gene and 90-gene panels). These were not extracted in this research stage.

2. **Exact TDDN component sequences**: The 12 ssDNA strands (T0a-d, T1a-d, T2a-d) that assemble to form the three TDDN monomers are in Supplementary Data 1 but were not extracted.

3. **Exact fluorophore-labeled detection oligo sequences**: The Cy5, Cy3, and FAM-labeled detection oligos that bind to T2 sticky ends are in Supplementary Data 1 but were not extracted.

4. **Genome build and transcriptome version**: The paper does not explicitly report which human or mouse genome build or transcriptome annotation version was used for probe design. This affects benchmark assessment.

5. **Primary probe 5'→3' orientation**: The schematic in Figure 2c shows target-binding on one end and readout sequence on the other, but the exact 5'→3' orientation was not verified against actual sequences.

6. **Disulfide bond position**: The primary probe contains a disulfide bond for TCEP cleavage, but its exact position within the 65 nt probe (e.g., between target and readout regions, or at a specific terminus) is not explicitly reported in the examined sources.

7. **TDDN sticky end sequences**: The exact sequences of the sticky ends on T0, T1, and T2 monomers are not reported in the main text, only in Supplementary Data 1.

8. **Codebook documentation**: For the F^N encoding system, the specific barcode assignments for each target RNA are not documented in a centralized codebook in the main text.

9. **TDDN-FISH code repository structure**: The Zenodo repositories (S02) were identified but the actual code was not downloaded or examined in detail. The repositories provide MATLAB scripts, ImageJ/TrackMate integration, and Cellpose segmentation, but specific dependencies, installation requirements, and documentation completeness are unknown.

10. **Later publications**: Whether the TDDN-FISH method has been applied by the same group or independent groups in subsequent publications is unknown.

11. **Independent validation**: Whether TDDN-FISH has been independently validated or applied by research groups other than the original developers is unknown.

12. **Official probe-design software**: The paper reports using OligoArray 2.0 (a general oligonucleotide design tool) for target sequences and MERFISH-derived readout sequences, but no dedicated TDDN-FISH probe-design software was identified.

## Research Stage Summary

**Status**: Research complete (research stage only; no independent review performed)

**Coverage**:
- Primary publication (S01) extensively examined
- Official code repository (S02) identified
- Supplementary Data 1 identified but not extracted
- No additional primary sources identified

**Key limitations**:
- Exact probe sequences not extracted from Supplementary Data 1
- Genome build and transcriptome version not reported
- Probe architecture inferred from schematics and text descriptions, not from exact sequences
