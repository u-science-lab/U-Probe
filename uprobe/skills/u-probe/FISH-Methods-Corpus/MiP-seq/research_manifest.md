# Research Manifest

Audit only; exclude from ordinary retrieval. This document records what was actually inspected during MiP-seq corpus research.

## Inspected Sources

| Source ID | Source type | Location | Files/tables/figures examined | Information extracted | Limitations | Unresolved questions | Inspection date |
|-----------|------------|----------|-----------------------------|---------------------|-------------|---------------------|----------------|
| S01 | Peer-reviewed article | Nature Biomedical Engineering website | Abstract, main text (partial), Figure 1-6, Extended Data Figs. 1-10, Methods section, Data availability | Official method name, full description, dual-barcode system, pairwise sequencing workflow, probe design rules, protocol parameters, species tested, multi-omics applications, software reference (InnerEye), dataset references | PMID not yet available; supplementary tables not directly accessed; exact probe sequences not verified | Complete barcode library sequences, anchor primer sequences, official design software name and parameters, genome build versions | 2026-10-15 |
| S02 | Preprint | bioRxiv website | Abstract, main text, supplementary figures S1-S19, Materials and Methods | Detailed probe design protocol, arm length specifications, initiator primer design, sequencing workflow, efficiency comparisons | Preprint status; publication may supersede | Whether preprint differences affect technical details | 2026-10-15 |
| S03 | Code repository | GitHub | README, installation instructions | InnerEye package description, conda installation workflow, SimpleElastix dependency | Repository may contain more detailed documentation not inspected | Complete API documentation, example usage, test data | 2026-10-15 |
| S04 | Dataset | Zenodo | Record page, file listing | Dataset availability, file sizes, version information | Large dataset (15.2 GB); actual file contents not inspected | Codebook linking barcodes to genes, decoded data format | 2026-10-15 |
| S05 | Dataset | Zenodo | Record page | 100 genes HDB dataset availability, file format (HDF5) | Dataset contents not directly inspected | Complete gene list, codebook mapping | 2026-10-15 |
| S06 | Dataset | GEO | Record page | RNA-seq data availability for validation | Dataset not directly accessed | Specific SNP positions, expression values | 2026-10-15 |
| S07 | Blog | RNA-Seq Blog | Summary article | Discovery purposes only | Non-peer-reviewed summary | Not applicable | 2026-10-15 |

## Unresolved Research Questions

### Probe Architecture

- **Complete probe sequences**: The exact nucleotide sequences of padlock probes, initiator primers, anchor primers, and query probes are contained in Supplementary Table S1 but were not directly accessed during this research phase.
- **Barcode library sequences**: The fixed dual-barcode library sequences used for encoding are not reported in the main text or supplementary figures.
- **Barcode length**: The exact length of barcode regions B1 and B2 is not explicitly stated.
- **Padlock orientation**: The 5'→3' orientation of the padlock probe arms relative to the barcode regions requires verification from supplementary figures.
- **Complete design parameters**: Full probe design filtering criteria (GC content, Tm, specificity, off-target alignment) are not provided.

### Reference Information

- **Genome build**: The reference genome version used for probe design is not explicitly stated in the examined sources.
- **Transcriptome annotation**: The annotation version for gene target selection is not reported.
- **Species-specific variations**: Whether probe design rules differ between mouse, rat, zebrafish, and human targets is not documented.

### Software and Analysis

- **Official design software**: A "self-developed probe searching program" is mentioned but the specific software name, version, and parameters are not provided.
- **Codebook format**: The exact format and contents of the gene codebook (mapping barcodes to genes) is not described.
- **Decoding algorithm**: Specific parameters for barcode matching and quality filtering are not detailed.

### Experimental Parameters

- **Complete reagent specifications**: Some buffer pH values, detailed washing volumes, and optimization parameters are not explicitly stated.
- **Probe validation criteria**: Specific thresholds for probe acceptance/rejection during design are not provided.
- **Efficiency benchmarks**: Direct comparison with other methods is limited to HCR 3.0 in the preprint.

### Publication Status

- **PMID**: The Nature Biomedical Engineering publication does not yet have a PubMed ID at the time of this inspection.
- **Supplementary file numbering**: Exact correspondence between supplementary tables in the preprint (Tables 1-4) and published paper requires verification.
