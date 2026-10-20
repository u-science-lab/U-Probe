# Probe Architecture

## Overview

MiP-seq employs a dual-component probe system comprising: (1) a padlock probe containing dual barcodes, and (2) an initiator primer (IP) that enables target-dependent rolling-circle amplification. The system also utilizes sequencing anchor primers and fluorescent query probes for readout. Probe architecture is represented in 5'→3' orientation for each molecular species.

## Canonical Architecture

### Padlock Probe (DNA oligo, linear during synthesis, circularized upon ligation)

Schematic 5'→3':

```
5' — [Target-binding arm 1 (12-16 bp)] — [Anchor primer complementary region] — [Barcode B1] — [Barcode B2] — [Initiator primer complementary region] — [Target-binding arm 2 (12-16 bp)] — 3'
```

Note: Exact orientation and component order requires verification from source figures. The above represents the functional architecture as described.

### Initiator Primer (IP, DNA oligo)

Schematic 5'→3':

```
5' — [Target-binding region (12-16 bp)] — [Padlock probe complementary region] — 3'
```

The initiator primer binds both the target RNA and the padlock probe, forming a three-way junction that increases specificity.

### Anchor Primers (for sequencing-by-ligation)

Three anchor primers are used: anchor-1, anchor-2, and anchor-3.

These primers bind adjacent to the barcode region on the circularized padlock probe and provide binding sites for fluorescent query probes.

### Query Probes (fluorescent readout oligos)

Four query probes are used, each conjugated to a different fluorophore:
- seq-1: Alexa Fluor 488
- seq-2: Alexa Fluor 546
- seq-3: Alexa Fluor 594
- seq-4: Alexa Fluor 647

Query probes interrogate paired bases in the dual barcodes during sequencing-by-ligation.

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|-----------|--------------|----------|----------|------------|--------|----------|----------------|----------|-------------|-----------------|------------------|
| Target-binding arm 1 | Not explicitly named | Hybridize to target sequence | 5' end of padlock | 5'→3' (complimentary to target) | 12-16 bp | Not provided | variable | yes | Complementary to target region | reported | S01, Fig. 1B |
| Target-binding arm 2 | Not explicitly named | Hybridize to target sequence | 3' end of padlock | 3'→5' (complimentary to target) | 12-16 bp | Not provided | variable | yes | Complementary to target region | reported | S01, Fig. 1B |
| Barcode B1 | Barcode 1, B1 | Encode gene identity (first position) | Central region of padlock | Not reported | Not reported | Not provided | barcode_library | yes | Fixed per gene; from barcode library | reported | S01, Fig. 1B |
| Barcode B2 | Barcode 2, B2 | Encode gene identity (second position) | Central region of padlock | Not reported | Not reported | Not provided | barcode_library | yes | Fixed per gene; from barcode library | reported | S01, Fig. 1B |
| Anchor primer binding region | Not explicitly named | Provide binding site for anchor primers | Adjacent to barcodes | Not reported | Not reported | Not provided | fixed | yes | Constant across all probes | reported | S01, Fig. 1B |
| Initiator primer binding region | Not explicitly named | Complementary to initiator primer | Central region of padlock | Not reported | Not reported | Not provided | fixed | conditional | Required only when initiator primer is used | reported | S01, Fig. 1B |
| IP target-binding region | Initiator primer (IP), "target arm" | Hybridize to target RNA | 5' end of IP | 5'→3' (complimentary to target) | 12-16 bp | Not provided | variable | conditional | Complementary to target; 3' end of targeting region | reported | S01, Fig. 1B, Materials and Methods |
| IP padlock-binding region | Initiator primer (IP) | Complementary to padlock probe | 3' end of IP | 3'→5' (complimentary to padlock) | Not reported | Not provided | fixed | conditional | Constant across all probes | reported | S01, Fig. 1B |
| Anchor primer | Anchor primer | Sequencing primer for SBL | Adjacent to barcode on RCP | Not reported | Not reported | Not provided | fixed | yes | Three variants: anchor-1, anchor-2, anchor-3 | reported | S01, Materials and Methods |
| Query probe | Query probe, seq-1/2/3/4 | Fluorescent readout | Binds to anchor-barcode complex | Not reported | Not reported | Not provided | fixed | yes | Four probes with different fluorophores | reported | S01, Materials and Methods |
| Detection probe | cy3-p4 | RCA efficiency testing | Binds to RCP | Not reported | Not reported | Not provided | fixed | no | Complementary to padlock backbone | reported | S01, Materials and Methods |

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|-------|---------------------|---------------------------|-----------------|-----------------|
| Target type | mRNA (primary), genomic DNA, antibody-conjugated oligonucleotides | V01, all applications | reported | S01, Fig. 1A |
| Target region | Exon region; preferably near 5' end of transcript | V01, RNA detection | reported | S01, Materials and Methods |
| Probe length | Not explicitly reported for full probe; arms are 12-16 bp | V01 | reported | S01, Materials and Methods |
| GC constraints | Not reported | — | unavailable | — |
| Tm constraints | Not explicitly reported | — | unavailable | — |
| Specificity filtering | Not reported | — | unavailable | — |
| Repeat masking | Not reported | — | unavailable | — |
| Exon selection | Preferentially select exon close to 5' end | V01, RNA detection | reported | S01, Materials and Methods |
| Isoform handling | Not reported | — | unavailable | — |
| Sequence complexity | Not reported | — | unavailable | — |
| Spacing | 2 target sites per gene | V01, all experiments | reported | S01, Materials and Methods |
| Overlap | Not applicable | — | not_applicable | — |
| Number of probes per target | 2 target sites (probe pairs) per gene for standard detection; 2-8 probe pairs tested | V01, efficiency comparisons | reported | S01, Fig. S3 |

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|-------|---------------------|---------------------------|-----------------|-----------------|
| Barcode length | Not explicitly reported | V01 | unavailable | — |
| Number of barcodes | 2 per probe (dual barcode: B1 and B2) | V01 | reported | S01, Fig. 1B |
| Codebook | 10 dual-base combinations per round; 10N total codes with N rounds | V01 | reported | S01, Fig. 1C |
| Hamming distance | Not reported | — | unavailable | — |
| Error correction | Pairwise sequencing provides base-pairing validation; dual-channel intensity correlation | V01 | inferred | S01, Fig. S2 |
| Round assignment | Genes assigned to barcode combinations across multiple rounds | V01 | reported | S01, Fig. 1B-VI |
| Channel assignment | Paired channels: AF488+AF546, AF488+AF594, AF488+AF647, AF546+AF594, AF546+AF647, AF594+AF647 | V01 | reported | S01, Fig. 1C |
| Fixed vs experiment-specific barcode sets | Experiment-specific assignment of genes to barcode combinations | V01 | reported | S01 |

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|-------|---------------------|---------------------------|-----------------|-----------------|
| Readout sequence | Dual barcode region interrogated by pairwise sequencing | V01 | reported | S01, Fig. 1B |
| Readout probe | Fluorescent query probes (seq-1 through seq-4) | V01 | reported | S01, Materials and Methods |
| Fluorophore | Alexa Fluor 488, Alexa Fluor 546, Alexa Fluor 594, Alexa Fluor 647 | V01 | reported | S01, Materials and Methods |
| Imaging round | Multiple rounds (typically 2-3 for 100+ genes) | V01 | reported | S01, Fig. 3 |
| Sequencing primer | Anchor primers (anchor-1, anchor-2, anchor-3) | V01 | reported | S01, Materials and Methods |
| Decoding role | Query probes ligate to anchor primers adjacent to barcodes; fluorescence indicates base identity | V01 | reported | S01, Fig. 1B-V |

## Primer / Linker / Adapter Components

See Probe Components table above. Key sequences:

- **Anchor primers**: Fixed sequences; three variants used for sequencing-by-ligation
- **Query probes**: Fixed fluorophore conjugations; four variants (seq-1 through seq-4)
- **Initiator primer binding site**: Fixed region on padlock probe
- **Detection probe (cy3-p4)**: Fixed sequence for RCA efficiency testing; conjugated to Cy3 fluorophore

Exact sequences for anchor primers, query probes, and detection probes are not provided in the examined sources.

## Amplification-Related Components

| Component | Original term | Function | Position | Evidence status | Source / locator |
|-----------|--------------|----------|----------|-----------------|------------------|
| Rolling-circle product (RCP) | Rolling-circle product, RCP | Contains tandem repeats of padlock sequence for signal amplification | Product of RCA | reported | S01, Fig. 1B-IV |
| Phi29 polymerase binding site | Not explicitly named | Region enabling Phi29 extension | Padlock backbone | inferred | S01, Fig. 1B |
| RCA primer | Padlock probe circularized form | Serves as primer for Phi29 extension | 3' end of circularized padlock | inferred | S01, Fig. 1B |

## Chemical Modifications

| Modification | Original term | Position | Function | Evidence status | Source / locator |
|-------------|--------------|----------|----------|-----------------|------------------|
| 5' Phosphate | 5' phosphorylation | 5' end of padlock probe | Required for ligation; prepared by T4 PNK | reported | S01, Materials and Methods |
| NHS-PEG4-Azide | NHS-PEG4-Azide | Antibody conjugation | Crosslinker for attaching DBCO-oligo to antibody | reported | S01, Materials and Methods |
| DBCO-modified oligonucleotide | DBCO-oligo | Antibody conjugation | Azide-reactive group for click chemistry conjugation | reported | S01, Materials and Methods |

## Probe-Design Workflow

1. **Target selection**: Use self-developed probe searching program to screen target sequences.
2. **Region selection**: Preferentially select exon close to 5' end of transcript.
3. **Probe pair design**: Select 2 target sites (probe pairs) for each gene.
4. **Padlock arm design**: Both ends of padlock probe are complementary to target sequence (12-16 bp each).
5. **Barcode assignment**: Assign unique barcode combination (B1+B2) from barcode library to each gene.
6. **Initiator primer design**: 5' end of initiator primer is complementary to 3' end of targeting sequence (12-16 bp); 3' end is complementary to padlock probe.
7. **Quality control**: Probes synthesized and validated (sequences listed in Supplementary Table S1 per S01).

## Official Probe-Design Software

Not yet examined. The preprint mentions a "self-developed probe searching program" but specific software name, repository, and design parameters were not located in the examined sources.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---------------|------------------|----------------|---------|---------------------|-------|--------|--------------|--------------|--------------|-------------------|-----------------|------------------|
| Not reported | Not reported | Not confirmed | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | unavailable | — |

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|------------|-------------------|---------------|---------|--------------------------|------------------|
| D01 | S01 Supplementary Table S1 | V01 | Probe sequences (padlock probes and initiator primers) for multiple genes | Available in Supplementary Data of S01 (Nature Biomedical Engineering) | S01, Supplementary Data |
| D02 | Zenodo 10837475 | V01 | Source image data; not probe sequences | Available at https://zenodo.org/records/10837475 | S01, Data availability |

Note: Supplementary probe sequence tables were reported as available with the publication. The exact supplementary file number and format requires verification.

## Version Differences

No version differences were identified. Only V01 has been documented.

## Unresolved Information

- **Exact padlock probe sequence**: Complete nucleotide sequences not provided in main text; Supplementary Table S1 contains sequences but was not directly accessed.
- **Barcode library sequences**: Fixed barcode sequences not provided.
- **Anchor primer sequences**: Fixed sequences not provided.
- **Query probe sequences**: Fixed sequences not provided.
- **Initiator primer padlock-binding region length**: Specific length not reported.
- **Complete target sites**: Full list of target coordinates/sequences not provided.
- **Genome/annotation versions**: Reference genome build not explicitly stated.
- **Official design software**: Specific software name and parameters not located.
- **Probe arm orientation**: Exact 5'→3' orientation of padlock probe arms requires verification from figures.
- **Barcode length**: Exact length of B1 and B2 barcodes not reported.
- **Probe validation criteria**: Specific filtering parameters not provided.
