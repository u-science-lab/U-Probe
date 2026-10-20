# Probe Architecture

## Overview

MiP-seq employs a dual-component probe system comprising: (1) a padlock probe containing dual barcodes, and (2) an initiator primer (IP) that enables target-dependent rolling-circle amplification. The system also utilizes sequencing anchor primers and fluorescent query probes for readout. Probe architecture is represented in 5'→3' orientation for each molecular species.

## Canonical Architecture

### Padlock Probe (DNA oligo, linear during synthesis, circularized upon ligation)

Schematic 5'→3':

```
5'P — [Target-binding arm 1 (13)] — CC — [Barcode B1 (3)] — TGCGTCTATTT — GT — [Barcode B2 (3)] — TAGTGGAGCCT — [Target-binding arm 2 (13)] — 3'
```

58 nt; lengths in nt. Example: Th-1-PLP-1 = GGTGGATTTTGGC-CC-GAG-TGCGTCTATTT-GT-CAT-TAGTGGAGCCT-CGGGTCTCTAAGT. All 434 58-nt padlocks in sheet 'Fig 3' follow this order [reported sequences, derived segmentation; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel); V01]. TGCGTCTATTT and TAGTGGAGCCT are the sequences of anchors Anchor-L-1 and Anchor-R-1 (sheet 'In Situ Sequencing probe'), so the anchors hybridize to the rolling-circle product next to B1 and B2; TAGTGGAGCCT is also the site of the initiator primer's 3′ tail. The Methods give 12–16 bp arms (S02, Methods, MiP-Seq probe design (3)); the published panel uses 13 + 13. Arm 1 (5′ end, 5′-phosphate) pairs with the 5′-side RNA segment and arm 2 with the segment immediately 3′ of it, so the ligation junction lies between the two arm sites [derived: Th-1 and Th-2 arms map contiguously to mouse Th RefSeq NM_009377.2, nt 359–384 and 478–503]. Other sheets of Supplementary Data 1 (e.g. 'Fig 1', 'Fig 2', 'Fig 4', 'Fig 5', 'Fig6') also use backbones read by fluorescent detection probes instead of pairwise sequencing.

### Initiator Primer (IP, DNA oligo)

Schematic 5'→3':

```
5' — [Target-binding region (12-16 bp)] — [short spacer] — AGGCTCCACTA — 3'
```

Example: Th-1-RCA-2 = GCCTCTGGGCAGG-C-AGGCTCCACTA (25 nt); the 432 primers of sheet 'Fig 3' are 25 nt (247) or 26 nt (185) and all end in AGGCTCCACTA [reported sequences, derived counts; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel); V01]. The 3′ tail AGGCTCCACTA is complementary to TAGTGGAGCCT on the padlock, so the primer's 3′ end sits on the circularized padlock and primes rolling-circle amplification ("RCA initiator primer probe for target-dependent RCA"; "the 3' end is complementary to the padlock probe": S02, Methods, MiP-Seq probe design, opening sentence and item 4). Its RNA-binding part binds 1 nt 3′ of the padlock arm sites in Th-1/Th-2 [derived, NM_009377.2]. The initiator primer binds both the target RNA and the padlock probe, forming a three-way junction that increases specificity.

### Anchor Primers (for sequencing-by-ligation)

The Methods name three anchor primers (anchor-1, anchor-2, anchor-3; S02, Methods, MiP-Seq probe design (6)); the published list has each as a left (B1) and right (B2) anchor, all 5′-phosphorylated [reported; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise'; V01]:

- Anchor-L-1 /5Phos/TGCGTCTATTT; Anchor-L-2 /5Phos/NTGCGTCTATT; Anchor-L-3 /5Phos/NNTGCGTCTAT
- Anchor-R-1 /5Phos/TAGTGGAGCCT; Anchor-R-2 /5Phos/NTAGTGGAGCC; Anchor-R-3 /5Phos/NNTAGTGGAGC

They carry the padlock's own sequence, so they hybridize to the rolling-circle product (the complement of the circle) next to B1 (L) and B2 (R). Each shift by N moves the read one barcode base further from the anchor: anchors -1, -2 and -3 read barcode bases 3, 2 and 1 [derived from the sequences]. A paired-end variant uses one middle anchor, Anchor-middle-1 /5Phos/GCGUCUAUUUAGUGGAGCCU (as published), with 5′-phosphorylated, 3′-dye queries 2merge2-R-F1–F4 (/5Phos/ANNNNN/3Alex488/, /5Phos/GNNNNN/3Alex546/, /5Phos/CNNNNN/3Alex594/, /5Phos/TNNNNN/3Alex647/) [reported; same sheet, section 'Paired-end'].

### Query Probes (fluorescent readout oligos)

Four query probes are used, each conjugated to a different fluorophore (seq-1 to seq-4 in S02, Methods, MiP-Seq probe design (6)); the published sequences are 5′-dye 6-mers whose 3′ base is the interrogated base [reported; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise'; V01]:
- 2merge1-L-F1 /5Alex488/NNNNNA: Alexa Fluor 488 = A
- 2merge1-L-F2 /5Alex546/NNNNNG: Alexa Fluor 546 = G
- 2merge1-L-F3 /5Alex594/NNNNNC: Alexa Fluor 594 = C
- 2merge1-L-F4 /5Alex647/NNNNNT: Alexa Fluor 647 = T

The query's 3′ end is ligated to the 5′-phosphate of an anchor, so in each round the same query set reads one B1 base (on an L anchor) and one B2 base (on an R anchor) of every RCP [derived from the sequences]. The same sheet lists a swapped set (section 'A-C base fluorophores swap Pairwise Sequencing': /5Alex594/NNNNNA, /5Alex488/NNNNNC).

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|-----------|--------------|----------|----------|------------|--------|----------|----------------|----------|-------------|-----------------|------------------|
| Target-binding arm 1 | Not explicitly named | Hybridize to target sequence | 5' end of padlock (5'-phosphate), nt 1-13 | Antisense (complementary to target); pairs with the 5'-side RNA segment | 12-16 bp (Methods); 13 nt in sheet 'Fig 3' | Per target; e.g. GGTGGATTTTGGC (Th-1-PLP-1) | variable | yes | Complementary to target region | reported (Methods, sequences); derived (position) | S02, Methods, MiP-Seq probe design (3); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| Target-binding arm 2 | Not explicitly named | Hybridize to target sequence | 3' end of padlock, nt 46-58 | Antisense (complementary to target); pairs with the RNA segment immediately 3' of arm 1's | 12-16 bp (Methods); 13 nt in sheet 'Fig 3' | Per target; e.g. CGGGTCTCTAAGT (Th-1-PLP-1) | variable | yes | Complementary to target region | reported (Methods, sequences); derived (position) | S02, Methods, MiP-Seq probe design (3); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| Barcode B1 | Barcode 1, B1 | Encode gene identity (first position) | nt 16-18, between CC and TGCGTCTATTT | Padlock sense | 3 nt | Per gene; e.g. GAG (Th) | barcode_library | yes | Fixed per gene | derived (from published padlock sequences) | S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| Barcode B2 | Barcode 2, B2 | Encode gene identity (second position) | nt 32-34, between GT and TAGTGGAGCCT | Padlock sense | 3 nt | Per gene; e.g. CAT (Th) | barcode_library | yes | Fixed per gene | derived (from published padlock sequences) | S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| Anchor primer binding region | Not explicitly named | Its sequence is the anchor sequence; anchors hybridize to the RCP copy of it | 3' of B1 (TGCGTCTATTT, nt 19-29) and 3' of B2 (TAGTGGAGCCT, nt 35-45) | Padlock sense | 11 nt each | TGCGTCTATTT; TAGTGGAGCCT | fixed | yes | Constant across the pairwise-sequencing padlocks | derived (from published sequences) | S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise' |
| Initiator primer binding region | Not explicitly named | Pairs with the initiator primer's 3' tail AGGCTCCACTA | nt 35-45, between B2 and arm 2 (same as the B2 anchor sequence) | Padlock sense | 11 nt | TAGTGGAGCCT | fixed | conditional | Required only when initiator primer is used | derived (from published sequences) | S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| IP target-binding region | Initiator primer (IP), "target arm" | Hybridize to target RNA | 5' end of IP | Antisense (complementary to target) | 12-16 bp | Per target; e.g. GCCTCTGGGCAGG (Th-1-RCA-2) | variable | conditional | Complementary to target; 3' end of targeting region | reported | S02, Methods, MiP-Seq probe design (4); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| IP padlock-binding region | Initiator primer (IP) | Pairs with TAGTGGAGCCT on the padlock; its 3' end primes RCA | 3' end of IP, after a short spacer (1 nt in Th-1-RCA-2) | Complementary to padlock | 11 nt | AGGCTCCACTA | fixed | conditional | Constant across all probes (all 432 primers in sheet 'Fig 3') | reported (Methods, sequences); derived (count) | S02, Methods, MiP-Seq probe design (4); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |
| Anchor primer | Anchor primer; Anchor-L-1/2/3, Anchor-R-1/2/3 | Sequencing primer for SBL | Hybridizes to the RCP next to B1 (L) or B2 (R) | Padlock sense; 5' phosphate | 11 nt (N-shifted forms include 1-2 N) | /5Phos/TGCGTCTATTT, /5Phos/NTGCGTCTATT, /5Phos/NNTGCGTCTAT; /5Phos/TAGTGGAGCCT, /5Phos/NTAGTGGAGCC, /5Phos/NNTAGTGGAGC | fixed | yes | Three positions (anchor-1, -2, -3) for each barcode | reported | S02, Methods, MiP-Seq probe design (6); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise' |
| Query probe | Query probe, seq-1/2/3/4; 2merge1-L-F1 to F4 | Fluorescent readout | Ligated to the 5' phosphate of an anchor on the RCP | 5'-dye; 3' base is read | 6 nt | /5Alex488/NNNNNA, /5Alex546/NNNNNG, /5Alex594/NNNNNC, /5Alex647/NNNNNT | fixed | yes | Four probes with different fluorophores, one per base | reported | S02, Methods, MiP-Seq probe design (6); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise' |
| Detection probe | cy3-p4; CY3-P4 | RCA efficiency testing | Binds to RCP | 5'-Cy3 | 20 nt | /5CY3/TGCGTCTATTTAGTGGAGCC | fixed | no | Complementary to the RCP | reported | S02, Methods, MiP-Seq probe design (5); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Fluorescent probe for detection' |

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
| Barcode length | 3 nt each (B1 and B2) | V01, pairwise-sequencing padlocks | derived (from published padlock sequences; three anchor positions read 3 bases) | S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel); S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise' |
| Number of barcodes | 2 per probe (dual barcode: B1 and B2) | V01 | reported | S01, Fig. 1B |
| Codebook | 10 dual-base combinations per round; 10^N total codes with N rounds | V01 | reported | S01, Fig. 1C |
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

- **Anchor primers**: /5Phos/TGCGTCTATTT, /5Phos/NTGCGTCTATT, /5Phos/NNTGCGTCTAT (L, for B1); /5Phos/TAGTGGAGCCT, /5Phos/NTAGTGGAGCC, /5Phos/NNTAGTGGAGC (R, for B2)
- **Query probes**: /5Alex488/NNNNNA, /5Alex546/NNNNNG, /5Alex594/NNNNNC, /5Alex647/NNNNNT
- **Initiator primer binding site**: TAGTGGAGCCT on the padlock; the primer's 3′ tail is AGGCTCCACTA
- **Detection probe (CY3-P4)**: /5CY3/TGCGTCTATTTAGTGGAGCC, for RCA efficiency testing

Sequences as published [reported; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'In Situ Sequencing probe', section 'Pairwise' and section 'Fluorescent probe for detection'; V01]. The same sheet lists detection probes for the other backbones (488-1, CY3-1, CY7-1, CY5-1).

## Amplification-Related Components

| Component | Original term | Function | Position | Evidence status | Source / locator |
|-----------|--------------|----------|----------|-----------------|------------------|
| Rolling-circle product (RCP) | Rolling-circle product, RCP | Contains tandem repeats of padlock sequence for signal amplification | Product of RCA | reported | S01, Fig. 1B-IV |
| RCA template | Circularized padlock probe | Template copied by Phi29 | Ligated padlock circle | reported | S02, Fig. 1B legend (III-IV) |
| RCA primer | RCA initiator primer (IP) | Its 3′ tail AGGCTCCACTA pairs with TAGTGGAGCCT on the circle; Phi29 extends the primer's 3′ end around the circle | 3′ end of the IP, on the padlock between B2 and arm 2 | reported (Methods); derived (sequence pairing) | S02, Methods, MiP-Seq probe design, opening sentence and item 4; S01 Supplementary Data 1 (probe sequences, 41551_2024_1205_MOESM5_ESM.xlsx), sheet 'Fig 3' (217-gene panel) |

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
7. **Quality control**: Probes synthesized and validated (sequences listed in S01 Supplementary Data 1; Table S1 in the preprint S02).

## Official Probe-Design Software

Not yet examined. The preprint mentions a "self-developed probe searching program" but specific software name, repository, and design parameters were not located in the examined sources.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---------------|------------------|----------------|---------|---------------------|-------|--------|--------------|--------------|--------------|-------------------|-----------------|------------------|
| Not reported | Not reported | Not confirmed | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | Not reported | unavailable | — |

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|------------|-------------------|---------------|---------|--------------------------|------------------|
| D01 | S01 Supplementary Data 1 (41551_2024_1205_MOESM5_ESM.xlsx) | V01 | Probe sequences: padlock probes, initiator primers, anchors, query and detection probes; sheets 'In Situ Sequencing probe', 'Fig 1', 'Fig 2', 'Fig 3', 'Fig 4', 'Fig 5', 'Fig6' | Downloaded from the publisher's Supplementary information, 2026-10-20 | S01, Supplementary information |
| D02 | Zenodo record 10837475 (version 4; concept DOI 10.5281/zenodo.8369209) | V01 | Source image data; not probe sequences | Available at https://zenodo.org/records/10837475 (checked 2026-10-20) | S01, Data availability |

Note: the publisher lists the probe file as "Supplementary Data" (probe sequences), followed by Supplementary Data 2 (efficiency comparison) and 3 (statistics, cluster gene lists); the preprint calls it Table S1.

## Version Differences

No version differences were identified. Only V01 has been documented.

## Unresolved Information

- **Complete target sites**: Full list of target coordinates/sequences not provided.
- **Genome/annotation versions**: Reference genome build not explicitly stated.
- **Official design software**: Specific software name and parameters not located.
- **Probe validation criteria**: Specific filtering parameters not provided.
