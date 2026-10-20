# Probe Architecture

## Overview

V01 has reverse-padlock DNA, linear priming DNA and encoding DNA sub-libraries. V02 changes target-binding lengths and coding/readout. The target RNA is not the RCA template: the circularized reverse-padlock DNA is. Reported: S01, Figure 1A and STAR Methods, Probe design.

## Canonical Architecture

Schematic sequences below are **not exact nucleotide strings**. S01, Probe design, explicitly gives V01 working-probe order 5′→3′:

```text
Reverse padlock before circularization (80 nt):
5′-phosphate — [common 20 nt] — [RNA-binding 30 nt] — [common 30 nt] — 3′
Linear priming probe used in RCA (60 nt):
5′ — [common 20 nt] — [RNA-binding 30 nt] — [padlock-binding primer 10 nt] — 3′
Ordered unamplified priming template (107 nt):
5′ — [above 60 nt] — [common manufacturing extension 47 nt] — 3′
Working encoding probe (150 nt):
5′ — [common 20] — [readout-binding 20] — [readout-binding 20]
   — [amplicon-binding 30] — [readout-binding 20] — [readout-binding 20] — [common 20] — 3′
```

The 47-nt extension contains a 27-nt reverse-complement T7-promoter region and a 20-nt PCR-priming region. The 30-nt encoding recognition sequence is the **same strand sequence as the padlock RNA-binding region** and recognizes its complementary RCA product, not RNA directly. The padlock RNA-binding segment is the reverse complement of target RNA. The padlock and priming RNA-binding sites are adjacent on RNA, without overlap, with no gap or a 1-nt gap. Reported: S01, Probe design, V01.

**Synthesis orientation matters:** S02 lists ordered unamplified templates. S01 explicitly states that the ordered encoding-template sequences are reverse-complemented because amplification reverse-complements them. Do not copy an S02 Encoding column as the final working encoding strand. The 107-nt priming template is also not the 60-nt oligo used during RCA. DNA is circular after padlock ligation; the circle has no free 5′/3′ ends. S01, Probe design and Probe synthesis.

V02 schematic, with uncertain unreported component ordering retained:

```text
Reverse padlock:
5′ — [common end] — [RNA-binding 40 nt: gRNA spacer + upstream U6-derived constant sequence]
   — [remaining common regions; see exact ordered template in S02] — 3′
Shared linear priming DNA (55 nt in S02; derived length):
5′ — [gRNA scaffold-binding 45 nt] — [padlock-binding 10 nt] — 3′
Encoding DNA:
5′ — [common/readout components; final order unresolved here]
   — [amplicon-binding 40 nt] — [remaining readout/common components] — 3′
```

V02 recognition roles are reported by S01, Probe design; the 55-nt shared primer and its 45+10 segmentation are derived/inferred from S02, Perturb-RAEFISH probe templates D3 and the S01 roles. Ordered V02 padlock templates are 100 nt and encoding templates 160 nt (derived character counts, S02 B3:C1760); they must not be conflated with V01 lengths. The V02 padlock has extra common sequence relative to simply replacing 30 nt with 40 nt in an 80-nt V01 oligo. Final strand/order of all V02 encoding components is unresolved; no exact working sequence is reconstructed.

## Probe Components

Sequence representations in the following table are schematic unless explicitly identified as exact. Original terms are from S01. Fixed means constant within the reported library, not across all future RAEFISH implementations.

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region_padlock | RNA-binding region | Recognize RNA | Internal to padlock | Reverse complement of RNA, oligo 5′→3′ | 30 nt | Variable region; S02 Padlock columns | variable | yes | Adjacent to priming target | reported | S01 Probe design; S02 human D/mouse D:E; V01 |
| target_binding_region_priming | RNA-binding region | Recognize adjacent RNA | Between common segments | DNA complementary to RNA | 30 nt | Variable region; S02 Priming columns | variable | yes | No overlap; 0/1-nt gap | reported | S01 Probe design; V01 |
| primer_rca | common region | Prime RCA by binding padlock | Priming oligo 3′ end | 5′→3′ synthesis from primer | 10 nt | Schematic; exact shared sequence shown below | fixed | yes | Must recognize padlock | reported | S01 Probe design; V01 |
| other_functional_sequence_padlock5 | common region | Invariant synthesis/ligation end | Padlock 5′ | 5′→3′ | 20 nt | Exact library sequence ACTCTCACTTGGGTTCCAAT | fixed | yes | Generated with phosphate-bearing RT primer | reported | S01 Probe design; S02 Library amplification oligos D4; V01/V02 |
| other_functional_sequence_padlock3 | common region | Invariant padlock end | Padlock 3′ | 5′→3′ | 30 nt | Schematic; not all roles assigned | fixed | yes | Splint joins padlock ends | reported | S01 Probe design; V01 |
| other_functional_sequence_encoding_target | sequence same as the RNA-binding region | Bind complementary RCA repeat | Central encoding segment | Working strand same as padlock binding region | 30 nt | Variable region | variable | yes | S02 ordered Encoding is reverse-complemented | reported | S01 Probe design; V01 |
| readout_sequence_1_to_4 | readout probe binding region | Present four positive code bits | Two each side of amplicon-binding region | Working encoding 5′→3′ | 20 nt each | Barcode-library sequences; exact templates in S02 | barcode_library | yes | Four sites implement weight 4 | reported | S01 Probe design; V01 |
| promoter_manufacturing | reverse complement of a T7 promoter | Probe pool manufacture | Within 47-nt priming-template extension | Source describes reverse complement | 27 nt | Schematic, not a new working-probe segment | fixed | conditional | Ordered template only | reported | S01 Probe design; V01 |
| primer_manufacturing | PCR priming region | Probe pool manufacture | Within extension | Source template orientation | 20 nt | See S02 amplification sheet | fixed | conditional | Not retained as a 47-nt extension on RCA primer | reported | S01 Probe design; V01 |
| target_binding_region_gRNA_padlock | RNA-binding regions | Recognize spacer and upstream constant RNA | Padlock internal region | RNA-complementary DNA | 40 nt | Exact ordered examples in S02; case preserved | variable | yes | 20-nt spacer + 20-nt U6-derived constant | reported | S01 Probe design; V02 |
| target_binding_region_gRNA_priming | RNA-binding region | Bind shared scaffold | Shared priming oligo | Complementary DNA | 45 nt | Exact D3 retained in S02 | fixed | yes | Shared across gRNA set | reported | S01 Probe design; S02 Perturb-RAEFISH templates D3; V02 |
| other_functional_sequence_gRNA_encoding | rolling circle amplicon binding regions | Bind gRNA-derived RCA repeat | Encoding recognition region | Final strand needs separate synthesis reconciliation | 40 nt | Variable region | variable | yes | Corresponds to padlock recognition sequence | reported | S01 Probe design; V02 |
| ligation_sequence_splint | splint oligo | Bring common padlock ends together | Separate DNA oligo | Exact order string 5′→3′; modified 3′ end | Source calls 16 nt; residue convention unresolved | Exact GAGAGTGGGTATA*T*C/3InvdT/ | fixed | yes | Two 6-nt binding regions; blocked against RCA priming | reported | S01 A549 protocol; V01/V02 |
| other_functional_sequence_toehold | toehold oligo | Remove splint | Separate modified oligo | Exact order string 5′→3′ | 15 printed bases (derived count) | Exact GA+TAT+ACC+CAC+TC+TC | fixed | yes | LNA notation + preserved | reported | S01 A549 protocol; V01/V02 |
| adapter_readout | readout probes with adapter sequences; adapter oligos | Connect encoding sites to common dyes | Separate oligos | Order strings 5′→3′; pairing assessed separately | Adapter 1 has 60 nt (derived); no universal length assumed | Exact example below; library in S02 | barcode_library | yes | Panel-specific readout identity | reported | S01 Probe design/readout methods; S02 readout sheet; V01 |
| fluorophore_common | Alexa Fluor 750; ATTO 647N-labeled common readout probes | Fluorescent detection | 5′ end of common oligos in order strings | Chemical group, not strand | Not applicable | Exact modification notation below | fixed | yes | Bind adapters in V01 | reported | S01 readout methods; S02 readout B2:B3; V01 |
| chemical_modification_padlock | phosphate-labeled primer | Enable padlock ligation | 5′ end | Chemical group | Not applicable | /5Phos/ | fixed | yes | Incorporated by manufacturing RT primer | reported | S01 Probe synthesis; S02 amplification D4/D13; V01/V02 |

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | Endogenous protein-coding RNA and lncRNA | V01 | reported | S01 Probe design |
| Target region | Adjacent 30-nt recognition sites | V01 | reported | S01 Probe design |
| Probe length | Padlock 80 nt; working priming 60 nt; ordered priming 107 nt; encoding 150 nt | V01 | reported | S01 Probe design |
| GC constraints | Numeric acceptance bounds not reported in examined design section | V01/V02 | unavailable | S01 Probe design; S02 workbook |
| Tm constraints | Numeric acceptance bounds not reported; Figure S1 compares Tm distributions rather than prescribing thresholds | V01/V02 | unavailable | S01 Probe design and Figure S1 legend |
| Specificity filtering | Exact alignment database, tool and thresholds unresolved | V01/V02 | unavailable | S01 Probe design; S03 archive inventory |
| Repeat masking | Not reported in examined design section | V01/V02 | unavailable | S01 Probe design |
| Exon selection | Not reported in examined design section | V01/V02 | unavailable | S01 Probe design |
| Isoform handling | No reproducible isoform-selection rule established | V01 | unavailable | S01 Probe design; S02 template sheets |
| Sequence complexity | No numerical filter established | V01/V02 | unavailable | S01 Probe design |
| Spacing | No gap or 1-nt gap on RNA | V01 | reported | S01 Probe design |
| Overlap | No overlapping recognition sites | V01 | reported | S01 Probe design |
| Number of probes per target | Human one set; mouse most two sets, with identical duplication when only one acceptable pair exists | V01 | reported | S01 Probe design |
| Reference annotation | GENCODE v46 human; vM35 mouse | V01 | reported | S01 Probe design |
| Excluded targets | 40 human and 26 mouse high-expression housekeeping genes excluded to reduce overcrowding | V01 | reported | S01 Probe design |
| gRNA recognition | Padlock 40 nt; shared primer recognition 45 nt | V02 | reported | S01 Probe design |

The human and mouse template sheets use numeric values under the header “Ensembl_GeneID” (e.g. 419 for DPM1, 1 for Gnai3). These are **not resolved stable Ensembl accessions**. Do not fabricate an ENSG/ENSMUSG prefix or use them as transcript coordinates (S02 A3:B3, reported table content).

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | 94 abstract binary bits; four physical 20-nt readout-binding sites | V01 | reported | S01 Probe design |
| Number of barcodes | 23,312 human/21,955 mouse target genes; decoding compares against 25,281 possible codes | V01; target vs decoder scope differ | reported | S01 Probe design and Image analysis |
| Codebook | MHD4_Code columns map row identities to 94-bit strings | V01 | reported | S02 human/mouse column C |
| Hamming distance | Minimum 4; Hamming weight 4 | V01 | reported | S01 Probe design |
| Error correction | Designed to tolerate single-bit errors; actual thresholds differ by sample | V01 | reported | S01 Results and Image analysis |
| Round assignment | 47 physical dual-color rounds; full ordered bit-to-round mapping requires explicit reconciliation | V01 | reported | S01 Sequential secondary probe hybridization and imaging; S02 readout sheet |
| Channel assignment | 647/638-nm and 750/749-nm channels | V01 | reported | S01 readout methods |
| Fixed vs experiment-specific barcode sets | Published assignments are panel-specific; no universal gene code asserted | V01 | reported | S02 human/mouse codebooks |
| Variant codebook | 14-bit, weight 4, distance 2; three gRNAs for one perturbation gene share a code | V02 | reported | S01 Probe design; S02 Perturb-RAEFISH templates E |

Distance-2 V02 codes must not be described as having V01's single-error-correction property. The 94-bit codeword is not a 94-nt contiguous molecular barcode.

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | Four 20-nt binding segments per working encoding probe | V01 | reported | S01 Probe design |
| Readout probe | 94 adapters plus two common dye-labeled oligos, rather than 94 individually dye-labeled oligos | V01 | reported | S01 Probe design/readout methods; S02 94-bit RAEFISH readout oligos |
| Fluorophore | ATTO 647N and Alexa Fluor 750; table order notation /5ATTO647NN/ and /5Alex750N/ | V01 | reported | S01 readout methods; S02 B2:B3 |
| Imaging round | 47 dual-color rounds | V01 | reported | S01 readout methods |
| Sequencing primer | Not applicable to sequential FISH imaging; manufacturing primers are separate | V01/V02 | reported | S01 workflow/Probe synthesis |
| Decoding role | Presence/absence pattern across bit measurements | V01 | reported | S01 Results, Figure 1 |
| Variant readout | 14 directly dye-labeled oligos, seven dual-color rounds | V02 | reported | S01 readout methods; S02 Perturb-RAEFISH readout oligos |

Exact published order strings (S02, “94-bit RAEFISH readout oligos”, B2:B6; V01, reported; DNA, 5′→3′ as synthesis order notation):

```text
ATTO647 dye labeled probe: /5ATTO647NN/AGGGACGGTTCCAATCGGATC
Alex750 dye labeled probe: /5Alex750N/GTCAAGGACGGCAGTGCAAA
647 blocking oligo: CGCAACGCTTGGGACGGTTCCAATCGGATC
750 blocking oligo: ACAACCGGCAGTCAAGGACGGCAGTGCAAA
Adapter 1: GACTTTGTGGGTACGTTGATGATCCGATTGGAACCGTCCCGATCCGATTGGAACCGTCCC
```

Adapter 1 is one exact library member, not a universal adapter. Its first 20 nt correspond to a library recognition sequence; the remaining sequence contains repeated segments partially complementary to the common dye oligo (inferred by sequence comparison of S02 B2 and B6). Do not infer one dye molecule per adapter from an abstract one-bit assignment.

## Primer / Linker / Adapter Components

Exact phosphate-bearing padlock RT primer: `/5Phos/ACTCTCACTTGGGTTCCAAT` (S02, “ Library amplification oligos”, D4 and D13; V01/V02; reported synthesis order string, DNA 5′→3′). V01 priming-library RT primer: `GTGGTAAAGCTCCGCGGCTT` (D7). V01 encoding-library RT primer: `TTGTCGCACGTTCGGTGTCG` (D10). These are manufacturing oligos, not RNA sequencing primers. The 10-nt RCA-priming terminus `ATTGGAACCC` is derived from S02 Priming columns, positions 51–60 (1-based), supported by S01's 60-nt working structure. No extra linker is invented.

## Amplification-Related Components

The reverse-padlock DNA is the circle template, and the adjacent linear probe supplies the RCA primer. Splint-assisted ligation joins the common ends away from the RNA-binding segment. Exact source strings from S01, A549 protocol/Key resources table, V01 (also reused in mouse/V02 protocols):

```text
Splint: GAGAGTGGGTATA*T*C/3InvdT/
Toehold: GA+TAT+ACC+CAC+TC+TC
```

These are published modified-oligo strings, not reconstructed sequences. S01 calls the splint “16-nt”, describes two 6-nt binding segments plus a 3-nt overhang, and prints 15 A/C/G/T letters followed by an inverted dT. Preserve all three descriptions: 15 printed standard bases + modified dT may explain 16 residues, but the counting convention is not explicitly resolved. Never silently change the sequence or label the entire oligo unmodified 16-nt DNA. The toehold removes splint; it is not an RCA primer.

## Chemical Modifications

Padlock library synthesis uses a 5′ phosphate-bearing RT primer (S01 Probe synthesis; S02 amplification sheet). S01 reports phosphorothioate linkage notation `*` and a 3′ inverted dT on the splint to prevent unwanted priming, and `+` for locked nucleic acids in the toehold. S02 explicitly places dyes at 5′ ends of the displayed readout-order strings. These are chemical groups, not nucleotide barcode positions. Reported, V01/V02 as scoped above.

## Probe-Design Workflow

Reported V01 inputs: human GENCODE v46 and mouse vM35 transcript sets. Select adjacent recognition sites without overlap/with 0–1-nt gap; assign panel-specific weight-4 HD4 codewords; construct padlock/priming/encoding templates; reverse-complement ordered encoding templates to account for manufacturing; order pools and manufacture working probes. Human has one set/target; mouse generally two with duplicates if only one acceptable pair. S01, Probe design. The exact site-selection/filtering algorithm, isoform rules and reference FASTA/GTF accession checksums were not established in examined sources. The GRCh38.p14 reference in S01 bulk RNA-seq analysis belongs to that comparator workflow; it is not proof of the reference used for every probe.

## Official Probe-Design Software

The S03 author archive contains four inspected GoodGeneList MAT files under Design240427/Design240516. They contain numeric ENSG/ENST/LongestENST fields, transcript Sequence strings and C/L binding sequences. Comparison by symbol plus both binding sequences matches 23,303 distinct published human sequence-pair records (derived: evidence/design_resource_inspection.json); this does not establish versioned stable transcript accessions, a reference build or a probe-selection algorithm. The human workbook has 23,312 non-control rows but 23,300 unique symbols; mouse has 21,955 rows but 21,954 unique symbols. Duplicate symbols must not be collapsed during extraction. The sampled mouse decoder RAEFISH/runM2_2_HPC.m sets RoundsOfHybs=94 and reshapes vectors with that size. This supports a 94-dimensional implementation for this archive example (inferred from code), while the publication's 96-dimensional wording remains an explicit source discrepancy.

No official probe-design program was established by the examined resources. This is unresolved (`null`), not a claim that none exists. The paper endorses **analysis** code at S03/S04. S03's archive directory contains Decoding and SingleCellAnalysis code; it does not show a probe-design module. ProbeDealer is a separate link in the laboratory navigation, not evidence that it generated these libraries.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RAEFISH analysis scripts | Yale DataAndCodes.zip | Author-linked analysis resource; not established as probe-design software | No release tag; observed Last-Modified 2025-05-17 | MATLAB (one script inspected) | Images/codebooks and processed data as indicated by file paths | Decoding and single-cell analyses | No probe-design algorithm identified | Sample script uses knnsearch/pdist2/bwconncomp/regionprops; complete dependency list unresolved | HTTP range access verified 2026-10-15 | No maintenance conclusion | inferred | S03 archive directory and sampled mouse runM2_2_HPC.m; V01/V02 |

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | Final article Table S1, publisher mmc1.xlsx | V01 human | 23,352 template/codebook rows, including 40 names ending NegativeControl; 23,312 other rows (derived counts) | Retrieved 2026-10-15 | S02 Human RAEFISH probe templates A:F |
| D02 | Same workbook | V01 mouse | 21,981 template/codebook rows, including 26 NegativeControl names; 21,955 other rows; two sets per row | Retrieved 2026-10-15 | S02 Mouse RAEFISH probe templates A:I |
| D03 | Same workbook | V02 | 1,758 padlock/encoding rows, shared priming oligo, 14-bit assignments and 14 labeled readouts; separate gRNA construction sheet | Retrieved 2026-10-15 | S02 Perturb-RAEFISH sheets |
| D04 | Mendeley 8kbv637pxh v2; author-linked ZIP | V01/V02; exact snapshot correspondence unresolved | Processed data, codebook file metadata, four human design-associated MAT files and sampled mouse decoder | Metadata, ZIP directory, four MAT files and sampled decoder inspected 2026-10-15 | S03/S04 |

The original workbook was retained as Table S1 (evidence/publisher_mmc1.xlsx); an inventory (evidence/workbook_inventory.json) records sheet headers and sample rows. Neither file is included in this repository (see README.md). Row counts are derived, not reported gene counts. “NegativeControl” names are preserved; no assumption is made that their sequences are the excluded high-expression gene's genuine targeting probes. Benchmark limitations are in benchmarks.md.

## Version Differences

V01 human/mouse: 30-nt recognition regions, 94-bit weight-4 HD4, 47 dual-color rounds. V02: 40-nt padlock recognition, 45-nt shared scaffold recognition, 40-nt encoding recognition, 14-bit weight-4 HD2, seven dual-color rounds. V02 templates are not assigned V01 final-strand geometry by assumption. S01 Probe design; S02 variant sheets. Final Table S1 consolidates resources that had different supplementary-table numbers in preprint S06.

## Unresolved Information

Exact site-selection algorithm, GC/Tm acceptance cutoffs, repeat/off-target handling, versioned transcript accessions/reference coordinates, numeric gene-ID accession reconstruction, final V02 encoding strand assembly, full bit/channel/round mapping verification, splint residue-count convention, and the endogenous decoding 94/96-dimensional wording remain open. S01/S02/S03 delimit the examined evidence. No complete reconstructed working-probe library or strict design ground truth is asserted.
