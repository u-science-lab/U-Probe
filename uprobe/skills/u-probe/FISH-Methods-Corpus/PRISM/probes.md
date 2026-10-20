# Probe Architecture

## Overview

PRISM combines RNA-targeting DNA padlocks, separate RCA primers and separate imaging probes. Four physical barcode segments encode an abstract four-channel intensity codeword. V01 denotes the 30/31-plex configuration; V02 denotes 64-plex. These internal labels do not establish distinct publication releases. [reported; S01 Fig. 1; S02]

## Canonical Architecture

Schematic, not a nucleotide sequence. The ordering below is **inferred** from exact backbone matches in S02 and the assembly convention in S07; 5′ identity is anchored by the padlock “5′ P” modification column. S02 does not explicitly label every displayed sequence 5′→3′. RNA pairing is antiparallel; original-reference alignments were not reconstructed.

```text
V01 listed mouse padlock, synthesized DNA strand (inferred 5′→3′):
5′-P — [20-nt terminal target arm] — [TC; role unresolved] —
[20-nt G] — [20-nt R] — [20-nt B] — [20-nt Y] — [20-nt terminal target arm] — 3′

V02 listed mouse-embryo padlock (inferred 5′→3′):
5′-P — [20-nt terminal target arm] — [20-nt G] — [20-nt R] —
[20-nt B] — [20-nt Y] — [20-nt terminal target arm] — 3′

Imaging probe (separate DNA molecule):
5′ — [20-nt sequence identical to a padlock barcode segment] — 3′
Label position varies; see component/readout tables.

RCA primer (separate DNA molecule):
5′ — [20-nt primer; two asterisk-marked linkages near 3′ end] — 3′
```

Ligation closes the padlock ends to form a circular DNA template; the circle has no free 5′/3′ termini. RCA produces repeated complementary copies that imaging probes bind. [reported; S01 Fig. 1b–c; S05 thin-tissue methods]

Derived length check: every complete listed backbone was located as an exact substring of its panel padlocks. V01 mouse oligos = 20+82+20 = 122 nt; V02 = 20+80+20 = 120 nt. HCC HBX (S02 HCC worksheet row 89) has 14+82+14 = 110 nt; its residual flanks are inferred as arms, not a reported universal 14-nt design. Do not generalize the preprint 40-nt target rule to HBX. All other 90 HCC records are 122 nt. [derived; S02; evidence/sequence_inventory.json]

The G-R-B-Y order is molecular order; the workbook four-digit codeword is not that order. Exact sequences remain in the original workbook, including modifications and leading non-nucleotide characters.

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| padlock_arm_5prime | binding target; padlock probe | RNA recognition; terminal ligation partner | 5-prime end | 5′→3′ oligo; antiparallel RNA pairing | 20 nt in mouse records; HBX exception 14 nt | Variable; exact strings in S02 | variable | yes | 40-nt target split into halves in preprint; final segmentation derived | inferred | S05, Probe and barcode design (V01 precursor); S02 padlock sheets (V01/V02) |
| padlock_arm_3prime | binding target; padlock probe | RNA recognition; terminal ligation partner | 3-prime end | 5′→3′ oligo; antiparallel RNA pairing | 20 nt in mouse records; HBX exception 14 nt | Variable; exact strings in S02 | variable | yes | Same scope as 5-prime arm | inferred | S05, Probe and barcode design; S02 padlock sheets (V01/V02) |
| other_functional_sequence_leader | No explicit term in source | Unresolved | Before G segment in listed V01 backbone | Listed-strand orientation | 2 nt | TC; exact shared prefix of listed V01 backbones | version_specific | unclear | No functional rule established | derived | S02, Imaging probe（31 panel）, I3:I33 (V01) |
| barcode_segment_g | Segment G | Channel-4 identity/intensity selection | First of four listed core segments | Padlock-sense sequence; imaging probe has same letters | 20 nt | Library sequences; see exact-sequence table | barcode_library | yes | V01 0G/1G; V02 0G/1G/2G | inferred | S02, imaging sheets (V01/V02); S01 Fig. 1c |
| barcode_segment_r | Segment R | Cy5 level selection | After G | As above | 20 nt | Library; not one universal sequence | barcode_library | yes | 0–4 V01; 0–5 V02 | inferred | S02, imaging sheets (V01/V02) |
| barcode_segment_b | Segment B | AF488 level selection | After R | As above | 20 nt | Library; not one universal sequence | barcode_library | yes | 0–4 V01; 0–5 V02 | inferred | S02, imaging sheets (V01/V02) |
| barcode_segment_y | Segment Y | Texas Red level selection | After B | As above | 20 nt | Library; not one universal sequence | barcode_library | yes | 0–4 V01; 0–5 V02 | inferred | S02, imaging sheets (V01/V02) |
| readout_probe | Decode probe; imaging probe | Competitive staining of RCA product | Separate oligo | Binds complementary RCA copy, not same-sense padlock | 20 nt | Exact strings in S02, column C | barcode_library | yes | Labeled/unlabeled pairs determine intensity | reported | S01 Fig. 1c; S02 imaging sheets (V01/V02) |
| amplification_primer | RCA primer | Initiate RCA | Separate oligo; complementary to R segment | 5′→3′ primer orientation supported by terminal modification notation | 20 nt excluding asterisks | Exact strings below / in S02 | barcode_library | yes | Five V01; twelve V02 entries | inferred | S02, RCA-primer sheets (V01/V02) |
| chemical_modification_padlock | Modification: 5′ P | 5-prime phosphorylation | 5-prime padlock terminus | 5-prime | Not applicable | Non-nucleotide group | fixed | yes | Recorded for listed padlocks | reported | S02 padlock sheets, Modification column (V01/V02) |
| chemical_modification_primer | phosphorothioate modification | Modified terminal primer linkages | Asterisk-marked positions near 3-prime end | Preserve source notation | Two asterisk-marked linkages per listed primer | See exact primer notation | fixed | yes | No substitution of asterisks by nucleotide letters | inferred | S02 RCA-primer sheets, sequence column and footnote (V01/V02) |
| fluorophore | Labeled | Fluorescence readout | Cy5/TxRed/AF488 at 3′; Cy3 at 5′; V02 2G has both ends | Separate from nucleotide order | Not applicable | Cy5; Texas Red; AlexaFluor 488; Cy3 | version_specific | conditional | Only labeled member carries dye; competitors unlabeled | reported | S02 imaging sheets, D3:D19 V01; D39:D41 V02 |

Component mappings and positional segmentation are inferred where the workbook supplies only full oligos; length arithmetic is derived. The original sequence strings, dye labels, phosphorylation and phosphorothioate legend are reported.

| Instance | Representation | Published sequence | Modification / qualification | Source / locator |
|---|---|---|---|---|
| V01 4R imaging probe | DNA; exact published instance | CGCAGCAGATAAATGAGCCA | 3′ Cy5; unlabeled counterpart also used | S02, Imaging probe（31 panel）, C3:D3 |
| V01 0G imaging probe | DNA; exact published instance | ACACGACGCTCTTCCGATCT | 5′ Cy3 label column applies to G series; theoretical 0G is unlabeled in pooling | S02, same sheet, C18:D19 and A39:C40 |
| V01 1G imaging probe | DNA; exact published instance | CTGCTGAACCGCTCTTCCGA | 5′ Cy3 | S02, same sheet, C19 and merged D18:D19 |
| V01 RCA primer 1 | DNA; exact published notation | TGGCTCATTTATCTGCTG*C*G | Asterisks are phosphorothioate linkages | S02, RCA Primer (31 panel), B2 and footnote |

Example strings are not universal constants. No reverse complements are presented as exact published sequences. Primer 1 is reverse-complementary to 4R after removing modification symbols [derived pairing check], consistent with RCA on the padlock template.

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | mRNA; direct RNA-templated padlock ligation | V01 precursor | reported | S05, Probe and barcode design / thin tissue |
| Target region | 40-nt window; final HBX has a shorter residual target region | V01 precursor; final V01 HBX exception | reported | S05, Probe and barcode design; see derived HBX analysis below |
| Probe length | Mouse V01: 122 nt; mouse V02: 120 nt; HCC: 90 × 122 nt and one 110-nt HBX record | Final listed oligos only | derived | S02 padlock sheets; evidence/sequence_inventory.json |
| GC constraints | No explicit historical GC-percent threshold in inspected preprint | V01 precursor | unavailable | S05, Probe and barcode design |
| Tm constraints | 65–75 °C; calculation model not specified | V01 precursor | reported | S05, Probe and barcode design |
| Specificity filtering | No BLAST matches to other genes in NCBI database; database release/thresholds unresolved | V01 precursor | reported | S05, Probe and barcode design |
| Repeat masking | Not reported in inspected design paragraph | V01 precursor | unavailable | S05, Probe and barcode design |
| Exon selection | Not reported as an explicit rule | V01 precursor | unavailable | S05, Probe and barcode design |
| Isoform handling | Cover all transcript variants if possible | V01 precursor | reported | S05, Probe and barcode design |
| Sequence complexity | Absence of 5 consecutive G; no significant DNA secondary structure (ΔG > −10; units not specified there) | V01 precursor | reported | S05, Probe and barcode design |
| Spacing | Not reported in historical paragraph; current software has spacing selection | Historical vs current software separated | unavailable | S05, Probe and barcode design; S07, README Strategies |
| Overlap | No explicit historical overlap threshold recovered | V01 precursor | unavailable | S05, Probe and barcode design |
| Number of probes per target | Preprint: 2 mouse targets, 3 HCC targets. Final lists: 1 per mouse gene; HCC generally 3, HBX 1 | Source-scoped discrepancy; do not merge | reported | S05, Probe and barcode design; S02 padlock sheets (counts derived) |

Current software defaults are not historical PRISM requirements: S07 configs/config_single_sequence.yaml sets a 40-nt search length, min_tm 50.0, max_tm 65.0, max_tm_diff 10.0, min_g_content 0.3, max_g_content 0.7 and final_probes_per_gene 3. The parameter is named G content, not GC content. config_template.yaml instead has max_tm 70.0. These documented configuration differences are preserved; no single default set is assigned to the paper. [reported; S07, current d7326be]

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | Four 20-nt segments (80 nt core); listed V01 backbones contain an additional 2-nt prefix | V01/V02 | reported | S05, Probe and barcode design; prefix derived from S02 |
| Number of barcodes | 30 or 31 including Ch4-only; expanded to 64 | V01/V02 | reported | S01 Fig. 1e / Extended Data Fig. 7 |
| Codebook | V01 H-column four-symbol codewords mapped to I-column physical backbones; gene assignment by exact backbone matching | V01 | derived | S02, Imaging probe（31 panel）, H3:I33; padlock sheets |
| Hamming distance | Not reported as a design criterion; geometric separation is used | V01/V02 | unavailable | S01 Fig. 1e; S04 Supplementary Fig. 2 |
| Error correction | Radius-vector distinguishability and cluster assignment; no independently established minimum-distance error-correcting code | V01/V02 | reported | S01 Fig. 1e–f |
| Round assignment | One pooled staining/imaging round | V01/V02 | reported | S01 Fig. 1c–d |
| Channel assignment | Workbook digits follow Y, B, R, G (TxRed, AF488, Cy5, Cy3), inferred from pure-channel records; code CSV column order differs | V01; inspect V02 separately | inferred | S02 H3:H7 + segment records; S06 configs/codebooks/prism30.csv |
| Fixed vs experiment-specific barcode sets | Published libraries are finite; genes and padlocks assigned within sample panels | V01/V02 | reported | S02 imaging and padlock sheets |

V01: first three intensity symbols sum to 4, producing 15 combinations, then Ch4 0/1 gives 30 plus optional 0001 [reported; S04 Fig. 2; S02 H3:H33]. V02: first three levels run 0–5 with sum 5, expanded with three fourth-channel levels and a fourth-channel-only code [S01 Extended Data Fig. 7; S02]. The V02 workbook labels high-G codes with a final digit **5** (rows 46–67), although S01 calls the high level **2**; preserve these source representations rather than rewrite codes.

V02 row 71 explains separate R/B/Y sequence families for high-Cy3 barcodes 43–63 because a reported FRET-like effect shifts color-space positions. This is a library/readout change, not simply another abstract intensity value. [reported; S02, Imaging probe (64 panel), rows 3–41, 71]

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | Four selected padlock-sense library segments; RCA copies are complementary | V01/V02 | reported | S01 Fig. 1a–c |
| Readout probe | 20-nt imaging probes; 17 V01 identities; 38 V02 identities including two cells with leading U+FEFF | V01/V02 | derived | S02 imaging sheets C3:C19; C3:C41 (non-empty sequence entries) |
| Fluorophore | 3′ Cy5, 3′ Texas Red, 3′ AlexaFluor 488; 5′ Cy3; V02 2G: 5′ and 3′ Cy3 | V01/V02 | reported | S02 imaging sheets label columns |
| Imaging round | 1; validation experiments may strip/re-stain to check decoding | V01/V02 | reported | S01 Fig. 1c; Extended Data Fig. 5/9/10 |
| Sequencing primer | Not applicable to fluorescence intensity readout; RCA primers are not sequencing primers | V01/V02 | reported | S01 Fig. 1b–c |
| Decoding role | Competitive occupancy sets color-intensity levels, classified in color space | V01/V02 | reported | S01 Fig. 1c,f |

The final 31-panel workbook records theoretical 3R labeled:unlabeled = 3:1 but experimental pooling = 9.6 µL:2.4 µL. These are different scoped records; the theoretical quarter-level model must not replace the empirical recipe. [reported; S02 same sheet rows 25 and 53]

## Primer / Linker / Adapter Components

The four barcode segments also provide binding features for imaging; R segments pair with separate RCA primers. Do not invent an additional universal adapter or linker. The two-nucleotide V01 prefix has unresolved function. [inferred/derived; S02 imaging/RCA-primer sheets; V01/V02]

S02 RCA Primer (31 panel) lists five 20-nt primers and RCA primer (64) lists twelve; all have two asterisks defined as phosphorothioate modifications. Each is prepared at 10 µM before pooling according to the footnote. Exact preparation volumes and final pooled per-primer concentrations are not resolved by that statement alone. [reported/derived; S02 primer sheets]

## Amplification-Related Components

The preprint uses Phi29 RCA and aminoallyl-dUTP, followed by BS(PEG)9 fixation [reported; S05 thin/thick tissue; V01 precursor]. S02 provides final-library RCA primers [V01/V02]. Selective amplification mixes active and “irreplicable” padlocks to reduce abundant-target signals [reported; S05 Probe and barcode design; S03 dilution table]. Exact inhibitory-oligo construction is not resolved here; no blocking modification is invented.

## Chemical Modifications

See distinct component rows: padlock 5′ P; asterisk-defined primer phosphorothioates; imaging-probe dye positions, including dual-end Cy3 on V02 2G. Aminoallyl-dUTP belongs to amplification chemistry rather than the synthesized padlock sequence. [reported; S02 label/modification columns and primer footnotes; S05 RCA recipe]

## Probe-Design Workflow

Historical workflow: select transcript windows → assess Tm, G runs, specificity and secondary structure → split into terminal arms → assign four-segment barcode → order padlock/imaging probes → pool with selective dilution where needed. Reference release, exact selection algorithm and original design software revision were not recovered. [reported; S05 Probe and barcode design; V01 precursor]

Current software: gene list → NCBI/Ensembl or supplied transcript/isoform data → candidate search → thermal filtering → BLAST specificity → scoring → regional ranking/spacing → selected sites → arm_5prime + backbone + arm_3prime assembly. The source reverses/complements target windows before arm assignment. This supports interpretation of arm orientation, but does not prove the original experiment ran this revision. [reported/inferred; S07 README, search_strategies.py SingleSequenceStrategy, probe_assembly.py _assemble_probe_sequence; current software]

## Official Probe-Design Software

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| probe-designer | https://github.com/HuangLab-PKU/probe-designer | Author-linked by S06 | d7326be; pyproject 0.2.0, CHANGELOG includes 0.2.1 | Python ≥3.10 | Gene list; mRNA/isoforms; configurable reference; backbone/assignment tables | Binding-site JSON; scored/selected sites; assembled XLSX/FASTA/JSON | Search, thermal filter, BLAST, scoring, regional selection, backbone assembly | Biopython/pandas/numpy/PyYAML/Typer; optional ViennaRNA; reference and BLAST setup | Public source downloaded 2026-10-15; not executed | 2026 changelog present; operational support not assessed | reported | S07 README; pyproject.toml; CHANGELOG.md; probe_assembly.py (current software) |
| PRISM-Code | https://github.com/HuangLab-PKU/PRISM-Code | Explicit S01 Code availability link | ab7b72d | Python | Stitched 2D/3D images and configuration | Spots, gene calls, segmentation outputs | Post-stitching processing; not historical probe-design code | See pyproject and docs; upstream spatial_img_core described as not public | Public source downloaded 2026-10-15 | No maintenance guarantee inferred | reported | S06 README; configs/codebooks (current software) |

S07 README/pyproject identify 0.2.0, while CHANGELOG contains 0.2.1 dated 2026-04-21. Use commit d7326be to identify the inspected code. Optional iLock code in probe_assembly.py is not evidence that PRISM uses iLock. PRISM-Analysis is downstream analysis, not probe-design software [S08 README].

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | 2025 Supplementary Table 2; SHA-256 in sources.md | V01 | 30 brain + 30 embryo + 91 HCC padlocks; 31 backbones; 17 imaging sequences; 5 RCA primers | Downloaded 2026-10-15 | S02, first five sheets |
| D02 | Same workbook | V02 | 64 embryo padlocks; 64 backbones; 38 imaging identities; 12 primers; pooling recipes | Downloaded 2026-10-15 | S02, 64-panel sheets |
| D03 | Zenodo 12755414 and paper-linked raw-image records | V01/V02; run-specific | Image/spot/analysis resources; not a probe-design truth set | 12755414 metadata/readme inspected 2026-10-15 | S09; S01 Data availability |

The workbook also supplies smFISH validation probes; these are comparison reagents, not additional PRISM target designs. Worksheet names, exact coordinates, merged dye-label cells and raw strings are preserved in evidence/supp4.xlsx. No design benchmark is Grade A solely because this table exists.

## Version Differences

V01 30/31 and V02 64 differ in barcode libraries, primer set, empirical pooling, and the V01 82-nt versus V02 80-nt listed backbones. These lengths are derived, not a reported universal synthesis-length rule. V02 uses separate high-green readout families. Preprint-to-final changes in mouse padlock number remain unresolved. [S01 Extended Data Fig. 7; S02; S05]

## Unresolved Information

Original transcript accessions/versions, genome/annotation release, coordinate convention, exact historical software/settings, unamplifiable competitor construction and the final-paper full Methods remain unresolved. The G-R-B-Y schematic is sequence-based interpretation; direct original-target alignments were not performed. V02 1G/2G cells have leading U+FEFF: preserve raw data and document normalization before computational use. The 2024 preprint and final data must remain separately scoped.
