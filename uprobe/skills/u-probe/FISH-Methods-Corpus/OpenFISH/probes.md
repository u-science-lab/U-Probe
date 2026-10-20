# Probe Architecture

## Overview

Highest-priority record; architecture is resolved at inspected-code level, with experimental chemistry and complete synthesized probe lists still unresolved. All implementation claims are scoped to S04 commit be490eb19503c742230665da1e665c6d6c3d97c8 and associated V01; they are not automatically settings used for all figures in S01.

## Canonical Architecture

Implementation schematic, not an exact published oligo specification. Define B as the 30-nt bridge library entry and r=reverse_complement(B). Python slices below are zero-based and end-exclusive. A and B in “probe A/B” are molecular names; the symbol B below denotes the bridge sequence only.

```text
Probe A: 5′ — [bridge[0:14], 14 nt] — [ta, 2 nt] — [target-binding P1, 30 nt] — 3′
Probe B: 5′ — [target-binding P2, 30 nt] — [ta, 2 nt] — [bridge[-14:], 14 nt] — 3′
Padlock precursor:
5′ — [r[-7:], 7 nt] — [ROx, 20 nt] — [connector, 2 nt] — [ROy, 20 nt] — [r[0:23], 23 nt] — 3′
```

Evidence: inferred from S04 utils.py findAllCandidates/get_full_probes1/2 and 03_Padlock_join cell 3 PadJoin.get_whole_pad/TEST_ORDER; V01 implementation snapshot. ROx/ROy may swap; connector is selected from ta/at/tt/aa. Derived precursor length: 7+20+2+20+23=72 nt; A/B each 14+2+30=46 nt. These are code-defined molecule lengths, not verified synthesis specifications.

In minus/mRNA mode, P1 is reverse complement of the first target region and P2 of the downstream region after a 2-nt gap. The padlock terminal segments are complementary to bridge-derived regions, not direct RNA recognition arms. Circle formation connects precursor ends; actual ligation chemistry, phosphate status, initiating oligo and mature RCA attachment are unresolved. Under ordinary template copying the RCA strand complements the padlock RO strings; that pairing relationship is inferred, while the actual labeled synthesis strands require supplementary verification.

**Exact source instances:** SEQUENCES/Bridge_sequences.csv row PadLock_001 is `TTTCGAGAGCCTTTAAGATAATGTGCCACA` (30 nt, reported file string; library instance, not universal). The `UsingPad.fa` PadLock_001 record is 144 nt and consists of two identical 72-nt units (derived by direct string comparison). All 110 entries have that duplication. The file is used for inter-padlock alignment (05_Probe_design cell 10), so it must not be transcribed as a 144-nt synthesized padlock. No unit was silently converted into an exact synthesized oligo.

03_Padlock_join cell 9 supplies exact RO example strings `GTGCCAACTCGTGCGTACTC` and `AACCGCTTACAACCGCTTAG`; these are notebook examples, not evidence of the experimental readout library or universal fluorophore conjugates. Source notation and case are preserved in local evidence. No reverse-complement nucleotide string has been asserted as a published exact sequence.

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region_a | Probe1 / probe A | Recognizes first selected RNA region | A nt 17–46 | Antisense in minus mode | 30 nt | Variable region; exact strings depend on target | variable | yes | Reverse complement of selected region | inferred | S04, V01 implementation snapshot, utils.py findAllCandidates; get_full_probes1 |
| target_binding_region_b | Probe2 / probe B | Recognizes second selected RNA region | B nt 1–30 | Antisense in minus mode | 30 nt | Variable region | variable | yes | Second region after 2-nt target gap | inferred | S04, V01 implementation snapshot, utils.py findAllCandidates; get_full_probes2 |
| other_functional_sequence_a | bridge_seq | Padlock-pairing overhang | A nt 1–14 | Same sequence as B[0:14], not its complement | 14 nt | Library substring; no universal sequence | barcode_library | yes | First 14 nt of selected 30-nt bridge | inferred | S04, V01 implementation snapshot, utils.py get_full_probes1 |
| other_functional_sequence_b | bridge_seq | Padlock-pairing overhang | B nt 33–46 | Same sequence as B[-14:] | 14 nt | Library substring | barcode_library | yes | Last 14 nt of selected bridge | inferred | S04, V01 implementation snapshot, utils.py get_full_probes2 |
| linker_a | No explicit term in source | Connects A overhang and recognition region | A nt 15–16 | 5′→3′ within oligo A | 2 nt | Exact code literal ta; connector role inferred | fixed | yes | Literal in inspected implementation only | inferred | S04, V01 implementation snapshot, utils.py get_full_probes1 |
| linker_b | No explicit term in source | Connects B recognition region and overhang | B nt 31–32 | 5′→3′ within oligo B | 2 nt | Exact code literal ta; connector role inferred | fixed | yes | Literal in inspected implementation only | inferred | S04, V01 implementation snapshot, utils.py get_full_probes2 |
| other_functional_sequence_pad_5 | reverse_random_bridge | Padlock terminal bridge-complementary segment | Padlock 5′ end | r[-7:], r=reverse complement of bridge | 7 nt | Library-derived substring | barcode_library | yes | Last seven bases of r | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3 get_whole_pad |
| other_functional_sequence_pad_3 | reverse_random_bridge | Padlock terminal bridge-complementary segment | Padlock 3′ end | r[0:23] | 23 nt | Library-derived substring | barcode_library | yes | First 23 bases of r | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3 get_whole_pad |
| readout_sequence_1 | RO1 | Presents one RO identity in padlock | Internal, order may swap with RO2 | Inserted as written; mature RCA strand complementary | 20 nt | Library-selected; see exact example below | barcode_library | yes | Joint folding screen | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3 TEST_ORDER |
| readout_sequence_2 | RO2 | Presents second RO identity in padlock | Internal, order may swap with RO1 | Inserted as written; mature RCA strand complementary | 20 nt | Library-selected | barcode_library | yes | Joint folding screen | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3 TEST_ORDER |
| linker_pad | No explicit term in source | Connects RO sequences | Between RO1 and RO2 | As selected in TEST_ORDER | 2 nt | Exact alternatives ta, at, tt, aa in code | variable | yes | Select order/connector passing screens | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3 TEST_ORDER |
| readout_probe | readout (RO) probe | Fluorescent interrogation of RCA product | Separate molecule | Synthesis strand/conjugation unresolved experimentally | 20 nt | RO library; exact conjugated oligo unresolved | barcode_library | yes | See RO library rules | reported | S04, V01 implementation snapshot, README Jupyter Notebooks Explanation; 01_RO_design cell 1 |

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | mRNA (minus), cDNA (plus); notebook also documents supplied genomic/artificial inputs | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 05_Probe_design cells 2,10 |
| Target region | Sliding 62-nt windows; two 30-nt recognition regions and a 2-nt intervening gap | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, utils.py findAllCandidates |
| Probe length | 30-nt target-binding region per A/B; 46 nt per assembled A/B (14+2+30) | V01 implementation snapshot | derived | S04, V01 implementation snapshot, utils.py get_full_probes1/2; sums shown here |
| GC constraints | Notebook 20–80%; Para constructor default 30–70%. Preserve both scopes | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 05_Probe_design cell 10; generator.py Para.__init__ |
| Tm constraints | Notebook sCelsius=47, nCelsius=37; target basicFilter rejects Tm below sCelsius. Model uses 30% formamide, 0.390 M monovalent and 0.0 M divalent ions in notebook comments | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 05_Probe_design cell 10; filter.py basicFilter |
| Specificity filtering | blastn-short/RefSeq RNA, taxonomy 10090, target-name classification and Tm screening; blacklist and inter-padlock checks. Notebook P1/P2 evalues 100/5; constructor defaults 100/50 | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, generator.py gene_generator, Para; cmdtools.py; filter.py |
| Repeat masking | No dedicated genome repeat-mask step established in examined design code; do not equate homopolymer filtering with repeat masking | V01 implementation snapshot | unavailable | S04, V01 implementation snapshot, generator.py; filter.py basicFilter |
| Exon selection | gtype may be gene, CDS or exon. Default gene. If unavailable, get_seqs falls back to whole record | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, generator_multi.py get_seqs; 05_Probe_design cell 10 |
| Isoform handling | Queries multiple NM_, XM_, NR_ records; identical target-region candidates are collapsed and counted. Frozen transcript accession set unresolved | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, generator_multi.py get_seqs; utils.py findAllCandidates |
| Sequence complexity | Reject N, Shannon entropy <1, GGGGG and CCCCC | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, filter.py basicFilter |
| Spacing | Within pair: 2-nt target gap. Between retained pairs: unresolved; manual selection | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, utils.py findAllCandidates; 05_Probe_design cell 14 |
| Overlap | Sliding candidate generation permits overlapping candidates; no universal final overlap rule established | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, utils.py findAllCandidates; 05_Probe_design cells 13–14 |
| Number of probes per target | User manually selects output indices; published experiment-specific pair counts unresolved | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 05_Probe_design cells 13–14 |

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | No single nucleotide barcode length: code is round plus color pair. Padlock has two 20-nt RO segments | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, 03_Padlock_join cells 3–6; 04_Panel_design cell 24 |
| Number of barcodes | Five colors give ten distinct pairs per round; notebook example is 110 genes in 11 rounds, not a universal limit | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 03_Padlock_join cell 4; 04_Panel_design cells 13,24 |
| Codebook | Columns Genes/RO1/RO2/Round in design output; demo uses gene/RO1/RO2/Round. Actual 20-gene two-round demo available as D02 | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 04_Panel_design cell 26; S05 CodingBook.CSV |
| Hamming distance | No explicitly reported error-correcting distance established. For five-bit constant-weight-two pair codes, minimum Hamming distance is 2 (derived: pairs sharing one channel differ at two positions) | V01 implementation snapshot | derived | S04, V01 implementation snapshot, 03_Padlock_join cell 6 combinations([0,1,2,3,4],2) |
| Error correction | Do not claim molecular error correction. Decoder classification and quality filtering differ from error-correcting code design | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, openDecode/decoding.py Decoding |
| Round assignment | GeneticChoose minimizes within-group Spearman-correlation objective after MAGIC imputation; user defines available slots | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, 04_Panel_design cells 18–26 |
| Channel assignment | Five-channel pair example AF488/AF546/AF594/AF647/AF750. Demo names last two Cy5/Cy7; labels must match image metadata and codebook | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 04_Panel_design cell 24; 06_Decode cell 2; S05 CodingBook.CSV |
| Fixed vs experiment-specific barcode sets | Round/gene assignment is panel-specific; bridge/RO libraries reusable. Decoder also supports single-color entries, without establishing a distinct chemistry version | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, 04_Panel_design; 06_Decode cell 3 |

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | 20 nt. RO-design example uses seqwalk max_orthogonality with four-base alphabet; three-base design is discussed as an alternative | V01 implementation snapshot | reported | S04, V01 implementation snapshot, 01_RO_design cells 1,4 |
| Readout probe | Separate RO molecules; exact synthesis modifications unresolved. RO strings are inserted directly in padlock by code | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, 03_Padlock_join cell 3; README |
| Fluorophore | Color labels given in codebooks; dye chemistry, supplier and terminal placement require experimental source inspection | V01 implementation snapshot | unavailable | S04, V01 implementation snapshot, 04_Panel_design cell 24; S05 CodingBook.CSV |
| Imaging round | One round contributes a color-pair identity. Demo has R1 and R2 | V01 implementation snapshot | reported | S04, V01 implementation snapshot, S05 CodingBook.CSV; 06_Decode cell 4 |
| Sequencing primer | Not established; do not rename RO probes sequencing primers | V01 implementation snapshot | unavailable | S04, V01 implementation snapshot, S04 notebooks 01–06 |
| Decoding role | Codebook maps round and channel combination to gene; spatial feature classification implemented | V01 implementation snapshot | inferred | S04, V01 implementation snapshot, openDecode/decoding.py Decoding; 06_Decode cell 3 |

## Primer / Linker / Adapter Components

See linker_a, linker_b and linker_pad above. The 30-nt bridge is a design-library sequence split across target probes and complemented in padlock termini; it is not established as a third separately supplied bridge oligo. No universal PCR adapter, promoter or sequencing primer is established (inferred/unavailable; S04 utils.py and 03_Padlock_join cell 3; V01 implementation).

## Amplification-Related Components

RCA is reported in S04 README. Padlock template and RO segments are evidenced; the exact amplification primer molecule and its 3′ initiation site remain unresolved. Do not assign primer activity to probe A or B solely from their sequence order (S04 01_RO_design cell 1; 03_Padlock_join; V01).

## Chemical Modifications

Not resolved from inspected sequence/code resources: 5′ phosphorylation, fluorophore attachment, blocking groups, hydrogel anchors and cleavage chemistry. Plain FASTA strings do not establish unmodified oligos or negate modifications (S04 SEQUENCES/UsingPad.fa; V01).

## Probe-Design Workflow

S04, V01 implementation: (1) design 20-nt ROs using sequence composition, repeat/folding, off-target and orthogonality screens (01); (2) generate 30-nt bridges and screen against reference RNA, ROs and each other (02); (3) assemble bridge-complementary padlock ends around two RO identities and folding-selected connector/order (03); (4) use scGIST or a supplied panel, MAGIC and a genetic algorithm to assign rounds/colors (04); (5) supply target CSV, bridge/padlock files and RefSeq RNA BLAST database; generate, filter, score and manually choose A/B pairs, then export ordering table (05). Notebook model temperatures are computational inputs. Reference release, random seeds and final manual selections must be frozen separately for reproducibility. Source code was inspected, not run.

## Official Probe-Design Software

Author provenance: S03 README identifies the study and links S04 as its step-by-step computational workflow; S07/S05 deposits identify Xinyang Li and link the associated study/workflow. These links support official author status, beyond a name match. Repositories were publicly cloned on 2026-10-15; software was not installed or run.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| openDesign and notebooks 01–05 | S04; companion S03 | Author-provided | S04 be490eb19503c742230665da1e665c6d6c3d97c8; S03 01ea661ef43e1c6c7bc19a3a6e3c214a420a262c | Python/Jupyter | Target CSV, target records or sequences, RefSeq RNA database, bridges, padlocks, optional background list | Candidate pair table, manually selected ordering CSV; RO/bridge/padlock libraries and panel codebook | Candidate scan, GC/Tm/complexity filters, BLAST, blacklist, folding, manual selection | Biopython, pandas, numpy, BLAST+, seqwalk, NUPACK, primer3-py, pandarallel; panel scGIST/MAGIC/Scanpy; see env files | Source available 2026-10-15 | No claim about continued maintenance; no pinned package release established | inferred | S04 V01 README; 01–05 notebooks; openDesign; env1/2_requirement.txt |

NUPACK is separately installed per README. Target-design code has a primer3 fallback for binding-efficiency calculation; this does not make the RO/padlock notebooks NUPACK-independent. The notebook and constructor differ in GC and P2 BLAST evalue defaults; do not merge them. generator.py docstring labels ion inputs mM while notebook and computation conversions indicate M; preserve notebook units and review this ambiguity before execution.

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | S04 commit be490eb19503c742230665da1e665c6d6c3d97c8; companion S03 | V01 implementation snapshot | 110 bridge sequences; 110 doubled padlock alignment references; blacklist and tRNA FASTA; single-target Gad2 input template; not a complete synthesized gene-probe table | Retrieved 2026-10-15 | S04 SEQUENCES and Example_ProbeDesign.csv |
| D02 | 10.5281/zenodo.20390362 | V01 protocol demo | 20-gene R1/R2 color codebook; imaging/decoded/metabolomics files listed, not a full probe list | Metadata and CodingBook.CSV retrieved 2026-10-15 | S05 record metadata/files and CodingBook.CSV |

D03/D04 are data-analysis resources, listed in benchmarks.md. Checksums for local source files are in evidence/protocol-sha256.json. Complete biological target-binding and chemically annotated synthesis tables were not located in the inspected repository trees; publisher supplements remain inaccessible. This is a partial-resource availability finding, not evidence that those tables do not exist.

## Version Differences

V01 is the original 2026 publication. S02 preprint-to-publication changes remain unresolved. S03 and S04 are different code snapshots; do not assume their defaults or filename conventions are identical. OpenDecode Version 2.0 is a software label (S04 06_Decode cell 0).

## Unresolved Information

Published supplementary probe and RO tables; experimentally used probe counts, positions and source transcript versions; true synthesized padlock strand/modifications; ligation/RCA initiation details; readout conjugation/removal; exact reference releases; final selected target probes. The filename mm39_tRNA.fa establishes a blacklist resource label only, not the assembly/annotation of every biological probe. Source-access gaps must remain explicit (S01–S04).
