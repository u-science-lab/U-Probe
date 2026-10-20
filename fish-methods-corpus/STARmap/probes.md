# Probe Architecture

## Overview

STARmap uses the SNAIL primer/padlock pair to recognize RNA. The primer splints padlock ligation and primes Phi29 RCA; the circularized padlock is the RCA template. V01 reports a five-base gene identifier with 1,024 possible nucleotide barcodes (S01 Results, SNAIL and SEDAL). S03 is a later developer protocol explicitly describing STARmap and STARmap PLUS probe architecture. Its detailed rules are not proof of the parameters used for every original V01 panel. Evidence below is reported unless explicitly marked inferred, derived or unavailable.

## Canonical Architecture

### Primer and padlock (linear precursors, 5′ → 3′)

```text
Primer: 5′ — [19–25 nt RNA-targeting region] — [12 nt padlock-pairing region] — 3′
Padlock: 5′ phosphate — [6 nt primer-pairing end] — [internal target/barcode/readout regions] — [6 nt primer-pairing end] — 3′ OH
```

Reported component inventory: S03, Padlock and primer probe design, pp. 8–9 (V01/V02 architecture). The padlock internal regions include a 19–25 nt target region, 1–5 nt barcode and 12 nt reading-probe-hybridizing region. The prose groups both terminal 6-nt segments together rather than giving an unambiguous continuous sequence map; internal order and the precise relationship between the internal barcode and barcode-containing terminal pairing regions are left unresolved here. S03 explicitly permits a barcode within the paired terminal regions when needed. Do not assert that the primer's 12 nt must complement only the separate internal barcode, or that the terminal arms exclude barcode information. Primer order is consistent with the example oligos in S03 Table 2 (inferred from sequence structure). For V01, S02 (Fig. S1A legend, p. 10) reports the 19–25 nt RNA-binding segment on each probe and the 6-nt primer–padlock pairing on each arm; it does not give the 12-nt or barcode segment lengths. This is a schematic, not a synthesis sequence.

The two RNA-hybridizing sites have a recommended **1–2 nt interspace** in S03 design (0–2 nt in the V01 design, S02 p. 3). They are distinct from the padlock ends that pair with the primer. Exact target-strand placement requires the selected oligo/reference pair; no assembled RNA/padlock linear diagram is asserted. Ligation closes the padlock's 3′ OH and 5′ phosphate. RCA produces tandem complements of that circle (S03 Overview).

### SEDAL oligos (5′ → 3′)

```text
Reading: 5′ phosphate — [cycle-specific sequence; mixed bases where specified] — [3′ inverted dT in S03 Table 4] — 3′
Decoding: 5′ Alexa dye — [six mixed positions] — [specified dinucleotide] — 3′ OH
```

Reported: S01 Fig. 1E (5′ phosphate and progressively degenerate reading positions); S03 Tables 4–5 (V02 protocol library). The decoding-probe 3′ end ligates to the reading-probe 5′ phosphate. Reading probes are not a fixed N11 oligo: Table 4 lists SeqA1–A6 and SeqB1–B6, while Step 61 names SeqA1–6 and SeqB1–4. Later SeqA probes gain mixed positions at the **5′ side**. The source overview describes 11-nt reading probes; Table 4 source notation, mixed-base annotations and inverted-dT modification should be consulted for each oligo rather than imposing a universal count. Tables 4–5 call the random positions ribonucleotides while also directing standard equal-ratio A/T/C/G mixed-base ordering; this source wording is retained as a chemistry ambiguity, not silently resolved. Exact published library strings are linked by table, not reconstructed here.

## Probe Components

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region.primer | gene-targeting region | Recognizes RNA | Primer target segment | 5′ → 3′ | 19–25 nt | variable_region; exact per-gene string not transcribed | variable | yes | 1–2 nt target-site interspace; see design table | reported | S03 Padlock and primer probe design, pp. 8–9; V01/V02 architecture |
| other_functional_sequence.primer_pairing | 12-nt sequence | Splints padlock ends; primes RCA | Primer 3′ region (sequence interpretation) | 5′ → 3′ | 12 nt | schematic_sequence; barcode-containing pairing segment | sample_specific | yes | Paired to padlock ends; barcode location not forced to internal barcode | reported | S03 design section, p. 8; Table 2; V01/V02 architecture |
| target_binding_region.padlock | gene-targeting region | Recognizes RNA | Internal padlock region | 5′ → 3′ | 19–25 nt | variable_region | variable | yes | See design table | reported | S03 design section, pp. 8–9; V01/V02 architecture |
| other_functional_sequence.padlock_5_end | 6 nt at the 5′ end | Pairs to primer; ligation substrate | Padlock 5′ end | 5′ → 3′ | 6 nt | schematic_sequence; exact bases depend on design | sample_specific | yes | Complementary to primer; barcode may be present if needed | reported | S03 design section, p. 8; V01/V02 architecture |
| other_functional_sequence.padlock_3_end | 6 nt at the 3′ end | Pairs to primer; ligation substrate | Padlock 3′ end | 5′ → 3′ | 6 nt | schematic_sequence | sample_specific | yes | Complementary to primer; not an RNA-targeting arm | reported | S03 design section, p. 8; V01/V02 architecture |
| barcode.padlock | gene-specific barcode | Gene identity | Padlock; internal/pairing instances must not be conflated | 5′ → 3′ | 5 nt in S01; 1–5 nt segments in S03 architecture; longer total barcodes discussed in S03 | barcode assignment; exact panel strings unexamined | barcode_library | yes | Avoid dinucleotide repeats and repeated-color patterns where possible (S03) | reported | S01 Results; S03 design section, pp. 8–9; V01/V02 |
| readout_sequence.padlock | reading probe-hybridizing region | Template recognition site for readout | Internal padlock; read as complement in amplicon | 5′ → 3′ | 12 nt | schematic_sequence; do not assert one universal sequence across SeqA/SeqB configurations | version_specific | yes | Consult selected library and S03 Table 2 | reported | S03 design section, p. 8; Tables 2, 4; V01/V02 architecture |
| chemical_modification.padlock_5_phosphate | 5′ phosphate | Enables padlock circularization | Padlock 5′ end | 5′ end | Not applicable | Not a nucleotide string | fixed | yes | Retain /5Phos/ source notation for exact oligos | reported | S03 design section; Table 2 note; V01/V02 |
| primer.reading | reading probe | Sets SEDAL position | Amplicon, not unamplified RNA | 5′ → 3′ | 11 nt in overview; consult Table 4 per oligo | degenerate_sequence library; not fixed across cycles | version_specific | yes | Cycle-specific; 5′ phosphate; 3′ inverted dT in Table 4 | reported | S01 Fig. 1E; S03 Overview, Table 4; V01 conceptual/V02 library |
| readout_probe.decoding | fluorescent decoding probe | Reports dinucleotide; ligates to reading 5′ phosphate | Amplicon; decoding 3′ OH at junction | 5′ → 3′ | 8 nt | degenerate_sequence: six mixed positions + specified dinucleotide | barcode_library | yes | 16 oligos; equal A/T/C/G mixed bases; four dye groups | reported | S03 Table 5; S01 Fig. 1E; V01 conceptual/V02 library |
| fluorophore.decoding | Alexa Fluor | Color readout | Decoding-probe 5′ end | 5′ attachment | Not applicable | Not a nucleotide string | barcode_library | yes | 488, 546, 594 or 647, assigned by Table 5 | reported | S03 Table 5; V02 protocol |
| chemical_modification.amplicon_amine | 5-(3-aminoallyl)-dUTP | Installs amines during RCA | RCA amplicon | Not applicable | Not applicable | Modified nucleotide; not a separate oligo | version_specific | yes | 20 µM with 250 µM dNTP mixture | reported | S03 Step 24; V02 protocol |
| chemical_modification.acrylic_NHS | acrylic acid N-hydroxysuccinimide esters | Hydrogel coupling | Amplicon before polymerization | Not applicable | Not applicable | Not an oligo | version_specific | yes | Numerical conditions not in examined S01 main text | reported | S01 Results, HTC; V01 |
| chemical_modification.MA_NHS | methacrylic acid N-hydroxysuccinimide esters | Amplicon/antibody hydrogel coupling | Before polymerization | Not applicable | Not applicable | Not an oligo | version_specific | yes | 25 mM in 100 mM bicarbonate, room temperature, 1 h | reported | S03 Steps 31–32; V02 protocol |

Original terms are retained beside functional identifiers. Primer/padlock pairing regions are not classified as RNA-targeting padlock arms. Component lengths are not summed into universal full-oligo lengths because barcode overlap and longer library configurations require sequence-level verification.

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | RNA recognition without reverse transcription | V01 | reported | S01 Results, SNAIL; V01 |
| Target region | Coding region of the shortest isoform for genes with several isoforms; non-coding RNAs used whole ("only the shortest isoforms were considered and the coding regions were used except for non-coding RNAs") | V01 | reported | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Probe length | 40–46 nt hybridization sequence per primer/padlock pair, designed with Picky 2.2, split into halves of 20–25 nt; the Fig. S1A legend gives 19–25 nt per probe; "All the probes were under 60 nt" | V01 | reported | S02 Materials and Methods, "STARmap probe design", p. 3; V01; S02 Fig. S1A legend, p. 10; V01 |
| GC constraints | None stated in the V01 design section | V01 | unavailable | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Tm constraints | Each probe's RNA-binding segment is designed to a Tm of 60 °C; the 40–46 nt sequence is split for "the best match of melting temperature (Tm) between the two halves"; the 6-nt primer–padlock pairing has a Tm below room temperature, negligible during DNA–RNA hybridization at 40 °C | V01 | reported | S02 Fig. S1A legend, p. 10; V01; S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Specificity filtering | No transcriptome off-target, BLAST or secondary-structure screen stated in the V01 design section | V01 | unavailable | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Isoform handling | Shortest isoform only | V01 | reported | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Spacing | 0–2 nt gap between the two halves, i.e. between the primer and padlock target sites | V01 | reported | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Number of probes per target | "4 sequences for each gene were designed" (4 primer/padlock pairs per gene) | V01 | reported | S02 Materials and Methods, "STARmap probe design", p. 3; V01 |
| Target region | Isoforms extracted from target-species transcriptome and consensus coding sequence databases; no universal ORF/UTR placement asserted | V02 design recommendations | reported | S03 design section, pp. 8–9; V02 design recommendations |
| Probe length | 40–46 nt candidate hybridization sequence split into two 19–25 nt segments; preserve these source ranges without assuming all endpoint combinations coexist | V02 | reported | S03 design section, pp. 8–9; V02 |
| GC constraints | Recommended 30–70% | V02 | reported | S03 design section, p. 9; V02 |
| Tm constraints | Recommended 56–70°C at 5 nM oligo and 300 mM salt; source gives salt without specifying monovalent composition | V02 | reported | S03 design section, p. 9; V02 |
| Specificity filtering | NUPACK for strong secondary structures; BLAST each probe against minus strand of the whole probe library for heterodimers; typically 10–20% removed. This BLAST step is not a reported transcriptome off-target screen | V02 | reported | S03 design section, p. 9; V02 |
| Repeat masking | Not reported in the examined S03 design section; probe-library heterodimer filtering is not repeat masking | V02 | unavailable | S03 design section, p. 9; V02 |
| Exon selection | Not reported as an exon-selection algorithm in examined S03; isoform/CCDS extraction alone does not establish one | V02 | unavailable | S03 design section, p. 9; V02 |
| Isoform handling | Extract isoforms from transcriptome and consensus coding sequence databases; reference releases and isoform prioritization unresolved | V02 | reported | S03 design section, p. 9; V02 |
| Sequence complexity | Strong secondary structures filtered with NUPACK; no separate numerical low-complexity cutoff reported in examined section | V02 | reported | S03 design section, p. 9; V02 |
| Spacing | 1–2 nt interspace between paired RNA-targeting regions | V02 | reported | S03 design section, p. 9; V02 |
| Overlap | Dedupe from BBMap removes overlaps and containments from candidate pool | V02 | reported | S03 design section, p. 9; V02 |
| Number of probes per target | Typically 3–6 primer/padlock pairs per transcript suggested for budget; ideally design >6 pairs before filtering. Two molecular species per pair does not mean only one pair per gene | V02 | reported | S03 design section, p. 9; V02 |

V01 rows come from the original supplementary materials (S02 PDF, read 2026-10-20); they name Picky 2.2 but no reference build, annotation release or other software settings. V02 rows are the later S03 recommendations, which do not identify the reference build, annotation release or exact software settings used for the original 160-/1,020-gene panels. The sources differ on the half length (20–25 nt in S02 Methods, 19–25 nt in the S02 Fig. S1A legend and S03) and the gap (0–2 nt in S02, 1–2 nt in S03).

## Barcode System

V01 uses a five-base gene-unique padlock identifier with 1,024 possible sequences (reported: S01 Results). This is a physical sequence distinct from the abstract barcode–gene assignment and round color code. S03 design section, p. 9, gives 4-nt (256) and 5-nt (1,024) examples and explicitly discusses **barcodes longer than 5 nt**, read in two parts using orthogonal SeqA/SeqB libraries. Thus 1–5 nt component descriptions are not a global V02 barcode-length ceiling. The theoretical count is not the experimentally demonstrated panel size or a filtered usable codebook size.

S03 recommends avoiding dinucleotide repeats and barcodes producing the same color throughout rounds where possible. Gene uniqueness is scoped to a panel; no prohibition on reuse across different panels is established. Minimum Hamming distance is not reported in the examined design section. Exact original gene assignments remain unexamined in S02; do not assume workbook completeness. Two-base decoding rejects inconsistent color strings rather than guaranteeing detection of every error (S01 SEDAL Results; S03 Overview).

## Readout System

Reported: S03 Overview, Tables 4–5 and Steps 50–61 (V02 protocol); S01 Fig. 1E (V01 concept). Readout occurs on RCA amplicons complementary to the padlock, not directly on RNA. Cycle-dependent reading probes set the position; 16 degenerate decoding oligos report dinucleotides through four colors. Overlapping interrogations shift by one base, followed by a quality-control round for the approximately 1,000-gene configuration.

| Specified decoding-oligo 3′ dinucleotide | 5′ dye | Source / locator |
|---|---|---|
| AA, CC, GG, TT | Alexa Fluor 488 | S03 Table 5; V02 library |
| AC, CA, TG, GT | Alexa Fluor 546 | S03 Table 5; V02 library |
| AG, GA, CT, TC | Alexa Fluor 594 | S03 Table 5; V02 library |
| AT, TA, GC, CG | Alexa Fluor 647 | S03 Table 5; V02 library |

The table specifies **oligo-side dinucleotides**; it is not a direct 5′→3′ RNA or padlock barcode table. Strand conversion and round ordering require the actual library and codebook. Table 4 explicitly gives phosphorylated cycle-specific reading sequences with 3′ inverted dT; those modified/degenerate source strings must not be described as a single fixed exact nucleotide sequence. Per-cycle mixture concentrations are in protocols.md (S03 Step 50).

## Primer / Linker / Adapter Components

The primer probe splints padlock ligation and primes RCA; the padlock is the RCA template. The S03 primer/padlock architecture does not specify a separate linker or adapter oligo. Exact primer-pairing segments depend on the selected design and should not be equated automatically with the internal barcode (reported inventory: S03 design section and Table 2; V01/V02 architecture).

## Amplification-Related Components

S03 Step 24 (V02 protocol) specifies 250 µM dNTP mixture and 20 µM aminoallyl-dUTP; retain the source concentration convention rather than silently interpreting total versus per-base concentrations. The primer extends on the circularized padlock template. The resulting amines support hydrogel coupling; numerical procedure conditions belong in protocols.md.

## Chemical Modifications

Reported modifications: padlock 5′ phosphate (S03 Table 2), reading-probe 5′ phosphate and 3′ inverted dT (Table 4), decoding-probe 5′ Alexa dyes (Table 5), and aminoallyl-dUTP incorporation (Step 24). V01 acrylic-NHS chemistry is reported by S01 HTC Results, but numerical V01 conditions are not in the examined main text; S02 was read for probe design only, not for these conditions. S11 Methods also reports acrylic NHS, at 20 mM in PBST for 1 h (V02 application), distinct from later S03 MA-NHS. S03 Steps 31–32 use MA-NHS at 25 mM in 100 mM bicarbonate, room temperature, 1 h, for the later V02 protocol. Nissl and DAPI are segmentation stains, not covalent oligo modifications. Protein comapping uses antibody incubation before MA-NHS and secondary antibody after sequencing (S03 Steps 28–30, 62–63).

## Probe-Design Workflow

S03 Padlock and primer probe design, pp. 8–9 (reported; later V02 recommendations):

1. Extract gene isoforms from species transcriptome and consensus coding sequence databases; exact release unresolved here.
2. Use Picky 2.2 for 40–46 nt hybridization candidates with recommended GC 30–70% and Tm 56–70°C at 5 nM oligo/300 mM salt.
3. Split candidates into 19–25 nt target segments separated by 1–2 nt and match segment Tm.
4. Remove overlaps/containments with Dedupe from BBMap.
5. Assign gene barcodes and assemble the paired probes; use separate SeqA/SeqB reading parts for >5 nt barcodes as described by S03.
6. Remove strong secondary structures with NUPACK and probe heterodimers by BLAST against the probe-library minus strand; typically 10–20% of probes are removed.

S03 does not establish that S05 implements this entire stack, or that the later stack generated the original 2018 panels. For V01, S02 ("STARmap probe design", p. 3) reports a shorter procedure: shortest isoform, coding region (non-coding RNAs whole); Picky 2.2 for 40–46 nt hybridization sequences, four per gene; each split into halves of 20–25 nt with a 0–2 nt gap and the best Tm match between halves; all probes under 60 nt, ordered from IDT. Exact original panel input/output files (S02 Tables S1–S2) and reference releases remain unexamined.

## Official Probe-Design Software

The **developer-described pipeline** is the S03 multi-tool stack below. Use of third-party tools by the developers does not make those tools STARmap-authored software. S05 is a related SNAIL tool hosted by Bava's lab; canonical implementation status and equivalence to S03 are not established. Repository popularity is not evidence of methodological authority or maintenance.

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Picky | S03 and S02 cited tool | Developer-used external tool | 2.2 | Not examined | Target isoform sequences | Hybridization candidates | GC/Tm selection (V02); 40–46 nt pair sequences, four per gene (V01) | Not examined | Named in protocol | Not examined | reported | S03 design section, p. 9; V02; S02 STARmap probe design, p. 3; V01 |
| Dedupe / BBMap | S03 cited tool | Developer-used external tool | Not specified in examined section | Not examined | Candidate sequences | Deduplicated pool | Overlap/containment removal | Not examined | Named in protocol | Not examined | reported | S03 design section, p. 9; V02 |
| NUPACK | S03 cited tool | Developer-used external tool | Not specified in examined section | Not examined | Full probes | Filtered candidates | Strong secondary structures | Not examined | Named in protocol | Not examined | reported | S03 design section, p. 9; V02 |
| BLAST | S03 cited tool | Developer-used external tool | Not specified in examined section | Not examined | Probes vs probe-library minus strand | Heterodimer-filtered library | Probe–probe interactions | Not examined | Named in protocol | Not examined | reported | S03 design section, p. 9; V02 |
| snail-probe-designer | https://github.com/bavalab/snail-probe-designer | Related Bava-lab SNAIL utility; original-panel use and S03 equivalence unverified | Repository lists v0.0.1–v0.0.3 executables; no tested release | Python (repository directory and README) | Detailed input format not examined | SNAIL probe designs; detailed output format not examined | README cites Wang 2018; rules not independently inspected in code | Not examined | Public page inspected 2026-10-15 | Unverified | reported | S05 README and root file listing; V01-related, not attributed to original panel generation |

S04, S06, S07 and S08 are analysis resources, not probe-design implementations. Their source records remain in sources.md.

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | aat5691_wang_sm_table-s1.xlsx | V01 | Original supplementary workbook candidate; exact sheets, sequence completeness and assignments unexamined | Referenced by S01/S02; direct content not retrieved | S02 Table S1; S01 supplementary section |
| D02 | aat5691_wang_sm_table-s2.xlsx | V01 | Original supplementary workbook candidate; exact contents unexamined; S02 p. 3 describes Table S2 as the SEDAL probes; file sizes are not interpreted | Referenced by S01/S02; direct content not retrieved | S02 Table S2; S01 supplementary section |
| D03 | ClusterMap/datasets/RawImage | V01-related | Image-resource candidate; full inventory and any accompanying codebooks unexamined | Prior existence record; not re-inventoried in Fix | S08 repository; see benchmarks.md |
| D04 | SCP1375 | V02 | STARmap PLUS data deposit; complete file inventory unexamined | Reported by S11 Data availability | S09; S11 Data availability |
| D05 | 10.5281/zenodo.7332091 | V02 | Companion deposit; complete inventory unexamined | Reported by S11 Data availability | S10; S11 Data availability |

Additional official resources in S03: Supplementary Tables 1–2 contain preselected human/mouse target regions (reported by design section; files not inspected); **Table 2** publishes example STARmap primer/padlock oligos for human ACTB, MALAT1, vtRNA1-1 and GAPDH; Tables 4–5 publish the reading/decoding libraries (table text inspected during Fix). These are later-protocol resources, not verified full replacements for original D01/D02. No benchmark grade is assigned merely from their existence.

## Version Differences

V01 has a five-base identifier and acrylic-NHS HTC (S01 Results). V02 STARmap PLUS integrates protein detection and improved probe design (S03 Overview). S11 Methods uses acrylic NHS, whereas later S03 uses MA-NHS; the chemistry difference belongs to the source-specific procedure rather than defining all V02 experiments. S11 maps 2,766 genes in eight sequencing rounds and uses PI segmentation; later S03 uses DAPI. Its design discussion permits longer barcodes and orthogonal reading libraries; 1–5 nt components are not a method-wide ceiling. Protein-retaining experiments skip proteinase K (Steps 42–43). RIBOmap and TEMPOmap have different tri-probe designs and remain sibling methods. S03 explicitly describes the di-probe architecture for both STARmap and STARmap PLUS, but later numerical protocol conditions are not retroactively assigned to the 2018 panels.

## Unresolved Information

Original S02 workbook contents (Tables S1–S2), the parts of the S02 supplementary procedure beyond probe design, reference assemblies/annotation releases, coordinates (including origin, interval convention and strand), original panel software settings, and complete dataset inventories remain unexamined. Exact placement/overlap of the internal barcode and barcode-containing terminal pairing segments is not resolved from the inspected S03 prose. S03 Table 4 and Step 61 differ in the listed SeqB range, and its overview length is not imposed on every listed oligo. Mixed-base chemistry wording in Tables 4–5 requires care. Explicit repeat-masking/exon-selection algorithms and a minimum barcode Hamming distance were not found in the examined S03 design section. Software maintenance and S05 equivalence to the later pipeline remain unverified. These gaps do not license invented sequences or Grade A/B ground truth.
