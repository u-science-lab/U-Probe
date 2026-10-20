# Probe Architecture

## Overview

V01 untargeted ExSeq uses modified random-octamer primers, cDNA circularization and RCA. V02/V03 original targeted experiments and V04 generalized targeted design use DNA padlocks ligated directly on RNA. V02 brain encoding, V03 cancer encoding and V04 recommended libraries are distinct scopes [S02, supplementary pp. 5-8, 18-23, 29-31; S12, abstract; reported]. Internal version grouping is described in description.md.

## Canonical Architecture

### Targeted padlock, V04 default implementation

S04 `Targeted_ExSeq_ProbeGeneration_Assembly.m` lines 463-482 establishes the synthesized order below; S12 Fig. 1A is the conceptual schematic. This is a **schematic sequence**, not a published full nucleotide string:

```text
5'-phosphate — [5' terminal homology arm, 16 nt] — [AAA, 3 nt]
 — [SOLiD barcode, 7 nt] — [backbone, 22 nt]
 — [Illumina barcode, 7 nt] — [AAA, 3 nt] — [3' terminal homology arm, 16 nt] — 3'-OH
```

The implementation names the terminal segments “5' homology arm” and “3' homology arm”. For a selected RNA segment written 5'→3', let P be its 32-nt reverse complement written 5'→3'. The probe starts with P[17:32] and ends with P[1:16], using one-based closed intervals. Thus the probe 5' terminal arm pairs with the RNA's upstream half, and the probe 3' terminal arm pairs with the downstream half; each duplex is antiparallel. The arms bring the probe 3'-OH and 5'-phosphate next to each other over an RNA-splinted **nick**, with no intervening RNA base. Strand/interval interpretation is **inferred** from S04 homology export and assembly lines 463-482; the two 16-nt terminal portions are **reported** in S02 supplementary p. 18 (V02/V03).

Default length is **derived**: 16+3+7+22+7+3+16=74 nt; the phosphate is not a nucleotide. S12 abstract p. 5 reports this default, and all 41 nonempty D01 workbook sequence cells have 74 nucleotides after removal of `/5Phos/` (**derived**, S04 `ExampleOutputAllProbes.xlsx`, Sheet1 column D). Original targeted backbone/spacers and experiment-specific barcode implementation are described in S02 supplementary pp. 18-20, 29 and Fig. S12; the V04 codebook must not be backfilled into them.

### Targeted circle, RCA product and separate primers

SplintR joins the terminal DNA ends to form a circle; the circle has no linear 5'/3' endpoints. Phi29 produces a linear concatamer of the reverse complement of the circle, with the starting phase determined by priming. No fixed copy count is established here [S02, supplementary pp. 20-22; S13, Steps 9-12; reported mechanism, inferred strand representation; V02-V04].

The default padlock backbone is 5'-`TCTCGGGAACGCTGAAGACGGC`-3'. The separate RCA primer is 5'-`TCTTCAGCGTTCCCGA*G*A`-3', **18 nt**, with `*` denoting phosphorothioate bonds. The SBS sequencing primer is 5'-`TCTCGGGAACGCTGAAGACGGC`-3', **22 nt**, identical to the padlock backbone but binding its complement on the amplicon [S13 Step 5; S14 materials/Step 10; S04 assembly lines 83-91; reported, V04]. Reverse-complementing backbone bases 1-18 yields the RCA primer letters (**derived**, ignoring modifications); it is not the reverse complement of the last 18 bases. The separate Cy3 detection oligo also binds the amplicon complement of the backbone [S13, Step 6; reported molecule, inferred pairing].

### Untargeted oligonucleotides and cDNA, V01

```text
5'-phosphate — [sample-specific constant region] — [random octamer NNNNNNNN] — 3'
RT → enzymatically processed single-stranded cDNA → cDNA circle → RCA concatamer
```

This is a schematic; the degenerate primer strings below are separately reported reagent notations. DNase I treatment precedes RT. RNase/RNase H and Endonuclease V processing precede CircLigase II circularization; inosine inclusion supports cDNA fragmentation. No target-specific gene panel is required, but functional oligonucleotides are required [S02, supplementary pp. 5-6; reported, V01].

## Probe Components

DNA strings are written 5'→3' in source synthesis notation. `/5Phos/`, `/5Cy3/` and `*` denote modifications, not bases. Fixed classes below refer to the stated implementation, not every possible ExSeq design. Chemical groups have no nucleotide length. The QC detection oligo is categorized as other_functional_sequence because it is an independent fluorescent hybridization reagent, without an anchoring or primer-initiation role. Suffixes identify physical instances without introducing new functional categories.

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| padlock_arm_5 | 5' homology arm | RNA recognition, upstream RNA half | Probe 5' terminus | DNA 5'→3'; antiparallel RNA | 16 nt default | Variable region P[17:32], defined above | variable | yes | Split reverse-complemented 32-mer | inferred | S04 assembly lines 463-482; V04 |
| spacer_5 | SPACER_SEQUENCE | Separates arm and barcode | After 5' arm | DNA 5'→3' | 3 nt default | Fixed-region exact default `AAA` | fixed | yes | Configurable parameter | reported | S04 assembly lines 93-97, 482; V04 |
| barcode_solid | SOLiD barcode | Readout encoding | Between first spacer and backbone | DNA 5'→3' | 7 nt default | Barcode-library region; conversion depends on specification | barcode_library | yes | See Barcode System | reported | S12 Fig. 1, Step 12.2; S04 assembly line 482; V04 |
| amplification_sequence_backbone | RCA_PRIMER_BINDING_SITE | Circle binding site for RCA primer; amplicon complement is SBS binding site | Between barcodes | DNA 5'→3' | 22 nt default | Fixed-region exact default `TCTCGGGAACGCTGAAGACGGC` | fixed | yes | Configurable backbone | reported | S04 assembly lines 83-91, 482; V04 |
| barcode_illumina | Illumina barcode | SBS encoding | Between backbone and second spacer | DNA 5'→3' | 7 nt default | Barcode-library region; default direct map 0/1/2/3→A/C/G/T on padlock | barcode_library | yes | Experiment-specific gene assignment | reported | S04 assembly Illumina mapping and line 482; V04 |
| spacer_3 | SPACER_SEQUENCE | Separates barcode and arm | Before 3' arm | DNA 5'→3' | 3 nt default | Fixed-region exact default `AAA` | fixed | yes | Same configured spacer | reported | S04 assembly lines 93-97, 482; V04 |
| padlock_arm_3 | 3' homology arm | RNA recognition, downstream RNA half | Probe 3' terminus | DNA 5'→3'; antiparallel RNA | 16 nt default | Variable region P[1:16], defined above | variable | yes | Split reverse-complemented 32-mer | inferred | S04 assembly lines 463-482; V04 |
| primer_rca | RCA primer | Phi29 initiation on padlock | Separate reagent, binds circle backbone | DNA 5'→3' | 18 nt | Exact sequence `TCT TCA GCG TTC CCG A*G*A` | fixed | yes | Phosphorothioate bonds as shown | reported | S13 Step 5; V04 |
| primer_sbs | sequencing primer | SBS initiation on RCA product | Separate reagent | DNA 5'→3' | 22 nt | Exact sequence `TCTCGGGAACGCTGAAGACGGC` | fixed | conditional | Illumina readout | reported | S14 materials/Step 10; V04 |
| other_functional_sequence_detection | amplicon detection oligonucleotide | Fluorescent QC hybridization, not an anchor or synthesis primer | Separate reagent | DNA 5'→3' | 18 nt | Exact notation `/5Cy3/TCTCGGGAACGCTGAAGA` | fixed | conditional | Optional QC; dye may be changed | reported | S13 Step 6; V04 |
| chemical_modification_phosphate | 5' phosphorylation | Enables circle closure | Padlock 5' end | 5' chemical group | Not applicable (group) | `/5Phos/`, not nucleotide sequence | fixed | yes | Included in assembly | reported | S12 Step 13.1; S04 assembly line 482; V04 |
| chemical_modification_ps | phosphorothioate backbone modification | Protects RCA primer from Phi29 exonuclease | Final two RCA-primer bonds | Bond positions in 5'→3' oligo | Not applicable (bonds) | `A*G*A` notation in primer_rca | fixed | yes | Preserve source modification notation | reported | S13 Steps 5, 10; V04 |
| chemical_modification_aadutp | aminoallyl dUTP (aadUTP) | Amines for amplicon cross-linking | RCA product | Not a linear sequence position | Not applicable (modified nucleotide) | No separate oligo sequence | fixed | yes | 40 µM added to 250 µM each dNTP; no incorporation fraction specified | reported | S13 Step 11.3; V04 |
| fluorophore_detection | /5Cy3/ | QC fluorescence | Detection oligo 5' end | 5' group | Not applicable (dye) | Cy3; Alexa 546 or another suitable dye permitted | version_specific | conditional | Dye substitution allowed | reported | S13 Step 6; V04 |
| primer_rt_15 | random octamer reverse transcription primer | RT on RNA, with constant region | Separate RT reagent | DNA 5'→3' | 31 nt (derived count) | Degenerate notation `/5Phos/ACTTCAGCTGCCCCGGGTGAAGANNNNNNNN` | sample_specific | conditional | 15-µm hippocampal slice | reported | S02 supplementary p. 5; V01; length derived from string |
| primer_rt_50 | random octamer reverse transcription primer | RT on RNA, with shorter constant region | Separate RT reagent | DNA 5'→3' | 26 nt (derived count) | Degenerate notation `/5Phos/TCTCGGGAACGCTGAAGANNNNNNNN` | sample_specific | conditional | Culture and 50-µm hippocampal slice | reported | S02 supplementary p. 5; V01; length derived from string |
| primer_rca_15 | rolling circle amplification primer | RCA on cDNA circle | Separate reagent | DNA 5'→3' | 23 nt (derived count) | Exact notation `TCTTCACCCGGGGCAGCTGAA*G*T` | sample_specific | conditional | 15-µm hippocampal slice | reported | S02 supplementary p. 6; V01; length derived from string |
| primer_rca_50 | rolling circle amplification primer | RCA on cDNA circle | Separate reagent | DNA 5'→3' | 18 nt | Exact notation `TCTTCAGCGTTCCCGA*G*A` | sample_specific | conditional | Culture and 50-µm hippocampal slice | reported | S02 supplementary p. 6; V01 |

Required=conditional in untargeted rows means required within the named sample branch. Random N positions are not an exact synthesized single-sequence molecule. SOLiD offset sequencing primers are a distinct reagent family; original sequences route to S02 Table S1 (not yet directly transcribed here), and generalized specification routes to S12 Step 12.2. No fictitious 1-nt ligation component, promoter or cleavage-site absence record is added.

### Worked-example sequences, V04 (D01)

| Gene | RefSeq input accession | Minimum Tm (°C) | Logical barcode | Nonempty full probes |
|---|---|---|---|---|
| AKT1 | NM_001014432.1 | 55 | 3212010 | 16 |
| CDK7 | NM_001799.4 | 53 | 2102033 | 13 |
| KRT17 | NM_000422.3 | 52 | 0221001 | 12 |

Inputs are **reported**, S04 `example/example_input/example_genelist.xlsx`, Sheet1 A1:D3. Counts are **derived**, `ExampleOutputAllProbes.xlsx`, Sheet1 columns A/D (41 nonempty sequence cells; blank output padding is excluded). Output column B names explicitly contain `design_embedded_len4`. The AKT1 logical-to-padlock Illumina map is **derived**: 3212010→`TGCGACA`, using 0/1/2/3→A/C/G/T; the output D1 contains that exact seven-base region immediately after the backbone. This is not a reverse-complement operation [S04 `construct_Illumina_nucleotide_from_logical.m`, assembly mapping, workbook D1].

Full synthesized sequence examples remain in the pinned workbook rather than being reconstructed here. SHA-256: `4e3de1bdcd60e9991a5a31e35a47386dca0c8bf0e0d5e4ade4a7a05396f43071`. The homology FASTA sequences are already reverse complements of target segments; they must not be mistaken for target-strand RNA.

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | RNA; targeted padlocks recognize selected transcripts | V02-V04 | reported | S02 supplementary p. 18; S12 abstract |
| Target region | Contiguous 32-mer; two adjacent 16-nt arms | V02/V03 and V04 default | reported | S02 pp. 18-19; S12 Step 9.2 |
| Probe length | 32-nt homology; 74-nt total default, configurable | V04 | reported | S12 abstract p. 5, Steps 9.2/13.2 |
| GC constraints | 40-65% per arm | V02/V03 original rule; V04 default | reported | S02 p. 19; S12 Step 9.2 |
| Tm constraints | Visual cortex threshold initially 55 °C, lowered to minimum 52 °C if fewer than six probes/gene; hippocampus threshold 50 °C | V02 | reported | S02 supplementary p. 19; ambiguity below |
| Tm constraints | Per-gene minimum; example 55/53/52 °C; maximum arm separation 8 °C | V04 | reported | S12 Step 9.2; S04 input workbook A1:D3 |
| Tm implementation | Both mean DNA-DNA Tm estimates over oligoprop Tm(2:6) must exceed threshold; salt 0.3 M and primer 100 nM; arm difference ≤8 °C default | V04 code | inferred | S04 screen_homology_physical.m |
| Specificity filtering | Original rule excludes non-target transcript hits >12 nt spanning the junction by at least three nt on each side | V02 original brain rule, inherited with stated exceptions by V03 | reported | S02 supplementary pp. 19, 29 |
| Specificity filtering | E-value cutoff 0.05; candidate exclusion requires nonexempt Plus/Plus hit with ≥13 identities and query-index span reaching at least four positions each side of junction | V04 code default | inferred | S04 screen_homology_blast.m lines 79-118; homology parameters |
| Repeat masking | Reject ≥5 identical consecutive bases and sequences lacking any one of A/C/G/T; not a genome repeat-mask annotation rule | V04 | inferred | S04 screen_homology_complexity.m |
| Exon selection | No mandatory CDS-only filter in inspected loop; search spans transcript, including UTRs | V04 | inferred | S04 Homology.m, search initialization and gb_sequences comments |
| Isoform handling | Select a transcript sharing the most exons across variants (often shortest, not always); custom GenBank input supported | V04 recommendation | reported | S12 Step 2 |
| Sequence complexity | Hairpin base 7, loop 6; dimer length 10 used in oligoprop screening | V04 default | inferred | S04 screen_homology_physical.m |
| Spacing | Successful original search advances to start 5 nt past prior window end; failed window advances 1 nt | V02 original rule | reported | S02 supplementary p. 19 |
| Spacing | SPACER_LENGTH=0 by default; acceptance increments tail by HOMOLOGY_LENGTH+SPACER_LENGTH | V04 code | inferred | S04 Homology.m, acceptance loop |
| Overlap | Successive accepted candidate regions are nonoverlapping with the default zero extra gap; junction-overlap cutoff concerns BLAST span, not probe spacing | V04 code | inferred | S04 Homology.m acceptance loop; screen_homology_blast.m |
| Number of probes per target | Maximum 8 visual cortex, maximum 16 hippocampus; random choice among passing regions | V02 original | reported | S02 supplementary p. 20 |
| Number of probes per target | 9-16 per selected gene; three too-short targets excluded from 300-gene starting panel | V03 original | reported | S02 supplementary p. 29 |
| Number of probes per target | MAX_NUMBER_OF_PROBES=16 default; assembly takes passing regions from transcript 5' end toward 3' end | V04 code | inferred | S04 Assembly.m lines 446-460 |

The original supplementary text says exclusion for “melting temperature for both 16-mer halves > threshold” (S02 p. 19), then describes lowering to a minimum Tm of 52 °C when needed. The inequality is internally ambiguous; it is preserved here rather than silently reversed. V04 code explicitly accepts both means **above** the threshold. That implementation does not establish the historical original code's exact behavior.

BLAST query-span coverage and identity count are separate code tests; the code does not check identity individually at all eight flanking positions. Hits with a matching accession or the substring `(` followed by the target gene symbol are exempted. Prefix artifacts such as CD4/CD44 are therefore possible; optional ValidateHomologyRegions re-BLAST/manual inspection is documented [S04 screen_homology_blast.m; S12 Step 10; inferred code interpretation, reported validation caveat; V04].

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | Four logical rounds, each physical barcode region 7 nt | V02 | reported | S02 supplementary pp. 19-20 |
| Codebook | R4=(R1+R2+R3) mod 4; hippocampus expression-weighted color balancing, cortex unweighted | V02 | reported | S02 p. 19 |
| Error correction | Original brain decoding retained perfect four-base matches; checksum is not described as one-error correction | V02 | reported | S02 p. 23 |
| Codebook / Hamming distance | DNABarcodes length 7, distance 3; one-substitution correction | V03 | reported | S02 pp. 29, 31 |
| Number of barcodes | 326 length-7 entries, priority counts 47/11/258/10; minimum pairwise distance 3 | V04 supplied library | derived | S04 utilities/prioritizedBarcodeLibrary.csv; counts and pairwise symbol comparisons |
| Codebook | 58 priority-1/2 entries have a distance-2 four-symbol prefix | V04 supplied subset | reported | S12 Step 4 |
| Hamming distance / error correction | Distance 3 supports one symbol substitution correction or two-substitution detection; distance 2 supports one-substitution detection (distinct decoding uses) | V04 design | reported | S12 Step 3 |
| Round assignment | Four SBS rounds for original brain; V04 full library seven SBS rounds, embedded prefix four; SOLiD readout uses offset-primer/dibase specification rather than one physical base per ligation | V02 / V04 separately | reported | S02 p. 23; S12 Step 12; S14 Step 11 |
| Channel assignment | Default padlock map 0/1/2/3→A/C/G/T; corresponding template bases T/G/C/A; physical channel mapping must use the experiment/chemistry | V04 | reported | S04 Assembly.m Illumina mapping; S12 Fig. 1B |
| Fixed vs experiment-specific sets | Gene-to-logical-code assignment is experiment-specific; all padlocks for a gene share its code | V02-V04 | reported | S02 pp. 19-20; S12 Steps 6-7 |

V03 uses `create.dnabarcodes(7, dist=3, heuristic="ashlock", cores=12, filter.gc=FALSE, population=500, iterations=500)`; it states the SOLiD and Illumina barcodes have the same sequence, with one interpreted in color space and the other in base space [S02, supplementary p. 29; reported]. V04 has `embedded_len4` and alternative consecutive/sequential specifications (S12 Step 12.2; S04 Assembly.m). These are not interchangeable with the original brain implementation, whose SOLiD barcode has a documented color-map reversal for one ligation due to a design error that did not affect the acquired Illumina data [S02 p. 20; reported].

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | Barcode regions internal to targeted padlock, read on complementary RCA product | V02-V04 | reported | S02 Fig. S12; S12 Fig. 1A |
| Readout probe | SBS primer plus reversible-terminator incorporation; separate QC detection oligo | V04 Illumina | reported | S13 Step 6; S14 Steps 10-14 |
| Fluorophore | Four-channel SBS imaging; detection oligo 5'-Cy3 or substitute dye | V04 | reported | S14 Step 13.2; S13 Step 6 |
| Imaging round | Four original brain SBS rounds; seven for generalized full-length barcode, four for prefix subset | V02 / V04 | reported | S02 p. 23; S12 Steps 3-4; S14 Step 11 |
| Sequencing primer | 22-nt SBS oligo identical to padlock backbone; binds the RCA complement | V04 | reported | S14 materials/Step 10; S04 Assembly.m backbone |
| Decoding role | Original brain exact four-base matching; cancer permits one mismatch in seven-symbol barcode | V02 / V03 | reported | S02 pp. 23, 31 |
| SOLiD readout | Five offset sequencing primers and dibase ligation; dedicated original experimental sections available | V01 / V03 | reported | S02 pp. 7-8, 30-31; Table S1 |

## Primer / Linker / Adapter Components

See the separate primer, detection and spacer records above. Two AAA spacers are part of the default padlock; the detection oligo is a fluorescent hybridization reagent, not an anchor or primer. The RCA primer binds the circle, whereas the SBS primer and QC detection oligo bind its complementary amplicon. SOLiD primers are defined by offset specification, rather than assumed identical to the full SBS primer [S13 Steps 5-6; S14 Step 10; S12 Step 12.2; S04 assembly line 482; reported reagents and inferred pairing; V04].

## Amplification-Related Components

SplintR circularization precedes Phi29 RCA; BS(PEG)9 subsequently cross-links aminoallyl-containing amplicons. Generalized RCA solution contains 250 µM **each** dNTP plus 40 µM aadUTP. No measured substitution/incorporation fraction follows from these concentrations [S13 Steps 9-12; reported, V04]. V01 has sample-specific RT/RCA reagents and processing [S02 pp. 5-6; reported].

## Chemical Modifications

Padlock 5'-phosphate, RCA-primer phosphorothioate bonds, aadUTP and the QC fluorophore are recorded above [S12 Step 13.1; S13 Steps 5-6, 11-12; reported, V04]. BS(PEG)9 reacts with amines to cross-link amplicons; this is not a separate nucleotide segment. Gel passivation is tissue/hydrogel chemistry and belongs in protocols.md [S11 Step 18; S13 Step 12; reported, V04].

## Probe-Design Workflow

For V04: select transcript accessions/GenBank sequences; obtain a species-specific RefSeq mRNA BLAST database; assign logical barcodes; screen transcript windows for complexity, arm GC/Tm, secondary structure and off-target alignment; optionally validate hits manually; then convert barcodes and assemble phosphorylated padlocks for ordering/pooling [S12 Steps 1-17; reported; S04 scripts, inferred execution behavior]. The workflow does not establish the historical database release used in original experiments.

In the pinned code, FASTA headers are `probe_<n>_index_<i>_MinTm_<threshold>`. `i` is the one-based start index in the input transcript's forward sequence; default target interval is [i,i+31], closed, on that transcript strand. Exported DNA is its reverse complement. This is **inferred** from MATLAB indexing, `idx=1`, the `TempSeq` slice and `probeindex` expression in S04 Homology.m. Example AKT1 FASTA begins `probe_1_index_264_MinTm_55`, corresponding to target positions 264-295 by that interpretation. These are transcript positions, not genome coordinates; accession versions route through the input GenBank/genelist. Assembly selects the first passing regions up to the configured maximum, whereas original V02 reports random selection [S04 assembly lines 446-460; S02 p. 20].

## Official Probe-Design Software

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TargetedExSeqProbeGen | https://github.com/anusinha/TargetedExSeqProbeGen | Explicitly linked by developers' protocol | 2fe7748277df9d9f92ea4fe6a1d05e98796c527a | MATLAB; barcode-generation R script | Genelist, GenBank transcript sequences, species-specific RefSeq BLAST database | Homology FASTAs; full-probe, plate and pooling spreadsheets | Window screening and barcode conversion/assembly as above | MATLAB Bioinformatics Toolbox, legacy BLAST interfaces; R/DNABarcodes for barcode generation | Public revision retrieved 2026-10-15 | Unresolved; accessible revision alone does not establish maintenance | reported | S12 abstract/Steps 1, 5; S04 file tree and scripts; V04 |
| ExSeqProcessing (analysis, not design) | https://github.com/dgoodwin208/ExSeqProcessing | Linked by original publication | Archived v1.0.0; wiki inspected separately | MATLAB pipeline | Multiround images and sequence dictionary | Spatial base-called reads | Setup, color correction, normalization, registration, puncta extraction, base calling | See repository documentation; full environment not reconstructed | S07 release metadata and S10 wiki available 2026-10-15 | Unresolved | reported | S01 Data and materials availability; S07 release; S10 Stage Overview; V01-V03 analysis |

DNABarcodes and NCBI BLAST are third-party dependencies, not official ExSeq design products. S12 Steps 1 and 5 document their use; no independent maintenance claim is made here.

## Public Probe Sets

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | S04 pinned example inputs/outputs | V04 | Full sequences for 16 AKT1, 13 CDK7, 12 KRT17 probes; logical assignments and reverse-complement homology FASTAs | Retrieved 2026-10-15 | S04 example input/output files, Sheet1; reported output; counts derived |
| D02 | Original Tables S9/S10 | V02 and V03 | Genes/accessions/barcodes and probe sequences; original methods explicitly route both brain panels and cancer here | Public links and captions verified; workbook cells not yet examined, 2026-10-15 | S02 supplementary pp. 20, 29, 90; reported resource scope |
| D03 | Original Tables S7/S8 | V02/V03 context | Technology and yield comparisons, not synthesized probe lists | Captions examined 2026-10-15 | S02 supplementary pp. 86-89; reported |

Associated non-probe resources: D04/S06 raw Illumina sequencing archive; D05/S07 software release; D06/S09 visualization software; D07/S08 HTAPP extracted reads/analysis files. They are separately typed in benchmarks.md and are not probe-sequence ground truth. S02 captions identify S11 as ExFISH comparison data and S12-S14 as read-count tables (cortex, hippocampus, cancer respectively), not probe lists [S02 p. 90; reported].

## Version Differences

V01: sample-specific modified random-octamer RT and cDNA/RCA, SOLiD readout. V02: original brain padlocks, four-symbol checksum and four SBS rounds, experiment-specific probe-count/Tm rules. V03: original cancer padlocks, length-7 distance-3 encoding, 9-16 probes/gene and SOLiD ligation. V04: generalized protocol/code defaults, 326-code library with 58-entry prefix subset, configurable design and readout specifications [S02 pp. 5-8, 18-23, 29-31; S12 abstract/Steps 3-4, 12-13; reported; internal grouping inferred].

## Unresolved Information

- Workbook cells for original S9/S10 are not yet examined; public existence does not establish every field's completeness.
- Historical off-target database releases and complete original probe-selection reproducibility remain unresolved. GRCm38/mm10 and UCSC RefSeq downloaded 2018-08-14 are reported for V01 ex situ annotation (S02 p. 11), not established as targeted design inputs.
- The original Tm inequality ambiguity is preserved above. No universal exact RCA copy count or incorporation fraction is asserted.
- Current maintenance commitment of the design and analysis projects remains unresolved; no abandonment is inferred from age or social counts.
