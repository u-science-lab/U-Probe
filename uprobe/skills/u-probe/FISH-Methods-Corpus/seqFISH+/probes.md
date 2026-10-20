# Probe Architecture

## Overview

Scope: V01, original RNA seqFISH+ (2019). Primary ssDNA probes recognize mRNA through a 28-nt region and present four readout-binding overhang sites corresponding to barcoding rounds I–IV. A gene is targeted by 24 primary probes; fluorescent 15-nt DNA readouts bind the overhangs. The four-symbol barcode is an abstract temporal codeword, not a contiguous four-base sequence. [S01, Fig. 1a; Methods, Primary-probe design, Readout-probe design and seqFISH+ encoding strategy; reported]

No exact target-specific primary probe or readout nucleotide string is transcribed here: none was resolved in the examined V01 resources. The codebook in S02 specifies gene-to-pseudocolour assignments, not oligo sequences. [S02, all sheets; S03, all workbook headers/content types; S06, recursive tree and sourcedata.zip inventory; unavailable within this scope]

## Canonical Architecture

### Mature primary probe: linear ssDNA, hybridized to RNA

Schematic only; the explicit 5′→3′ record must retain unresolved order:

```text
5′ — [Unresolved ordered arrangement of: one 28-nt RNA-binding region,
      four 15-nt readout-binding sites for rounds I–IV,
      and unresolved residual/linking sequence] — 3′
```

Reported: the primary contains RNA-binding sequence and barcode-bearing overhangs; the readout-binding sequences are reverse complements of the designed readout probes. Inferred from Fig. 1a: the four sites are depicted as two flanks around the target-binding segment. The figure does not mark strand polarity or establish a synthesized 5′→3′ order; do not assign I–II to the 5′ flank solely from the drawing. Binding to RNA and readouts involves complementary, antiparallel strands; no complete strand-resolved oligo map was recovered. [S01, Fig. 1a and Methods, Readout-probe design; V01]

The final manufactured primary is reported as approximately 93 nt. The sum 28 + 4 × 15 = 88 nt is derived from the component lengths; this does not identify the remaining sequence or establish five T spacers. Do not turn the approximate total into an exact 93-nt design specification. [S01, Methods, Primary-probe construction and Fig. 1a; reported total, derived component sum; V01]

### Fluorescent readout: separate linear DNA oligo

```text
5′ — [15-nt readout DNA sequence] — 3′
       dye coupled through the precursor's 5′ amine
```

Readouts are ordered as 5′-amine-modified oligos and conjugated to Alexa Fluor 647, Cy3B or Alexa Fluor 488 using the reported synthesis procedure. The dye is a chemical conjugate, not a nucleotide block. Its binding site on the primary is the reverse complement; exact strings are unresolved. [S01, Methods, Readout-probe design/synthesis; reported; V01]

### Anchoring reagent and manufacturing intermediates

LNA polyT15/polyT30 oligos bind poly(A) RNA and participate in hydrogel anchoring; they are separate reagents, not primary-probe tails. S01 Extended Data Fig. 1 shows acrydite functionalization; the inspected S05 mirror explicitly identifies Acrydite/polyT LNA probes, with official-body equivalence unresolved. Their modification placement and complete synthesis notation are not resolved here. [S01, Extended Data Fig. 1 and Primary-probe construction; S05 mirror, Primary probes hybridization; reported (S01); inferred (mirror-only identity); V01]

Manufacturing includes PCR products, an IVT RNA intermediate and RT-derived ssDNA followed by primer removal. Do not equate array-pool oligos, primers, RNA intermediates and final hybridization probes. Their exact 5′→3′ structures remain unavailable. [S01, Methods, Primary-probe construction; reported workflow/unavailable exact structures; V01]

## Probe Components

All rows concern V01. Sequence entries distinguish variable regions or schematics from exact strings. Functional identifiers implement the corpus terminology without replacing author terms.

| Component | Original term | Function | Position | Orientation | Length | Sequence | Sequence class | Required | Design rule | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|
| target_binding_region | gene-specific primary probes; target sequences | Direct RNA recognition | Internal placement inferred from Fig. 1a; exact index unresolved | Complementary to RNA; complete synthesized strand map unresolved | 28 nt | Variable region; exact strings unavailable | variable | yes | Coding exons first; UTR exons if needed | reported (roles/lengths); inferred (position) | S01, Primary-probe design and Fig. 1a; V01 |
| readout_sequence_I–IV | overhang sequences (I–IV) | Present four temporal-codeword readout sites | Flanking placement inferred from Fig. 1a; exact 5′→3′ ordering unresolved | Reverse complement of corresponding readout oligo | 15 nt per site | Library-selected variable regions; strings unavailable | barcode_library | yes | Assigned by gene codeword; do not concatenate integer symbols as DNA | reported (roles/lengths); inferred (position) | S01, Fig. 1a, Readout-probe design and encoding strategy; V01 |
| readout_sequence_detector | readout probes | Bind primary overhangs | Separate oligo | 5′→3′ DNA; complementary binding partner | 15 nt | Library-selected variable region; strings unavailable | barcode_library | yes | Random candidates; GC 40–60%; specificity/cross-hybridization filters | reported | S01, Readout-probe design/synthesis; V01 |
| chemical_modification_readout | 5′-amine modified | Dye conjugation | Readout 5′ terminus | Chemical group, not strand sequence | Not applicable | Chemical modification; no nucleotide string | sample_specific | yes | Couple dye as described below | reported | S01, Readout-probe synthesis; V01 |
| fluorophore_readout | Alexa Fluor 647; Cy3B; Alexa Fluor 488 | Fluorescent detection | Conjugated through readout 5′ amine precursor | Not applicable to dye | Not applicable | Non-nucleotide conjugate | sample_specific | yes | Gene barcode stays in one channel | reported | S01, Readout-probe synthesis and encoding strategy; V01 |
| primer_forward | forward primer containing a uracil nucleotide | RT during manufacture; removable primer | Manufacturing intermediate; removed before use | Exact sequence/strand map unresolved | Not reported | Exact sequence unavailable | unresolved | conditional | USER cleavage next to uracil | reported | S01, Primary-probe construction; V01 |
| primer_reverse | reverse complement of T7 sequences | Primer in manufacture | Manufacturing construct; final retained extent unresolved | Relation to T7 stated; exact oligo string unresolved | Not reported | Do not substitute a generic T7 sequence | unresolved | conditional | Use source-specific primer; sequence not recovered | reported | S01, Primary-probe construction; V01 |
| anchor_sequence_polyT | LNA polyT15 and LNA polyT30 | RNA/hydrogel anchoring reagent | Separate oligos | Bind RNA poly(A); modification orientation unresolved | 15 or 30 nt, as named | Homopolymer class; complete modified oligo notation unresolved | fixed | conditional | S01 mixes both; S05 mirror permits one or both | reported (S01); inferred (mirror-only option) | S01, Primary-probe construction and Extended Data Fig. 1; S05, Primary probes hybridization; V01 |

## Target-Binding Design

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Target type | mRNA | V01, 10,000-gene RNA experiment | reported | S01, abstract and encoding strategy |
| Target region | Coding exons first; coding plus UTR exons when insufficient | V01 design | reported | S01, Primary-probe design |
| Probe length | 28-nt target-binding region; final primary approximately 93 nt | V01; different molecular lengths | reported | S01, Primary-probe design/construction |
| GC constraints | 45–65% | V01 target sequences | reported | S01, Primary-probe design |
| Tm constraints | No numeric requirement recovered | V01 examined methods/protocol | unavailable | S01, Primary-probe design; S05 |
| Specificity filtering | Local BLAST versus mouse transcriptome; other-gene 15-nt matches treated as off targets; reject expected summed off-target expression exceeding 10,000 FPKM using ENCODE RNA-seq-derived table | V01; do not replace expression threshold with automatic rejection of every 15-nt hit | reported | S01, Primary-probe design |
| Cross-hybridization | Probe-versus-probe BLAST; matches ≥17 nt removed by dropping matched probe from larger probe set | V01 | reported | S01, Primary-probe design |
| Repeat masking | UCSC masked genome and annotation used; exact masking release/options unspecified | V01 | reported | S01, Primary-probe design |
| Exon selection | Coding exons preferentially; UTR fallback | V01 | reported | S01, Primary-probe design |
| Isoform handling | No consensus percentage, transcript-selection algorithm or isoform release recovered | V01 examined resources | unavailable | S01, Primary-probe design; S02–S04; S06 inventory |
| Sequence complexity | Reject ≥5 consecutive identical bases | V01 | reported | S01, Primary-probe design |
| Spacing | No numeric spacing rule recovered | V01 | unavailable | S01, Primary-probe design; S05 |
| Overlap | Explicit overlap-selection rule not recovered | V01 | unavailable | S01, Primary-probe design; S05 |
| Number of probes per target | Genes with <24 accepted candidates dropped; 24 primaries per gene used | V01 10,000-gene experiment | reported | S01, Primary-probe design and encoding strategy |

The reference is mouse and UCSC-derived, but neither assembly nor annotation release is specified in these examined sections. ENCODE inputs, BLAST version/options and exact candidate enumeration/tie-breaking are unresolved. Do not import spacing, isoform-consensus or primer defaults from third-party software. [S01, Primary-probe design; S11, pipeline defaults; unavailable V01 implementation details]

## Barcode System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Barcode length | Four base-20 symbols, one per barcoding round; not four nt | V01 | reported | S01, encoding strategy; S02, columns C–F |
| Number of barcodes | 10,000 assigned genes; stated capacity 8,000 per channel ×3 =24,000 with redundancy | V01 | reported | S01, encoding strategy; S02 |
| Codebook | Gene names in column A; blank B; round I–IV in C–F; symbols 1–20 | V01 published panel | reported | S02, three sheets, row 1 and data rows |
| Channel membership | 647 nm: 3,334; 561 nm: 3,333; 488 nm: 3,333 | V01 actual table | derived | S02, count nonempty gene rows in each sheet; evidence/inspection.json |
| Hamming distance | Minimum symbol Hamming distance 2 within each published channel | V01 published subset, not an author-quoted constant | derived | S02, C–F; unique three-position projections and repeated two-position projections; evidence/inspect_resources.py |
| Error correction | Fourth round provides redundancy; calling uses distance-one matching with ambiguity rejection and ≥3/4 seed agreement | V01 original decoder | reported | S01, encoding strategy and Barcode calling |
| Round assignment | 20 serial readout hybridizations per barcoding round ×4 =80; hybridization 81 repeats 1 | V01 | reported | S01, Fig. 1 and seqFISH+ imaging |
| Channel assignment | Each gene decoded within one channel; Alexa Fluor 647/Cy3B/Alexa Fluor 488 | V01 | reported | S01, encoding strategy and Readout-probe synthesis; S02 sheet names |
| Fixed vs experiment-specific barcode sets | Published assignments fixed for this panel; not an invariant library for all seqFISH+ studies | V01 | reported | S02, gene-specific assignments; S01, experiment scope |

Derivation of minimum distance: within each sheet all four-position codewords and every three-position projection are unique, so no pair differs in zero or one symbols. Every two-position projection has 400 unique values for >3,300 rows, so some pair agrees at two positions. Together these establish minimum symbol distance 2. This is not binary one-hot distance and does not guarantee correction of every substitution. [S02; derived; V01]

## Readout System

| Field | Reported value / rule | Version / experiment scope | Evidence status | Source / locator |
|---|---|---|---|---|
| Readout sequence | 15 nt; random A/T/G/C candidates; GC 40–60%; BLAST versus mouse transcriptome; remove pairs with ten contiguously matching bases | V01 | reported | S01, Readout-probe design |
| Readout probe | Fluorescent DNA oligo binds its reverse-complement site on primary | V01 | reported | S01, Readout-probe design/synthesis |
| Fluorophore | Alexa Fluor 647, Cy3B, Alexa Fluor 488; 5′ amine precursors | V01 | reported | S01, Readout-probe synthesis |
| Imaging round | Three distinct readouts per serial hybridization; 80 encoding hybridizations; repeat first in 81 | V01 | reported | S01, seqFISH+ imaging |
| Readout library size | 240 unique readouts stated in the S05 mirror (official-body equivalence unresolved); independently 20 ×4 ×3 =240 slots | V01 | inferred (mirror statement); derived (slot count) | S05, Readouts hybridization; S01, imaging/encoding strategy |
| Sequencing primer | Not applicable: hybridization imaging does not sequence nucleotides | V01 | reported | S01, workflow and seqFISH+ imaging |
| Decoding role | Pseudocolour within each round forms one symbol; no complete nucleotide-to-slot mapping recovered | V01 | reported | S01, Fig. 1; S02; exact mapping unavailable in S02/S06 |

## Primer / Linker / Adapter Components

The manufacturing forward primer includes uracil; USER treatment removes the primer. The reverse primer is described relative to T7 sequence. Exact primer strings, primer lengths, linker strings and the full array-pool construct were not recovered. No sequence is inferred from common T7 promoters or from later DNA seqFISH+ primer lists. [S01, Primary-probe construction; unavailable sequence details; V01]

## Amplification-Related Components

PCR/IVT/RT produce probes outside the specimen. No padlock, ligation junction, RCA primer or HCR initiator belongs to the reported V01 RNA workflow. These headings are not evidence that such components exist. [S01, Primary-probe construction and experimental methods; reported workflow; V01]

## Chemical Modifications

Readouts use 5′-amine precursors for dye coupling. Separate LNA polyT anchoring reagents and post-hybridization Label-IT (cells) or label-X followed by re-embedding (tissue) support retention. Primary probes are not assumed to have an ordered 5′ acrydite solely because the workflow uses acrydite anchoring oligos. [S01, Readout-probe synthesis, cell/tissue methods, Extended Data Fig. 1; reported; V01]

## Probe-Design Workflow

Reported V01 order: retrieve mouse gene sequences using UCSC masked genome/annotation; enumerate 28-nt exon candidates, coding first and UTR fallback; apply GC/homopolymer filters; assess transcriptome off targets and expression-weighted burden; remove probe–probe matches; retain genes with at least 24 candidates and use 24 probes/gene. Independently generate 15-nt readout candidates, filter GC and sequence cross-reactivity, and append their reverse complements to primary probes according to the gene codeword. Physical append order is unresolved. [S01, Primary-probe design, Readout-probe design and encoding strategy; reported; V01]

## Official Probe-Design Software

No positive evidence of a publicly released, original-author V01 probe-design pipeline was found in the examined paper, protocol, CaiGroup repository README/tree or source-data archive. Availability remains unresolved (`null`), not proven false. The linked official software is an image-processing/decoding pipeline. [S01, Code availability; S05; S06; unavailable public design implementation within examined scope; V01]

| Software name | Source/repository | Official status | Version | Programming language | Input | Output | Design logic | Dependencies | Availability | Maintenance status | Evidence status | Source / locator |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| seqFISH-PLUS | https://github.com/CaiGroup/seqFISH-PLUS | Official, linked by original paper; analysis software | Inspected tree ee6c416006fe09262526eb364c6f66edd35adc88; no release assigned | MATLAB | Images, PSF, ROIs, thresholds | Decoded gene calls and locations | Not a documented oligo-design implementation | MATLAB; radialcenter.m; ImageJ ROIs | Public 2026-10-15 | Maintenance not established from this inspection | reported | S01, Code availability; S06, README/recursive tree; V01 |
| SeqFISH+ Probe Designer, Oligo Designer Toolsuite | https://oligo-designer-toolsuite.readthedocs.io/en/latest/_pipelines/seqfishplus_probe_designer.html | Third-party implementation; not original-author software | Mutable latest documentation; package release not pinned | Python | Target/reference FASTA and gene identifiers | Designed probes, readouts and codebook | Separate implementation with configurable defaults | See project documentation; not tested here | Documentation public 2026-10-15 | Not assessed | reported | S11, Pipeline Description and references; contemporary third-party implementation, not V01 evidence |

## Public Probe Sets

`public_probe_sets_available: true` means a public barcode-assignment resource is verified, not that complete synthesized probe sequences are available.

| Dataset ID | Release / accession | Method version | Contents | Availability / access date | Source / locator |
|---|---|---|---|---|---|
| D01 | Supplementary Table 1; MOESM2 | V01 | 10,000 gene names, channel membership and four-symbol codewords; no target-binding or full oligo strings | Downloaded and all sheets inspected, 2026-10-15 | S02, three Codebook sheets |
| D02 | CaiGroup/seqFISH-PLUS sourcedata.zip; Zenodo 2669683 v1 | V01 | Processed counts/locations, segmentation resources; not a sequence-level probe library | Repository archive downloaded; Zenodo landing/file metadata inspected, 2026-10-15 | S06, archive inventory; S07, Description and Files |

Supplementary Tables 2–4 contain cell clusters, RNA localization and ligand/receptor analyses, respectively, rather than probe sequences. DNA seqFISH+ tissue probe files are separate extension resources and do not establish V01 RNA sequence availability. [S03, all workbook sheets; S10, README/probes listing; reported resource scope]

## Version Differences

V01 cell/tissue preparations share the published 10,000-gene RNA panel. The smFISH comparison uses two readout sites per primary, unlike the canonical four-round panel. DNA seqFISH+ targets genomic loci and requires separate architecture/protocol treatment; its sequences are not substituted here. [S01, main text and Methods, smFISH; S09, Fig. 1; reported]

## Unresolved Information

Exact primary/readout/primer strings; full component polarity/order; residual sequence/spacers underlying the approximately 93-nt total; array-oligo layout; exact UCSC reference release and ENCODE inputs; BLAST settings; overlap/spacing, isoform rules and tie-breaking; an official public design implementation; and complete original synthesized probe sets remain unresolved in S01–S08. A current third-party design recipe cannot fill these gaps. See research_manifest.md for access and inspection boundaries.
