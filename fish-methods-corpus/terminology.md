# Terminology Policy

## Author terminology and functional categories
Preserve the original authors' exact terms alongside every standardized functional label. Categories support retrieval; they do not establish equivalence between method-specific molecules. Possible source terms below are illustrative vocabulary, not claims that any listed method uses them. A mapping needs source context, method version, and a locator.

Acceptable: standard category target_binding_region; original term “homology arm”. Replacing “homology arm” throughout without preserving the original term is unacceptable. If no original term is given, record “No explicit term in source”; do not invent a quotation. Preserve spelling and capitalization in Original term. One component may have multiple evidenced roles; identify physical instances so length and membership are not double-counted.

## Method classification vocabulary

These controlled labels are governed by `corpus_schema.md` (`Method classification` section). This file defines only their prose definitions and cross-references; do not redefine or extend the controlled set here. Record labels via the `classification.probe_architecture` and `classification.readout_strategy` keys in `metadata.yaml`. Plain-language rationale belongs in `description.md` under `## Classification`.

### Primary classification: probe_architecture

| Label | Functional meaning | Typical original terms (illustrative) |
|---|---|---|
| linear | Contiguous oligo(s) carrying recognition and encoding/readout features inline; no mandatory in situ circularization of a recognition oligo. | primary probe; encoding probe; oligonucleotide probe set. |
| padlock | A recognition oligo whose 5' and 3' ends are brought together on a target or splint and ligated to form a circular template (typically amplified by RCA). | padlock probe; reverse padlock; SNAIL probe pair with padlock circularization. |
| crispr | Sequence specificity conferred by a Cas-family protein programmed by a guide RNA; readout is attached to the protein or to the guide. | dCas9/sgRNA; crRNA:tracrRNA; fRNP. |
| nanostructure | DNA or other nucleic-acid nanostructures used as signal amplifiers attached to a primary recognition element. | TDDN (tetrahedral DNA dendrimer); branched DNA amplifier. |

These labels describe probe architecture, not imaging modality or chemistry of readout.

### Secondary classification: readout_strategy

| Label | Functional meaning | Typical original terms (illustrative) |
|---|---|---|
| sequential_fish | Identity is decoded from a temporal pattern of fluorescent-oligo hybridization, imaging, and stripping cycles. | sequential FISH; MERFISH; seqFISH; SEDAL rounds. |
| in_situ_sequencing | Identity is decoded by an in situ sequencing reaction (sequencing-by-ligation, SOLiD, SBS) on amplified templates. | SEDAL sequencing-by-ligation; SOLiD; in situ sequencing. |
| fluorescent_barcode | Identity is encoded by intensity ratios, color combinations, or co-localization patterns of fluorophores within one imaging round. | intensity barcode; PRISM color-intensity codeword. |
| live_imaging | The assay is designed for live-cell imaging, without requiring fixation or destructive in situ chemistry beyond reagent delivery. | live-cell FISH; live-cell CRISPR imaging. |

These labels describe readout mechanisms, not imaging modalities in general. A method that uses confocal microscopy to read a sequential FISH signal is still `sequential_fish`. A method that adds fixed-cell compatibility to an otherwise `in_situ_sequencing` workflow does not automatically gain `live_imaging`.

### Distinguishing readout from imaging mode

Do not infer `live_imaging` from "uses fluorescence microscopy". A method that is fixed-endpoint and image-based is not `live_imaging` unless it is documented to support live-cell time-lapse acquisition. Conversely, a live-cell method that uses a CRISPR-based label and one round of imaging is `crispr` + `live_imaging` + (typically) `fluorescent_barcode`, not `sequential_fish`.

## Functional vocabulary

| Standard corpus concept | Possible source-specific terms | Functional meaning | Usage rule | Common ambiguity |
|---|---|---|---|---|
| probe | probe; oligonucleotide | A molecular reagent used for recognition, encoding, detection, or another documented role. | Specify molecular species and role; retain author term. | May mean one oligo, a complex, or a whole set. |
| primary_probe | primary probe | An author-designated initial probe layer. | Use only with source definition; record what it binds. | Does not necessarily bind the biological target directly. |
| encoding_probe | encoding probe | A probe participating in assignment or presentation of identity information. | Describe physical components and encoding role separately. | Not universally synonymous with primary probe. |
| readout_probe | readout probe | A physical reagent that detects or interrogates a readout feature. | Identify binding partner and signal role. | Different from the readout sequence or process. |
| detection_probe | detection probe | A probe used in a documented detection step. | Retain author label and explain detection mechanism. | May overlap functionally with readout or reporter probe. |
| reporter_probe | reporter probe | A probe carrying or enabling a reporting signal. | State signal carrier and binding partner. | Not necessarily directly labeled. |
| target_binding_region | target-binding region; homology arm | Probe segment that hybridizes directly to the stated target. | Identify target molecule, position, and orientation. | Not the target-side sequence; may comprise several arms. |
| target_binding_arm | target-binding arm; homology arm | A discrete target-hybridizing segment described as an arm. | Preserve arm identity and relationship to other segments. | An arm is not automatically a padlock arm. |
| padlock_arm | padlock arm; homology arm | A target-binding arm of a documented padlock architecture. | Use only when padlock context is supported. | Do not infer topology or ligation from the word arm. |
| barcode | barcode; barcode sequence | An identifier encoded by a physical sequence or a documented abstract pattern. | State sequence versus abstract pattern and what identity it encodes. | Can refer to sample, target, molecule, spatial identity, or codeword. |
| codeword | codeword; barcode | An abstract assignment of symbols or states used to represent identity. | Record symbol alphabet and ordering, rounds/channels where reported. | Not necessarily a contiguous nucleotide sequence. |
| readout_sequence | readout sequence; readout handle | A sequence used as a recognition feature during readout. | Specify which molecule carries it and its binding partner. | May be conflated with barcode or readout probe. |
| primer | primer | A molecular reagent or sequence used to initiate a documented synthesis reaction. | Record role, molecule, binding site, and direction. | The oligo, its binding site, and its functional region differ. |
| sequencing_primer | sequencing primer | A primer used in a sequencing reaction. | Specify reaction and bound template. | Not every readout probe is a sequencing primer. |
| amplification_primer | amplification primer | A primer used in an amplification reaction. | Identify amplification reaction and binding site. | A primer-binding sequence is not the primer molecule. |
| linker | linker | A connector with a reported linking role. | State physical nature and position. | May be nucleotide or non-nucleotide; not automatically inert. |
| spacer | spacer | A segment or group with a reported spacing role. | Preserve chemistry, length, and role if reported. | Not universally interchangeable with linker. |
| adapter | adapter; adaptor | A component connecting a molecule to another process or recognition system. | State its documented interface or function. | May contain primers or other roles, not one universal sequence type. |
| anchor_sequence | anchor sequence; anchor | A sequence with a documented anchoring or reference-binding role. | Specify what is anchored and to what. | Anchor may also refer to a reagent or non-sequence feature. |
| ligation_site | ligation site; ligation junction | The documented position of a joining reaction. | Record partners, orientation, and junction context. | A site is not necessarily a distinct oligo component. |
| ligation_sequence | ligation sequence | A sequence participating in a documented ligation role. | Distinguish the sequence from the junction and ligase. | Do not infer a conserved motif. |
| amplification_sequence | amplification sequence; amplification handle | A sequence with a documented amplification-related role. | State whether template, binding site, or another function. | Does not by itself identify amplification chemistry. |
| promoter | promoter | A sequence documented to direct transcription initiation. | Record recognized system and orientation if reported. | Do not classify a primer site as a promoter by appearance. |
| fluorophore | fluorophore; dye | A fluorescent signal-bearing chemical entity. | Record chemical identity and attachment position when reported. | A channel label is not an exact dye identity. |
| chemical_modification | chemical modification | A documented chemical alteration or conjugated group. | Preserve exact notation and location. | Do not treat a modification as a nucleotide by default. |
| cleavage_site | cleavage site | A site or group participating in documented cleavage. | Record chemistry and position only as evidenced. | Not automatically a sequence motif. |
| other_functional_sequence | author-specific functional sequence | An evidenced sequence role outside established categories. | Preserve author term and define role; justify new categories. | Do not use as a substitute for unresolved function. |
| target_sequence | target sequence | The sequence of the biological or intermediate molecule being recognized. | Identify molecule, reference, strand, and scope. | Can mean target-side sequence or a selected design input; distinguish both. |
| target_region | target region | A region selected or recognized on a target molecule. | Record reference coordinates/conventions if available. | Different from the probe's target-binding segment. |
| probe_set | probe set; probe library | A documented collection of probe molecules or designs. | State membership, target scope, release, and completeness. | A subset or example is not automatically the full set. |
| codebook | codebook | A mapping between identities and codewords. | Preserve version, alphabet, order, and assignments. | A barcode sequence list alone may not define assignments. |
| imaging_round | imaging round; imaging cycle | A documented iteration of an imaging/readout procedure. | Use the author's round boundaries and numbering. | Not necessarily a sequencing cycle or one exposure. |
| sequencing_round | sequencing round; sequencing cycle | An author-defined sequencing iteration. | Preserve original term and chemistry-specific granularity. | Round and cycle are not universally equivalent. |

## Ambiguous terms requiring qualification
- “barcode”: specify whether a nucleotide sequence, abstract codeword, or assignment; record identity scope and ordering.
- “readout”: distinguish sequence, physical probe, signal acquisition process, and decoding observation.
- “primer”: distinguish the oligo from its binding site and distinguish sequencing from amplification roles.
- “probe”: identify molecular species, functional layer, target/binding partner, and whether the term denotes one molecule or a collection.

These are contextual warnings, not assertions that all methods contain those entities.

## Sequence vocabulary

| Standard corpus concept | Possible source-specific terms | Functional meaning | Usage rule | Common ambiguity |
|---|---|---|---|---|
| exact_sequence | sequence; oligo sequence | A fully specified nucleotide string explicitly provided in a textual/table/file source. | Preserve molecule type, orientation, modifications, notation, source locator, and version. | A diagram or inferred string is not an exact published sequence. |
| sequence_example | example sequence | A source-provided instance illustrating a design. | Identify as an example and its experiment scope, even if its letters are exact. | One example does not define a universal or complete probe set. |
| consensus | consensus sequence | A source-defined summary across sequences. | Preserve notation and derivation as reported. | A consensus is not necessarily an individual molecule. |
| motif | motif | A reported pattern or functional sequence feature. | Record constraints and context. | A motif does not determine a full oligo. |
| schematic_sequence | schematic; sequence diagram | A symbolic representation of component organization. | Use labeled placeholders and mark as schematic. | A schematic sequence is not an exact nucleotide sequence. |
| degenerate_sequence | degenerate sequence | A sequence with symbols representing multiple nucleotide choices. | Preserve source symbols and document their defined alphabet. | Do not expand choices into an asserted exact molecule. |
| reverse_complement | reverse complement | A strand transformation of an identified sequence. | Retain the original, identify input alphabet, and label computed output derived with the operation. | Complement and reverse complement differ; do not assume the displayed strand. |
| variable_region | variable region | A region whose sequence varies within a stated design scope. | Specify the scope and selection rule if reported. | Variable does not mean random or unknown. |
| fixed_region | fixed region; constant region | A region constant within a specified version/library/experiment. | State the scope of constancy. | Fixed in one set does not mean fixed across versions. |

## Representation and evidence
Sequence representation and sequence class are independent. probes.md uses sequence classes fixed, variable, barcode_library, sample_specific, version_specific, unresolved; none establishes exactness. A sequence example may be exact at the instance level while sample_specific at the class level.

Evidence states are reported, derived, inferred, unavailable, exactly as defined in corpus_schema.md. Functional labels and source terms do not establish evidence on their own. Pending inspection is Not yet examined.; missing reported information and genuine inapplicability remain distinct. No fabricated nucleotide strings are permitted as examples.

