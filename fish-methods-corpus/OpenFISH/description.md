# Method Overview

## Purpose

Map RNA in tissue and optionally relate the resulting cell identities to metabolomic signals (reported; S01 abstract; V01).

## Core Principle

The author implementation separates target-dependent probe A/B design from reusable padlock/readout design. Target-probe overhangs are selected from a bridge library; padlocks contain two RO sequences and bridge-complementary terminal regions (inferred from S04, openDesign/utils.py get_full_probes1/2 and 03_Padlock_join.ipynb PadJoin). V01 denotes the original 2026 publication; code observations refer specifically to the inspected author implementation snapshot, not necessarily every experiment in that paper.

## Overall Workflow

Design RO and bridge libraries; assemble padlocks; select a gene panel and assign round/color combinations; design and manually select target probes; acquire images; register, detect spots, segment and decode; optionally register MALDI data and transfer ion signals to cells (reported in S04 README notebook descriptions; inferred implementation order from notebooks 01–09). Wet-laboratory details remain incomplete in this draft.

## Classification

Controlled two-level classification (see `corpus_schema.md`, "Method classification"; `terminology.md`, "Method classification vocabulary"):

- `probe_architecture`: `padlock`.
- `readout_strategy`: `sequential_fish`.

Structural evidence: a target-binding probe pair (two 30-nt regions separated by a 2-nt target gap) brings a reusable padlock together; the padlock carries two RO sequences and bridge-complementary terminal regions, is circularized, and serves as the RCA template (see `probes.md` Canonical Architecture; Probe Components include `padlock` species and an RCA step). The padlock is the recognition/encoding molecule; the per-target probe pair supplies specificity. Identity is decoded across multiple imaging rounds, with each gene tagged with one round identity and a pair of five color channels (ten combinations per round); readout is `sequential_fish` because temporal rounds are essential to the codeword. Source: S04 notebooks 01_RO_design, 03_Padlock_join, 04_Panel_design, 06_Decode; S01 abstract (V01).

## Target Recognition Strategy

For the mRNA setting, the code reverse-complements two target regions; each is 30 nt, with a 2-nt gap on the target. The padlock is not itself the pair of direct RNA-binding arms (inferred; S04 utils.py findAllCandidates and get_full_probes1/2).

## Encoding Strategy

The notebook assigns genes to one round and a pair of five color channels, with ten combinations per round. This is an abstract imaging code, physically implemented through RO sequences in padlocks; it is not a universal contiguous gene barcode (reported/inferred; S04 03_Padlock_join cells 4–6, 04_Panel_design cells 17–26).

## Amplification Strategy

Rolling circle amplification is explicitly mentioned in the author README and RO-design notebook. Polymerase identity, actual amplification conditions, initiation and anchoring chemistry are unresolved here (reported/unavailable within inspected scope; S04 README and 01_RO_design cell 1).

## Imaging / Readout Strategy

The author decoder accepts projected CZI images, round/channel codebooks and calibration parameters. A 20X setting is supplied as an example; this is not an experimentally universal acquisition specification (reported; S04 06_Decode cells 0–4).

## Decoding Strategy

OpenDecode includes preprocessing, registration, spot detection, segmentation, decoding and matrix generation. The notebook calls Spotiflow and Cellpose. The decoding implementation uses spatial/channel features and a random-forest classifier (reported workflow; inferred algorithm; S04 06_Decode cells 0–6 and openDecode/decoding.py Decoding).

## Major Versions or Variants

V01 is the 2026 publication. S02 is a 2025 preprint with a different title. No separate experimental version is established from that title difference. The OpenDecode Version 2.0 label is software nomenclature, not evidence of an OpenFISH V2 chemistry (S04 06_Decode cell 0).

## Strengths

The inspected code separates reusable RO/padlock design from target sequence design and exposes editable parameters and ordering outputs (inferred; S04 notebooks 01–05).

## Limitations

The notebook retrieves live NCBI records unless sequences are supplied, leaves the RefSeq database release unspecified, uses stochastic selection and requires manual probe selection. Exact published panel reproduction cannot be assumed (inferred; S04 03–05, generator_multi.py get_seqs).

## Unresolved or Version-Specific Points

Unresolved: the publisher full methods and supplementary files could not be retrieved. The inspected author repositories do not supply a complete wet-laboratory protocol. This is an access/extraction gap, not evidence that the paper omits the information (S01; S04 README opening notice; V01). Full published probe lists, chemical modifications, ligation/RCA mechanism details and preprint-to-publication changes require further source access.
