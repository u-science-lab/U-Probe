# Benchmark Resources

## Benchmark Summary

Public candidates are evidenced, but none is accepted as Grade A or as strict target-selection ground truth. Highest scoped grade is **B for D03 synthetic amplifier/readout assembly only**; it does not describe RNA targeting or end-to-end probe design. D01/D02/D04 remain C under currently inspected evidence; D05 is unassessed. Assessment date: 2026-10-15. Grades are reasoned research assessments, not validation outcomes.

- A — Fully reproducible: Enough information exists to reconstruct the relevant probe-design output and directly compare another probe-design pipeline against it.
- B — Mostly reproducible: Most relevant information is available, but minor information remains missing.
- C — Architecture-level only: Useful for understanding probe structure or design principles, but unsuitable as strict probe-design ground truth.
- D — Insufficient information: The design cannot be reliably reconstructed or compared.

A published sequence list is not automatically Grade A. Reference identity, target regions, strand, design parameters and historical output scope matter. Gene-expression matrices, measured RNA locations, and RNA-seq data are not synthesized-probe ground truth.

## Candidate Datasets

### D01 — Original 140-/1001-gene design resources

| Field | Value |
|---|---|
| Dataset ID | D01 |
| Dataset name | Original 140-/1001-gene design resources |
| Method version | V01 |
| Species | Homo sapiens |
| Cell / tissue / sample | IMR90 human fibroblasts |
| Genome build | Exact design FASTA/build not independently established; do not infer from year |
| Transcriptome / annotation version | GENCODE v18 reported; Cufflinks v2.1 used for isoform abundance |
| Target genes | 130 biological RNAs + 10 controls; 985 + 16 controls. Gene lists in Tables S1/S3 not parsed |
| Target regions | Most abundant isoform; 30-nt sites; 1-kb computational fragments |
| Genomic/transcript coordinates | Assembly FASTA checksum, transcript versions and coordinate origin/interval conventions not recovered |
| Target-binding sequences | Reported within Table S5 template designs; not directly parsed |
| Full probe sequences | Table S5 is explicitly template sequence, not automatically final synthesized ssDNA |
| Barcode assignments | Tables S1/S3 identified; two shuffled MHD4 codebooks; binary files inaccessible |
| Readout sequences | 16 exact 30-nt readout-probe sequences inspected in S01 |
| Primer/linker sequences | Indexed primers and common RT primer inspected in S01; member-specific junctions unresolved |
| Design software | OligoArray2.1 + BLAST+ + custom scripts; exact historical environment unresolved |
| Design parameters | S01 thermodynamic, homology and probe-number rules recorded in probes.md |
| Source | S01 Probe design and Probe synthesis; S02 Tables S1/S3/S5 |
| Availability | Article inspected; supplemental endpoints returned challenge/error pages on 2026-10-15 |
| Missing information | Unparsed designs/codebooks; reference snapshot; full final-probe boundaries; historical run settings |
| Benchmark grade | C |
| Grade justification | Current inspected evidence supports architecture and rules, not reconstruction of the published design output. This grade describes the usable evidence in this draft; it does not assert that the published supplements are inherently incomplete. |
| Comparison scope | Architecture/rule comparison only; no strict U-Probe target-selection ground truth |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource existence; inferred benchmark assessment |

### D02 — High-throughput 130-RNA encoding library

| Field | Value |
|---|---|
| Dataset ID | D02 |
| Dataset name | High-throughput 130-RNA encoding library |
| Method version | V02 |
| Species | Homo sapiens |
| Cell / tissue / sample | U-2 OS |
| Genome build | hg38 / GRCh38, explicitly reported |
| Transcriptome / annotation version | cDNA annotation release unresolved; release-84 ncRNA URL explicitly given, not generalized to all inputs |
| Target genes | 130 RNAs, including 85 from prior library and 45 sampled in expression range; 10 blank codewords |
| Target regions | 30-nt sites; nonoverlapping selection from each isoform |
| Genomic/transcript coordinates | Target-coordinate conventions not recovered from inaccessible workbook |
| Target-binding sequences | Included in encoding probes per SI; not independently parsed |
| Full probe sequences | Dataset S1 reports ordered template oligos; final product termini unverified |
| Barcode assignments | 16-bit MHD4 weight 4; gene-to-word table not parsed |
| Readout sequences | S04 Table S1 contains 16 exact 20-nt probes with dye assignments |
| Primer/linker sequences | 20-nt primer design and A spacers reported; exact per-member strings not parsed |
| Design software | Author MATLAB pipeline; S06 pinned example inspected, not run |
| Design parameters | GC 43–63%; Tm 66–76 °C; specificity indices 75–100%; ncRNA screen |
| Source | S04 SI Encoding Probe Design/Table S1; S13 |
| Availability | Article and code accessible; pnas.1612826113.sd01.xlsx binary not parsed 2026-10-15 |
| Missing information | Complete workbook; target accessions and coordinates; historical seeds/environment |
| Benchmark grade | C |
| Grade justification | Sequence-level target-selection reproduction has not been established. Architecture, reference build and design constraints are useful but not strict ground truth. |
| Comparison scope | Architecture/design-rule comparisons; not full output matching |
| Assessment date | 2026-10-15 |
| Evidence status | reported resource and parameter facts; inferred grade |

### D03 — 16-bit 5×5 bDNA amplifier and readout library

| Field | Value |
|---|---|
| Dataset ID | D03 |
| Dataset name | 16-bit 5×5 bDNA amplifier and readout library |
| Method version | V05 |
| Species | Homo sapiens experimental context; amplifier strings are synthetic |
| Cell / tissue / sample | U-2 OS experiments; 130-RNA library and reduced-probe tests |
| Genome build | Not required for the restricted synthetic-amplifier assembly comparison; underlying RNA design reference not established here |
| Transcriptome / annotation version | Not required for synthetic sequence pairing; not established for underlying target library |
| Target genes | Underlying 130-RNA panel reused; this workbook is a bit/amplifier library, not a target sequence set |
| Target regions | Amplifier handles and repeat-binding regions; not RNA-targeting sites |
| Genomic/transcript coordinates | Not applicable to synthetic sequences. No genomic interval convention is implied |
| Target-binding sequences | No RNA-target-binding sequence library in S12 |
| Full probe sequences | 16 primary + 16 secondary amplifier strings, original and final readout columns parsed; each amplifier 130 nt after whitespace removal |
| Barcode assignments | Bits 1–16 present; gene-to-codeword mapping is not contained in this table |
| Readout sequences | 16 exact final readout strings and 16 original column strings parsed |
| Primer/linker sequences | 2-nt spacers within amplifiers; dye/linker chemical metadata not fully specified by workbook |
| Design software | No runnable official bDNA-specific design program established in examined S08/S12; do not confuse S06 RNA design code with a demonstrated bDNA generator |
| Design parameters | Three-letter design, five repeats, empirical amplifier screening reported |
| Source | S08 Fig. 1 and Amplifier staining; S12 Sheet1 rows 2–18 and SI caption |
| Availability | Publisher workbook downloaded and parsed 2026-10-15; SHA-256 in evidence/file_manifest.json |
| Missing information | No historical random seed/software release; final readout dye/linker mapping incomplete at workbook level; no RNA target design output |
| Benchmark grade | B |
| Grade justification | Mostly reproducible ONLY for synthetic amplifier/readout assembly: all 16 sequence sets and bit identities are present; 3′ binding placement is reported and sequence pairing checks agree. Missing historical generation seed and label metadata are minor for this sequence-pairing task. This is not Grade B for target selection, binding specificity or full wet-lab reproduction. |
| Comparison scope | Restricted handle/repeat/readout assembly and sequence-pairing comparison; not U-Probe target-selection ground truth |
| Assessment date | 2026-10-15 |
| Evidence status | reported strings/placement; derived lengths and pair checks; inferred scope-limited grade |

### D04 — 10,050-gene MERFISH codebook and encoding resource

| Field | Value |
|---|---|
| Dataset ID | D04 |
| Dataset name | 10,050-gene MERFISH codebook and encoding resource |
| Method version | V06 |
| Species | Homo sapiens |
| Cell / tissue / sample | U-2 OS |
| Genome build | Not independently established from inspected main-text design account |
| Transcriptome / annotation version | Not independently established; supplementary design reference requires inspection |
| Target genes | 10,050 genes; 9,050 nonoverlap + 1,000 overlap subset |
| Target regions | 30-nt targeting regions; up to 48 probes for nonoverlap design; 20-nt overlap allowed in subset |
| Genomic/transcript coordinates | Not established; do not borrow V02 coordinates/reference |
| Target-binding sequences | Encoding resource associated with Dataset S2 cited in design account; binary not parsed |
| Full probe sequences | Member-specific synthesized/template boundaries unresolved |
| Barcode assignments | 69-bit HD4 weight 4; 10,050 assigned of 12,903 words; 2,853 blanks. Dataset S1 cited |
| Readout sequences | 69 readout identities; full table not parsed |
| Primer/linker sequences | Not independently established |
| Design software | Specific design script/release not independently inspected |
| Design parameters | Three readouts per probe; 23 three-color readout rounds; other design details require SI |
| Source | S09 Results first subsection and Supplementary Material |
| Availability | Main article/inventory inspected 2026-10-15; sequence CSV and SI not parsed |
| Missing information | Binary designs/codebook; full reference; exact oligos; coordinate conventions; complete implementation |
| Benchmark grade | C |
| Grade justification | Inspected evidence establishes architecture and coding/overlap strategy, but not reproduction of the library. Full published resources may enable a later higher grade after inspection. |
| Comparison scope | Architecture, overlap and code-capacity comparisons only |
| Assessment date | 2026-10-15 |
| Evidence status | reported main-text facts; inferred grade |

### D05 — Author-hosted example data/design inputs

| Field | Value |
|---|---|
| Dataset ID | D05 |
| Dataset name | Author-hosted example data/design inputs |
| Method version | V01/V02 software examples |
| Species | Human example context reported; archive members not inspected |
| Cell / tissue / sample | Examples from original/fast-design workflows |
| Genome build | Script names Homo_sapiens.GRCh38.ncrna.fa; exact member bytes not checked |
| Transcriptome / annotation version | Transcript FASTA and isoform FPKM tracking named in example code; releases unresolved |
| Target genes | Archive codebook.csv named; contents not inspected |
| Target regions | Target design examples, not a verified complete historical library |
| Genomic/transcript coordinates | Unexamined |
| Target-binding sequences | Not inspected |
| Full probe sequences | Not inspected |
| Barcode assignments | Author page reports original MHD2/MHD4 barcode sets; design example expects codebook.csv |
| Readout sequences | Design example expects readouts.fasta; binary not accessed |
| Primer/linker sequences | Primer output is generated by example code; archive contents unexamined |
| Design software | S06 library_design_example.m |
| Design parameters | Example parameters are in pinned source; execution not performed |
| Source | S07 Additional Data and Resources; S06 README and example script |
| Availability | Author page accessible 2026-10-15; archive redirected to S3 and access returned 403 |
| Missing information | Archive contents and integrity; exact reference bytes; historical output equivalence |
| Benchmark grade | Unassessed |
| Grade justification | No grade assigned because archive contents were not examined. Resource discovery is not evidence of reproducibility. |
| Comparison scope | Potential software workflow replay after archive access is resolved |
| Assessment date | 2026-10-15 |
| Evidence status | reported linked resource; unavailable inspection of binary contents |
