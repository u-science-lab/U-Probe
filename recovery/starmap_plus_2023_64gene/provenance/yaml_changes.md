# Design configuration changes

The original specification was a source-grounded architecture description, not an executable U-Probe protocol. The measured run used a frozen executable configuration. The portable case preserves the design choices; its runner was reorganized after the measured run and was not subjected to another full candidate-generation run during curation.

| Original specification | Executable choice | Basis and consequence |
|---|---|---|
| Standalone meta/probe_design fields | Native genome/extracts/probes/encoding configuration | Use existing UProbeAPI extraction, layout and construction functions without changing core algorithms. |
| Reference database unspecified | GRCm39 genomic FASTA plus matching NCBI RefSeq RS_2024_02 GTF | Available reference; independent checks confirmed published arms are compatible. Historical database hashes are recorded in run_metadata.yaml. |
| Transcript selection unspecified | All 308 isoforms, spliced separately | Include 157 NM, 140 XM, 4 NR and 7 XR accessions. Do not silently restrict to coding or canonical transcripts. All full-length native extractions match independent GTF splicing and prebuilt transcript FASTA. |
| Both binding arms observed at 19-25 nt | Enumerate each arm independently at 19-25 nt | Do not constrain combined length to the narrower values observed in the finite published set. |
| Binding-site gap unresolved | RNA order padlock_site, gap 1-2 nt, primer_site | One-time geometry audit found 62 reference pairs with gap 1 and 194 with gap 2. This is a panel-wide inferred architecture, not a parameter explicitly established from paper text or per-probe design input. |
| Window and sampling unresolved | Total 39-52 nt, step 1, all 98 length layouts | Exhaustive site enumeration avoids stride or top-k omissions. |
| Quality thresholds unresolved | Uppercase ACGT and legal architecture only; no GC, Tm or off-target filters | Original numerical computational thresholds were not established. The result is a permissive architecture-constrained candidate-coverage upper bound. |
| Candidate export unresolved | Keep all candidates; deduplicate identical gene/window/layout across positions and isoforms | Preserve every distinct eligible sequence/layout for recovery evaluation. No final ranking or four-pair limit was applied. |
| Functional structure specified | Preserve the fixed backbone, connection order and gene-specific codebook | Padlock: A{code}TA{rc(padlock_site)}AATTATTAC{code}CATACACTAAAGATA. Primer: {rc(primer_site)}TA{rc(code)}TTATCTT. Record padlock 5-prime phosphate separately from DNA. |

## Freeze and execution revisions

The measured design was frozen at 2026-10-31T03:28:07Z. Its SHA256 is `4bff8b65d58e2feb5d218f064b7279ea14d568f13c483b962cf4d72a64bcd608`. Historical genome-path configuration SHA256 is `cf76dcffcdba4e87f3cd651e7521962e2e387d60b9b5710195a5f26fae509fe8`. These identify the measured historical files, not a promise that portable path/comment changes have identical bytes.

Resource revisions changed batching, worker-loss detection, parallel construction and lossless compression only. No frozen lengths, gap, transcript policy, filters, codebook or backbone changed. The initial cgroup OOM and subsequent long-transcript execution are recorded in run_metadata.yaml. Npy decompressed records were byte-identical across all three execution strategies; an independent window-count oracle confirmed all 64 per-gene totals and 44,271,276 total pairs.

Published binding sequences and full oligos were reserved for integrity auditing and evaluation. Design input contained gene identifiers, the public functional codebook/backbone, genome and annotation; it did not contain reference binding strings or reference coordinates. The reported 100% SRR/TRR/MTR/CPSR describes full-sequence recovery from the exhaustive independent candidate collection. Final selected-set metrics are not available.

The current portable design SHA256 is `db3d3555784580ac79d71d3082774ba492b926a28526e41b26d554f479aa7b81`. The only configuration edits during curation changed the protocol display name to `STARmap_PLUS_2023_64gene_unfiltered_candidate_recovery` and renamed the execution annotation key to `recovery_execution`; generation parameters, genes, barcodes and grammar remain unchanged.
