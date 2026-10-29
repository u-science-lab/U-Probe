# Web presets (reference copies)

The web interface loads built-in method configurations from the `probes/` directory in the server's data directory via `/workflow/builtin_probes`. These files are reference copies; the interface does not read this directory. They are not library templates and are not included in the wheel.

The eight additional configurations are scoped research templates, not complete reproductions of their published methods. Consult their comments and the corresponding `fish-methods-corpus` entries for version, missing design rules and evidence limitations. Supply panel-matched external readout, barcode, bridge and primer sequences where required; a generated sequence is not automatically a validated codebook. Shared library components must be consistent within the relevant library.

| File | SHA-256 of the copy |
| --- | --- |
| `ExSeq_web.yaml` | `4e0326c2ffdabb31fa22fbc5ebbdc589eb712a627f4260ea5484790c37cd8378` |
| `MERFISH_web.yaml` | `8b2d1697bf63b50c55f346c4be9755639d378f8a3b780be8c6fcff11b409febf` |
| `MiP-seq_web.yaml` | `3a954b2f6ba810a913042ea90dcaf10516bde59bc3f4f87607ad975a00041daa` |
| `Oligo-LiveFISH_web.yaml` | `7d370bd9e5803e57d74810c467c0e60aedc2b139186d277ba499cf9ef0b96a0d` |
| `OpenFISH_web.yaml` | `cf3ae978c06d9b61b7fa9494731c1fce9083a27a6804d7a82b3269032c0c7e6d` |
| `PRISM_web.yaml` | `86b11af1c39ed95b206b313d13d9031b6fc5785d883f1f668fb90918f95ef72d` |
| `RAEFISH_web.yaml` | `b0170a3a3c11ef3d0e2d55727aeed8a63fc174c6e065a47a74914b15d5e9ec80` |
| `STARmap_web.yaml` | `d5c94a3aafa8e437aa2af245af9c12c0d71d68c6eba1c43f60e8505394c36684` |
| `TDDN-FISH_web.yaml` | `19b178747ae44b5583dee938a4e2a39a71a8f304b9d6c0c6075d87e367bdf5c2` |
| `seqFISH+_web.yaml` | `422bff872b8d577f4a927c0cf1cd1fcca948c673c96ce1dfa03922a3766afa48` |

The server files are the source of truth. Compare each copy with `<data_dir>/probes/<method>.yaml` and update its hash when the server file changes.

## Differences from the library templates (U-Probe 0.6.0)

As described in issue #5 and checked against the files.

| | Web preset | Library template |
| --- | --- | --- |
| **MiP-seq** target | 40-nt exon windows | mRNA, `transcripts: mane`, explicit 13 + 13 + 1 + 13 nt layout |
| MiP-seq construction | `BC1 + A + BC2` between the padlock arms; primer tail `rc(BC2)[:-2]` | B1/B2 barcodes in a fixed backbone; fixed RCA-primer tail `AGGCTCCACTA` |
| MiP-seq specificity | Bowtie2 on the whole 40-nt window, `target_mapped_genes <= 4` | Bowtie 1, each oligo's site searched on its own, up to 2 substitutions, `target_mapped_genes <= 1` |
| MiP-seq filters | Tm 36-44 C for all three target segments, fold-score and terminal-GC filters | Tm >= 28 C per binding segment, target start within the first 1,000 nt, two sites per gene |
| **STARmap** barcodes | four fields `BC1`-`BC4`, backbone supplied through `encoding` | one `code` field in fixed padlock and primer backbones |
| STARmap layout | 40-46-nt CDS windows, gap 0-2 nt | 40-46-nt CDS windows, gap 1-2 nt |
| STARmap filters | whole-window GC 0.3-0.7, window Tm 56-70 C, both arm Tm >= 60 C | 56-70 C per binding segment under stated conditions |
| STARmap specificity | no mapped-gene filter | paired-oligo specificity filter, Tm-balanced split, four sites per gene |

The two sets are different designs, not two spellings of one. Unifying them, with the fixed sequences
(backbone, barcodes, primer tail) as replaceable modules and the published sequence as the default, is
tracked in issue #5.
