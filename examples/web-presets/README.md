# Web presets (reference copies)

The web interface loads its built-in MiP-seq and STARmap presets (`/workflow/builtin_probes`) from the
`probes/` directory of the server's data directory, not from this repository. The two files here are
reference copies of those presets, added in `b149db2` (2026-10-28) so that the configurations behind the
first version's RCA experiments can be read next to the library templates. **The web interface does not
read these files, and they are not library templates**: `uprobe templates`, `uprobe run -t ...` and the
Agent Skill use `uprobe/templates/mipseq.yaml` and `starmap.yaml`, which carry hybridization conditions,
literature sources, example inputs and limitations. This directory is not part of the wheel.

| File | SHA-256 of the copy |
| --- | --- |
| `MiP-seq_web.yaml` | `3a954b2f6ba810a913042ea90dcaf10516bde59bc3f4f87607ad975a00041daa` |
| `STARmap_web.yaml` | `d5c94a3aafa8e437aa2af245af9c12c0d71d68c6eba1c43f60e8505394c36684` |

The server's own files are the source of truth. Compare `sha256sum <data_dir>/probes/<method>.yaml` with the
table before relying on a copy, and update the copy and the hash together when the server's file changes.

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
