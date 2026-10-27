# Targets, regions and isoforms

Read this when you decide *where* on a gene the probes go: mature mRNA vs CDS vs UTR vs intron,
which transcript (MANE Select, Ensembl canonical, a listed isoform), isoform-specific designs, exon
junction windows, and short targets. U-Probe options: `extracts.target_region.source` and `transcripts`
(`references/uprobe/yaml-reference.md`).

## Region choices

| Region | U-Probe | Binds | Use when | Caveats |
| --- | --- | --- | --- | --- |
| Mature mRNA (5' UTR + CDS + 3' UTR) | `source: mRNA` | spliced RNA (and the exonic part of pre-mRNA) | counting mRNA molecules; most room for probes | 3' UTR length varies (alternative polyadenylation) |
| Coding sequence | `source: CDS` (`exon` is the same) | CDS of coding genes; all exons of non-coding genes | reproducing a protocol that tiled the CDS (Raj et al. 2008 smFISH; seqFISH+, STARmap templates) | fewer sites; CDS is shared most with paralogs and processed pseudogenes |
| UTRs | `source: UTR` | 5'/3' UTR features | separating paralogs that share CDS; adding sites | no windows across junctions; genes without UTR features give none |
| Introns / pre-mRNA | `source: genome` with the intron coordinates, or a literal sequence | nascent transcripts at the transcription site | transcription activity, nascent RNA, splicing kinetics | orientation and specificity need care (below) |
| Genomic locus (DNA FISH) | `source: genome` | denatured DNA, either strand | chromatin tracing, Oligopaint | repeats; see `references/concepts/specificity.md` |

Evidence for these uses:
- smFISH on the coding region of one reference mRNA: Raj A et al. 2008 Nat Methods 5:877,
  doi:10.1038/nmeth.1253. U-Probe's smFISH template uses the whole mRNA instead, because paralogs and
  pseudogenes share coding sequence more than UTRs (README, measured on Ensembl 113: the reference mRNA
  gave six housekeeping genes 119 probes, the CDS of all isoforms 11, without an expression table).
- Isotype-specific actin 3' UTRs ("isotype specific but ... conserved in evolution"), used as
  isotype-specific probes: Ponte P, Gunning P, Blau H, Kedes L 1983 Mol Cell Biol 3:1783,
  doi:10.1128/mcb.3.10.1783-1791.1983.
- The opposite trade-off: an smFISH tissue protocol asks for about 1,054 nt of target for 48 x 20-nt
  probes and advises avoiding UTRs (Lyubimova A et al. 2013 Nat Protoc 8:1743, doi:10.1038/nprot.2013.109);
  UTRs vary between isoforms and can carry repeats. Check `_off_targets` when you include them.
- Regions shared by all spliced isoforms: seqFISH designs used "consensus regions of all spliced
  isoforms" (Eng CHL et al. 2017 Nat Methods 14:1153, doi:10.1038/nmeth.4500), U-Probe's seqFISH+
  template approximates them with the reference transcript's CDS.
- Alternative polyadenylation: ~54% of human and ~32% of mouse genes have alternative poly(A) sites
  (Tian B et al. 2005 Nucleic Acids Res 33:201, doi:10.1093/nar/gki158); proliferating cells shorten 3'
  UTRs (Sandberg R et al. 2008 Science 320:1643, doi:10.1126/science.1155390). Probes only in the distal
  3' UTR miss short-UTR isoforms.
- Intron probes label nascent RNA at transcription sites: "the intron spot marked the site of
  transcriptionally active genes" (about 16 x 20-nt probes in the first intron; Levesque MJ & Raj A 2013
  Nat Methods 10:246, doi:10.1038/nmeth.2372); exon + intron probe sets separate co- and
  post-transcriptional splicing (Vargas DY et al. 2011 Cell 147:1054, doi:10.1016/j.cell.2011.10.024);
  intron seqFISH used 25 probes on 5' introns (Shah S et al. 2018 Cell 174:363,
  doi:10.1016/j.cell.2018.05.035).

## Which transcript

Most human multi-exon genes are alternatively spliced (~92-95%; Wang ET et al. 2008 Nature 456:470,
doi:10.1038/nature07509; Pan Q et al. 2008 Nat Genet 40:1413, doi:10.1038/ng.259). The transcript
set decides which windows exist.

| `transcripts` | Windows from | Pick it when |
| --- | --- | --- |
| `mane` | MANE Select (human); Ensembl canonical where there is none (other species, most ncRNAs) | default for human mRNA counting |
| `canonical` | GTF tag `Ensembl_canonical` (one per gene, all species in Ensembl) | non-human, or to match Ensembl's display |
| `[ENST...]` | the listed transcripts (one per target needs to be listed) | a specific isoform, or a transcript the user names |
| `all` (default) | union of every annotated transcript's exons | the most windows; input for the isoform analysis below (not a "shared by all isoforms" filter) |

- MANE Select: "a representative transcript for each human protein-coding gene", matched end to end
  between RefSeq and Ensembl/GENCODE (UTRs included) and aligned perfectly to GRCh38 (Morales J et al.
  2022 Nature 604:310, doi:10.1038/s41586-022-04558-8; v0.95 covered 18,584 genes, 97%; MANE v1.0 has
  19,062 MANE Select transcripts per NCBI's release file).
- Ensembl canonical: a single representative transcript per locus, replaced by MANE Select where
  available (Yates AD et al. 2020 Nucleic Acids Res 48:D682, doi:10.1093/nar/gkz966). Exists for all
  Ensembl species and non-coding genes.
- `all` includes retained-intron, NMD and fragmentary transcripts: windows can lie in sequence the
  reference mRNA does not have. README (Ensembl 113): with `all`, only 16-38% of the final probes of six
  housekeeping genes lay on the reference mRNA.
- `summary.json -> targets.<t>.transcripts` (`chosen_by`: `MANE_Select`, `Ensembl_canonical`, `listed`,
  `all`) records what was used; report it.
- The reference transcript and the method's windows come first; coverage of other isoforms only breaks
  ties between windows that both satisfy the method. Do not shorten probes, move off the reference or add
  filters to reach isoforms the user did not ask for.

### Reference transcript with ~0 TPM

A transcript-level expression table can give the reference transcript (MANE Select or Ensembl canonical)
about 0 TPM while the gene is well expressed. That is weak evidence, not a reason to leave it: transcript
TPM of overlapping isoforms is an estimate. Salmon, for instance, infers abundances "given the collection of
ambiguously mapping sequenced fragments" with an EM algorithm (Patro R et al. 2017 Nat Methods 14:417,
doi:10.1038/nmeth.4197; kallisto and RSEM do the same), so reads that the isoforms share are split between
them by the model, not observed. MANE Select is by design one representative transcript per protein-coding
gene, an exact match between the Ensembl/GENCODE and RefSeq annotations (Morales J et al. 2022 Nature 604:310,
doi:10.1038/s41586-022-04558-8).
- Keep the units on sequence of the reference transcript. If the table suggests that another isoform carries
  the signal, choose windows the reference shares with the expressed isoforms: run on the reference, then
  check which units also lie on the expressed isoform with `check_probes.py --transcript TARGET=ID
  --transcripts transcript.fa` (units not on it fail).
- If the shared sequence holds fewer units than requested, keep the units it does hold and report the
  shortfall. Add units on other isoforms only as labelled extras, never in place of a unit that satisfies
  both; do not switch the whole design to another transcript.
- Do not add coverage filters (for example a minimum share of the expressed isoforms' TPM) or arm-Tm floors
  the method does not publish to force the choice (SKILL.md rule 10).
- For the off-target rule the gene's total still counts: `min_fraction_of_target` uses the sum over all the
  target's transcripts, so a ~0 reference transcript does not make the relative threshold zero.

## Isoform-specific designs

A probe set is isoform-specific only if its windows lie in sequence that the other isoforms lack:
an isoform-unique exon (or exon part), or an exon-exon junction only that isoform has.
- Unique exons: tile them with a full probe set if long enough. Waks Z, Klein AM, Silver PA 2011 Mol
  Syst Biol 7:506, doi:10.1038/msb.2011.32 measured alternative splicing in single cells with one smFISH
  set per alternative terminal exon (tens to 96 x 20-nt probes; isoforms had to differ by > 800 bp). A
  constitutive-exon set in a second colour gives total mRNA.
- U-Probe does **not** score isoform specificity: `mapped_genes` counts genes, so other isoforms of the
  target gene are always "the target". Find unique sites yourself:
  1. Run with `transcripts: all`, the off-target filter, and `remove_overlap` (it merges the transcript
     lists of identical sites into `transcript_names`).
  2. Keep rows whose `transcript_names` contain only the wanted isoform(s). Check against the GTF that
     no other transcript of the gene contains the window (sites in overlapping but different exon
     entries are listed separately before `remove_overlap`).
- `transcripts: [ENST_wanted]` alone is **not** isoform-specific: it tiles the whole isoform, shared
  exons included.
- Sites shared by all major isoforms (to count every isoform equally): the same `transcript_names`
  analysis, keeping rows that list every isoform you care about.

## Repeats inside the target

A target with tandem or internal repeats (tandem-repeat genes, repeated domains) contains the same
window at several positions, and one oligo binds every copy. U-Probe treats it as one oligo:
- Identical windows of a target count once: `remove_overlap` keeps such an oligo only if none of its sites
  is within `location_interval` of a kept one and then occupies all of them, and `equal_space` ranks it once.
  So spacing holds between oligos over all their sites, and the count is of distinct oligos.
- The kept row gets `n_sites` and `oligo_sites` (`start-end;...` along the target), and `summary.json` a
  `multi_site_oligos` warning listing each such oligo of the final table. `scripts/check_probes.py` reports
  them, counts distinct oligos (`--min-probes`) and checks `--min-gap` over all their sites (from
  `oligo_sites`, or from the transcript with `--transcripts`).
- A multi-site oligo adds signal per probe (it binds several copies); report them, and keep the others
  as the evenly spread set. A shortfall in distinct oligos is reported before any change of selection
  (`references/uprobe/troubleshooting.md`, "Repeats inside the target").

## Exon-exon junction windows

When `transcripts` selects exactly one transcript per target, `mRNA`/`CDS`/`exon` tile its spliced
sequence, so windows also cross junctions (`exon_name` lists the joined exons, `;`-separated). Such
windows bind only the spliced RNA, not pre-mRNA or genomic DNA. README: ACTB's MANE mRNA gives 1,793
20-nt windows instead of 1,698 within exons (19 more per junction).
- A junction offers only (window length − 1) windows, overlapping each other: enough for one probe or
  one probe unit, not for a many-probe set. Junction-specific detection therefore suits single-probe-signal
  methods (padlock/RCA, a split pair spanning the junction), not smFISH-type pools.
- Selecting only junction windows: there is no numeric column for it; filter the TSV afterwards
  (`exon_name` contains `;`).
- Off-target checks: `mapped_genes` searches spliced transcripts and handles junction windows;
  `mapped_sites` and `kmer_count` search the genome and do not (`references/uprobe/troubleshooting.md`).
- With `transcripts: all` or several transcripts per target, windows stay inside single exons.

## Introns and pre-mRNA in U-Probe

- `source: genome` gives the reference **+ strand**. The probe must be the reverse complement of the
  RNA: for a gene on the + strand use `rc(target_region)`; for a gene on the − strand the RNA is
  `rc(target_region)`, so the probe is `target_region`. Getting this wrong labels nothing.
  Alternative: supply the intron as a literal target in RNA sense and use `rc(target_region)`.
- `mapped_genes` searches only mature transcripts (`transcript.fa`); off-targets in other introns are
  invisible to it. Use genome-based checks for intronic windows: `mapped_sites` (Bowtie2, windows >= 22
  nt) and `kmer_count` (Jellyfish), plus repeat filtering; introns are repeat-rich (`references/concepts/specificity.md`).
- Nascent signal is a few molecules at one site: use many probes over the 5' part of long introns, or
  amplification.

## Short targets

Capacity bound: with window length w and gap g (`location_interval`), a target of usable length L holds
at most about (L + g) / (w + g) non-overlapping windows, before filters remove some.

| Design | w + g | Probes wanted | Usable sequence needed |
| --- | --- | --- | --- |
| smFISH 20 nt, gap 2 (Stellaris default) | 22 | 48 | ~1,050 nt |
| MERFISH 30 nt, gap 0 | 30 | 92 | ~2,760 nt |
| HCR v3 pair 52 nt, gap 3 | 55 | 20 pairs | ~1,100 nt |
| seqFISH+ 28 nt, gap 0 | 28 | 24 | ~670 nt |

Filters (specificity, Tm, GC) typically remove a sizeable share of positions, so a ~1-kb mRNA rarely
yields 48 smFISH probes. For short or crowded targets, in this order:
1. Widen the region: `mRNA` (UTRs) instead of `CDS`; a single reference transcript tiles across junctions.
2. Smaller `location_interval` down to the method's smallest published value. Report the shortfall of the
   published-length set first: a `length` range is not a rescue, because it changes the length of most of the
   set (`references/uprobe/troubleshooting.md`, "Length ranges are not a way to fill the count"); take
   another length only at sites where the standard one fails, and only inside the method's published range.
3. Rescue specificity losses legitimately: expression table, `exempt_copies` (`references/concepts/specificity.md`).
4. Accept fewer probes if the method tolerates it (state the expected dimmer spots), or suggest an
   amplified method (HCR, padlock/RCA; `references/techniques/`) for very short RNAs: switch only if the
   user agrees. Never fill the count by loosening a published rule (SKILL.md rule 10).

## Common goals -> region

| Goal | Region and transcripts |
| --- | --- |
| Count mRNA molecules of a gene | `mRNA` + `mane` (human) / `canonical` (other species); `CDS` if reproducing a CDS-based protocol |
| Detect all isoforms equally | windows in exons shared by the major isoforms (`all` + `transcript_names` analysis) |
| One isoform | its unique exons, or its unique junction with a single-probe-signal method |
| Transcription sites / nascent RNA | introns (5' introns of long genes), with an exon set in a second colour |
| Paralog-rich family | `mRNA` (UTRs diverge most); check `_off_targets` and `_copies` |
| Non-coding RNA | `mRNA` (CDS falls back to exons) with `canonical`/`mane` |
| Reporter, transgene, virus, fusion junction | literal sequence target (`name: SEQUENCE`), RNA sense |
| Genomic locus | `source: genome`, region coordinates |
