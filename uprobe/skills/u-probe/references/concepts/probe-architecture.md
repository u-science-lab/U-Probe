# Probe architecture and orientation

Read this when you write or check a `probes:` block: which part binds the target, which parts are
flanks, which strand each part must be, and how split pairs and padlocks are assembled. Ends with a
glossary. U-Probe syntax: `references/uprobe/yaml-reference.md`.

## Parts of a FISH probe

| Part | Role | Typical length | Example methods |
| --- | --- | --- | --- |
| Binding region (homology, arm, target-binding sequence) | hybridizes to the target RNA, cDNA or DNA | 17-22 nt (smFISH), 25-30 (HCR halves, seqFISH+, MERFISH), 2 x 13-25 (padlock arms), 35-41 (OligoMiner) | all |
| Readout / overhang / encoding sequence | binding site for a labelled secondary (readout) oligo; encodes the gene's barcode | 15-20 nt each, 1-4 per probe | MERFISH, seqFISH+, smiFISH (FLAP) |
| Initiator (split or whole) | triggers hybridization chain reaction amplification | 2 x ~18 nt halves for split initiators | HCR v3 |
| Backbone (padlock) | joins the two arms of a padlock; carries the gene ID/barcode and an anchor (common detection site) | e.g. 40 nt in HybISS (20-nt ID + 20-nt anchor) | padlock ISS, HybISS, STARmap, MiP-seq |
| Primer site | PCR amplification of oligo pools, or the RCA primer of a two-oligo design | ~20 nt | Oligopaint, MERFISH libraries; SNAIL/MiP-seq primer |
| Spacer | a few nt that separate functional parts | 1-3 nt | HCR (`AA`/`TA`), MERFISH (`A`) |

The fluorophore of a directly labelled oligo is a synthesis modification (e.g. a 3' amine coupled to a
dye, Raj et al. 2008), not part of the designed sequence.

Sources: Raj A et al. 2008 Nat Methods 5:877, doi:10.1038/nmeth.1253 (smFISH, 48 x 20-nt singly
labelled oligos); Chen KH et al. 2015 Science 348:aaa6090, doi:10.1126/science.aaa6090 and Moffitt JR
et al. 2016 PNAS 113:11046, doi:10.1073/pnas.1612826113 (MERFISH encoding probes: 30-nt target region,
readout sequences, PCR primers); Eng CHL et al. 2019 Nature 568:235, doi:10.1038/s41586-019-1049-y
(seqFISH+); Choi HMT et al. 2018 Development 145:dev165753, doi:10.1242/dev.165753 (HCR v3 split
initiators); Nilsson M et al. 1994 Science 265:2085, doi:10.1126/science.7522346 (padlock probes);
Beliveau BJ et al. 2012 PNAS 109:21301, doi:10.1073/pnas.1213818110 (Oligopaints); Tsanov N et al.
2016 Nucleic Acids Res 44:e165, doi:10.1093/nar/gkw784 (smiFISH FLAP).

## Orientation rules

Nucleic acid duplexes are antiparallel: a probe binds the sequence that is its **reverse complement**.
Write every sequence 5'->3'; U-Probe templates concatenate parts 5'->3' from left to right.

What U-Probe's `target_region` is:
- RNA sources (`mRNA`, `exon`, `CDS`, `UTR`): the **RNA sense sequence** of the window, 5'->3', as the
  mature RNA reads (genes on the - strand are reverse-complemented). Equal to the cDNA/mRNA in Ensembl.
- `source: genome`: the reference **+ strand** of the window.
- A literal sequence target: exactly the sequence you supplied.

| Target molecule | It reads (5'->3') | Probe binding part | U-Probe `expr` |
| --- | --- | --- | --- |
| mRNA / pre-mRNA / ncRNA (smFISH, MERFISH, seqFISH+, HCR, direct-RNA padlocks) | `target_region` | reverse complement | `rc(target_region)` (or `rc(target_parts['...'])`) |
| cDNA made in situ by reverse transcription (padlock ISS, HybISS) | `rc(target_region)` | RNA-sense sequence | `target_region[...]` |
| Denatured genomic DNA (Oligopaint, DNA FISH) | both strands | either; one strand per probe | `target_region` (binds the - strand) or `rc(target_region)` (binds the + strand) |
| Readout probe R (as ordered, labelled) binding the primary probe | primary carries `rc(R)` | - | `rc(encoding[target]['R1'])` when `encoding` stores the readout probe's own sequence |
| PCR primers of an oligo pool | forward primer as is at 5'; reverse primer as its reverse complement at 3' | - | template constants |

Checks:
- A probe for an RNA must **not** contain the RNA's own sequence; if `binding == target_region`, it
  binds only the antisense strand.
- In a DNA pool that uses both strands, a probe on one strand and a probe on the other strand from
  overlapping windows are complementary to each other and can pair in solution; keep overlapping
  windows on one strand (`remove_overlap` acts on positions, not strands).
- Store in `encoding` whichever sequence the source paper lists, and say so in a comment: MERFISH and
  seqFISH+ templates store the readout probe and take `rc(...)`; HCR templates store the initiator as it
  appears in the probe.

For off-target searches the orientation matters too: search `target_region` with `strand: sense`
(sites where the window itself occurs in a transcript = RNAs the probe can bind); search a probe's
binding part with `strand: antisense` (`references/concepts/specificity.md`).

## Multi-part probes

### Split pairs (HCR v3, RNAscope-like ZZ pairs)

Two probes bind adjacent sites; signal requires both. Choi et al. 2018 split the HCR initiator between
two probes (2 x 25-nt binding regions, 2-nt gap): HCR is "triggered if, and only if, both ... bind
specifically to their adjacent binding sites". RNAscope's "double Z" probes (18-25 nt each, 10-20
pairs per target) must "hybridize to the target sequence in tandem" for amplification (Wang F et al.
2012 J Mol Diagn 14:22, doi:10.1016/j.jmoldx.2011.08.002).

Geometry: model the whole site as one window `left + gap + right` (RNA sense). Each half-initiator must
sit at the probe end that faces the gap:

```yaml
extracts: {target_region: {source: mRNA, transcripts: mane, length: 52,
  layout: {template: "{left}{gap}{right}", parts: {left: {length: 25}, gap: {length: 2}, right: {length: 25}}}}}
probes:
  probe_odd:  {template: "{initiator}{spacer}{binding}", parts: {initiator: {expr: "encoding[target]['init_odd']"},
               spacer: {expr: "encoding[target]['spacer_odd']"}, binding: {expr: "rc(target_parts['left'])"}}}
  probe_even: {template: "{binding}{spacer}{initiator}", parts: {binding: {expr: "rc(target_parts['right'])"},
               spacer: {expr: "encoding[target]['spacer_even']"}, initiator: {expr: "encoding[target]['init_even']"}}}
```

`rc(left)` has its 5' end over the gap, so the odd initiator goes 5'; `rc(right)` has its 3' end over
the gap, so the even initiator goes 3'. Run specificity on the whole 52-nt window: the pair is the unit
that produces signal.

### Padlock probes

A padlock is one linear oligo whose two terminal arms hybridize side by side on the target; a ligase
joins the ends only if they are juxtaposed and matched at the junction, and the circle is copied by
rolling-circle amplification (RCA) into one bright product per ligation event (Nilsson et al. 1994;
Larsson C et al. 2004 Nat Methods 1:227, doi:10.1038/nmeth723).

Rule: for the strand the padlock binds, with 5' half H5 and 3' half H3 of its site,
**padlock = rc(H5) + backbone + rc(H3)**; the 5' arm pairs with the 5' half, and the two ends meet at the
H5|H3 junction.
- Direct RNA (MiP-seq template): site = `target_region` = A B -> `rc(A) + backbone + rc(B)`.
- cDNA (padlock ISS, HybISS; Ke R et al. 2013 Nat Methods 10:857, doi:10.1038/nmeth.2563; Gyllborg D et al.
  2020 Nucleic Acids Res 48:e112, doi:10.1093/nar/gkaa792): site = `rc(target_region)`, so with
  `target_region` = W1 W2 the padlock is `W2 + backbone + W1`
  (`arm_5: target_region[15:30]`, `arm_3: target_region[0:15]` for 2 x 15 nt).
- Check: joining the 3' arm to the 5' arm (arm_3 + arm_5) gives the reverse complement of the bound
  strand's site; for cDNA padlocks, that is `target_region` itself.

### Padlock + primer and other two-oligo units

STARmap SNAIL (Wang X et al. 2018 Science 361:eaat5691, doi:10.1126/science.aat5691) and MiP-seq bind a
padlock and a primer next to each other on the RNA; the primer holds the padlock's ends for ligation
(a splint) and primes RCA, so a padlock is circularized and amplified only when both bind the same RNA
molecule, which excludes signal from single misbound probes (the double-hybridization requirement). In U-Probe make one window
with a `layout` (`padlock_site`, `gap`, `primer_site`), bind each oligo to its part with `rc(...)`, and
run the off-target search on the whole window.

### Several probes per window

`probes:` may define any number of oligos per window, and parts may reference other probes
(`rc(circle_probe.part2.barcode2)`), e.g. an amplification oligo that must bind a barcode carried by a
circle probe. Expressions form a DAG; cycles are rejected.

## Glossary

- **Window / target region**: the target stretch one probe (or probe unit) binds; `target_region`.
- **Probe set / pool**: all probes of one target; a spot needs many of them in multi-probe methods.
- **Primary (encoding) probe**: binds the cellular target and carries readout sites.
- **Readout (secondary) probe**: labelled oligo that binds a readout site; changed between rounds.
- **Initiator**: HCR trigger sequence; **amplifier**: the HCR hairpins it opens.
- **Padlock**: linear probe circularized by ligation on the target; **backbone**: its non-binding middle;
  **anchor**: a common detection site in the backbone; **ID/barcode**: gene-specific backbone sequence.
- **RCA / RCP**: rolling-circle amplification / its product, a long concatemer of the circle's
  complement that collapses into one bright, countable spot.
- **Splint**: an oligo (or the target) that holds two ends together for ligation.
- **Tiling, spacing**: placing windows along a target; the gap between adjacent windows.
- **Tm**: temperature at which half of the duplexes are formed (`references/concepts/hybridization-and-tm.md`).
- **Off-target**: any other molecule the probe can bind (`references/concepts/specificity.md`).
- **Paralog / pseudogene / retrocopy**: duplicated gene / non-functional copy / copy made by reverse
  transcription of an mRNA (`references/concepts/specificity.md`).
- **MANE Select, Ensembl canonical, isoform, junction window**: `references/concepts/targets-and-isoforms.md`.
