# DNA Oligopaint 2020: chrII sequence recovery

[Woglar et al., PLOS Biology (2020)](https://doi.org/10.1371/journal.pbio.3000817):
[S2 Data](https://doi.org/10.1371/journal.pbio.3000817.s009) contains 9,984 complete
120-nt oligos, representing 9,975 unique target/sequence records across 100 bins.
[S1 Table](https://doi.org/10.1371/journal.pbio.3000817.s010) supplies readout sequences.
The observed 5′→3′ structure is index5 (20 nt), fiducial (20), barcode (20), genomic
binding (40), index3 (20). Bins 053/054 use the S2 replacements Read053B188/Read054B189;
the remaining barcodes are reverse complements of S1 readout binding segments.

The reference is UCSC ce11 / WBcel235 (WS245),
[whole-genome FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz),
compressed SHA-256 `b5f80e37d22eb8f1cce1efff09a7e2e2ac2b6c73650f674c41c0b2e281c97241`.
Other checksums and dependency versions are in [sources.json](report/sources.json)
and [provenance.json](report/provenance.json).
Published inclusive bin boundaries cover chrII:11,500,001–12,500,000:
99 intervals of 10,001 bp followed by 9,901 bp. The paper names the nominal
region as 11,500,001–12,500,001; configuration follows the explicit S2 bin metadata.

## Results

Metrics use the unchanged [common evaluator](../common/sequence_recovery.py),
which deduplicates identical target/role/full-DNA records. CPSR is N/A.

| Stage | Candidates | Exact matches / 9,975 | SRR | TRR | MTR | Candidates per bin |
|---|---:|---:|---:|---:|---:|---:|
| candidate | 1,972,876 | 9,879 | 99.037594% | 100% | 99.037893% | 17,924–19,924 |
| filtered | 1,303,762 | 9,676 | 97.002506% | 100% | 97.002080% | 7,304–17,628 |
| non_overlap | 19,344 | 135 | 1.353383% | 31% | 1.353451% | 130–235 |

All stages retain candidates in 100/100 bins; non-overlap matches published oligos
in 31/100 bins. Candidate recovery covers 9,888 original rows including duplicates.
See [per-bin counts](report/bin_counts.csv), [summary](report/summary.csv),
[per-target metrics](report/per_gene.csv), [sequence matches](report/sequence_matches.csv)
and [pair matches](report/pair_matches.csv).

Native genome extraction and YAML/DAG assembly enumerate both directions of
996,100 ACGT windows. Identical full oligos within a bin are collapsed while
preserving independent sites and native IDs. No reference binding DNA or inferred
site coordinates enter generation.

Quality filtering reads the saved candidates: one native Bowtie 1 call (`-v 0 -a`,
both strands) requires one physical exact site across all seven ce11 records;
native GC filtering uses 30%–70%. Native primer3/SantaLucia reports Tm under
the explicit benchmark conditions in YAML without imposing a threshold: no
reproducible paper design cutoff is given. These are benchmark settings.

Native `remove_overlap` uses ascending stable order of genomic start, end, strand
(+ before −), full DNA and native ID, with zero gap. Both strands occupy the same
physical interval. No published matches or per-bin count influence selection.
Its recovery measures selection retention, not quality.
[Filtering details](report/filtering.json) give retention fractions, Tm distribution
and stepwise losses.

## Reproduction

Use an installed U-Probe environment with Bowtie 1.3.1 on `PATH`; dependency
versions are recorded in provenance. Run from the repository root into a fresh
validation directory. Raw data, candidates and indexes remain outside version control.

```bash
export PYTHONPATH=.
m=recovery/oligopaint_2020_chrII_1mb
raw=../raw_data
run=../validation
mkdir -p "$raw" "$run/reference" "$run/index"
curl -fL 'https://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.3000817.s009&type=supplementary' -o "$raw/s009.fa"
curl -fL 'https://journals.plos.org/plosbiology/article/file?id=10.1371/journal.pbio.3000817.s010&type=supplementary' -o "$raw/s010.xlsx"
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz -o "$raw/ce11.fa.gz"
gzip -dc "$raw/ce11.fa.gz" > "$run/reference/ce11.fa"
bowtie-build --threads 4 "$run/reference/ce11.fa" "$run/index/ce11"
python "$m/prepare_inputs.py" --raw "$raw" --output "$run/inputs" --config "$m/design.yaml"
python "$m/generate_full.py" --config "$m/design.yaml" --fasta "$run/reference/ce11.fa" --output "$run/candidate"
python "$m/evaluate_stages.py" --reference "$run/inputs/reference_oligos.csv" --generated "$run/candidate/candidate.csv.gz" --stage candidate --report "$run/report"
cp "$run/inputs/sources.json" "$run/candidate/generation.json" "$run/report/"
python "$m/filter_candidates.py" --candidate "$run/candidate/candidate.csv.gz" --reference "$run/inputs/reference_oligos.csv" --config "$m/design.yaml" --index "$run/index/ce11" --output "$run/filtering" --report "$run/report" --threads 4
python "$m/evaluate_stages.py" --diagnose --reference "$run/inputs/reference_oligos.csv" --fasta "$run/reference/ce11.fa" --config "$m/design.yaml" --exact-counts "$run/filtering/exact_counts.csv.gz" --report "$run/report"
```

## Limitations

The 96 unrecovered unique binding sequences have no exact site on ce11 chrII;
the underlying assembly, sample or sequence differences are unresolved. No second
assembly was tested. [Unrecovered records](report/unrecovered.json) preserve the evidence.
Of the recovered unique references, quality filtering removes 185 for multiple
genomic exact sites and 18 for GC. Genome exact-site uniqueness is not a complete
hybridization specificity prediction. Paper-specific ranking, final-count selection,
secondary structure and experimental performance are not evaluated. Novel candidates
are not treated as false positives. Final Recovery: **Not evaluated**.
