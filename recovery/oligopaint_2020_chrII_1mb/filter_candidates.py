"""Filter the saved Stage 1 pool once; never extracts or constructs probes."""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from uprobe.core.attributes import add_attributes
from uprobe.core.attributes._attributes import _bowtie_alignments
from uprobe.core.process import filter_table, remove_overlap

from evaluate_stages import evaluate_stage
from generate_full import FIELDS


def rc(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def canonical(seq):
    return min(seq, rc(seq))


def rows(path):
    with gzip.open(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle)


def specificity(candidate, index, output, threads):
    # Search each binding / reverse-complement equivalence class once, across
    # all bins. The existing aligner emits every exact alignment on both strands.
    sequences = {}
    for row in rows(candidate):
        seq = canonical(row['binding'])
        if seq not in sequences:
            sequences[seq] = str(len(sequences))
    queries = {identifier: seq for seq, identifier in sequences.items()}
    counts = collections.Counter()
    palindrome_sites = {}
    for query, strand, chrom, offset, mismatches in _bowtie_alignments(
            str(output), queries, 'exact', str(index), 0, threads, 'both'):
        assert mismatches == 0 and chrom in {'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX', 'chrM'}
        # A palindromic 40-mer aligns twice at the same physical locus. Count it once.
        if queries[query] == rc(queries[query]):
            sites = palindrome_sites.setdefault(query, set())
            if (chrom, offset) in sites:
                continue
            sites.add((chrom, offset))
        counts[query] += 1
    assert all(counts[q] >= 1 for q in queries), 'FASTA-derived query lacks its exact reference site'
    with gzip.open(output / 'exact_counts.csv.gz', 'wt', newline='') as handle:
        writer = csv.writer(handle); writer.writerow(['binding_canonical', 'exact_sites'])
        writer.writerows((seq, counts[q]) for seq, q in sequences.items())
    return {seq: counts[q] for seq, q in sequences.items()}


def distributions(values):
    array = np.asarray(values, dtype=float)
    return dict(zip(['min', 'p05', 'p25', 'median', 'p75', 'p95', 'max'],
                    map(float, np.quantile(array, [0, .05, .25, .5, .75, .95, 1])))) if len(array) else None


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    cfg = yaml.safe_load(args.config.read_text())
    evidence = {'status': 'running', 'specificity': {'aligner': 'Bowtie 1', 'mismatches': 0,
                'all_alignments': True, 'strand': 'both', 'unique_rule': 'one physical whole-genome exact site',
                'query_pool': 'all Stage 1 binding sequences, canonical reverse-complement deduplication',
                'aligner_invocations': 1}, 'gc_fraction_range_inclusive': [0.30, 0.70],
                'tm': {'method': 'U-Probe add_attributes / primer3 SantaLucia',
                       'conditions': cfg['hybridization'], 'filter': None,
                       'reason': 'Article gives experimental conditions but no reproducible design Tm threshold; benchmark reporting conditions only'},
                'non_overlap_order': 'genomic start, end, + before -, full DNA, native ID; ascending stable order',
                'non_overlap_gap_nt': 0, 'final_recovery': 'Not evaluated',
                'stage2_extraction_calls': 0, 'stage2_assembly_calls': 0}
    (args.report / 'filtering.json').write_text(json.dumps(evidence, indent=2))
    try:
        print('specificity started', flush=True)
        hits = specificity(args.candidate, args.index, args.output, args.threads)
        evidence['specificity']['unique_queries'] = len(hits)
        evidence['specificity']['query_site_count_histogram'] = dict(collections.Counter(hits.values()))
        evidence['specificity']['version'] = subprocess.check_output(['bowtie', '--version'], text=True).splitlines()[0].split('/')[-1]
        print('specificity completed', len(hits), flush=True)
        bins, tm_values, totals = [], [], collections.Counter()
        extra = ['exact_sites', 'gc', 'tm']
        with gzip.open(args.output / 'filtered.csv.gz', 'wt', newline='', compresslevel=1) as qf, \
             gzip.open(args.output / 'non_overlap.csv.gz', 'wt', newline='', compresslevel=1) as nf:
            quality_writer = csv.DictWriter(qf, FIELDS + extra, extrasaction='ignore')
            selected_writer = csv.DictWriter(nf, FIELDS + extra, extrasaction='ignore')
            quality_writer.writeheader(); selected_writer.writeheader()
            for gene, group in itertools.groupby(rows(args.candidate), lambda row: row['gene']):
                frame = pd.DataFrame(list(group))
                for col in ['start', 'end']:
                    frame[col] = frame[col].astype(int)
                before = len(frame)
                frame['exact_sites'] = frame.binding.map(lambda s: hits[canonical(s)])
                frame = filter_table(frame, {'specificity': cfg['post_process']['filters']['unique_exact']})
                specific = len(frame)
                frame = add_attributes(frame.copy(), {'attributes': {'gc': cfg['attributes']['gc']}}, {})
                frame = filter_table(frame, {'gc': cfg['post_process']['filters']['gc_range']})
                frame = add_attributes(frame.copy(), {'attributes': {'tm': cfg['attributes']['tm']},
                                                      'hybridization': cfg['hybridization']}, {})
                quality = len(frame)
                tm_values.extend(frame.tm.tolist())
                frame['stage'] = 'filtered'
                quality_writer.writerows(frame.to_dict('records'))
                frame['target'] = gene
                frame['target_region'] = frame.binding
                frame['_strand_order'] = frame.strand.map({'+': 0, '-': 1})
                frame = frame.sort_values(['start', 'end', '_strand_order', 'full_sequence', 'native_probe_id'], kind='stable')
                selected = remove_overlap(frame, cfg, location_interval=0)
                assert all(int(b) > int(a) for a, b in zip(selected.end.iloc[:-1], selected.start.iloc[1:])), gene
                selected['stage'] = 'non_overlap'
                selected_writer.writerows(selected.to_dict('records'))
                detail = dict(bin=gene, candidate=before, specificity=specific, filtered=quality,
                              non_overlap=len(selected), specificity_removed=before-specific,
                              gc_removed=specific-quality, tm_removed=0, overlap_removed=quality-len(selected))
                bins.append(detail); totals.update({k: v for k, v in detail.items() if k != 'bin'})
                if len(bins) % 10 == 0:
                    print('filtered bins', len(bins), flush=True)
        assert len(bins) == 100
        evidence['tm']['quality_distribution_celsius'] = distributions(tm_values)
        evidence['counts'] = dict(totals)
        evidence['candidate_retention'] = {stage: totals[stage] / totals['candidate']
                                            for stage in ['specificity', 'filtered', 'non_overlap']}
        evidence['metrics'] = {}
        for stage in ['filtered', 'non_overlap']:
            evidence['metrics'][stage] = evaluate_stage(args.reference, args.output / (stage + '.csv.gz'), stage, args.report)
        # Audit whole-region non-overlap, including adjacent bins and both DNA directions.
        occupied = sorted((int(r['start']), int(r['end'])) for r in rows(args.output / 'non_overlap.csv.gz'))
        assert all(b[0] > a[1] for a, b in zip(occupied, occupied[1:]))
        per_gene = list(csv.DictReader((args.report / 'per_gene.csv').open()))
        recoveries = {(r['stage'], r['gene']): r for r in per_gene if r['oligo_type'] == 'all'}
        for detail in bins:
            for stage in ['candidate', 'filtered', 'non_overlap']:
                r = recoveries[(stage, detail['bin'])]
                detail['reference_count'] = int(r['reference_sequences'])
                detail[stage + '_matched'] = int(r['matched_sequences'])
        with (args.report / 'bin_counts.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, list(bins[0])); writer.writeheader(); writer.writerows(bins)
        evidence['per_bin_candidate_counts'] = {stage: {'covered_bins': sum(b[stage] > 0 for b in bins),
            'min': min(b[stage] for b in bins), 'max': max(b[stage] for b in bins)}
            for stage in ['candidate', 'filtered', 'non_overlap']}
        evidence['candidate_sha256_after_filtering'] = hashlib.sha256(args.candidate.read_bytes()).hexdigest()
        evidence['status'] = 'completed'
    except Exception as error:
        evidence['status'] = 'incomplete'
        evidence['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        (args.report / 'filtering.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'counts': evidence['counts'], 'metrics': evidence['metrics']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ['candidate', 'reference', 'config', 'index', 'output', 'report']:
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--threads', type=int, default=4)
    run(parser.parse_args())
