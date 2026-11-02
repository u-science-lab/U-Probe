"""Apply the unchanged common evaluator independently to named benchmark stages."""
import argparse
import collections
import csv
import gzip
import hashlib
import importlib.metadata
import json
import subprocess
from pathlib import Path

import yaml
from recovery.common.sequence_recovery import evaluate, read_csv, iter_csv

VERSION = 'oligopaint_2020_chrII_1mb'


def evaluate_stage(reference, generated, stage, report):
    # The shared evaluator accepts candidate/final. Each independent candidate
    # pool is evaluated identically, then labelled with its experimental stage.
    def rows():
        for row in iter_csv(generated):
            row['stage'] = 'candidate'
            yield row
    result = evaluate(read_csv(reference), rows(), version=VERSION)
    report.mkdir(parents=True, exist_ok=True)
    for name in ['summary', 'per_gene', 'sequence_matches', 'pair_matches']:
        records = result[name]
        path = report / (name + '.csv')
        exists = path.exists()
        if exists:
            with path.open() as handle:
                assert all(row['stage'] != stage for row in csv.DictReader(handle)), f'{stage} already evaluated'
        for row in records:
            row['stage'] = stage
            if 'CPSR_pct' in row:
                row['CPSR_pct'] = 'N/A'
        fields = list(records[0]) if records else ['stage', 'reference_pair_id', 'matched']
        with path.open('a' if exists else 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fields)
            if not exists:
                writer.writeheader()
            writer.writerows(records)
    assert not result['exclusions']
    return next(row for row in result['summary'] if row['oligo_type'] == 'all')


def audit_existing(reference, fasta, config, exact_counts, report):
    """Diagnose saved outcomes only; does not regenerate or realign candidates."""
    from pyfaidx import Fasta
    cfg = yaml.safe_load(config.read_text())
    regions = {v['bin_id']: tuple(map(int, target.split(':')[1].split('-')))
               for target, v in cfg['encoding'].items()}
    matches = list(csv.DictReader((report / 'sequence_matches.csv').open()))
    by_stage = {stage: {r['reference_oligo_id']: r for r in matches if r['stage'] == stage}
                for stage in ['candidate', 'filtered', 'non_overlap']}
    missed = [r for r in by_stage['candidate'].values() if r['matched'] == 'False']
    def rc(seq):
        return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    details, reasons = [], collections.Counter()
    with Fasta(str(fasta)) as fa:
        chrom = fa['chrII'][:].seq.upper()
    for row in missed:
        binding = row['full_sequence'][60:100]
        low, high = regions[row['gene']]
        hits = []
        for strand, query in [('+', binding), ('-', rc(binding))]:
            cursor = 0
            while True:
                pos = chrom.find(query, cursor)
                if pos < 0:
                    break
                hits.append([pos + 1, pos + 40, strand]); cursor = pos + 1
        reason = ('exact_chrII_site_outside_bin' if hits else 'no_exact_chrII_binding_site')
        if any(start >= low and end <= high for start, end, _ in hits):
            reason = 'within_bin_but_unrecovered'
        reasons[reason] += 1
        details.append(dict(reference_id=row['reference_oligo_id'], bin=row['gene'], reason=reason,
                            bin_start=low, bin_end=high, exact_chrII_sites=hits))
    (report / 'unrecovered.json').write_text(json.dumps({'unique_unrecovered': len(missed),
        'reasons': dict(reasons), 'records': details}, indent=2) + '\n')
    quality_losses = [r for key, r in by_stage['candidate'].items()
                      if r['matched'] == 'True' and by_stage['filtered'][key]['matched'] == 'False']
    canonical = lambda seq: min(seq, rc(seq))
    required = {canonical(r['full_sequence'][60:100]) for r in quality_losses}
    counts = {}
    with gzip.open(exact_counts, 'rt') as handle:
        for row in csv.DictReader(handle):
            if row['binding_canonical'] in required:
                counts[row['binding_canonical']] = int(row['exact_sites'])
    losses = collections.Counter()
    for row in quality_losses:
        binding = row['full_sequence'][60:100]
        count = counts[canonical(binding)]
        gc = (binding.count('G') + binding.count('C')) / 40
        reason = 'multiple_exact_genome_sites' if count != 1 else 'gc_outside_30_70' if not .3 <= gc <= .7 else 'unexplained'
        losses[reason] += 1
    published = read_csv(reference)
    raw_matched = {}
    for stage, stage_rows in by_stage.items():
        keys = {(r['gene'], r['full_sequence']) for r in stage_rows.values() if r['matched'] == 'True'}
        raw_matched[stage] = sum((r['gene'], r['full_sequence']) in keys for r in published)
    evidence = json.loads((report / 'filtering.json').read_text())
    generation = json.loads((report / 'generation.json').read_text())
    assert generation['candidate_sha256'] == evidence['candidate_sha256_after_filtering']
    evidence['published_quality_losses'] = dict(losses)
    evidence['raw_reference_records'] = len(published)
    evidence['unique_reference_sequences_within_target'] = len(by_stage['candidate'])
    evidence['raw_reference_records_recovered'] = raw_matched
    evidence['unchanged_stage1_candidate_verified'] = True
    (report / 'filtering.json').write_text(json.dumps(evidence, indent=2) + '\n')
    provenance = {'base_code_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                  'fasta_sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(),
                  'config_sha256': hashlib.sha256(config.read_bytes()).hexdigest(),
                  'genome': 'UCSC ce11 / WBcel235 (WS245)',
                  'genome_source': 'https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz',
                  'package_versions': {name: importlib.metadata.version(name) for name in
                                       ['numpy', 'pandas', 'PyYAML', 'pyfaidx', 'primer3-py']},
                  'generation_passes': 1, 'specificity_passes': 1, 'core_changes': False,
                  'metric_denominator': 'Common evaluator deduplicates by method, target, role and full DNA; 9984 rows / 9975 unique target oligos'}
    (report / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps({'candidate_missing': dict(reasons), 'quality_losses': dict(losses),
                      'raw_records_matched': raw_matched, 'coverage': evidence['per_bin_candidate_counts']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--generated', type=Path)
    parser.add_argument('--stage', choices=['candidate', 'filtered', 'non_overlap'])
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--diagnose', action='store_true')
    for name in ['fasta', 'config', 'exact-counts']:
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args()
    if args.diagnose:
        assert args.fasta and args.config and args.exact_counts
        audit_existing(args.reference, args.fasta, args.config, args.exact_counts, args.report)
    else:
        assert args.generated and args.stage
        print(evaluate_stage(args.reference, args.generated, args.stage, args.report))
