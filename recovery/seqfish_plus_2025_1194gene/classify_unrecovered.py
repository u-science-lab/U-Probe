#!/usr/bin/env python3
'''Classify unmatched reference probe bodies using existing native candidate shards.'''
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--matches', type=Path, required=True)
p.add_argument('--run-dir', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
with a.matches.open() as f:
    rows = [r for r in csv.DictReader(f) if r['matched'].lower() == 'false']
by_gene = defaultdict(list)
for r in rows:
    by_gene[r['gene']].append(r)
details = {d['gene'].casefold(): d for d in json.loads((a.run_dir/'run_details.json').read_text())['genes']}
counts = Counter()
for gene, records in by_gene.items():
    detail = details[gene.casefold()]
    wanted = {r['full_sequence'][30:65] for r in records}
    present = set()
    with gzip.open(a.run_dir/'candidates_by_gene'/(detail['gene']+'.csv.gz'), 'rt') as f:
        for candidate in csv.DictReader(f):
            binding = candidate['full_sequence'][30:65]
            if binding in wanted:
                present.add(binding)
    for r in records:
        reason = ('reference_annotation_missing' if not detail['transcripts'] else
                  'readout_arrangement_or_sequence_mismatch' if r['full_sequence'][30:65] in present else
                  'binding_not_in_annotated_rna_candidates')
        r['reason'] = reason
        counts[reason] += 1
a.output.parent.mkdir(parents=True, exist_ok=True)
with a.output.open('w', newline='') as f:
    fields = ['gene', 'reference_oligo_id', 'reason']
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    w.writerows(rows)
print(json.dumps(dict(counts), sort_keys=True))
