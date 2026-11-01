#!/usr/bin/env python3
"""Classify shared-evaluator misses using independently generated window shards."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matches', type=Path, required=True)
    p.add_argument('--run-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    with args.matches.open(newline='') as f:
        misses = [r for r in csv.DictReader(f) if r['stage'] == 'candidate' and r['matched'] == 'False']
    coverage = {r['gene']: r for r in json.loads((args.run_dir/'coverage.json').read_text())}
    output = []
    for gene in sorted({r['gene'] for r in misses}):
        records = [r for r in misses if r['gene'] == gene]
        needed = {}
        for r in records:
            seq = r['full_sequence']
            needed[r['reference_oligo_id']] = rc(seq[:13])+rc(seq[-13:]) if r['oligo_type'] == 'padlock' else rc(seq[:13])
        present = set()
        shard = args.run_dir/'candidates_by_gene'/(gene+'.windows.txt.gz')
        if shard.exists():
            with gzip.open(shard, 'rt') as f:
                for line in f:
                    window = line.strip()
                    for r in records:
                        key = r['reference_oligo_id']
                        region = window[:26] if r['oligo_type'] == 'padlock' else window[-13:]
                        if region == needed[key]:
                            present.add(key)
        for r in records:
            key = r['reference_oligo_id']
            if not coverage.get(gene, {}).get('transcripts', 0):
                reason = 'transcript_not_covered'
                # Post-run annotation evidence; no alias substitution during generation.
                evidence = ('Input mm10 GTF annotates Ccn3 as Nov (ENSMUSG00000037362.8); '
                            'gene-name alias mapping was not applied, so no Ccn3 transcript was extracted'
                            if gene == 'Ccn3' else 'No same-name annotated transcript in the input GTF')
            elif key not in present:
                reason = 'target_site_difference'
                evidence = 'Published binding region absent from legal same-gene generated windows'
            elif r['oligo_type'] == 'rca' and r['full_sequence'][13:-11] not in {'A','C','G','T','AC','CC','GC','TC'}:
                reason = 'rca_spacer_difference'
                evidence = 'Binding region present; published spacer outside generated variants'
            else:
                reason = 'probe_structure_difference'
                evidence = 'Binding region present; complete oligo absent'
            output.append(dict(gene=gene, oligo_type=r['oligo_type'], reference_oligo_id=key,
                               reference_pair_id=r['reference_pair_id'], reason=reason, evidence=evidence))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene','oligo_type','reference_oligo_id','reference_pair_id','reason','evidence'])
        writer.writeheader()
        writer.writerows(output)
    print(json.dumps({reason: sum(r['reason'] == reason for r in output) for reason in sorted({r['reason'] for r in output})}))


if __name__ == '__main__':
    main()
