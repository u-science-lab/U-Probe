#!/usr/bin/env python3
"""Evaluate STARmap PLUS 2023 64-gene sequence recovery from generated long CSV."""
import argparse
import json
import sys
from pathlib import Path

METHOD_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(METHOD_ROOT.parent / 'common'))
from sequence_recovery import evaluate, iter_csv, read_csv, write_csv

VERSION = 'STARmap_PLUS_2023_64gene'


def evaluate_files(reference_oligos, reference_pairs, generated, output):
    refs, pairs = read_csv(reference_oligos), read_csv(reference_pairs)
    if any(r.get('oligo_type') not in ('padlock', 'primer') for r in refs):
        raise ValueError('STARmap references must contain padlock and primer roles')
    result = evaluate(refs, iter_csv(generated), pairs, VERSION)
    # Keep the archived method report schema; generic results use member_matches.
    method_pairs = []
    for row in result['pair_matches']:
        members = json.loads(row['member_matches'])
        if set(members) != {'padlock', 'primer'}:
            raise ValueError('STARmap reference group must contain padlock and primer')
        converted = {k: v for k, v in row.items() if k not in ('member_matches', 'all_members_individually_matched', 'matched', 'generated_pair_ids')}
        converted.update(padlock_matched=members['padlock'], primer_matched=members['primer'],
                         all_members_individually_matched=row['all_members_individually_matched'],
                         matched=row['matched'], generated_pair_ids=row['generated_pair_ids'])
        method_pairs.append(converted)
    result['pair_matches'] = method_pairs
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name, records in result.items():
        if name == 'exclusions':
            (output/'reference_exclusions.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
        else:
            write_csv(output/(name+'.csv'), records, list(records[0]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-oligos', type=Path, default=METHOD_ROOT/'reference/reference_oligos.csv')
    parser.add_argument('--reference-pairs', type=Path, default=METHOD_ROOT/'reference/reference_pairs.csv')
    parser.add_argument('--generated', type=Path, required=True, help='Generated long CSV or .csv.gz')
    parser.add_argument('--output', type=Path, default=METHOD_ROOT/'report')
    args = parser.parse_args()
    result = evaluate_files(args.reference_oligos, args.reference_pairs, args.generated, args.output)
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    main()
