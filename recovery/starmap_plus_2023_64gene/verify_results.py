#!/usr/bin/env python3
"""Check compact STARmap recovery reports; optionally independently rescan generated CSV."""
import argparse
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

METHOD_ROOT = Path(__file__).resolve().parent
VERSION = 'STARmap_PLUS_2023_64gene'


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        yield from csv.DictReader(handle)


def dna(value):
    value = re.sub(r'\s+', '', value).upper()
    value = value.replace('/5PHOS/', '')
    if not value or re.fullmatch('[ACGT]+', value) is None:
        raise ValueError('Invalid DNA or unsupported modification token')
    return value


def require(condition, message):
    if not condition:
        raise ValueError(message)


def boolean(value):
    require(value in ('True', 'False'), 'Expected an explicit True/False report value')
    return value == 'True'


def ids(value):
    values = value.split(';') if value else []
    require(len(values) == len(set(values)) and all(values), 'Duplicate or empty linked ID')
    return set(values)


def number(row, name, expected):
    require(name in row and row[name] != '', 'Missing numeric field: '+name)
    require(math.isclose(float(row[name]), expected, rel_tol=0, abs_tol=1e-9),
            'Incorrect '+name+': '+str(row))


def verify(reference_oligos, reference_pairs, report, generated=None):
    report = Path(report)
    reference = {}
    ref_links = defaultdict(dict)
    original = list(rows(reference_oligos))
    for row in original:
        role, gene = row['oligo_type'], row['gene'].casefold()
        require(role in ('padlock', 'primer') and bool(gene), 'Invalid reference gene/role')
        seq = dna(row['full_sequence'])
        if row.get('dna_sequence'):
            require(dna(row['dna_sequence']) == seq, 'Reference DNA fields disagree')
        key = (gene, role, seq)
        reference.setdefault(key, row)
        group = (gene, row['pair_id'])
        require(role not in ref_links[group] or ref_links[group][role] == key,
                'Conflicting member in reference group')
        ref_links[group][role] = key
    wide = list(rows(reference_pairs))
    seen = set()
    for row in wide:
        group = (row['gene'].casefold(), row['pair_id'])
        require(group not in seen and set(ref_links.get(group, {})) == {'padlock', 'primer'},
                'Invalid or duplicate reference pair')
        seen.add(group)
        for role in ('padlock', 'primer'):
            require(dna(row.get(role+'_full') or row[role+'_dna']) == ref_links[group][role][2],
                    'Wide/long reference sequence mismatch')
    require(seen == set(ref_links), 'Wide/long reference pair sets differ')
    summaries = list(rows(report/'summary.csv'))
    stages = {r['stage'] for r in summaries}
    require(bool(stages) and stages <= {'candidate', 'final'}, 'Invalid report stage')
    seq_rows = list(rows(report/'sequence_matches.csv'))
    pair_rows = list(rows(report/'pair_matches.csv'))
    gene_rows = list(rows(report/'per_gene.csv'))
    matched = defaultdict(set)
    evidence = {}
    physical = {}
    for row in seq_rows:
        stage = row['stage']
        require(stage in stages and row['method_version'] == VERSION, 'Sequence version/stage mismatch')
        key = (row['gene'].casefold(), row['oligo_type'], dna(row['full_sequence']))
        require(key in reference and row['full_sequence'] == key[2], 'Reported DNA/gene/role is not a reference')
        ref = reference[key]
        require(row['reference_oligo_id'] == ref.get('oligo_id', ref.get('probe_id', '')) and
                row['reference_pair_id'] == ref['pair_id'], 'Reference identifier mismatch')
        require((stage, key) not in evidence, 'Duplicate sequence report record')
        probe_ids, pair_ids = ids(row['generated_probe_ids']), ids(row['generated_pair_ids'])
        require(probe_ids == {p+'_'+key[1] for p in pair_ids}, 'Generated probe/pair/role links disagree')
        require(boolean(row['matched']) == bool(pair_ids), 'Sequence matched flag disagrees with evidence')
        evidence[(stage, key)] = (probe_ids, pair_ids)
        if pair_ids:
            matched[stage].add(key)
        for pair in pair_ids:
            physical_key = (stage, key[0], pair, key[1])
            require(physical_key not in physical or physical[physical_key] == key[2],
                    'One generated pair/role claims different full DNA sequences')
            physical[physical_key] = key[2]
    require(set(evidence) == {(stage, key) for stage in stages for key in reference}, 'Missing sequence report rows')
    recovered_pairs = defaultdict(set)
    pair_seen = set()
    for row in pair_rows:
        stage, group = row['stage'], (row['gene'].casefold(), row['reference_pair_id'])
        require(stage in stages and group in ref_links and (stage, group) not in pair_seen,
                'Unknown or duplicate pair report record')
        require(row['method_version'] == VERSION and row['required_roles'] == 'padlock;primer' and
                boolean(row['cpsr_applicable']), 'Incorrect pair method/required roles')
        pair_seen.add((stage, group))
        members = ref_links[group]
        pad = evidence[(stage, members['padlock'])][1]
        pri = evidence[(stage, members['primer'])][1]
        linked = pad & pri
        require(ids(row['generated_pair_ids']) == linked, 'Reported complete pair links are not member intersection')
        require(boolean(row['padlock_matched']) == bool(pad) and boolean(row['primer_matched']) == bool(pri) and
                boolean(row['all_members_individually_matched']) == bool(pad and pri) and
                boolean(row['matched']) == bool(linked), 'Incorrect pair recovery flag')
        if linked:
            recovered_pairs[stage].add(group)
    require(pair_seen == {(stage, group) for stage in stages for group in ref_links}, 'Missing pair report rows')
    summary_seen, gene_seen, computed = set(), set(), []
    for stage in sorted(stages):
        for role in ('all', 'padlock', 'primer'):
            selected = {k for k in reference if role == 'all' or k[1] == role}
            totals, hits = Counter(k[0] for k in selected), Counter(k[0] for k in selected & matched[stage])
            expected = dict(reference_sequences=len(selected), matched_sequences=sum(hits.values()),
                SRR_pct=100*sum(hits.values())/len(selected), reference_targets=len(totals),
                recovered_targets=len(hits), TRR_pct=100*len(hits)/len(totals),
                MTR_pct=100*sum(hits[g]/n for g,n in totals.items())/len(totals))
            if role == 'all':
                expected.update(reference_pairs=len(ref_links), matched_pairs=len(recovered_pairs[stage]),
                                CPSR_pct=100*len(recovered_pairs[stage])/len(ref_links))
            found = [r for r in summaries if r['stage']==stage and r['oligo_type']==role]
            require(len(found)==1, 'Missing/duplicate summary role')
            row = found[0]
            if role != 'all':
                require(all(row[n]=='' for n in ('reference_pairs','matched_pairs','CPSR_pct')), 'Role-only CPSR must be inapplicable')
            summary_seen.add((stage, role))
            for name, value in expected.items():
                number(row, name, value)
            number(row, 'excluded_reference_records', 0)
            number(row, 'duplicate_reference_records', len(original)-len(reference))
            number(row, 'generated_unique_reference_hits', len(matched[stage]))
            for gene, total in totals.items():
                entries = [r for r in gene_rows if r['stage']==stage and r['oligo_type']==role and r['gene'].casefold()==gene]
                require(len(entries)==1, 'Missing/duplicate per-gene row')
                item = entries[0]
                gene_seen.add((stage, role, gene))
                number(item, 'reference_sequences', total)
                number(item, 'matched_sequences', hits[gene])
                number(item, 'SRR_pct', 100*hits[gene]/total)
                if role == 'all':
                    n = sum(g==gene for g,p in ref_links)
                    m = sum(g==gene for g,p in recovered_pairs[stage])
                    for name,value in dict(reference_pairs=n, matched_pairs=m, CPSR_pct=100*m/n).items():
                        number(item, name, value)
                else:
                    require(all(item[n]=='' for n in ('reference_pairs','matched_pairs','CPSR_pct')), 'Role-only per-gene CPSR must be inapplicable')
            computed.append(dict(stage=stage, oligo_type=role, **expected))
    require(len(summaries)==len(summary_seen) and len(gene_rows)==len(gene_seen), 'Unexpected extra summary/per-gene rows')
    for stage in stages:
        counts_reported = {int(r['generated_records']) for r in summaries if r['stage']==stage}
        require(len(counts_reported)==1 and next(iter(counts_reported))>=len(matched[stage]), 'Inconsistent generated row counts')
    counts = None
    if generated:
        counts, actual = Counter(), defaultdict(lambda: (set(), set()))
        current, members = None, {}
        def finish_pair():
            require(set(members)=={'padlock','primer'}, 'Generated pair must have two adjacent roles')
        for row in rows(generated):
            stage, gene, role = row.get('stage') or 'candidate', row['gene'].casefold(), row['oligo_type'].casefold()
            require(stage in stages and role in ('padlock','primer'), 'Unexpected generated stage/role')
            require((row.get('method_version') or VERSION)==VERSION, 'Generated method version mismatch')
            pair, probe = row['pair_id'], row.get('probe_id') or row.get('oligo_id')
            require(bool(pair) and probe==pair+'_'+role, 'Invalid generated probe/pair links')
            seq = dna(row['full_sequence'])
            key = (gene, role, seq)
            counts[stage] += 1
            group = (stage, gene, pair)
            if current != group:
                if current is not None:
                    finish_pair()
                current, members = group, {}
            require(role not in members, 'Duplicate generated role in one pair')
            members[role] = seq
            if key in reference:
                actual[(stage,key)][0].add(probe)
                actual[(stage,key)][1].add(pair)
        if current is not None:
            finish_pair()
        for key, value in evidence.items():
            require(actual[key]==value, 'Full generated scan disagrees with compact sequence evidence')
        for row in summaries:
            number(row, 'generated_records', counts[row['stage']])
    return dict(status='passed', scope='independent_full_generated_scan' if generated else 'compact_report_consistency',
        full_generated_scan_performed=bool(generated),
        limitation=None if generated else 'Reported hit IDs are evidence, not proof of the contents of an unprovided generated file.',
        reference_sequences=len(reference), reference_pairs=len(ref_links), checked_sequence_rows=len(seq_rows),
        checked_pair_rows=len(pair_rows), generated_records=dict(counts) if counts is not None else None,
        unverified_fields=[] if generated else ['generated_records'], recomputed_metrics=computed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-oligos', type=Path, default=METHOD_ROOT/'reference_oligos.csv')
    parser.add_argument('--reference-pairs', type=Path, default=METHOD_ROOT/'reference_pairs.csv')
    parser.add_argument('--report', type=Path, default=METHOD_ROOT/'report')
    parser.add_argument('--generated', type=Path, help='Optional complete long CSV/.csv.gz; expensive for the archived run')
    parser.add_argument('--output', type=Path, help='Optional verification JSON output')
    args = parser.parse_args()
    result = verify(args.reference_oligos, args.reference_pairs, args.report, args.generated)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
