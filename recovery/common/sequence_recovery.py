#!/usr/bin/env python3
"""Exact sequence recovery; stdlib only. Never treats unpublished output as false positives."""
import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path


def normalize_sequence(value):
    """Remove whitespace and explicit /modification/ tokens, never DNA letters."""
    value = re.sub(r'\s+', '', str(value or '')).upper()
    mods = re.findall(r'/([^/]+)/', value)
    if any(re.fullmatch('[ACGT]+', token) for token in mods):
        raise ValueError('DNA-only slash token cannot be interpreted as a chemical modification')
    dna = re.sub(r'/[^/]+/', '', value)
    if not dna or re.fullmatch('[ACGT]+', dna) is None:
        raise ValueError('Missing/invalid full DNA sequence (only A,C,G,T accepted)')
    return dna, mods


def read_csv(path):
    return list(iter_csv(path))


def iter_csv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, mode='rt', newline='', encoding='utf-8-sig') as handle:
        yield from csv.DictReader(handle)


def field(row, *names, default=''):
    return next((str(row[n]).strip() for n in names if row.get(n) not in (None, '')), default)


def canonical(row, default_version):
    gene = field(row, 'target_id', 'gene')
    role = field(row, 'oligo_type').casefold()
    if not gene or not role:
        raise ValueError('Missing target or oligo_type')
    # full_sequence is authoritative; dna_sequence is a consistency check.
    dna, mods = normalize_sequence(field(row, 'full_sequence', 'dna_sequence'))
    if row.get('dna_sequence'):
        declared, _ = normalize_sequence(row['dna_sequence'])
        if declared != dna:
            raise ValueError('dna_sequence disagrees with full_sequence')
    version = field(row, 'method_version', default=default_version)
    return (version, gene.casefold(), role, dna), gene, mods


def write_csv(path, rows, columns):
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def evaluate(reference_rows, generated_rows, pair_rows=None, version='unspecified'):
    references = {}
    excluded = []
    ref_groups = defaultdict(dict)
    duplicates = 0
    for index, row in enumerate(reference_rows, 2):
        try:
            key, gene, mods = canonical(row, version)
        except ValueError as error:
            excluded.append({'source': 'reference', 'csv_row': index, 'reason': str(error)})
            continue
        pair = field(row, 'probe_set_id', 'pair_id')
        if key in references:
            duplicates += 1
        else:
            references[key] = dict(row=row, gene=gene, mods=mods)
        if pair:
            group = (key[0], key[1], pair)
            if key[2] in ref_groups[group] and ref_groups[group][key[2]] != key:
                raise ValueError('Reference group has multiple different members for one role: ' + str(group))
            ref_groups[group][key[2]] = key
    if not references:
        raise ValueError('No evaluable reference records')
    # Wide pair table must agree exactly with long oligo table.
    if pair_rows is not None:
        seen_pairs = set()
        for row in pair_rows:
            gene = field(row, 'target_id', 'gene').casefold()
            pair = field(row, 'probe_set_id', 'pair_id')
            group = (field(row, 'method_version', default=version), gene, pair)
            if group in seen_pairs:
                raise ValueError('Duplicate reference pair ID: ' + str(group))
            seen_pairs.add(group)
            if group not in ref_groups:
                raise ValueError('Group not represented in reference oligo CSV: ' + str(group))
            for role in ref_groups[group]:
                dna, _ = normalize_sequence(field(row, role + '_full', role + '_dna'))
                if dna != ref_groups[group][role][3]:
                    raise ValueError('Wide and long reference sequences disagree: ' + str(group))
        if seen_pairs != set(ref_groups):
            raise ValueError('Wide and long reference pair sets disagree')
    generated = defaultdict(lambda: defaultdict(list))
    generated_groups = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(set))))
    stage_counts = defaultdict(int)
    for index, row in enumerate(generated_rows, 2):
        stage = field(row, 'stage', default='candidate').casefold()
        if stage not in ('candidate', 'final'):
            raise ValueError('Unknown generated stage: ' + stage)
        try:
            key, _, mods = canonical(row, version)
        except ValueError as error:
            # Invalid generated records are implementation failures, never silently excluded.
            raise ValueError('Generated CSV row %d: %s' % (index, error)) from error
        pair = field(row, 'probe_set_id', 'pair_id')
        probe = field(row, 'probe_id', 'oligo_id', default='csv_row_%d' % index)
        stage_counts[stage] += 1
        if key not in references:
            continue
        generated[stage][key].append((probe, pair, mods, field(row, 'modification_5p')))
        if pair:
            generated_groups[stage][(key[0], key[1])][pair][key[2]].add(key[3])
    # Empty generation is still a measured candidate set; final only if explicitly supplied.
    stages = ['candidate'] + (['final'] if 'final' in stage_counts else [])
    matches, pairs, per_gene, summary = [], [], [], []
    for stage in stages:
        matched = set(references).intersection(generated[stage])
        for key, ref in references.items():
            hits = generated[stage].get(key, [])
            matches.append(dict(stage=stage, method_version=key[0], gene=ref['gene'],
                oligo_type=key[2], reference_oligo_id=field(ref['row'], 'oligo_id', 'probe_id'),
                reference_pair_id=field(ref['row'], 'probe_set_id', 'pair_id'), full_sequence=key[3],
                reference_modifications=json.dumps(ref['mods']),
                reference_modification_5p=field(ref['row'], 'modification_5p'), matched=bool(hits),
                generated_probe_ids=';'.join(sorted({h[0] for h in hits})),
                generated_pair_ids=';'.join(sorted({h[1] for h in hits if h[1]})),
                generated_modification_5p=';'.join(sorted({h[3] for h in hits if h[3]})),
                generated_modifications=json.dumps(sorted({json.dumps(h[2]) for h in hits}))))
        for group, members in ref_groups.items():
            recovered_roles = {r: key in matched for r, key in members.items()}
            group_hits = []
            for generated_pair, generated_members in generated_groups[stage].get(group[:2], {}).items():
                if all(key[3] in generated_members.get(role, set()) for role, key in members.items()):
                    group_hits.append(generated_pair)
            pairs.append(dict(stage=stage, method_version=group[0], gene=references[next(iter(members.values()))]['gene'],
                reference_pair_id=group[2], required_roles=';'.join(sorted(members)),
                cpsr_applicable=len(members) > 1,
                member_matches=json.dumps(recovered_roles, sort_keys=True),
                all_members_individually_matched=all(recovered_roles.values()),
                matched=bool(group_hits), generated_pair_ids=';'.join(sorted(group_hits))))
        roles = ['all'] + sorted({key[2] for key in references})
        for role in roles:
            keys = [key for key in references if role == 'all' or key[2] == role]
            genes = sorted({key[1] for key in keys})
            counts = []
            for gene in genes:
                gene_keys = [key for key in keys if key[1] == gene]
                n = len(gene_keys)
                m = sum(key in matched for key in gene_keys)
                gene_pairs = [p for p in pairs if p['stage'] == stage and p['gene'].casefold() == gene and p['cpsr_applicable']]
                counts.append((n, m))
                per_gene.append(dict(stage=stage, oligo_type=role, gene=references[gene_keys[0]]['gene'],
                    reference_sequences=n, matched_sequences=m, SRR_pct=100*m/n,
                    reference_pairs=len(gene_pairs) if role == 'all' and gene_pairs else '',
                    matched_pairs=sum(p['matched'] for p in gene_pairs) if role == 'all' and gene_pairs else '',
                    CPSR_pct=100*sum(p['matched'] for p in gene_pairs)/len(gene_pairs) if role == 'all' and gene_pairs else ''))
            stage_pairs = [p for p in pairs if p['stage'] == stage and p['cpsr_applicable']]
            n, m = sum(c[0] for c in counts), sum(c[1] for c in counts)
            summary.append(dict(stage=stage, oligo_type=role, reference_sequences=n, matched_sequences=m,
                SRR_pct=100*m/n, reference_targets=len(genes), recovered_targets=sum(c[1] > 0 for c in counts),
                TRR_pct=100*sum(c[1] > 0 for c in counts)/len(genes),
                MTR_pct=100*sum(c[1]/c[0] for c in counts)/len(genes),
                reference_pairs=len(stage_pairs) if role == 'all' and stage_pairs else '',
                matched_pairs=sum(p['matched'] for p in stage_pairs) if role == 'all' and stage_pairs else '',
                CPSR_pct=100*sum(p['matched'] for p in stage_pairs)/len(stage_pairs) if role == 'all' and stage_pairs else '',
                generated_records=stage_counts[stage], generated_unique_reference_hits=len(generated[stage]),
                excluded_reference_records=len(excluded), duplicate_reference_records=duplicates))
    return dict(summary=summary, per_gene=per_gene, sequence_matches=matches, pair_matches=pairs,
                exclusions=excluded)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-oligos', required=True)
    parser.add_argument('--reference-pairs')
    parser.add_argument('--generated', required=True)
    parser.add_argument('--output', '--output-dir', dest='output_dir', required=True)
    parser.add_argument('--method-version', default='unspecified')
    args = parser.parse_args()
    result = evaluate(read_csv(args.reference_oligos), iter_csv(args.generated),
        read_csv(args.reference_pairs) if args.reference_pairs else None, args.method_version)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in result.items():
        if name == 'exclusions':
            (out / 'reference_exclusions.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
        else:
            columns = list(rows[0]) if rows else ['stage', 'reference_pair_id', 'matched']
            write_csv(out / (name + '.csv'), rows, columns)
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    main()
