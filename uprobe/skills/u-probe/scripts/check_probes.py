#!/usr/bin/env python3
"""Check a U-Probe probe table before delivery.

Reads the final table of `uprobe run` (<name>_<time>_probes.tsv) and checks that
- every probe sequence column holds only A, C, G and T (no N placeholders left);
- each binding part is the reverse complement of the target: rc(binding) lies in
  `target_region` (equal for one-part probes, a sub-part of split pairs and padlocks);
  parts given with --sense (padlock arms that bind cDNA) must lie in it unchanged;
- with --transcripts, every window lies on the transcript U-Probe designed on (ids
  from --summary, or --transcript TARGET=ID), so the probes are on the reference mRNA;
- with --gc, --tm, --min-gap or --min-probes, every probe or target meets those limits;
and summarises per target: probes vs requested, span along the transcript, spacing,
binding length and GC, Tm columns and off-target columns. Spacing uses target_start/target_end,
or sub_region (s-e along the region) for genomic regions (source: genome).

Standard library only. Exit status: 0 all checks passed, 1 a check failed, 2 unusable input.
Example:
  python check_probes.py results/run1/smfish_<time>_probes.tsv \\
      --summary results/run1/summary.json --transcripts genomes/<genome>/transcript.fa
"""
import argparse
import csv
import json
import re
import statistics
import sys

ACGT = re.compile(r'^[ACGT]+$')
LETTERS = re.compile(r'^[A-Za-z]+$')
TM_COLUMN = re.compile(r'(^|_)tm\d*$')
VERSION = re.compile(r'\.\d+$')
SUB_REGION = re.compile(r'^(\d+)-(\d+)$')  # window along a genomic region (source: genome), 1-based


def rc(seq):
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def gc(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0.0


def number(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return None if x != x else x  # NaN -> None


def stats(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return {'min': round(min(values), 3), 'median': round(statistics.median(values), 3), 'max': round(max(values), 3)}


def read_table(path):
    limit = sys.maxsize
    while True:  # off-target lists can be long
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10
    with open(path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
        return reader.fieldnames or [], list(reader)


def record_ids(header):
    """Ids a FASTA header stands for: its first word, the first '|' field (GENCODE) and, for
    U-Probe's transcript.fa records <gene_id>_<transcript_id>, the transcript id; unversioned."""
    word = (header.split() or [''])[0].split('|')[0]
    ids = {word} | ({word.split('_', 1)[1]} if '_' in word else set())
    return {VERSION.sub('', i) for i in ids}


def read_fasta(path, wanted=None, max_records=None):
    """{id: sequence} of the records whose ids are in `wanted`; without `wanted`, the first
    `max_records` records under their first word."""
    wanted = None if wanted is None else {VERSION.sub('', w) for w in wanted}
    found, name, chunks = {}, None, []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    found[name] = ''.join(chunks).upper()
                if wanted is None and max_records is not None and len(found) >= max_records:
                    return found
                hits = record_ids(line[1:].strip()) & wanted if wanted is not None else None
                name = (sorted(hits)[0] if hits else None) if wanted is not None else (line[1:].split() or ['?'])[0]
                chunks = []
            elif name:
                chunks.append(line.strip())
    if name:
        found[name] = ''.join(chunks).upper()
    return found


def window_position(row):
    """(start, end) of a window along its target, 1-based inclusive: target_start/target_end, or the
    sub_region "s-e" of a genomic-region window (source: genome); None when the row has neither."""
    start, end = number(row.get('target_start')), number(row.get('target_end'))
    if start is not None and end is not None:
        return start, end
    match = SUB_REGION.match((row.get('sub_region') or '').strip())
    return (int(match.group(1)), int(match.group(2))) if match else None


def gene_labels(text):
    """Gene names of an <attribute>_off_targets or _copies cell ('GENE:mismatches, GENE2')."""
    labels = []
    for item in (text or '').split(','):
        label, colon, count = item.strip().rpartition(':')
        item = label if colon and count.isdigit() else item.strip()
        if item:
            labels.append(item)
    return labels


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split('\n', 1)[1])
    p.add_argument('tsv', help='final probe table of uprobe run (*_probes.tsv)')
    p.add_argument('--summary', help='summary.json of the same run (transcript ids, requested probes)')
    p.add_argument('--transcripts', help="transcript FASTA, e.g. transcript.fa next to the genome FASTA")
    p.add_argument('--transcript', action='append', default=[], metavar='TARGET=ID',
                   help='transcript a target was designed on (overrides --summary); repeatable')
    p.add_argument('--binding', action='append', default=[], metavar='COLUMN',
                   help='part that must be rc(target) [default: every column named <probe>.binding]')
    p.add_argument('--sense', action='append', default=[], metavar='COLUMN',
                   help='part that carries the target sense sequence (cDNA padlock arms)')
    p.add_argument('--gc', nargs=2, type=float, metavar=('MIN', 'MAX'), help='GC fraction window of binding parts')
    p.add_argument('--tm', nargs=2, type=float, metavar=('MIN', 'MAX'), help='window for the Tm columns')
    p.add_argument('--tm-column', action='append', default=[], metavar='COLUMN',
                   help='Tm column(s) [default: numeric columns named *_tm or tm<N>]')
    p.add_argument('--min-gap', type=int, metavar='NT', help='minimum gap between neighbouring windows')
    p.add_argument('--min-probes', type=int, metavar='N', help='minimum probes per target')
    p.add_argument('--json', action='store_true', help='print the report as JSON')
    return p.parse_args(argv)


def check(args):
    columns, rows = read_table(args.tsv)
    if not rows:
        raise SystemExit(f'{args.tsv}: no probes in the table')
    for needed in ('target', 'target_region'):
        if needed not in columns:
            raise SystemExit(f'{args.tsv}: no {needed!r} column; is this a U-Probe *_probes.tsv?')
    parts = [(c, 'rc') for c in (args.binding or [c for c in columns if c.split('.')[-1] == 'binding'])]
    parts += [(c, 'sense') for c in args.sense]
    missing = [c for c, _ in parts if c not in columns]
    if missing or not parts:
        raise SystemExit(f"binding column(s) {missing or 'none found'}; name them with --binding/--sense "
                         f"(columns: {', '.join(columns)})")
    probes = {c.split('.')[0] for c, _ in parts}
    seq_columns = [c for c in columns if c.split('.')[0] in probes
                   and all(LETTERS.match(r[c]) for r in rows if r[c])]
    tm_columns = args.tm_column or [c for c in columns if TM_COLUMN.search(c)
                                    and all(number(r[c]) is not None for r in rows if r[c])]
    mapped = [c for c in columns if c.endswith('mapped_genes') and all(number(r[c]) is not None for r in rows if r[c])]

    summary = {}
    if args.summary:
        with open(args.summary, encoding='utf-8') as fh:
            summary = json.load(fh)
    chosen = {}
    for name, entry in (summary.get('targets') or {}).items():
        info = (entry or {}).get('transcripts') or {}
        if info.get('ids'):
            chosen[name] = {'ids': list(info['ids']), 'chosen_by': info.get('chosen_by')}
    for item in args.transcript:
        target, _, tid = item.partition('=')
        if not tid:
            raise SystemExit(f'--transcript {item!r}: expected TARGET=ID')
        chosen[target] = {'ids': [tid], 'chosen_by': 'given'}
    sequences = {}
    if args.transcripts:
        wanted = {i for c in chosen.values() for i in c['ids']}
        # Without ids only a one-record FASTA can be used (for a one-target table).
        sequences = read_fasta(args.transcripts, wanted) if wanted else read_fasta(args.transcripts, max_records=2)

    failures, warnings, report = [], [], {}

    def fail(row, message):
        failures.append(f"{row.get('probe_id') or row['target']}: {message}")

    for row in rows:
        window = row['target_region'].upper()
        for c in seq_columns:
            if row[c] and not ACGT.match(row[c]):
                at = next(i for i, letter in enumerate(row[c]) if letter not in 'ACGT')
                fail(row, f'{c} is not pure ACGT: {row[c][at]!r} at position {at + 1} of {len(row[c])}')
        for c, orientation in parts:
            seq = row[c].upper()
            expected = rc(seq) if orientation == 'rc' else seq
            if not seq:
                fail(row, f'{c} is empty')
            elif expected not in window:
                fail(row, f"{c} is not the {'reverse complement' if orientation == 'rc' else 'sense sequence'} "
                          f"of target_region" + (' (it is the sense sequence)' if seq in window else ''))

    targets = list(dict.fromkeys(r['target'] for r in rows))
    for target in targets:
        trows = [r for r in rows if r['target'] == target]
        entry = {'probes': len(trows)}
        requested = ((summary.get('targets') or {}).get(target) or {}).get('requested_probes')
        if requested:
            entry['requested'] = requested
            if len(trows) < requested:
                warnings.append(f'{target}: {len(trows)} probes, {requested} requested (report the shortfall)')
        if args.min_probes and len(trows) < args.min_probes:
            failures.append(f'{target}: {len(trows)} probes, fewer than --min-probes {args.min_probes}')
        # Where the windows lie on the transcript(s) the target was designed on.
        ids = (chosen.get(target) or {}).get('ids')
        tx = {key: sequences[key] for key in (VERSION.sub('', i) for i in ids or []) if key in sequences}
        if not ids and len(sequences) == 1 and len(targets) == 1:
            tx = dict(sequences)
        if args.transcripts and not tx:
            warnings.append(f'{target}: no transcript sequence to check against (give --summary or --transcript)')
        positions = []
        entry['chosen_by'] = (chosen.get(target) or {}).get('chosen_by')
        if entry['chosen_by'] == 'all':
            warnings.append(f'{target}: designed on all annotated transcripts (transcripts: all), not on one '
                            'reference transcript; intended?')
        if tx:
            entry['transcripts'] = {i: len(s) for i, s in tx.items()}
            for r in trows:
                hits = [(i, s.find(r['target_region'].upper())) for i, s in tx.items()]
                hits = [(i, pos) for i, pos in hits if pos >= 0]
                if not hits:
                    fail(r, f"window not on the designed transcript ({', '.join(tx)})")
                else:
                    positions.append((hits[0][0], hits[0][1] + 1, hits[0][1] + len(r['target_region'])))
            for i, length in entry['transcripts'].items():
                spans = [(s, e) for t, s, e in positions if t == i]
                if spans:
                    first, last = min(s for s, _ in spans), max(e for _, e in spans)
                    entry.setdefault('span', {})[i] = {'from': first, 'to': last, 'length': length,
                                                       'fraction': round((last - first + 1) / length, 3)}
        # Spacing along the target: gap = next start - previous end - 1 (1-based, inclusive).
        located = [window_position(r) for r in trows]
        if all(located):
            located = sorted(located)
        else:
            located = sorted((s, e) for _, s, e in positions)
        if args.min_gap is not None and len(trows) > 1 and len(located) < len(trows):
            failures.append(f'{target}: cannot check --min-gap: windows without a position (no '
                            'target_start/target_end or sub_region, and none found with --transcripts)')
        gaps = [int(b[0] - a[1] - 1) for a, b in zip(located, located[1:])]
        if gaps:
            entry['gap'] = stats(gaps)
            entry['overlaps'] = sum(g < 0 for g in gaps)
            if args.min_gap is not None and min(gaps) < args.min_gap:
                failures.append(f'{target}: smallest gap {min(gaps)} nt < --min-gap {args.min_gap}')
            elif entry['overlaps']:
                warnings.append(f"{target}: {entry['overlaps']} overlapping windows (intended?)")
        for c, _ in parts:
            seqs = [r[c].upper() for r in trows if r[c]]
            gcs = [gc(s) for s in seqs]
            entry.setdefault('parts', {})[c] = {'length': stats([len(s) for s in seqs]), 'gc': stats(gcs)}
            if args.gc:
                out = sum(not args.gc[0] <= g <= args.gc[1] for g in gcs)
                if out:
                    failures.append(f'{target}: {out} {c} outside GC {args.gc[0]}-{args.gc[1]}')
        for c in tm_columns:
            tms = [number(r[c]) for r in trows]
            entry.setdefault('tm', {})[c] = stats(tms)
            if args.tm:
                out = sum(t is None or not args.tm[0] <= t <= args.tm[1] for t in tms)
                if out:
                    failures.append(f'{target}: {out} {c} outside Tm {args.tm[0]}-{args.tm[1]}')
        for c in mapped:
            counts = stats([number(r[c]) for r in trows])
            hit = [r for r in trows if gene_labels(r.get(f'{c}_off_targets'))]
            entry.setdefault('off_targets', {})[c] = {'max_genes': counts and counts['max'],
                                                      'rows_with_off_targets': len(hit)}
            if hit:
                examples = '; '.join(f"{r.get('probe_id')}: {r[f'{c}_off_targets'][:60]}" for r in hit[:3])
                warnings.append(f'{target}: {len(hit)} probes list off-targets in {c}_off_targets ({examples})')
            copies = sorted({g for r in trows for g in gene_labels(r.get(f'{c}_copies'))})
            if copies:
                warnings.append(f"{target}: exempt copies in {c}_copies: {', '.join(copies[:10])} (name them in the report)")
        report[target] = entry
    return {'ok': not failures, 'table': args.tsv, 'binding_columns': [c for c, _ in parts],
            'failures': failures, 'warnings': warnings, 'targets': report}


def render(result):
    lines = [f"{result['table']}: {'OK' if result['ok'] else 'FAILED'} "
             f"({len(result['failures'])} failures, {len(result['warnings'])} warnings)"]
    for target, e in result['targets'].items():
        req = f" of {e['requested']} requested" if 'requested' in e else ''
        lines.append(f"{target}: {e['probes']} probes{req}" + (f" ({e['chosen_by']})" if e.get('chosen_by') else ''))
        for i, s in (e.get('span') or {}).items():
            lines.append(f"  on {i} ({s['length']} nt): nt {s['from']}-{s['to']} ({s['fraction']:.0%} of the transcript)")
        if 'gap' in e:
            lines.append(f"  gap between windows: {e['gap']['min']:g}-{e['gap']['max']:g} nt "
                         f"(median {e['gap']['median']:g}); overlaps: {e['overlaps']}")
        for c, p in (e.get('parts') or {}).items():
            lines.append(f"  {c}: length {p['length']['min']:g}-{p['length']['max']:g}, "
                         f"GC {p['gc']['min']:.2f}-{p['gc']['max']:.2f} (median {p['gc']['median']:.2f})")
        for c, t in (e.get('tm') or {}).items():
            if t:
                lines.append(f"  {c}: {t['min']:.1f}-{t['max']:.1f} (median {t['median']:.1f})")
        for c, o in (e.get('off_targets') or {}).items():
            lines.append(f"  {c}: max {o['max_genes']}, {o['rows_with_off_targets']} probes list off-target genes")
    lines += [f'FAIL {m}' for m in result['failures'][:50]]
    if len(result['failures']) > 50:
        lines.append(f"... {len(result['failures']) - 50} more failures (use --json)")
    lines += [f'WARN {m}' for m in result['warnings']]
    return '\n'.join(lines)


def main(argv=None):
    args = parse_args(argv)
    try:
        result = check(args)
    except SystemExit as e:
        if isinstance(e.code, str):
            print(f'check_probes: {e.code}', file=sys.stderr)
            return 2
        raise
    except (OSError, ValueError) as e:
        print(f'check_probes: {e}', file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2) if args.json else render(result))
    return 0 if result['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
