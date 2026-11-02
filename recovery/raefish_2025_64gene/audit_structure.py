"""Stream native RAEFISH groups and check structure independently of recovery."""
import argparse
import csv
import gzip
import json
from itertools import groupby
from pathlib import Path


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle)


def audit(generated, codebook, readouts, output):
    codes = {r['target_id']: r for r in rows(codebook)}
    adapters = {int(r['bit']): r['sequence'] for r in rows(readouts)}
    checked = 0
    gaps = {0: 0, 1: 0}
    directions = {'padlock_upstream': 0, 'padlock_downstream': 0}
    for (gene, sid), members in groupby(rows(generated), lambda r: (r['target_id'], r['probe_set_id'])):
        members = list(members)
        by_role = {r['oligo_type']: r for r in members}
        assert len(members) == 3 and set(by_role) == {'padlock', 'priming', 'encoding'}, sid
        assert all(r['method_version'] == 'raefish_2025_64gene' for r in members), sid
        origin = {(r['transcript_id'], r['target_start'], r['target_end'], r['gap'], r['target_region'], r['orientation']) for r in members}
        assert len(origin) == 1, sid
        tx, start, end, gap, target, orientation = next(iter(origin))
        start, end, gap = int(start), int(end), int(gap)
        assert gap in gaps and len(target) == 60 + gap and end-start+1 == len(target), sid
        assert tx and set(target) <= set('ACGT'), sid
        p, l, e = (by_role[r]['full_sequence'] for r in ('padlock', 'priming', 'encoding'))
        assert tuple(map(len, (p, l, e))) == (80, 107, 150), sid
        assert p[:20] == 'ACTCTCACTTGGGTTCCAAT' and p[50:] == 'TTCCTCTAACTCTCTAATCCACTAATACCC', sid
        assert l[:20] == 'GTGGTAAAGCTCCGCGGCTT' and l[50:60] == 'ATTGGAACCC', sid
        assert l[60:] == 'TATAGTGAGTCGTATTATCCGTACGGCGATTGGTCAATGCGGAACGA', sid
        assert l[50:60] == rc(p[10:20]), sid
        assert orientation in ('padlock_upstream', 'padlock_downstream'), sid
        assert sid == f'{gene}|{tx}|{start}:{end}|{orientation}', sid
        ps, ls = (start, start+30+gap) if orientation == 'padlock_upstream' else (start+30+gap, start)
        assert all(int(r['padlock_binding_start']) == ps and int(r['priming_binding_start']) == ls for r in members), sid
        pb, lb = (target[:30], target[30+gap:]) if orientation == 'padlock_upstream' else (target[30+gap:], target[:30])
        assert p[20:50] == rc(pb) and l[20:50] == rc(lb), sid
        code = codes[gene]
        bits = [i+1 for i, bit in enumerate(code['code']) if bit == '1']
        assert len(code['code']) == 94 and bits == [int(code[f'bit{i}']) for i in range(1,5)], gene
        assert all(code[f'R{i}'] == adapters[bit] for i,bit in enumerate(bits,1)), gene
        assert e[:20] == 'GCGTCGTTATGGTGCAACGT' and e[130:] == 'CGACACCGAACGTGCGACAA', sid
        assert e[60:90] == rc(p[20:50]), sid
        assert all(e[a:a+20] == code[f'R{i}'] for i,a in enumerate((20,40,90,110),1)), sid
        assert by_role['encoding']['functional_encoding'] == rc(e), sid
        checked += 1
        gaps[gap] += 1
        directions[orientation] += 1
    assert checked > 0
    result = dict(native_probe_sets_checked=checked, gap0_sets=gaps[0], gap1_sets=gaps[1],
                  orientation_counts=directions,
                  native_group_id_rule='target_id|transcript_id|start:end|orientation; checked against every origin',
                  all_checks_passed=True,
                  checks=['three roles from one native record', 'shared target/transcript/window origin',
                          '80/107/150-nt synthesis templates', 'adjacent nonoverlapping 30-nt RNA sites',
                          'RCA priming complement', 'encoding dependence on padlock binding',
                          '94-bit weight-4 readout assignment', 'synthesis-template orientation'],
                  note='Structural acceptance is separate from exact published-sequence recovery.')
    Path(output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generated', required=True)
    parser.add_argument('--codebook', required=True)
    parser.add_argument('--readouts', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    audit(args.generated, args.codebook, args.readouts, args.output)
