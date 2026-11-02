"""Parse the published split-FISH tables; audit pairing separately from generation."""
import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import openpyxl
import yaml
from pyfaidx import Fasta
from uprobe.core.gen.fun import read_gtf, GTF_RUN_FIELDS, extract_spliced


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def write(path, rows, fields):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--templates', type=Path, required=True)
    p.add_argument('--codebook', type=Path, required=True)
    p.add_argument('--bridges', type=Path, required=True)
    p.add_argument('--author-bridges', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    table = lambda path: list(openpyxl.load_workbook(path, read_only=True, data_only=True).active.values)
    raw = table(a.templates)[2:]
    codes = {g: [i+1 for i, b in enumerate(c.split()) if b == '1'] for g, c in table(a.codebook)[2:]}
    bridges = {}
    for bit, seq in table(a.bridges)[2:]:
        bridges[int(bit[1:])] = seq.split()[1].upper()
    author = [r[0].upper() for r in table(a.author_bridges)]
    assert len(author) == len(bridges) == 26
    assert all(bridges[i] == author[i-1] for i in bridges)
    groups = defaultdict(list)
    for excel_row, (gene, tx, seq) in enumerate(raw, 3):
        parts = seq.upper().split()
        dna = ''.join(parts)
        assert len(dna) == 88 and re.fullmatch('[ACGT]+', dna)
        assert parts[0] == 'AACGAACGGAGGGTCATTGG' and parts[-1] == 'CCCTAATGCGAGTAGAGCCT'
        if parts[1] == 'GGTAC':
            role, half, binding = 'p1', parts[2], parts[4]
            assert parts[3] == 'TAT' and parts[5] == 'GAATTC'
            bits = [b for b in codes[gene] if half == rc(bridges[b][:9])]
            assert len(bits) == 1
            bit = bits[0]
        else:
            role, half, binding, bit = 'p2', parts[4], parts[2], None
            assert parts[1] == 'GGTACC' and parts[3] == 'TAT' and parts[5] == 'AATTC'
            assert all(half == rc(bridges[b][9:]) for b in codes[gene])
        assert len(binding) == 25
        groups[gene].append(dict(gene=gene, transcript_id=tx, oligo_id=f'S1_row{excel_row}',
                                 oligo_type=role, full_sequence=dna, binding=binding, bit=bit))
    genes = set(groups)
    subset = a.output/'audit_annotation.gtf'
    with a.gtf.open() as src, subset.open('w') as dst:
        for line in src:
            if line.startswith('#'):
                continue
            m = re.search(r'gene_name "([^"]+)"', line)
            if m and m[1] in genes:
                dst.write(line)
    gtf = read_gtf(str(subset), extract_fields=GTF_RUN_FIELDS)
    oligos, pairs, audits, uncertain = [], [], [], []
    with Fasta(str(a.fasta)) as fa:
        for gene, rows in groups.items():
            assert len(rows) == 144 and len(codes[gene]) == 2
            txs = {r['transcript_id'] for r in rows}
            assert len(txs) == 1
            tx = next(iter(txs))
            ann = gtf[gtf.gene_name.eq(gene) & gtf.transcript_id.eq(tx)]
            sequence = ''
            if not ann.empty:
                _, segments = extract_spliced(ann, fa, gene, 'mRNA')
                sequence = ''.join(s for _, s in segments).upper()
            for block, bit in enumerate(codes[gene]):
                p1s, p2s = rows[block*72:block*72+36], rows[block*72+36:block*72+72]
                assert all(r['oligo_type'] == 'p1' and r['bit'] == bit for r in p1s)
                assert all(r['oligo_type'] == 'p2' for r in p2s)
                for i, (p1, p2) in enumerate(zip(p1s, p2s), 1):
                    pair = f'{gene}:B{bit}:{i}'
                    # Audit only: published arms never enter the generator or its config.
                    pattern = re.compile(re.escape(rc(p2['binding']))+'[ACGT]{2}'+re.escape(rc(p1['binding'])))
                    sites = [m.start()+1 for m in pattern.finditer(sequence)]
                    confirmed = bool(sites)
                    audit = dict(gene=gene, reference_pair_id=pair, transcript_id=tx, bit=bit,
                                 confirmed=confirmed, transcript_starts=sites,
                                 reason='' if confirmed else 'No exact 25+2+25 site in published M4 transcript')
                    audits.append(audit)
                    if confirmed:
                        pairs.append(dict(gene=gene, pair_id=pair, p1_full=p1['full_sequence'], p2_full=p2['full_sequence']))
                    else:
                        uncertain.append(audit)
                    for row in (p1, p2):
                        oligos.append({**row, 'pair_id': pair if confirmed else ''})
    write(a.output/'reference_oligos.csv', [{k:r[k] for k in ['gene','oligo_id','oligo_type','full_sequence','pair_id']} for r in oligos],
          ['gene','oligo_id','oligo_type','full_sequence','pair_id'])
    write(a.output/'reference_pairs.csv', pairs, ['gene','pair_id','p1_full','p2_full'])
    write(a.output/'codebook.csv', [dict(gene=g, bit=b, bridge=bridges[b]) for g in groups for b in codes[g]], ['gene','bit','bridge'])
    cfg = yaml.safe_load(Path(__file__).with_name('design.yaml').read_text())
    cfg['targets'] = list(groups)
    cfg['encoding'] = {g: {'bridges':[bridges[b] for b in codes[g]], 'bits':codes[g]} for g in groups}
    (a.output/'design.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    hashes = {key: hashlib.sha256(getattr(a,key).read_bytes()).hexdigest() for key in ['templates','codebook','bridges','author_bridges']}
    result = dict(genes=len(groups), oligos=len(oligos), proposed_pairs=len(audits), confirmed_pairs=len(pairs),
                  uncertain_pairs=uncertain, codebook_blank_controls=[g for g in codes if g not in genes],
                  bridge_core_matches_author=True, input_sha256=hashes,
                  pairing_rule='Per gene and ascending on-bit: 36 P1 followed by 36 P2; same ordinal, verified 25+2+25 on published transcript. Unverified pairs omitted from CPSR only.')
    (a.output/'reference_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    (a.output/'pair_audit.json').write_text(json.dumps(audits)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'uncertain_pairs'}), flush=True)
    print('uncertain_pairs',len(uncertain),flush=True)


if __name__ == '__main__':
    main()
