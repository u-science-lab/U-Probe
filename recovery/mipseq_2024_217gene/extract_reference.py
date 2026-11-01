#!/usr/bin/env python3
"""Extract only the published 217-gene Fig 3 panel; Python standard library."""
import argparse
import collections
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

SOURCE = 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-024-01205-7/MediaObjects/41551_2024_1205_MOESM5_ESM.xlsx'
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
PLP = re.compile(r'^([ACGT]{13})CC([ACGT]{3})TGCGTCTATTTGT([ACGT]{3})TAGTGGAGCCT([ACGT]{13})$')

def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('workbook', type=Path)
    ap.add_argument('--output-dir', type=Path, default=Path(__file__).parent)
    args = ap.parse_args()
    z = zipfile.ZipFile(args.workbook)
    strings = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS)]
    workbook = ET.fromstring(z.read('xl/workbook.xml'))
    sheet = next(s for s in workbook.find('m:sheets', NS) if s.attrib['name'] == 'Fig 3')
    rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
    rel = next(r for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels')) if r.attrib['Id'] == rid)
    rows = []
    for row in ET.fromstring(z.read('xl/' + rel.attrib['Target'])).findall('.//m:row', NS):
        cells = {}
        for cell in row.findall('m:c', NS):
            val = cell.find('m:v', NS)
            val = val.text if val is not None else ''
            cells[re.sub(r'\d', '', cell.attrib['r'])] = strings[int(val)] if cell.attrib.get('t') == 's' and val else val
        if cells.get('B', '').strip() in ('', 'Name'):
            continue
        rows.append((int(row.attrib['r']), cells['B'], cells['C']))
    assert len(rows) % 2 == 0
    oligos, pairs, codes, anomalies = [], [], [], []
    for pad, rca in zip(rows[::2], rows[1::2]):
        pr, pn, ps = pad
        rr, rn, rs = rca
        pm = PLP.fullmatch(ps)
        assert pm and re.fullmatch('[ACGT]{14,15}AGGCTCCACTA', rs), (pr, rr)
        prefix = pn.strip().rsplit('-', 2)[0]
        assert prefix == rn.strip().rsplit('-', 2)[0], (pr, rr)
        gene = re.sub(r'-[12]$', '', prefix)
        pair = prefix if prefix != gene else gene + '-unindexed'
        b1, b2 = pm.group(2, 3)
        code = dict(gene=gene, B1=b1, B2=b2)
        if code not in codes:
            codes.append(code)
        pairrow = dict(gene=gene, pair_id=pair, padlock_full=ps, rca_full=rs, padlock_source_row=pr, rca_source_row=rr)
        pairs.append(pairrow)
        for rownum, name, sequence, role in [(pr, pn, ps, 'padlock'), (rr, rn, rs, 'rca')]:
            note = ''
            if role == 'rca' and '-PLP-' in name:
                note = 'RCA role normalized: 25-nt RCA suffix, adjacent matching padlock prefix, no padlock backbone'
                anomalies.append(dict(source_row=rownum, original_name=name, normalization=note))
            if name != name.strip():
                anomalies.append(dict(source_row=rownum, original_name=name, normalization='Outer whitespace ignored for name parsing; original retained'))
            if gene == prefix:
                note = (note + '; ' if note else '') + 'Missing pair index; unique adjacent pair retained'
            oligos.append(dict(gene=gene, pair_id=pair, oligo_type=role, oligo_id=f'Fig3_row_{rownum}', full_sequence=sequence, dna_sequence=sequence, source_sheet='Fig 3', source_row=rownum, original_name=name, normalization=note))
    genes = list(dict.fromkeys(o['gene'] for o in oligos))
    assert len(genes) == 217 and len(pairs) == 434 and len(oligos) == 868
    assert len(set(p['pair_id'] for p in pairs)) == len(pairs)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    write(out / 'targets.csv', [dict(gene=g) for g in genes])
    write(out / 'reference_oligos.csv', oligos)
    write(out / 'reference_pairs.csv', pairs)
    write(out / 'codebook.csv', codes)
    audit = dict(source_url=SOURCE, workbook_sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest(), source_sheet='Fig 3', targets=len(genes), raw_oligos=len(rows), valid_oligos=len(oligos), valid_pairs=len(pairs), sequence_lengths=dict(collections.Counter(len(o['full_sequence']) for o in oligos)), rca_spacer_counts=dict(collections.Counter(o['full_sequence'][13:-11] for o in oligos if o['oligo_type']=='rca')), role_normalizations=[a for a in anomalies if 'RCA role' in a['normalization']], name_whitespace=[a for a in anomalies if 'whitespace' in a['normalization']], unindexed_pair=dict(gene='Syt1', source_rows=[91,92]), unresolved_reference_records=0, pair_evidence='Adjacent padlock/RCA rows share the name prefix; 58-nt padlock and 25/26-nt RCA structures validated; no sequence rewritten.')
    (out / 'reference_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({k:audit[k] for k in ['targets','raw_oligos','valid_oligos','valid_pairs','sequence_lengths','role_normalizations','unindexed_pair','unresolved_reference_records']}))

if __name__ == '__main__':
    main()
