"""Prepare a deterministic, annotation-stratified RAEFISH pilot from Table S1."""
import argparse
import csv
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
METHOD = 'raefish_2025_64gene'
SOURCES = {'supplement': 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867425010372-mmc1.xlsx',
           'publication': 'https://doi.org/10.1016/j.cell.2025.09.006',
           'annotation_checksums': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/MD5SUMS',
           'annotation': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz'}


def read_workbook(path):
    with zipfile.ZipFile(path) as z:
        strings = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS)]
        rels = {x.attrib['Id']: x.attrib['Target'] for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets = {}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets', NS):
            title = sheet.attrib['name']
            if title not in {'Human RAEFISH probe templates', '94-bit RAEFISH readout oligos'}:
                continue
            name = rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            root = ET.fromstring(z.read('xl/' + name))
            rows = []
            for row in root.findall('m:sheetData/m:row', NS):
                values = {}
                for cell in row:
                    value = cell.find('m:v', NS)
                    text = value.text if value is not None else ''
                    if cell.attrib.get('t') == 's':
                        text = strings[int(text)]
                    values[re.sub(r'\d', '', cell.attrib['r'])] = text
                rows.append(values)
            sheets[title] = rows
        return sheets


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(path, fields, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def prepare(workbook, gtf, out, audit_path):
    sheets = read_workbook(workbook)
    raw = sheets['Human RAEFISH probe templates'][2:]
    readouts = {int(r['A'].split()[1]): r['B'][:20]
                for r in sheets['94-bit RAEFISH readout oligos'] if r.get('A', '').startswith('Adapter ')}
    assert set(readouts) == set(range(1, 95))
    rows = []
    for row in raw:
        gid = 'ENSG' + str(int(row['A'])).zfill(11)
        bits = [i + 1 for i, c in enumerate(row['C']) if c == '1']
        assert len(row['C']) == 94 and set(row['C']) <= {'0', '1'} and len(bits) == 4
        padlock, priming, encoding = row['D'], row['E'], row['F']
        assert all(set(s) <= set('ACGT') for s in (padlock, priming, encoding))
        assert tuple(map(len, (padlock, priming, encoding))) == (80, 107, 150)
        assert padlock[:20] == 'ACTCTCACTTGGGTTCCAAT'
        assert padlock[50:] == 'TTCCTCTAACTCTCTAATCCACTAATACCC'
        assert priming[:20] == 'GTGGTAAAGCTCCGCGGCTT'
        assert priming[50:60] == 'ATTGGAACCC'
        assert priming[60:] == 'TATAGTGAGTCGTATTATCCGTACGGCGATTGGTCAATGCGGAACGA'
        assert encoding[:20] == 'GCGTCGTTATGGTGCAACGT'
        assert encoding[130:] == 'CGACACCGAACGTGCGACAA'
        assert encoding[60:90] == rc(padlock[20:50])
        assert all(encoding[start:start+20] == readouts[bit] for start, bit in zip((20,40,90,110), bits))
        rows.append(dict(target_id=gid, gene_symbol=row['B'], code=row['C'], bits=bits,
                         sequences=dict(padlock=padlock, priming=priming, encoding=encoding)))
    assert len(rows) == len({r['target_id'] for r in rows}) == 23352
    excluded = [r for r in rows if r['gene_symbol'].endswith('NegativeControl')]
    effective = [r for r in rows if not r['gene_symbol'].endswith('NegativeControl')]
    assert len(excluded) == 40 and len(effective) == 23312
    genes = {}
    gtf_text_hash = hashlib.sha256()
    opener = gzip.open if str(gtf).endswith('.gz') else open
    with opener(gtf, 'rt') as f:
        for line in f:
            gtf_text_hash.update(line.encode())
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            genes[attrs['gene_id'].split('.')[0]] = attrs
    missing = [r['target_id'] for r in effective if r['target_id'] not in genes]
    # Freeze before any sequence matching: first 32 IDs in each declared biotype.
    selected = []
    for biotype in ('protein_coding', 'lncRNA'):
        group = sorted((r for r in effective if genes.get(r['target_id'], {}).get('gene_type') == biotype), key=lambda r:r['target_id'])
        assert len(group) >= 32
        selected.extend(group[:32])
    selected.sort(key=lambda r:r['target_id'])
    out.mkdir(parents=True, exist_ok=True)
    targets = [dict(target_id=r['target_id'],gene_symbol=r['gene_symbol'],gene_type=genes[r['target_id']]['gene_type']) for r in selected]
    write_csv(out/'targets.csv', ['target_id','gene_symbol','gene_type'], targets)
    codes = []
    for r in selected:
        d = dict(target_id=r['target_id'],code=r['code'])
        d.update({f'bit{i}':bit for i,bit in enumerate(r['bits'],1)})
        d.update({f'R{i}':readouts[bit] for i,bit in enumerate(r['bits'],1)})
        codes.append(d)
    write_csv(out/'codebook.csv', ['target_id','code','bit1','bit2','bit3','bit4','R1','R2','R3','R4'], codes)
    write_csv(out/'readouts.csv', ['bit','sequence'], [dict(bit=b,sequence=s) for b,s in sorted(readouts.items())])
    refs = []
    for r in selected:
        for role, sequence in r['sequences'].items():
            refs.append(dict(method_version=METHOD,target_id=r['target_id'],probe_set_id=r['target_id']+'|published',oligo_type=role,full_sequence=sequence))
    write_csv(out/'reference_oligos.csv', ['method_version','target_id','probe_set_id','oligo_type','full_sequence'], refs)
    gtf_md5 = digest(gtf, 'md5')
    if str(gtf).endswith('.gz'):
        assert gtf_md5 == '68070e3c87f6a0e70af6e4ff4d55f61c', 'GENCODE v46 official MD5 mismatch'
    audit = dict(gtf_uncompressed_sha256=gtf_text_hash.hexdigest(),
                 gtf_md5=gtf_md5, gtf_official_compressed_md5='68070e3c87f6a0e70af6e4ff4d55f61c',
                 output_sha256={name:digest(out/name) for name in ('targets.csv','codebook.csv','readouts.csv','reference_oligos.csv')},
                 sources=SOURCES, workbook_sha256=digest(workbook),gtf_sha256=digest(gtf),
                 human_raw_records=len(rows),human_raw_unique_ids=len(rows),human_excluded_records=len(excluded),
                 exclusion_rule='gene_symbol ends with NegativeControl',
                 excluded_gene_symbols=[r['gene_symbol'] for r in excluded],human_effective_reference_records=len(effective),
                 human_effective_reference_oligos=len(effective)*3,
                 human_annotation_missing_ids=missing,
                 human_effective_biotypes=dict(Counter(genes.get(r['target_id'],{}).get('gene_type','missing') for r in effective)),
                 pilot_rule='Sort Ensembl IDs independently within GENCODE v46 protein_coding and lncRNA; take first 32 per stratum; combine and sort. No sequence matching is used.',
                 pilot_targets=targets,pilot_reference_oligos=len(refs),pilot_reference_probe_sets=len(selected),
                 integrity_checks='All Human rows have valid DNA, 80/107/150 lengths, 94-bit weight-4 codes, shared fixed regions, rc(padlock-binding) encoding middle, and adapter-prefix readouts in ascending bit order.',
                 annotation='GENCODE v46 / GRCh38',
                 template_encoding_rule='E5 + R1 + R2 + rc(padlock.RNA_binding) + R3 + R4 + E3; functional encoding is rc(template).',
                 target_id_rule='Numeric Table S1 Ensembl_GeneID padded to 11 digits after ENSG; annotation version suffix stripped.')
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k not in ('pilot_targets','excluded_gene_symbols')},indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workbook',type=Path,required=True)
    p.add_argument('--gtf',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    a=p.parse_args()
    prepare(a.workbook,a.gtf,a.output_dir,a.audit)
