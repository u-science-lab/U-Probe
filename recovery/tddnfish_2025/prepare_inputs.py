"""Extract the complete nine-gene HeLa panel; references never feed generation."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import openpyxl


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workbook', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    w = openpyxl.load_workbook(a.workbook, data_only=True)
    assert hashlib.sha256(a.workbook.read_bytes()).hexdigest() == '513a441145496471228bd2f4e2a67f2b59807b1502196396b43007f775ce6855', 'Unexpected source workbook'
    rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    td = dict(list(w['TDN_sequence'].values)[2:32])
    # Independent fixed readouts: first 30 bases of the published T0a strands.
    readouts = {rc(re.sub(r'\s', '', td['T0a-' + name])[:30]): name
                for name in ('hyb1-Cy5', 'hyb1-Cy3', 'hyb1-FAM',
                             'hyb2-Cy5', 'hyb2-Cy3', 'hyb2-FAM')}
    refs, codes = [], set()
    for i, row in enumerate(list(w['Figure2'].values)[81:135], 82):
        name, seq = row
        gene = {'TGFB 1': 'TGFB1', 'NF-KB': 'NFKB1'}.get(name, name)
        seq = re.sub(r'\s', '', seq).upper()
        assert len(seq) == 65 and re.fullmatch('[ACGT]+', seq)
        assert seq[:30] in readouts
        refs.append(dict(target_id=gene, oligo_type='primary',
                         oligo_id=f'Figure2:B{i}', full_sequence=seq,
                         published_name=name, species='Homo sapiens'))
        codes.add((gene, seq[:30], readouts[seq[:30]]))
    assert len(refs) == 54 and len({r['target_id'] for r in refs}) == 9
    assert len({(r['target_id'], r['full_sequence']) for r in refs}) == 54
    a.output.mkdir(parents=True, exist_ok=True)
    def write(name, rows):
        with (a.output / name).open('w', newline='') as f:
            writer = csv.DictWriter(f, list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    write('reference_oligos.csv', refs)
    write('codebook.csv', [dict(target_id=g, readout=s, channel=c)
                           for g, s, c in sorted(codes)])
    audit = dict(source='https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41467-025-64294-1/MediaObjects/41467_2025_64294_MOESM3_ESM.xlsx',
                 workbook_sha256=hashlib.sha256(a.workbook.read_bytes()).hexdigest(),
                 scope='Figure2 rows 82-135: complete two-round nine-gene HeLa panel',
                 reference_records=54, unique_reference_oligos=54, targets=9,
                 excluded_records=0, duplicate_records=0,
                 architecture_5prime_to_3prime='readout(30)-binding(35)',
                 normalizations={'TGFB 1': 'TGFB1', 'NF-KB': 'NFKB1'},
                 nfkb_evidence='All six binding sequences match NFKB1 ENST00000226574.9; mapping is evaluation-only, never candidate window selection.',
                 outside_scope='Other Figure1/2 panels, influenza Figure3, mouse Figure4/5; no enlargement of the complete HeLa panel',
                 candidate_inputs='Only gene identities and fixed external readouts; no published binding DNA, full primary DNA or inferred coordinates')
    (a.output / 'input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')


if __name__ == '__main__':
    main()
