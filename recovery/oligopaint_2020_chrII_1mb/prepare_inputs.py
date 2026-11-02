"""Prepare public metadata and evaluation-only references. No sequence mapping."""
import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import yaml

HASHES = {
    's009.fa': '2cc4b2327b2794dc14856df1882883497d0cc1aad4235a64fd023a4304300651',
    's010.xlsx': '21f0df561c5ebff501712aaf10307e4e6829ce90c4257427a1b6497e0a85c04d',
    'ce11.fa.gz': 'b5f80e37d22eb8f1cce1efff09a7e2e2ac2b6c73650f674c41c0b2e281c97241',
}
VERSION = 'oligopaint_2020_chrII_1mb'
# S2 explicitly labels these replacement barcodes Read053B188 / Read054B189.
# Only their fixed 20-nt barcode segments are used; binding DNA is never configuration.
BARCODE_OVERRIDES = {53: 'ACGACATGTTGCGACTGTCG', 54: 'TAGATGCGCACACCGCGATC'}


def rc(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def prepare(raw, output, config):
    output.mkdir(parents=True, exist_ok=True)
    for name, expected in HASHES.items():
        assert hashlib.sha256((raw / name).read_bytes()).hexdigest() == expected, name
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(raw / 's010.xlsx') as archive:
        strings = [''.join(x.itertext()) for x in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', ns)]
        rows = ET.fromstring(archive.read('xl/worksheets/sheet1.xml')).findall('.//s:row', ns)
        readouts = {}
        for row in rows[1:]:
            values = [strings[int(c.find('s:v', ns).text)] if c.attrib.get('t') == 's' else ''.join(c.itertext()) for c in row]
            number = re.search(r'readout(\d+)', values[0])
            assert number
            readouts[int(number[1])] = values[1][:20].upper()
    # Published bin boundaries: consecutive inclusive intervals from S2 headers.
    # Their formula is fixed metadata, never coordinates inferred from binding DNA.
    targets, encoding = [], {}
    for number in range(1, 101):
        start = 11500001 + (number - 1) * 10001
        end = min(start + 10000, 12500000)
        target = f'chrII:{start}-{end}'
        targets.append(target)
        encoding[target] = {'bin_id': f'bin{number:03d}', 'index5': 'ACGTCCGCCGCATCTACGAG',
                            'fiducial': 'CATCAACGCCACGATCAGCT',
                            'barcode': BARCODE_OVERRIDES.get(number, rc(readouts[number])),
                            'index3': 'CTGAACGGCGCACGTGTCTT'}
    cfg = {'name': VERSION, 'genome': 'ce11', 'targets': targets, 'encoding': encoding,
           'extracts': {'target_region': {'source': 'genome', 'length': 40, 'step': 1, 'masked': 'keep'}},
           'probes': {}}
    for role, binding in [('forward', 'target_region'), ('reverse', 'rc(target_region)')]:
        cfg['probes'][role] = {'template': '{index5}{fiducial}{barcode}{binding}{index3}',
            'parts': {name: {'expr': binding if name == 'binding' else f"encoding[target]['{name}']"}
                      for name in ['index5', 'fiducial', 'barcode', 'binding', 'index3']}}
    # Conditions are benchmark reporting conditions, not reconstructed paper thresholds.
    cfg['hybridization'] = {'na_mM': 50, 'mg_mM': 0, 'dntp_mM': 0.6, 'probe_nM': 50, 'formamide_pct': 0}
    cfg['attributes'] = {'gc': {'type': 'gc_content', 'target': 'binding'},
                         'tm': {'type': 'annealing_temperature', 'target': 'binding', 'tm_method': 'primer3'}}
    cfg['post_process'] = {'filters': {'unique_exact': {'condition': 'exact_sites == 1'},
                                      'gc_range': {'condition': 'gc >= 0.30 & gc <= 0.70'}},
                           'remove_overlap': {'location_interval': 0}}
    config.parent.mkdir(parents=True, exist_ok=True)
    text = yaml.safe_dump(cfg, sort_keys=False)
    config.write_text(text)
    refs = []
    for block in (raw / 's009.fa').read_text().split('>')[1:]:
        lines = block.splitlines()
        name, seq = lines[0], ''.join(lines[1:]).upper()
        assert re.fullmatch('[ACGT]{120}', seq)
        match = re.search(r'Read(\d+)[^_]*__boundary_10k_(\d+)_chr2:(\d+)-(\d+)', name)
        assert match and int(match[1]) == int(match[2])
        target = f'chrII:{match[3]}-{match[4]}'
        fixed = encoding[target]
        assert seq[:60] == fixed['index5'] + fixed['fiducial'] + fixed['barcode']
        assert seq[100:] == fixed['index3']
        refs.append({'method_version': VERSION, 'gene': fixed['bin_id'], 'oligo_type': 'primary',
                     'oligo_id': name, 'full_sequence': seq})
    assert len(refs) == 9984 and len({r['gene'] for r in refs}) == 100
    with (output / 'reference_oligos.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(refs[0]))
        writer.writeheader(); writer.writerows(refs)
    (output / 'sources.json').write_text(json.dumps({'sha256': HASHES, 'reference_count': len(refs),
        'bin_boundaries': 'Published inclusive S2 header boundaries; 99 x 10001 bp + 9901 bp',
        'barcode_exceptions': 'S2 Read053B188 / Read054B189 replace S1 readout053 / readout054; fixed barcode DNA only',
        'generator_inputs': 'FASTA, fixed public bin boundaries, S1 readout-derived barcodes and fixed index/fiducial DNA only'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.raw, args.output, args.config)
