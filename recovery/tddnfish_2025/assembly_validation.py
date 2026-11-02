"""Validate declared scaffold interfaces and native modification exports; no recovery metric."""
import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import openpyxl
import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.probe import DAG
from uprobe.core.validate import validate_protocol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', required=True)
    parser.add_argument('--config', default=str(Path(__file__).with_name('design.yaml')))
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    p = yaml.safe_load(Path(args.config).read_text())
    assert not [f for f in validate_protocol(p) if f.severity == 'error']
    workbook = openpyxl.load_workbook(args.workbook, data_only=True)
    td = {}
    for row in workbook['TDN_sequence'].values:
        cells = [str(v).strip() for v in row if v is not None]
        if len(cells) >= 2:
            td[cells[0]] = cells[1]
    reference = {f'T{g}{s}': td['T0a-hyb1-Cy5' if (g, s) == (0, 'a') else f'T{g}{s}']
                 for g in range(3) for s in 'abcd'}
    reference['Detection'] = td['Cy5-ssDNA']
    clean = lambda s: s.replace('-SS-', '').replace('Cy5-', '').replace(' ', '')
    rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    api = UProbeAPI(p, {}, out, require_genome=False)
    windows = api.generate_target_seqs()
    assert len(windows) == 1
    native = api.construct_probes(windows)
    row = native.iloc[0]
    matches = {n: row[n] == clean(s) for n, s in reference.items()}
    assert all(matches.values())
    # Check all six tetrahedron edges in each layer as sequence complements.
    core_edges = []
    for g in range(3):
        for i, left in enumerate('abcd'):
            for right in 'abcd'[i+1:]:
                x, y = row[f'T{g}{left}.core'], row[f'T{g}{right}.core']
                length = next(n for n in range(min(len(x), len(y)), 0, -1)
                              if any(rc(x[j:j+n]) in y for j in range(len(x)-n+1)))
                sites = [(j, y.find(rc(x[j:j+length]))) for j in range(len(x)-length+1)]
                sites = [(j, k) for j, k in sites if k >= 0]
                mismatches = min(sum(a != b for a, b in zip(rc(x[j:j+17]), y[k:k+17]))
                                 for j in range(len(x)-16) for k in range(len(y)-16))
                core_edges.append(dict(left=f'T{g}{left}', right=f'T{g}{right}',
                                       longest_exact_complement=length, offsets=sites,
                                       expected_edge_bases=17, exact_17nt_edge=length >= 17,
                                       best_17nt_mismatches=mismatches))
    variants = []
    for name in ('hyb1-Cy5', 'hyb1-Cy3', 'hyb1-FAM', 'hyb2-Cy5', 'hyb2-Cy3', 'hyb2-FAM'):
        published = td['T0a-' + name]
        q = copy.deepcopy(p)
        sticky = re.sub(r'\s', '', published)[:30]
        q['probes']['primary']['parts']['readout'] = {'template': rc(sticky), 'parts': {}}
        if '-SS-' not in published:
            q['probes']['T0a'].pop('modifications')
        color = name.split('-')[-1]
        q['probes']['Detection']['modifications']['5prime']['name'] = color
        built = UProbeAPI(q, {}, out, require_genome=False).construct_probes(windows).iloc[0]
        assert built['T0a'] == re.sub(r'\s', '', published).replace('-SS-', '')
        assert built['Detection'] == re.sub(r'^(Cy5|Cy3|FAM)-\s*', '', td[color+'-ssDNA'])
        assert built['Detection.modifications'][0]['name'] == color
        variants.append(dict(variant=name, exact=True, fluorophore=color,
                             disulfide_explicit='-SS-' in published))
    interfaces = [('primary.readout', 'T0a.sticky')]
    interfaces += [(f'T{g}{s}.sticky', f'T{g+1}a.sticky') for g in range(2) for s in 'bcd']
    interfaces += [(f'T2{s}.sticky', 'Detection.binding') for s in 'bcd']
    assert all(rc(row[a]) == row[b] for a, b in interfaces)
    dag = DAG()
    dag.from_config(p)
    edges = {n.name: [d.name for d in n.deps] for n in dag.nodes}
    mutations = []
    for source, target in [('primary.readout', 'T0a.sticky'), ('T0b.sticky', 'T1a.sticky'),
                           ('T1b.sticky', 'T2a.sticky'), ('T2b.sticky', 'Detection.binding')]:
        q = copy.deepcopy(p)
        probe, part = source.split('.')
        replacement = 'A' * len(row[source])
        q['probes'][probe]['parts'][part] = {'template': replacement, 'parts': {}}
        changed = UProbeAPI(q, {}, out, require_genome=False).construct_probes(windows).iloc[0]
        assert changed[target] == rc(replacement) and changed[target] != row[target]
        mutations.append(source + ' -> ' + target)
    bond = len(td['T0a-hyb1-Cy5'].split('-SS-')[0])
    expected = [{'probe': 'T0a', 'type': 'disulfide', 'position': 'internal', 'after_base': bond}]
    assert bond == 30 and row['T0a.modifications'] == expected
    detection_mod = row['Detection.modifications']
    assert detection_mod == [{'probe': 'Detection', 'type': 'fluorophore', 'position': '5prime', 'name': 'Cy5'}]
    api.post_process_probes(native, raw_csv=True)
    checked = []
    for key in ('tsv',):
        table = pd.read_csv(out / api.output_files[key], sep='\t')
        assert json.loads(table.loc[0, 'T0a.modifications']) == expected
        assert json.loads(table.loc[0, 'Detection.modifications']) == detection_mod
        assert table.loc[0, 'T0a'] == clean(reference['T0a'])
        checked.append(key)
    from uprobe.core.report.run_summary import save_tsv
    save_tsv(api._report_raw_data, out / 'raw-check.tsv')
    raw = pd.read_csv(out / 'raw-check.tsv', sep='\t')
    assert json.loads(raw.loc[0, 'T0a.modifications']) == expected
    assert json.loads(raw.loc[0, 'Detection.modifications']) == detection_mod
    checked.append('raw_tsv')
    for key in ('raw_xlsx', 'xlsx'):
        with ZipFile(out / api.output_files[key]) as archive:
            values = [''.join(e.itertext()) for e in ET.fromstring(archive.read('xl/sharedStrings.xml'))]
            assert any(v.startswith('[{') and json.loads(v) == expected for v in values)
            assert any(v.startswith('[{') and json.loads(v) == detection_mod for v in values)
        checked.append(key)
    result = {'scope': 'Assembly consistency; fixed bodies are declared inputs, not recovered sequences',
              'source_sha256': hashlib.sha256(Path(args.workbook).read_bytes()).hexdigest(),
              'target': 'synthetic 35-nt control; no published binding sequence supplied',
              'scaffold_matches': matches, 'scaffold_strands_exact': '12/12',
              'groups': p['meta']['groups'], 'groups_semantics': 'declarative metadata; dependencies enforced by native DAG',
              'fixed_inputs': p['meta']['fixed_inputs'], 'dag_edges': edges,
              'complementary_interfaces': len(interfaces), 'dependency_perturbations_passed': mutations,
              'tetrahedral_core_edges': core_edges, 'panel_variants': variants,
              'limitations': 'Native DAG assembles strand sequences and modification annotations; tetrahedron membership and 1:3:9:27 copy topology are metadata, not enforced molecular complexes. Two hyb2 T0a entries lack an explicit SS marker; none is inferred.',
              'chemical_modifications': {'disulfide_after_base': bond, 'detection': detection_mod, 'passed': True},
              'native_exports_passed': checked, 'CPSR': 'N/A', 'folding_or_experimental_performance': 'Not evaluated'}
    (out / 'assembly.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('scaffold_strands_exact', 'complementary_interfaces', 'native_exports_passed')}))


if __name__ == '__main__':
    main()
