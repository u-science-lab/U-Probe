"""Chemical annotations stay separate from native DNA and survive both exports."""
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from uprobe.core.api import UProbeAPI
from uprobe.core.gen.probe import construct_probes
from uprobe.core.modifications import ModificationError, normalize_modifications
from uprobe.core.validate import validate_protocol
from uprobe.templates import load_template


def protocol():
    return {'name': 'chemical', 'genome': 'test', 'targets': {'G': 'ACGTACGT'},
            'extracts': {'target_region': {'source': 'mRNA', 'length': 8, 'step': 8}},
            'encoding': {}, 'probes': {'p': {'template': '{a}{b}', 'parts': {
                'a': {'expr': 'rc(target_region)'}, 'b': {'template': 'ACGT', 'parts': {}}}}},
            'attributes': {}, 'post_process': {}, 'summary': {}}


def test_legacy_sequence_and_columns_unchanged(tmp_path):
    p = load_template('padlock_iss')
    p['encoding'] = {'G': {'ID': 'ACGT' * 5}}
    p['targets'] = {'G': 'ACGT' * 10}
    for probe in p['probes'].values():
        probe.pop('modifications', None)
    api = UProbeAPI(p, {}, tmp_path, require_genome=False)
    targets = api.generate_target_seqs()
    expected = construct_probes(p, api._probe_contexts(targets, list(range(len(targets)))))
    actual = api.construct_probes(targets)
    pd.testing.assert_frame_equal(actual, expected)
    assert not any(c.endswith('.modifications') for c in actual)


def test_padlock_terminal_validation_and_api(tmp_path):
    p = load_template('padlock_iss')
    p['encoding'] = {'G': {'ID': 'ACGT' * 5}}
    p['targets'] = {'G': 'ACGT' * 10}
    probe = next(iter(p['probes']))
    p['probes'][probe]['modifications'] = {'5prime': 'phosphate'}
    assert not [f for f in validate_protocol(p) if f.severity == 'error' and 'modifications' in f.path]
    api = UProbeAPI(p, {}, tmp_path, require_genome=False)
    row = api.construct_probes(api.generate_target_seqs()).iloc[0]
    assert row[probe + '.modifications'] == [{'probe': probe, 'type': 'phosphate', 'position': '5prime'}]
    assert set(row[probe]) <= set('ACGTN')


def test_internal_and_terminal_survive_raw_final_exports(tmp_path, monkeypatch):
    from uprobe.core import api as module
    p = protocol()
    p['probes']['p']['modifications'] = {'5prime': 'phosphate', 'internal': [
        {'type': 'disulfide', 'after_base': 4}, {'type': 'phosphorothioate', 'after_base': 8}]}
    p['probes']['label'] = {'template': '{body}', 'parts': {'body': {'expr': 'rc(p)'}}, 'modifications': {
        '3prime': {'type': 'fluorophore', 'name': 'Cy5'}}}
    api = UProbeAPI(p, {}, tmp_path, require_genome=False)
    native = api.construct_probes(pd.DataFrame({'target': ['G'], 'target_region': ['ACGTACGT']}))
    assert native.loc[0, 'p'] == 'ACGTACGTACGT'
    assert len(native.loc[0, 'p']) == 12
    assert native.loc[0, 'p.modifications'][1]['after_base'] == 4
    monkeypatch.setattr(module, 'add_attributes', lambda df, *args: df.copy())
    final = api.post_process_probes(native, raw_csv=True)
    assert final.loc[0, 'p.modifications'] == native.loc[0, 'p.modifications']
    from uprobe.core.report.run_summary import save_tsv
    raw_tsv = tmp_path / 'raw-check.tsv'
    save_tsv(native, raw_tsv)
    assert api._report_raw_data.loc[0, 'p.modifications'] == native.loc[0, 'p.modifications']
    assert json.loads(pd.read_csv(raw_tsv, sep='\t').loc[0, 'p.modifications']) == native.loc[0, 'p.modifications']
    for key in ('tsv', 'raw_xlsx', 'xlsx'):
        path = tmp_path / api.output_files[key]
        table = pd.read_excel(path) if key.endswith('xlsx') else pd.read_csv(path, sep='\t')
        for probe in ('p', 'label'):
            assert json.loads(table.loc[0, probe + '.modifications']) == native.loc[0, probe + '.modifications']
            assert table.loc[0, probe] == native.loc[0, probe]


@pytest.mark.parametrize('value', [
    {'5prime': 'unknown'}, {'middle': 'phosphate'},
    {'internal': [{'type': 'disulfide', 'after_base': x} for x in [0]]},
    {'internal': [{'type': 'disulfide', 'after_base': True}]},
    {'internal': [{'type': 'disulfide', 'after_base': 1.5}]},
    {'internal': [{'type': 'disulfide', 'after_base': 2}, {'type': 'disulfide', 'after_base': 2}]},
    {'internal': [{'type': 'disulfide', 'after_base': 2}, {'type': 'phosphorothioate', 'after_base': 2}]},
    {'5prime': 'disulfide'}, {'internal': [{'type': 'phosphate', 'after_base': 2}]},
])
def test_invalid_modifications_report_yaml_path(value):
    p = protocol()
    p['probes']['p']['modifications'] = value
    errors = [f for f in validate_protocol(p) if f.severity == 'error']
    assert any(f.path.startswith('probes.p.modifications') for f in errors)
    with pytest.raises(ModificationError, match='probes.p.modifications'):
        normalize_modifications(value, 'p')


def test_runtime_length_and_dna_symbols(tmp_path):
    p = protocol()
    p['probes']['p']['modifications'] = {'internal': [{'type': 'disulfide', 'after_base': 12}]}
    api = UProbeAPI(p, {}, tmp_path, require_genome=False)
    with pytest.raises(ValueError, match='probes.p.modifications.internal'):
        api.construct_probes(pd.DataFrame({'target': ['G'], 'target_region': ['ACGTACGT']}))
    p = protocol()
    p['probes']['p']['parts']['b']['template'] = 'AC-SS-GT'
    assert any(f.severity == 'error' and f.path.startswith('probes.p.parts.b') for f in validate_protocol(p))


def test_nested_modifications_are_rejected():
    p = protocol()
    p['probes']['p']['parts']['a']['modifications'] = {'5prime': 'phosphate'}
    assert any(f.severity == 'error' and f.path == 'probes.p.parts.a.modifications' for f in validate_protocol(p))


def test_duplicate_yaml_modifications_rejected(tmp_path):
    from uprobe.core.validate import validate_files
    text = yaml.safe_dump(protocol(), sort_keys=False).replace(
        '    parts:', '    modifications:\n      5prime: phosphate\n      5prime: hydroxyl\n    parts:', 1)
    path = tmp_path / 'duplicate.yaml'
    path.write_text(text)
    assert any(f.severity == 'error' and 'modifications' in f.path for f in validate_files(path))
    with pytest.raises(ValueError, match='duplicate'):
        UProbeAPI(path, {}, tmp_path, require_genome=False)


def test_reindexed_sequence_rechecks_modification_bonds(tmp_path):
    p = protocol()
    p['probes']['p']['parts']['a']['expr'] = "encoding[target]['pieces'][probe_index % 2]"
    p['encoding'] = {'G': {'pieces': ['ACGTACGT', 'AC']}}
    p['probes']['p']['modifications'] = {'internal': [{'type': 'disulfide', 'after_base': 8}]}
    api = UProbeAPI(p, {}, tmp_path, require_genome=False)
    target = pd.DataFrame({'target': ['G', 'G', 'G'], 'target_region': ['ACGTACGT'] * 3,
                           'start': [0, 8, 16], 'end': [8, 16, 24]})
    # Both retained provisional candidates have the long sequence; renumbering
    # the second retained candidate produces a shorter final strand.
    from uprobe.core.gen.probe import construct_probes
    provisional = construct_probes(p, api._probe_contexts(target.iloc[[0, 2]], [0, 2]))
    provisional = pd.concat([target.iloc[[0, 2]].reset_index(drop=True), provisional.reset_index(drop=True)], axis=1)
    provisional = api._annotate_modifications(provisional)
    with pytest.raises(ValueError, match='probes.p.modifications.internal'):
        api._number_final_probes(provisional)


@pytest.mark.parametrize('name', ['ExSeq_web.yaml', 'STARmap_web.yaml'])
def test_existing_modified_presets(name):
    path = Path(__file__).parents[1] / 'examples' / 'web-presets' / name
    if not path.exists():
        pytest.skip('preset not distributed')
    p = yaml.safe_load(path.read_text())
    assert not [f for f in validate_protocol(p) if f.severity == 'error' and 'modifications' in f.path]
