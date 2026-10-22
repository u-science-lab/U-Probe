import json
import os

import pytest
import yaml

from uprobe.http.utils.paths import get_probes_dir
from uprobe.http.utils.probe_library import load_builtin_probes, load_builtin_probe
from uprobe.http.utils.migrate_probes import migrate


@pytest.fixture
def library(tmp_path, monkeypatch):
    monkeypatch.setenv('UPROBE_SERVER_ROOT', str(tmp_path))
    (tmp_path / 'config.ini').write_text(f'[Paths]\ndata_dir = {tmp_path / "configured-data"}\n', encoding='utf8')
    root = get_probes_dir()
    root.mkdir()
    return root


def write(root, name, value='ACGT'):
    config = {'probes': {'pad': {'template': '{arm}', 'parts': {'arm': {'expr': repr(value)}}}}}
    (root / f'{name}.yaml').write_text(yaml.safe_dump(config), encoding='utf8')
    return config


def test_dynamic_add_edit_delete_and_single_bad_file(library):
    assert load_builtin_probes() == {}
    first = write(library, 'MiP-seq')
    assert load_builtin_probe('MiP-seq') == first
    changed = write(library, 'MiP-seq', 'TTTT')
    assert load_builtin_probe('MiP-seq') == changed
    write(library, 'ExSeq')
    assert set(load_builtin_probes()) == {'MiP-seq', 'ExSeq'}
    (library / 'ExSeq.yaml').write_text('probes: [broken', encoding='utf8')
    assert load_builtin_probes() == {'MiP-seq': changed}
    (library / 'MiP-seq.yaml').unlink()
    assert load_builtin_probes() == {}


def test_duplicate_keys_are_not_silently_overwritten(library):
    write(library, 'Good')
    (library / 'Bad.yaml').write_text('probes: {}\nprobes: {}\n', encoding='utf8')
    assert list(load_builtin_probes()) == ['Good']


def test_missing_directory_does_not_fall_back_to_json(library):
    (library.parent / 'probe.json').write_text('{"old": {"probes": {}}}', encoding='utf8')
    library.rmdir()
    with pytest.raises(FileNotFoundError):
        load_builtin_probes()


def test_config_directory_changes_are_dynamic(library, tmp_path):
    other = tmp_path / 'other/probes'
    other.mkdir(parents=True)
    expected = write(other, 'New')
    (tmp_path / 'config.ini').write_text(f'[Paths]\ndata_dir = {other.parent}\n', encoding='utf8')
    assert get_probes_dir() == other
    assert load_builtin_probes() == {'New': expected}


def test_migration_preserves_every_field_and_unfinished_methods(tmp_path):
    original = {'MiP-seq': {'probes': {}, 'extracts': {'target_region': {'length': [40, 46]}},
                            'references': ['original'], 'design_notes': ['unchanged']},
                'OpenFISH': {'probes': {'pad': {'template': '{missing}', 'parts': {}}}}}
    source = tmp_path / 'probe.json'
    source.write_text(json.dumps(original), encoding='utf8')
    target = tmp_path / 'probes'
    migrate(source, target)
    assert {p.stem: yaml.safe_load(p.read_text()) for p in target.glob('*.yaml')} == original
    assert json.loads(source.read_text()) == original
    with pytest.raises(FileExistsError):
        migrate(source, target)


def test_existing_http_contract(library):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from uprobe.http.routers.workflow import workflow
    expected = {'MiP-seq': write(library, 'MiP-seq')}
    app = FastAPI()
    app.include_router(workflow)
    client = TestClient(app)
    assert client.get('/workflow/builtin_probes').json() == expected
    write(library, 'MiP-seq', 'GGGG')
    assert client.get('/workflow/builtin_probes').json() == load_builtin_probes()


@pytest.mark.skipif(os.name == 'nt', reason='agent store requires Unix fcntl')
def test_ai_resolves_same_configured_directory(library, tmp_path, monkeypatch):
    from uprobe.http.agent.output_sandbox import resolve_agent_paths, apply_agent_runtime_env, probe_resource_context
    for key in ['UPROBE_PROJECT_ROOT', 'UPROBE_WORKSPACE_ROOT', 'UPROBE_DATA_DIR',
                'UPROBE_PROBES_DIR', 'UPROBE_GENOMES_PATH', 'UPROBE_OUTPUT_DIR', 'UPROBE_PROBE_JSON']:
        monkeypatch.delenv(key, raising=False)
    paths = resolve_agent_paths(entrypoint='web', output_dir=tmp_path / 'out', memory_dir=tmp_path / 'memory')
    assert paths.probes_dir == library
    assert paths.genomes_path == library.parent / 'genomes.yaml'
    apply_agent_runtime_env(paths)
    assert os.environ['UPROBE_PROBES_DIR'] == str(library)
    assert 'UPROBE_PROBE_JSON' not in os.environ
    assert str(library) in probe_resource_context()
