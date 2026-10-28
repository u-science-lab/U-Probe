import asyncio
from copy import deepcopy
from io import BytesIO
from types import SimpleNamespace
import pytest
import yaml
from fastapi import HTTPException, UploadFile
from uprobe.http.utils.target_inputs import normalize_target_inputs, normalize_target_name

@pytest.mark.parametrize('edge', [' ', '\t', '\r\n', '\u00a0', '\u3000', '\u200b', '\ufeff', '\u2060'])
def test_boundary_cleanup(edge):
    assert normalize_target_name(edge + 'Gapdh' + edge) == 'Gapdh'

@pytest.mark.parametrize('name', ['', ' \t\u200b', 'Gap dh', 'Gap\tdh', 'Gap\ndh', 'Gap\u200bdh', 'Gap\x00dh', 123])
def test_invalid_names(name):
    with pytest.raises(ValueError):
        normalize_target_name(name)

def test_preserve_sequences_and_sync_keys():
    config = {'targets': [{' Gapdh\t': 'ACTG'}], 'encoding': {' Gapdh\t': {'BC1': 'ACT'}}, 'post_process': {'equal_space': {' Gapdh\t': {'number_desired': 10}}}}
    original = deepcopy(config)
    result = normalize_target_inputs(config)
    assert result['targets'] == [{'Gapdh': 'ACTG'}]
    assert result['encoding'] == {'Gapdh': {'BC1': 'ACT'}}
    assert 'Gapdh' in result['post_process']['equal_space']
    assert config == original

@pytest.mark.parametrize('targets', [[], 'Gapdh', ['Gapdh', ' Gapdh\t'], ['Gapdh', {'Gapdh': 'ACTG'}], [None], [{}]])
def test_invalid_target_lists(targets):
    with pytest.raises(ValueError):
        normalize_target_inputs({'targets': targets})

def test_mapping_collision():
    with pytest.raises(ValueError, match='Duplicate'):
        normalize_target_inputs({'targets': ['Gapdh'], 'encoding': {'Gapdh': {}, 'Gapdh\t': {}}})

def test_submit_normalizes_before_persist(monkeypatch):
    import uprobe.http.routers.workflow as route
    original = b'"attributes":\n  "arm3_tm":\n    "target": "target_parts.arm3"\n    "type": "annealing_temperature"\n  "arm5_tm":\n    "target": "target_parts.arm5"\n    "type": "annealing_temperature"\n  "target_gc":\n    "target": "target_region"\n    "type": "gc_content"\n  "target_length":\n    "target": "target_region"\n    "type": "length"\n  "target_tm":\n    "target": "target_region"\n    "type": "annealing_temperature"\n"description": "Protocol for designing STARmap probes from species GRCm39"\n"encoding":\n  "Gapdh\\t":\n    "BC1": "GCTTATGATTGT"\n    "BC2": "CCGCTT"\n    "BC3": "CTCCCACTGACA"\n    "BC4": "CGATTG"\n"extracts":\n  "target_region":\n    "layout":\n      "parts":\n        "arm3":\n          "length":\n          - !!int "19"\n          - !!int "25"\n        "arm5":\n          "length":\n          - !!int "19"\n          - !!int "25"\n        "gap":\n          "length":\n          - !!int "0"\n          - !!int "2"\n      "template": "{arm5}{gap}{arm3}"\n    "length":\n    - !!int "40"\n    - !!int "46"\n    "source": "CDS"\n    "step": !!int "1"\n"genome": "GRCm39"\n"name": "starmap_514181"\n"post_process":\n  "filters":\n    "arm3_tm":\n      "condition": "arm3_tm >= 60"\n    "arm5_tm":\n      "condition": "arm5_tm >= 60"\n    "target_gc":\n      "condition": "target_gc >= 0.3 & target_gc <= 0.7"\n    "target_length":\n      "condition": "target_length >= 40 & target_length <= 46"\n    "target_tm":\n      "condition": "target_tm >= 56 & target_tm <= 70"\n"probes":\n  "pad_probe":\n    "modifications":\n      "3prime": "hydroxyl"\n      "5prime": "phosphate"\n    "parts":\n      "part1":\n        "expr": "encoding[target][\'BC2\']"\n      "part2":\n        "expr": "rc(target_parts[\'arm5\'])"\n      "part3":\n        "expr": "encoding[target][\'BC1\']"\n      "part4":\n        "expr": "encoding[target][\'BC3\']"\n      "part5":\n        "expr": "encoding[target][\'BC4\']"\n    "template": "{part1}{part2}{part3}{part4}{part5}"\n  "primer_probe":\n    "parts":\n      "arm3":\n        "expr": "rc(pad_probe.part1) + rc(pad_probe.part5)"\n      "arm5":\n        "expr": "rc(target_parts[\'arm3\'])"\n    "template": "{arm5}{arm3}"\n"summary":\n  "attributes":\n  - "target_length"\n  - "target_gc"\n  - "target_tm"\n  "report_name": "rna_report"\n"targets":\n- "Gapdh\\t"\n'
    saved = []
    monkeypatch.setattr(route, 'update_task_in_db', lambda user, task: saved.append(task))
    asyncio.run(route.submit_task(UploadFile(file=BytesIO(original), filename='protocol.yaml'), SimpleNamespace(username='test')))
    config = yaml.safe_load(saved[0].yaml_content)
    assert config['targets'] == ['Gapdh']
    assert list(config['encoding']) == ['Gapdh']

def test_submit_rejects_before_persist(monkeypatch):
    import uprobe.http.routers.workflow as route
    monkeypatch.setattr(route, 'update_task_in_db', lambda *args: pytest.fail('Must not persist invalid input'))
    with pytest.raises(HTTPException) as error:
        asyncio.run(route.submit_task(UploadFile(file=BytesIO(b'targets: ["Gap dh"]'), filename='protocol.yaml'), SimpleNamespace(username='test')))
    assert error.value.status_code == 422
