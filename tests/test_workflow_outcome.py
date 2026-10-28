import json

import pytest

from uprobe.http.utils.workflow_outcome import workflow_failure, stored_workflow_failure


@pytest.mark.parametrize('status', ['no_valid_targets', 'no_target_sequences', 'no_probes_constructed', 'error'])
def test_early_exit_is_failure_with_reason(status):
    message = workflow_failure({'status': status, 'message': 'Stopped before report.',
                                'inputs': {'genome': 'hg19', 'invalid_targets': ['Cd4']}})
    assert 'Stopped before report.' in message
    assert 'hg19' in message and 'Cd4' in message


@pytest.mark.parametrize('status', ['ok', 'partial', 'no_probes'])
def test_normal_results_and_raw_fallback_remain_successful(status):
    assert workflow_failure({'status': status}) is None


def test_historical_summary_and_missing_summary(tmp_path):
    assert stored_workflow_failure(tmp_path) is None
    path = tmp_path / 'summary.json'
    path.write_text('{bad json')
    assert stored_workflow_failure(tmp_path) is None
    path.write_text(json.dumps({'status': 'no_valid_targets', 'inputs': {'invalid_targets': ['Cd4']}}))
    assert 'including case' in stored_workflow_failure(tmp_path)


def test_task_list_repairs_historical_state_before_filtering(tmp_path, monkeypatch):
    import asyncio
    from types import SimpleNamespace
    from uprobe.http.routers import task as router

    directory = tmp_path / 't1'
    directory.mkdir()
    (directory / 'summary.json').write_text(json.dumps({'status': 'no_valid_targets'}))
    record = dict(id='t1', name='test', genome='hg19', parameters={'probe_length': 32},
                  status='completed', progress=100, result_url='empty.zip',
                  created_at='2026-10-27T10:00:00', updated_at='2026-10-27T10:00:00')
    saved = []
    monkeypatch.setattr(router, 'load_user_tasks', lambda _: [record])
    monkeypatch.setattr(router, 'get_results_dir', lambda: tmp_path)
    monkeypatch.setattr(router, 'update_task_in_db', lambda user, task: saved.append(task))
    tasks = asyncio.run(router.get_tasks(status_filter='failed', search=None,
                                        current_user=SimpleNamespace(username='tester')))
    assert len(tasks) == 1 and saved[0].status == 'failed'
    assert tasks[0].result_url is None and tasks[0].progress == 0
    assert 'exact gene names' in tasks[0].error_message


def test_runner_does_not_package_early_exit(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from uprobe.core import api
    from uprobe.http.utils import task_control, paths
    from uprobe.http.utils.uprobe_runner import run_uprobe_workflow

    monkeypatch.setattr(task_control, 'register_task_process', lambda *a, **k: None)
    monkeypatch.setattr(task_control, 'cleanup_task_process', lambda *a, **k: None)
    monkeypatch.setattr(paths, 'get_genomes_yaml', lambda: tmp_path / 'absent.yaml')
    monkeypatch.setattr(paths, 'get_user_genomes_yaml', lambda _: tmp_path / 'absent.yaml')
    monkeypatch.setenv('PATH', __import__('os').environ['PATH'])
    fake = SimpleNamespace(summary={'status': 'no_valid_targets'}, run_workflow=lambda **k: None)
    monkeypatch.setattr(api, 'UProbeAPI', lambda **k: fake)
    result = run_uprobe_workflow(protocol_yaml='name: test', username='tester', task_id='t1',
                                 output_dir=str(tmp_path), threads=1, raw_csv=False,
                                 continue_invalid_targets=True, log_path=str(tmp_path / 'run.log'))
    assert result['ok'] is False
    assert not list(tmp_path.glob('*.zip'))
