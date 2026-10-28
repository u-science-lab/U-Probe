import asyncio
from contextlib import asynccontextmanager
from datetime import datetime

import pytest


@pytest.mark.asyncio
async def test_task_validation_failure_keeps_reason_and_configuration(tmp_path, monkeypatch):
    from uprobe.http.routers import task
    from uprobe.http.utils import task_queue

    @asynccontextmanager
    async def slot():
        yield

    monkeypatch.setattr(task_queue, 'task_slot', slot)
    monkeypatch.setattr(task_queue, 'get_task_semaphore', lambda: asyncio.Semaphore(1))
    monkeypatch.setattr(task, 'get_tasks_dir', lambda: tmp_path / 'tasks')
    monkeypatch.setattr(task, 'get_results_dir', lambda: tmp_path / 'results')
    monkeypatch.setattr(task, 'get_process_pool', lambda: None)
    loop = asyncio.get_running_loop()
    future = loop.create_future()
    monkeypatch.setattr(loop, 'run_in_executor', lambda *args: future)
    now = datetime.now()
    record = task.TaskRead(id='phase-test', name='test', genome='mouse', parameters={'probe_length':32},
                           status='pending', phase='queued', progress=0, created_at=now, updated_at=now,
                           yaml_content='name: retained-configuration')
    task.update_task_in_db('tester', record)
    background = asyncio.create_task(task._run_uprobe_task('tester', record.id))
    await asyncio.sleep(0)
    current = task.find_task_by_id('tester', record.id)
    assert current.status == 'running' and current.phase == 'validating'
    future.set_result({'ok':False, 'validation_failed':True, 'error':'Missing canonical tag'})
    await background
    current = task.find_task_by_id('tester', record.id)
    assert current.status == 'failed' and current.phase == 'validation_failed'
    assert current.error_message == 'Missing canonical tag'
    assert current.yaml_content == record.yaml_content


def test_restart_marks_queued_work_failed_but_keeps_unstarted_tasks(tmp_path, monkeypatch):
    from uprobe.http.routers import task
    monkeypatch.setattr(task, 'get_tasks_dir', lambda: tmp_path)
    now = datetime.now()
    for identifier, phase in [('queued', 'queued'), ('manual', None)]:
        task.update_task_in_db('tester', task.TaskRead(
            id=identifier, name=identifier, genome='mouse', parameters={'probe_length':32},
            status='pending', phase=phase, progress=0, created_at=now, updated_at=now))
    task.reset_stuck_tasks_on_startup()
    queued = task.find_task_by_id('tester', 'queued')
    assert queued.status == 'failed' and queued.phase is None
    assert 'Server restarted' in queued.error_message
    assert task.find_task_by_id('tester', 'manual').status == 'pending'
