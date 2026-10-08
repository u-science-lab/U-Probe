import json
import os
import tempfile

import pytest

# The HTTP modules resolve data_dir at import time; keep them off the configured server paths.
os.environ.setdefault("UPROBE_SERVER_ROOT", tempfile.mkdtemp(prefix="uprobe-http-test-"))

from uprobe.http.routers import auth  # noqa: E402
from uprobe.http.utils import task_queue  # noqa: E402


@pytest.fixture
def users_db(tmp_path, monkeypatch):
    path = tmp_path / "users_db.json"
    monkeypatch.setattr(auth, "USERS_DB_FILE", path)
    monkeypatch.setattr(auth, "USERS_DB_LOCK", path.with_suffix(".lock"))
    return path


def test_no_default_root_password(users_db, monkeypatch):
    monkeypatch.delenv("UPROBE_ROOT_PASSWORD", raising=False)
    assert auth.load_users_db() == {}
    assert json.loads(users_db.read_text()) == {}


def test_root_seeded_from_environment(users_db, monkeypatch):
    monkeypatch.setenv("UPROBE_ROOT_PASSWORD", "s3cret-pass")
    db = auth.load_users_db()
    assert auth.pwd_context.verify("s3cret-pass", db["root"]["hashed_password"])
    assert not auth.pwd_context.verify("123456", db["root"]["hashed_password"])


def test_unreadable_users_db_is_not_replaced(users_db):
    users_db.write_text('{"alice": ')
    with pytest.raises(json.JSONDecodeError):
        auth.load_users_db()
    assert users_db.read_text() == '{"alice": '


def test_secret_key_is_persisted(tmp_path):
    path = tmp_path / ".secret_key"
    first = auth.load_or_create_secret_key(path)
    assert auth.load_or_create_secret_key(path) == first
    assert path.stat().st_mode & 0o777 == 0o600


def test_task_slots_are_shared_between_lock_holders(tmp_path, monkeypatch):
    # Separate open() calls stand in for separate uvicorn workers: flock conflicts across them.
    monkeypatch.setattr(task_queue, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(task_queue, "MAX_CONCURRENT_TASKS", 2)
    held = [task_queue._try_acquire_slot_file() for _ in range(2)]
    assert all(held)
    assert task_queue._try_acquire_slot_file() is None
    held[0].close()
    again = task_queue._try_acquire_slot_file()
    assert again is not None
    for handle in (held[1], again):
        handle.close()


@pytest.mark.asyncio
async def test_task_slot_waits_for_a_free_slot(tmp_path, monkeypatch):
    import asyncio

    monkeypatch.setattr(task_queue, "get_data_dir", lambda: tmp_path)
    monkeypatch.setattr(task_queue, "MAX_CONCURRENT_TASKS", 1)
    other_worker = task_queue._try_acquire_slot_file()
    entered = asyncio.Event()

    async def run():
        async with task_queue.task_slot(poll_seconds=0.01):
            entered.set()

    job = asyncio.create_task(run())
    await asyncio.sleep(0.1)
    assert not entered.is_set()
    other_worker.close()
    await asyncio.wait_for(job, timeout=2)
    assert entered.is_set()
