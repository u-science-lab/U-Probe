import asyncio
import contextlib
import fcntl
import multiprocessing
import os
import logging
from uprobe.http.utils.paths import get_config, get_data_dir

log = logging.getLogger(__name__)

def _load_queue_config():
    config = get_config()
    
    # Default threads per task
    threads = 10
    if config.has_section("TaskQueue") and config.has_option("TaskQueue", "task_threads"):
        threads = config.getint("TaskQueue", "task_threads")
    
    # Allow environment variable override
    threads = int(os.getenv("UPROBE_TASK_THREADS", str(threads)))
    
    total_cores = multiprocessing.cpu_count()
    
    # Calculate max concurrent tasks to ensure at least 1 task runs
    max_tasks = max(1, total_cores // threads)
    
    if config.has_section("TaskQueue") and config.has_option("TaskQueue", "max_concurrent_tasks"):
        val = config.get("TaskQueue", "max_concurrent_tasks").strip().lower()
        if val != "auto":
            try:
                configured_max_tasks = int(val)
                safe_max_tasks = max(1, total_cores // threads)
                if configured_max_tasks > safe_max_tasks:
                    log.warning(f"Configured max_concurrent_tasks ({configured_max_tasks}) exceeds safe limit based on CPU cores ({total_cores}) and threads per task ({threads}). Limiting to {safe_max_tasks}.")
                    max_tasks = safe_max_tasks
                else:
                    max_tasks = configured_max_tasks
            except ValueError:
                log.warning(f"Invalid max_concurrent_tasks value '{val}' in config, falling back to auto ({max_tasks})")
                
    return threads, max_tasks

TASK_THREADS, MAX_CONCURRENT_TASKS = _load_queue_config()

_task_semaphore = None

def get_task_semaphore() -> asyncio.Semaphore:
    """
    Get the global task concurrency semaphore.
    Must be initialized within an active asyncio event loop.
    """
    global _task_semaphore
    if _task_semaphore is None:
        log.info(f"Init task queue: total_cores={multiprocessing.cpu_count()}, threads_per_task={TASK_THREADS}, max_concurrent_tasks={MAX_CONCURRENT_TASKS}")
        _task_semaphore = asyncio.Semaphore(MAX_CONCURRENT_TASKS)
    return _task_semaphore


def _try_acquire_slot_file():
    """Return an flock-held slot file, or None when every slot is taken."""
    slot_dir = get_data_dir() / "task_slots"
    slot_dir.mkdir(parents=True, exist_ok=True)
    for i in range(MAX_CONCURRENT_TASKS):
        handle = open(slot_dir / f"slot-{i}.lock", "a")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            return handle
        except BlockingIOError:
            handle.close()
    return None


@contextlib.asynccontextmanager
async def task_slot(poll_seconds: float = 1.0):
    """
    Hold one of MAX_CONCURRENT_TASKS execution slots.

    The asyncio semaphore only limits a single server worker; the slot files
    extend the limit to all uvicorn workers sharing the same data_dir. flock
    locks are released when the holding process exits, so a crashed worker
    does not leak its slot.
    """
    async with get_task_semaphore():
        handle = _try_acquire_slot_file()
        while handle is None:
            await asyncio.sleep(poll_seconds)
            handle = _try_acquire_slot_file()
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()
