"""Translate core workflow summaries into actionable web task failures."""
import json
from pathlib import Path


def workflow_failure(summary):
    if not isinstance(summary, dict):
        return None
    status = summary.get('status')
    if status not in {'no_valid_targets', 'no_target_sequences', 'no_probes_constructed', 'error'}:
        return None
    from uprobe.http.utils.task_errors import PUBLIC_ERRORS
    return PUBLIC_ERRORS.get(status, PUBLIC_ERRORS['unknown'])


def stored_workflow_failure(directory: Path):
    try:
        return workflow_failure(json.loads((directory / 'summary.json').read_text(encoding='utf-8')))
    except (OSError, ValueError, TypeError):
        return None
