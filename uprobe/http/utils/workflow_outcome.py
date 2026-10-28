"""Translate core workflow summaries into actionable web task failures."""
import json
from pathlib import Path


def workflow_failure(summary):
    if not isinstance(summary, dict):
        return None
    status = summary.get('status')
    if status not in {'no_valid_targets', 'no_target_sequences', 'no_probes_constructed', 'error'}:
        return None
    message = summary.get('message') or f'Workflow stopped: {status}.'
    inputs = summary.get('inputs') or {}
    if inputs.get('genome'):
        message += f" Genome: {inputs['genome']}."
    invalid = inputs.get('invalid_targets') or []
    if invalid:
        message += ' Targets not found: ' + ', '.join(map(str, invalid)) + '.'
    advice = {
        'no_valid_targets': 'Check the selected species/genome and exact gene names (including case) or transcript IDs against its annotation.',
        'no_target_sequences': 'Check target coordinates, extraction source, region length and transcript selection.',
        'no_probes_constructed': 'Check probe expressions, required encoding keys and barcode sequences.',
    }
    if status in advice:
        message += ' ' + advice[status]
    return message


def stored_workflow_failure(directory: Path):
    try:
        return workflow_failure(json.loads((directory / 'summary.json').read_text(encoding='utf-8')))
    except (OSError, ValueError, TypeError):
        return None
