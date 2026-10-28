import subprocess

import pytest

from uprobe.http.utils.task_errors import (
    TaskInputError, task_error_message, validate_target_whitespace,
)


def test_names_target_and_exact_whitespace_position():
    with pytest.raises(TaskInputError, match='Gapdh.*2.*第 5') as caught:
        validate_target_whitespace({'targets': [{'Gapdh': 'ACGT A\nC'}]})
    assert 'whitespace' in task_error_message(caught.value)
    assert 'Gapdh' not in task_error_message(caught.value)


@pytest.mark.parametrize('targets', [
    ['Gapdh'], [{'Gapdh': 'ACGTN'}], [{'region': 'chr1:1-40'}],
])
def test_other_target_forms_remain_valid(targets):
    validate_target_whitespace({'targets': targets})


def test_exit_code_does_not_claim_input_is_cause():
    message = task_error_message(subprocess.CalledProcessError(134, 'bowtie2 -x /private/path'))
    assert 'Sequence alignment failed' in message
    assert '/private/path' not in message
    assert '134' not in message


def test_unknown_error_keeps_details_in_log_only():
    message = task_error_message(RuntimeError('/private/path secret'))
    assert 'RuntimeError' not in message
    assert 'secret' not in message


def test_bowtie_diagnostic_explains_fastq_failure():
    error = subprocess.CalledProcessError(134, 'bowtie2')
    error.tool_diagnostic = 'Error: Read 16 has more quality values than read characters.'
    assert 'lengths do not match' in task_error_message(error)


def test_legacy_error_privacy():
    from uprobe.http.utils.task_errors import public_task_error, PUBLIC_ERRORS
    for raw in ['Traceback SECRET /private/path', 'FileNotFoundError: SECRET', 'bowtie2 --password=SECRET', '目标 SECRET 的输入序列含有空白字符']:
        result = public_task_error(raw)
        assert result in PUBLIC_ERRORS.values()
        assert 'SECRET' not in result
        assert public_task_error(result) == result

def test_task_serialization_privacy():
    from uprobe.http.routers.task import TaskRead
    task = TaskRead(id='test', name='test', genome='test', parameters={'probe_length': 20},
                    status='failed', progress=0, created_at='2026-10-28T00:00:00',
                    updated_at='2026-10-28T00:00:00', error_message='Traceback SECRET')
    assert 'SECRET' not in task.model_dump()['error_message']
    task.error_message = 'FileNotFoundError: SECRET'
    assert 'SECRET' not in task.model_dump_json()

def test_workflow_privacy():
    from uprobe.http.utils.workflow_outcome import workflow_failure
    result = workflow_failure({'status': 'no_valid_targets', 'message': 'SECRET',
                               'inputs': {'genome': '/private', 'invalid_targets': ['SECRET']}})
    assert 'No valid targets' in result and 'SECRET' not in result and '/private' not in result

def test_all_public_errors_are_english_and_legacy_reasons_are_preserved():
    from uprobe.http.utils.task_errors import PUBLIC_ERRORS, LEGACY_ERRORS, public_task_error
    assert all(message.isascii() for message in PUBLIC_ERRORS.values())
    for legacy, key in LEGACY_ERRORS.items():
        assert public_task_error(legacy) == PUBLIC_ERRORS[key]
