"""User-facing task failures; full tracebacks remain in run.log."""
import re
import subprocess


class TaskInputError(ValueError):
    pass


def validate_target_whitespace(protocol):
    """Reject whitespace in explicit sequences without changing user input."""
    targets = protocol.get('targets', [])
    if isinstance(targets, dict):
        targets = [targets]
    if not isinstance(targets, list):
        return
    for entry in targets:
        if not isinstance(entry, dict):
            continue
        for target, sequence in entry.items():
            if not isinstance(sequence, str):
                continue
            # Coordinate expressions and other target forms are not sequences.
            if not re.fullmatch(r'[ACGTUNRYKMSWBDHVacgtunrykmswbdhv\s]+', sequence):
                continue
            positions = [i + 1 for i, c in enumerate(sequence) if c.isspace()]
            if positions:
                raise TaskInputError(
                    f'目标 {target} 的输入序列含有 {len(positions)} 个空白字符，'
                    f'首个位置为第 {positions[0]} 个字符（从 1 开始）。'
                    '请在 targets 中删除该序列的空格、换行或制表符后重新运行。'
                    '任务已在设计前停止，尚未生成探针。')



PUBLIC_ERRORS = {'unknown': 'Task failed. The cause could not be determined. Please contact the administrator.', 'whitespace': 'Input sequences contain whitespace. Remove spaces, line breaks or tabs and retry.', 'quality': 'Sequence alignment failed: sequence and quality lengths do not match. Check the input sequences.', 'alignment': 'Sequence alignment failed. Please contact the administrator.', 'external': 'An external computation failed. Please contact the administrator.', 'missing': 'A required file is missing. Please contact the administrator.', 'permission': 'The task cannot access a required resource. Please contact the administrator.', 'memory': 'Insufficient memory. Reduce the task size and retry.', 'timeout': 'The task timed out. Please retry later.', 'restart': 'The task was interrupted by a server restart. Please rerun it.', 'results': 'The task did not produce complete results. Check the design parameters and retry.', 'no_valid_targets': 'No valid targets were found. Check the species, gene names or transcript IDs.', 'no_target_sequences': 'No target sequences were extracted. Check coordinates, region lengths and transcript selection.', 'no_probes_constructed': 'No probes were generated. Check probe expressions, encoding and barcode settings.'}

LEGACY_ERRORS = {'任务执行失败，暂无法确定具体原因，请联系管理员。': 'unknown', '输入序列含空白字符，请删除空格、换行或制表符后重试。': 'whitespace', '序列比对失败：序列与质量值长度不一致，请检查输入序列。': 'quality', '序列比对程序执行失败，请联系管理员检查。': 'alignment', '外部计算程序执行失败，请联系管理员检查。': 'external', '所需文件缺失，请联系管理员检查。': 'missing', '任务无法访问所需资源，请联系管理员检查。': 'permission', '计算内存不足，请减少任务规模后重试。': 'memory', '任务执行超时，请稍后重试。': 'timeout', '服务重启导致任务中断，请重新运行。': 'restart', '任务未生成完整结果，请检查设计参数后重试。': 'results', '未找到有效靶标，请检查物种、基因名称或转录本编号。': 'no_valid_targets', '未提取到靶标序列，请检查坐标、区域长度和转录本选择。': 'no_target_sequences', '未生成探针，请检查探针表达式、编码和条形码配置。': 'no_probes_constructed'}

def public_task_error(value):
    """Only fixed public text is returned, including for legacy failures."""
    if value is None:
        return None
    text = str(value)
    if text in PUBLIC_ERRORS.values():
        return text
    if text in LEGACY_ERRORS:
        return PUBLIC_ERRORS[LEGACY_ERRORS[text]]
    lower = text.lower()
    rules = [
        ('whitespace', ('空白字符',)),
        ('quality', ('more quality values than read characters', 'fewer quality values than read characters', '序列与质量值长度不一致')),
        ('no_valid_targets', ('no_valid_targets', 'targets not found', 'no valid targets')),
        ('no_target_sequences', ('no_target_sequences', 'no target sequences')),
        ('no_probes_constructed', ('no_probes_constructed', 'no probes constructed')),
        ('restart', ('server restarted before this task finished',)),
        ('results', ('zip archive is missing', 'without result files')),
        ('memory', ('memoryerror', 'cannot allocate memory', 'out of memory')),
        ('permission', ('permissionerror', 'permission denied')),
        ('missing', ('filenotfounderror', 'no such file or directory')),
        ('timeout', ('timeouterror', 'timeoutexpired')),
        ('alignment', ('bowtie2',)),
        ('external', ('calledprocesserror', '外部计算程序')),
    ]
    for key, markers in rules:
        if any(marker in lower for marker in markers):
            return PUBLIC_ERRORS[key]
    return PUBLIC_ERRORS['unknown']

def task_error_message(error):
    if isinstance(error, TaskInputError):
        return PUBLIC_ERRORS['whitespace']
    if isinstance(error, subprocess.CalledProcessError):
        diagnostic = getattr(error, 'tool_diagnostic', '')
        if any(marker in diagnostic for marker in ('more quality values than read characters', 'fewer quality values than read characters')):
            return PUBLIC_ERRORS['quality']
        return PUBLIC_ERRORS['alignment' if 'bowtie2' in str(error.cmd) else 'external']
    for kind, key in ((FileNotFoundError, 'missing'), (PermissionError, 'permission'), (MemoryError, 'memory'), (TimeoutError, 'timeout'), (subprocess.TimeoutExpired, 'timeout')):
        if isinstance(error, kind):
            return PUBLIC_ERRORS[key]
    return PUBLIC_ERRORS['unknown']
