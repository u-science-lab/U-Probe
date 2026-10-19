"""The chat agent is the optional `web-agent` extra: the core, the CLI and the server work without it.

Each check runs in a subprocess, because importing `uprobe.http.server` reconfigures logging and
redirects stdout/stderr, and because the import blocker must be in place before anything imports
`uprobe`. With the blocker, every `import pantheon...` fails as if PantheonOS were not installed.
"""
import ast
import json
import os
import random
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / 'uprobe' / 'core'

BLOCK_PANTHEON = '''
import sys

class _BlockPantheon:
    def find_spec(self, name, path=None, target=None):
        if name == 'pantheon' or name.startswith('pantheon.'):
            raise ModuleNotFoundError(f'No module named {name!r} (blocked by the test)', name=name)
        return None

sys.meta_path.insert(0, _BlockPantheon())
'''

# Stand-ins for the PantheonOS names the agent imports, used only when PantheonOS is not installed.
FAKE_PANTHEON = '''
import sys, types

def _module(name, **attrs):
    module = types.ModuleType(name)
    module.__dict__.update(attrs)
    sys.modules[name] = module
    return module

try:
    import pantheon.chatroom
except ImportError:
    for name in ('pantheon', 'pantheon.factory', 'pantheon.utils'):
        _module(name, __path__=[])
    _module('pantheon.chatroom', ChatRoom=type('ChatRoom', (), {}))
    _module('pantheon.factory.template_io', UnifiedMarkdownParser=type('UnifiedMarkdownParser', (), {}))
    _module('pantheon.utils.vision', parse_image_mentions=lambda content: content)
    _module('pantheon.utils.log', set_level=lambda level: None)
'''

# The /agent API the web front end (uprobe-web-ui) calls; it must not change.
AGENT_API = {
    ('GET', '/agent/conversations'),
    ('POST', '/agent/conversations'),
    ('GET', '/agent/conversations/{conversation_id}'),
    ('PATCH', '/agent/conversations/{conversation_id}'),
    ('DELETE', '/agent/conversations/{conversation_id}'),
    ('POST', '/agent/conversations/{conversation_id}/message'),
    ('POST', '/agent/conversations/{conversation_id}/message/stream'),
    ('POST', '/agent/conversations/{conversation_id}/rewind'),
    ('POST', '/agent/conversations/{conversation_id}/rewind/stream'),
    ('POST', '/agent/conversations/{conversation_id}/stop'),
    ('POST', '/agent/conversations/{conversation_id}/clear'),
    ('POST', '/agent/conversations/{conversation_id}/upload'),
    ('DELETE', '/agent/conversations/{conversation_id}/upload/{attachment_id}'),
    ('GET', '/agent/conversations/{conversation_id}/uploads/{filename}'),
    ('GET', '/agent/conversations/{conversation_id}/files/{file_path}'),
}

# Starts the app (startup events included), calls it and prints one RESULT line: the server's
# logging setup sends everything else written to stdout through its log handlers.
SERVER_PROBE = '''
import json
from fastapi.testclient import TestClient
from uprobe.http.agent import web_agent_available
import uprobe.http.server as server

with TestClient(server.app) as client:
    paths = server.app.openapi()['paths']
    responses = {
        'root': client.get('/').status_code,
        'list': client.get('/agent/conversations'),
        'stream': client.post('/agent/conversations/c1/message/stream', json={'content': 'hi'}),
    }
result = {
    'available': web_agent_available(),
    'routes': sorted([method.upper(), path] for path, ops in paths.items() for method in ops),
    'root': responses['root'],
    'list': [responses['list'].status_code, responses['list'].json()],
    'stream': responses['stream'].status_code,
    'modules': sorted(m for m in sys.modules if m.startswith(('pantheon', 'uprobe.http.routers.agent'))),
}
print('RESULT ' + json.dumps(result), file=sys.__stdout__, flush=True)
'''


def run_python(code, tmp_path):
    env = dict(os.environ, UPROBE_SERVER_ROOT=str(tmp_path / 'server'))
    return subprocess.run([sys.executable, '-c', code], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=300)


def probe_result(process):
    lines = [line for line in process.stdout.splitlines() if line.startswith('RESULT ')]
    assert process.returncode == 0 and len(lines) == 1, process.stdout[-3000:] + process.stderr[-3000:]
    return json.loads(lines[0][len('RESULT '):])


def uprobe_cli(tmp_path, *args):
    """Run the `uprobe` console script (uprobe.core.cli:main) with PantheonOS blocked."""
    code = BLOCK_PANTHEON + textwrap.dedent(f'''
        sys.argv = ['uprobe', *{list(args)!r}]
        from uprobe.core.cli import main
        main()
    ''')
    return run_python(code, tmp_path)


def tiny_genome(tmp_path):
    rng = random.Random(7)
    chrom = ''.join(rng.choices('ACGT', k=3000))
    fasta = tmp_path / 'tiny.fa'
    fasta.write_text('>1\n' + ''.join(chrom[i:i + 80] + '\n' for i in range(0, len(chrom), 80)))
    attrs = 'gene_id "G1"; gene_name "G1"; transcript_id "G1.t1";'
    gtf = tmp_path / 'tiny.gtf'
    gtf.write_text(''.join(f'1\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}\n' for feature, start, end in (
        ('gene', 501, 1900), ('transcript', 501, 1900), ('exon', 501, 900), ('exon', 1501, 1900))))
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump({'tiny': {'fasta': str(fasta), 'gtf': str(gtf), 'align_index': []}}))
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({
        'name': 'no_agent', 'genome': 'tiny', 'targets': ['G1'], 'encoding': {},
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 40}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': {'target_gc': {'target': 'target_region', 'type': 'gc_content'}},
        'post_process': {'filters': {'target_gc': {'condition': 'target_gc >= 0'}}},
    }))
    return protocol, genomes


def test_the_blocker_hides_pantheon(tmp_path):
    process = run_python(BLOCK_PANTHEON + 'import pantheon\n', tmp_path)
    assert process.returncode != 0 and 'blocked by the test' in process.stderr


def test_import_uprobe_without_pantheon(tmp_path):
    code = BLOCK_PANTHEON + textwrap.dedent('''
        import pkgutil, importlib
        import uprobe, uprobe.core
        assert uprobe.UProbeAPI is uprobe.core.UProbeAPI
        for module in pkgutil.walk_packages(uprobe.core.__path__, 'uprobe.core.'):
            importlib.import_module(module.name)
        print(sorted(m for m in sys.modules if m.startswith(('pantheon', 'uprobe.http.agent', 'uprobe.http.routers'))))
    ''')
    process = run_python(code, tmp_path)
    assert process.returncode == 0, process.stderr[-3000:]
    assert process.stdout.strip() == '[]'


def test_uprobe_help_without_pantheon_lists_no_agent_command(tmp_path):
    process = uprobe_cli(tmp_path, '--help')
    assert process.returncode == 0, process.stderr[-3000:]
    commands = [line.split()[0] for line in process.stdout.split('Commands:')[1].splitlines() if line.strip()]
    assert {'run', 'validate', 'server'} <= set(commands) and 'agent' not in commands
    assert uprobe_cli(tmp_path, 'agent').returncode == 2


def test_uprobe_run_without_pantheon(tmp_path):
    protocol, genomes = tiny_genome(tmp_path)
    output = tmp_path / 'out'
    process = uprobe_cli(tmp_path, 'run', '-p', str(protocol), '-g', str(genomes), '-o', str(output),
                         '--threads', '1', '--json')
    assert process.returncode == 0, process.stderr[-3000:]
    summary = json.loads(process.stdout)
    assert summary['inputs']['valid_targets'] == ['G1'] and summary['counts']['final_probes'] > 0
    assert (output / 'summary.json').exists()


def test_server_starts_without_the_extra(tmp_path):
    process = run_python(BLOCK_PANTHEON + SERVER_PROBE, tmp_path)
    result = probe_result(process)
    assert result['available'] is False
    assert result['root'] == 200
    # Every other router is mounted; the stand-in for /agent stays out of the API docs.
    prefixes = {path.split('/')[1] for _, path in result['routes']}
    assert {'genome', 'workflow', 'task', 'auth', 'user', 'custom_probes'} <= prefixes
    assert 'agent' not in prefixes
    status, body = result['list']
    assert status == 503 and "uprobe[web-agent]" in body['detail']
    assert result['stream'] == 503
    assert result['modules'] == []
    warnings = [line for line in process.stdout.splitlines() if 'Web chat agent disabled' in line]
    assert len(warnings) == 1 and "uprobe[web-agent]" in warnings[0]


def test_chat_routes_are_mounted_with_the_extra(tmp_path):
    process = run_python(FAKE_PANTHEON + SERVER_PROBE, tmp_path)
    result = probe_result(process)
    assert result['available'] is True
    routes = {tuple(route) for route in result['routes']}
    assert {route for route in routes if route[1].startswith('/agent')} == AGENT_API
    # The real router answers: authentication first, not 503.
    assert result['list'][0] == 401
    assert 'uprobe.http.routers.agent' in result['modules']
    assert 'Web chat agent disabled' not in process.stdout


def core_imports(path):
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    package = '.'.join(path.relative_to(ROOT).with_suffix('').parts[:-1])
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.rsplit('.', node.level - 1)[0] if node.level > 1 else package
                module = f'{base}.{node.module}' if node.module else base
            else:
                module = node.module
            yield module
            yield from (f'{module}.{alias.name}' for alias in node.names)


def is_agent_module(name):
    agent = ('uprobe.http.agent', 'uprobe.http.routers.agent', 'uprobe.core.agent')
    return name.split('.')[0] == 'pantheon' or any(name == m or name.startswith(m + '.') for m in agent)


@pytest.mark.parametrize('path', sorted(CORE.rglob('*.py')), ids=lambda p: str(p.relative_to(CORE)))
def test_core_does_not_import_the_agent(path):
    assert [name for name in core_imports(path) if is_agent_module(name)] == []


def test_the_import_check_sees_agent_imports():
    for module in ('pantheon.chatroom', 'uprobe.http.agent', 'uprobe.http.agent.session_manager',
                   'uprobe.http.routers.agent'):
        assert is_agent_module(module)
    assert not is_agent_module('uprobe.http.agentx') and not is_agent_module('uprobe.core.api')


def test_no_agent_code_left_in_core_or_core_requirements():
    # A stale uprobe/core/agent/__pycache__ from an old checkout is no source.
    assert list((CORE / 'agent').rglob('*.py')) == []
    assert 'pantheon' not in (ROOT / 'requirements.txt').read_text().lower()
