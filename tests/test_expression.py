import copy
import io
import json
import random
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from uprobe.core.gen import probe as probe_module
from uprobe.core.gen.expression import (
    CONTEXT_NAMES, ExpressionError, check_expression, check_probe_expressions,
    evaluate_expression, expression_references,
)
from uprobe.core.gen.probe import construct_probes
from uprobe.core.gen.utils import parse_expression
from uprobe.core.utils import reverse_complement

HERE = Path(__file__).parent
ROOT = HERE.parent


def legacy_build(self, context):
    """ExprProbe.build before the safe evaluator (eval-based), used as reference."""
    eval_globals = {"rc": reverse_complement}
    eval_locals = context.copy()
    dep_name_map = {}
    for dep in self.deps:
        if not dep.done:
            raise RuntimeError(f"dependency {dep.name} was not built before {self.name}")
        if dep.result is None:
            raise RuntimeError(f"dependency {dep.name} has no result in memory")
        safe_dep_name = dep.name.replace('.', '_')
        dep_name_map[dep.name] = safe_dep_name
        eval_locals[safe_dep_name] = str(dep.result)
    modified_expr = self.expr
    for original_name in sorted(dep_name_map.keys(), key=len, reverse=True):
        modified_expr = modified_expr.replace(original_name, dep_name_map[original_name])
    self.result = str(eval(modified_expr, eval_globals, eval_locals))
    self.done = True


def random_seq(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def make_contexts(protocol, n=5, length=70, layout=None, seed=0):
    """Contexts shaped like UProbeAPI.construct_probes builds them."""
    rng = random.Random(seed)
    targets = []
    for target in protocol.get('targets') or ['G']:
        targets.append(next(iter(target)) if isinstance(target, dict) else target)
    contexts = []
    for i in range(n):
        target_region = random_seq(rng, length)
        target_parts = {}
        if layout:
            offset = 0
            for name, size in layout.items():
                target_parts[name] = target_region[offset:offset + size]
                offset += size
        contexts.append({
            'target_region': target_region,
            'target': targets[i % len(targets)],
            'encoding': protocol.get('encoding') or {},
            'target_parts': target_parts,
        })
    return contexts


# Probe structures from the README example and a MiP-seq-like structure whose
# expected sequences are asserted in tests/test_sampling.py.
README_PROTOCOL = {
    'targets': ['CD4'],
    'encoding': {'CD4': {'BC1': 'ACGAGCCTTCCA', 'BC2': 'CGGTAATGGACT'}},
    'probes': {
        'probe_1': {'template': '{part1}{part2}', 'parts': {
            'part1': {'expr': 'rc(target_region[0:20])'},
            'part2': {'template': 'CC{barcode1}TGCGTCTATTT{barcode2}TAGTGGAGCCT', 'parts': {
                'barcode1': {'expr': "encoding[target]['BC1']"},
                'barcode2': {'expr': "encoding[target]['BC2']"}}}}},
        'probe_2': {'template': '{part1}AGGCTCCACTA', 'parts': {
            'part1': {'expr': 'rc(target_region[-10:])'}}},
    },
}
MIP_LIKE_PROTOCOL = {
    'targets': ['G'],
    'encoding': {'G': {'BC1': 'ACGTAC', 'BC2': 'ATCGTA'}},
    'probes': {
        'mRNA': {'template': '{part1}{part2}', 'parts': {
            'part1': {'expr': 'target_region[:26]'},
            'part2': {'expr': 'target_region[27:40]'}}},
        'pad_probe': {'template': '{part1}{part2}{part3}', 'parts': {
            'part1': {'expr': 'rc(target_region[0:13])'},
            'part2': {'template': '{barcode1}A{barcode2}', 'parts': {
                'barcode1': {'expr': "encoding[target]['BC1']"},
                'barcode2': {'expr': "encoding[target]['BC2']"}}},
            'part3': {'expr': 'rc(target_region[13:26])'}}},
        'amp_probe': {'template': '{part1}{part2}{part3}', 'parts': {
            'part1': {'expr': 'rc(target_region[27:40])'},
            'part2': {'expr': 'target_region[26]'},
            'part3': {'expr': 'rc(pad_probe.part2.barcode2)[:-2]'}}},
    },
}
LAYOUT_PROBES = {'pad': {'template': '{left}{right}', 'parts': {
    'left': {'expr': "rc(target_parts['left'])"},
    'right': {'expr': "rc(target_parts['right'][-6:])"}}}}


def load_protocols():
    protocols = {}
    for name in ['DNA_format.yaml', 'RNA_format.yaml', 'double_hyb_rca.yaml', 'probes.yaml']:
        protocols[name] = (yaml.safe_load((HERE / 'data' / name).read_text()), None)
    example = json.loads((ROOT / 'examples' / 'probe.template.json').read_text())['ExampleProbe']
    example.update(targets=['G'], encoding={'G': {'BC1': 'ACGTACGTAC'}})
    protocols['probe.template.json'] = (example, None)
    protocols['README'] = (README_PROTOCOL, None)
    protocols['MiP-seq-like'] = (MIP_LIKE_PROTOCOL, None)
    protocols['layout'] = ({'probes': LAYOUT_PROBES}, {'left': 20, 'gap': 2, 'right': 18})
    protocols['sampling'] = ({'probes': {'pad': {'template': '{part1}', 'parts': {
        'part1': {'expr': 'rc(target_region)'}}}}}, None)
    return protocols


PROTOCOLS = load_protocols()


@pytest.mark.parametrize('name', list(PROTOCOLS))
def test_existing_protocols_build_identical_probes(monkeypatch, name):
    protocol, layout = PROTOCOLS[name]
    contexts = make_contexts(protocol, layout=layout)
    new = construct_probes(copy.deepcopy(protocol), copy.deepcopy(contexts))
    monkeypatch.setattr(probe_module.ExprProbe, 'build', legacy_build)
    old = construct_probes(copy.deepcopy(protocol), copy.deepcopy(contexts))
    pd.testing.assert_frame_equal(new, old)
    assert (new != '').any().any()


def test_mip_like_sequences():
    contexts = make_contexts(MIP_LIKE_PROTOCOL, length=40)
    probes = construct_probes(MIP_LIKE_PROTOCOL, contexts)
    for i, context in enumerate(contexts):
        target = context['target_region']
        assert probes.iloc[i]['mRNA'] == target[:26] + target[27:40]
        assert probes.iloc[i]['pad_probe'] == reverse_complement(target[:13]) + 'ACGTAC' + 'A' + 'ATCGTA' + reverse_complement(target[13:26])
        assert probes.iloc[i]['amp_probe'] == reverse_complement(target[27:40]) + target[26] + reverse_complement('ATCGTA')[:-2]


SURVEY_CONTEXT = {
    'target_region': 'ATGAAGGCCTGCCGGTTATGAAGGCCTGCCGGTTATGAAG',
    'target': 'G',
    'gene_id': 'G',
    'encoding': {'G': {'BC1': 'ACT', 'BC2': 'GGA', 'BC3': 'TTC', 'barcode1': 'AAAT', 'barcode2': 'ATTTG'}},
    'target_parts': {'left': 'ATGAAGGCCTGC', 'right': 'CGGTTATGAAGG'},
    'probe_1': 'ACGTACGTACGTAC',
    'circle_probe.part2.barcode2': 'ATTTATTGCGATTTGGC',
    'pad_probe.part2.barcode2': 'ATCGTA',
}
# Every expression found in the repository (protocols, examples, README, agent templates).
SURVEY = [
    "rc(target_region[0:20])", "rc(target_region[-10:])", "encoding[target]['BC1']",
    "encoding[target]['BC2']", "target_region[:]", "rc(target_region[0:13])",
    "rc(target_region[13:26])", "rc(target_region[-13:])", "target_region[26]",
    "rc(target_region[0:11])", "target_region[0:15]", "probe_1[0:10]",
    "encoding[gene_id]['BC2']", "encoding[gene_id]['BC3']", "encoding[gene_id]['barcode1']",
    "encoding[gene_id]['barcode2']", "rc(circle_probe.part2.barcode2)",
    "rc(pad_probe.part2.barcode2)[:-2]", "rc(target_region)", "rc(target_parts['left'])",
    "rc(target_parts['right'][-6:])", "'ACTACTACTA'", "'CATCATCATC'",
]


@pytest.mark.parametrize('expr', SURVEY)
def test_survey_expressions_match_eval(expr):
    namespace = dict(SURVEY_CONTEXT)
    legacy = expr
    for name in sorted((n for n in namespace if '.' in n), key=len, reverse=True):
        legacy = legacy.replace(name, name.replace('.', '_'))
        namespace[name.replace('.', '_')] = namespace[name]
    assert check_expression(expr, SURVEY_CONTEXT) == []
    assert evaluate_expression(expr, SURVEY_CONTEXT) == eval(legacy, {'rc': reverse_complement}, namespace)


def test_concatenation_and_index_arithmetic():
    ns = {'target_region': 'AACCGGTT', 'p.a': 'GG'}
    assert evaluate_expression("'TT' + rc(target_region[0:2]) + p.a", ns) == 'TTTTGG'
    assert evaluate_expression('target_region[1 + 1:-1]', ns) == 'CCGGT'
    with pytest.raises(ExpressionError, match='cannot add str and int'):
        evaluate_expression("target_region + 1", ns)


def test_references_exclude_functions():
    assert parse_expression('rc(circle_probe.part2.barcode2)[:-2]') == {'circle_probe.part2.barcode2'}
    assert expression_references("encoding[target]['BC1']") == {'encoding', 'target'}
    with pytest.raises(ExpressionError, match="open"):
        parse_expression("open('x')")


@pytest.mark.parametrize('expr,message', [
    ("__import__('os')", "call to '__import__()' is not allowed"),
    ("__import__('os').system('true')", "method or attribute call is not allowed"),
    ("().__class__.__bases__", "attribute access is not allowed"),
    ("open('x')", "call to 'open()' is not allowed"),
    ("lambda: 1", "lambda expression is not allowed"),
    ("[c for c in target_region]", "list comprehension is not allowed"),
    ("getattr(rc, 'x')", "call to 'getattr()' is not allowed"),
    ("'abc'.upper", "attribute access is not allowed"),
    ("'abc'.upper()", "method or attribute call is not allowed"),
    ("target_region.upper()", "method or attribute call is not allowed"),
    ("target_region.__class__", "dunder name is not allowed"),
    ("__builtins__", "dunder name is not allowed"),
    ("target_region * 1000", "operator '*' is not allowed"),
    ("target_region if target else ''", "conditional expression is not allowed"),
    ("f'{target}'", "f-string is not allowed"),
    ("(x := 1)", "assignment expression is not allowed"),
    ("{}", "dict literal is not allowed"),
    ("1.5", "literal 1.5 is not allowed"),
    ("rc(target_region, 1)", "rc() takes exactly one positional argument"),
    ("rc(seq=target_region)", "rc() takes exactly one positional argument"),
    ("import os", "invalid syntax"),
    ("target_region; open('x')", "invalid syntax"),
    ("a+" * 100000 + "a", "too long or too deeply nested"),
])
def test_rejected_expressions(expr, message):
    errors = check_expression(expr)
    assert errors and message in errors[0]
    with pytest.raises(ExpressionError, match='.') as exc:
        evaluate_expression(expr, {'target_region': 'ACGT', 'target': 'G'})
    assert message in str(exc.value)


def test_rejected_expression_is_never_executed(tmp_path):
    marker = tmp_path / 'pwned'
    expr = f"open({str(marker)!r}, 'w')"
    protocol = {'probes': {'p': {'template': '{a}', 'parts': {'a': {'expr': expr}}}}}
    with pytest.raises(ExpressionError, match="probe part 'p.a': call to 'open\\(\\)'"):
        construct_probes(protocol, make_contexts(protocol))
    assert not marker.exists()


def test_unknown_names():
    assert check_expression('gene_id', CONTEXT_NAMES) == ["unknown name 'gene_id'"]
    assert 'attribute access is not allowed' in check_expression('target_region.upper', CONTEXT_NAMES)[0]
    assert check_expression('target_region.upper') == []  # grammar only
    with pytest.raises(ExpressionError, match="unknown name 'target_region.upper'"):
        evaluate_expression('target_region.upper', {'target_region': 'ACGT'})


def test_check_probe_expressions_names_parts():
    probes = copy.deepcopy(MIP_LIKE_PROTOCOL['probes'])
    assert check_probe_expressions(probes) == []
    probes['amp_probe']['parts']['part3']['expr'] = 'rc(pad_probe.part9)'
    probes['mRNA']['parts']['part1']['expr'] = 'target_region.lower()'
    errors = check_probe_expressions(probes)
    assert len(errors) == 2
    assert errors[0].startswith("probe part 'mRNA.part1': method or attribute call")
    assert errors[1].startswith("probe part 'amp_probe.part3': unknown name 'pad_probe.part9'")
    assert check_probe_expressions(probes, check_names=False) == [errors[0]]
    assert check_probe_expressions(None) == []
    probes = PROTOCOLS['double_hyb_rca.yaml'][0]['probes']
    assert check_probe_expressions(probes) == []
    legacy = yaml.safe_load(yaml.safe_dump(probes).replace('encoding[target]', 'encoding[gene_id]'))
    errors = check_probe_expressions(legacy)
    assert len(errors) == 2 and all("unknown name 'gene_id'" in e for e in errors)


def test_runtime_errors_are_clear():
    ns = {'target_region': 'ACGT', 'target': 'G', 'encoding': {'G': {'BC1': 'AA'}}}
    with pytest.raises(ExpressionError, match=r"key 'BC9' not found in 'encoding\[target\]'"):
        evaluate_expression("encoding[target]['BC9']", ns)
    with pytest.raises(ExpressionError, match='index 10 out of range'):
        evaluate_expression('target_region[10]', ns)
    with pytest.raises(ExpressionError, match=r'rc\(\) expects a sequence string, got dict'):
        evaluate_expression('rc(encoding)', ns)
    with pytest.raises(ExpressionError, match='slice bounds must be integers'):
        evaluate_expression("target_region['a':]", ns)
    with pytest.raises(ExpressionError, match='expression must be a string'):
        evaluate_expression(5, ns)


@pytest.mark.asyncio
async def test_http_submit_rejects_unsafe_expressions(monkeypatch):
    from fastapi import HTTPException, UploadFile
    from uprobe.http.routers import workflow
    saved = []
    monkeypatch.setattr(workflow, 'update_task_in_db', lambda user, task: saved.append(task))

    async def submit(protocol):
        payload = yaml.safe_dump(protocol).encode()
        return await workflow.submit_task(UploadFile(filename='test.yaml', file=io.BytesIO(payload)),
                                          SimpleNamespace(username='test'))

    protocol = yaml.safe_load((HERE / 'data' / 'RNA_format.yaml').read_text())
    assert (await submit(protocol))['status'] == 'success'
    for expr, message in [("__import__('os').system('true')", 'method or attribute call is not allowed'),
                          ("encoding[gene_id]['BC1']", "unknown name 'gene_id'"),
                          ("rc(target_region", 'invalid syntax')]:
        bad = copy.deepcopy(protocol)
        bad['probes']['probe_2']['parts']['part1']['expr'] = expr
        with pytest.raises(HTTPException) as exc:
            await submit(bad)
        assert exc.value.status_code == 422
        assert "probe part 'probe_2.part1'" in exc.value.detail and message in exc.value.detail
    assert len(saved) == 1
