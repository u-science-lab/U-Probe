"""Safe evaluation of probe part expressions (``expr``).

Expressions are never passed to ``eval``. They are parsed with ``ast`` and only
this small Python-like language is accepted:

* string and integer literals: ``'ACGT'``, ``13``, ``-6``
* names from the probe context (``target_region``, ``target``, ``encoding``,
  ``target_parts``) and names of other probes or parts, including dotted part
  names such as ``circle_probe.part2.barcode2``
* subscripts with an index, key or slice: ``target_region[26]``,
  ``target_region[-13:]``, ``encoding[target]['BC1']``
* ``+`` to concatenate strings (or add integers, e.g. in slice bounds)
* calls to whitelisted functions: ``rc(seq)`` (reverse complement)

Everything else (attribute access, other calls and operators, comprehensions,
lambdas, ...) is rejected with an error naming the offending construct.
"""
import ast
import functools
import typing as T

from ..utils import reverse_complement


class ExpressionError(ValueError):
    """Raised when a probe expression is invalid or cannot be evaluated."""


# Functions callable from expressions; each takes one sequence string.
FUNCTIONS = {'rc': reverse_complement}
# Names provided to every expression by UProbeAPI.construct_probes.
CONTEXT_NAMES = ('target_region', 'target', 'encoding', 'target_parts')

_CONSTRUCTS = {
    ast.Lambda: 'lambda expression',
    ast.ListComp: 'list comprehension',
    ast.SetComp: 'set comprehension',
    ast.DictComp: 'dict comprehension',
    ast.GeneratorExp: 'generator expression',
    ast.List: 'list literal',
    ast.Tuple: 'tuple',
    ast.Set: 'set literal',
    ast.Dict: 'dict literal',
    ast.Compare: 'comparison',
    ast.BoolOp: 'boolean operator',
    ast.IfExp: 'conditional expression',
    ast.JoinedStr: 'f-string',
    ast.NamedExpr: 'assignment expression',
    ast.Starred: 'starred expression',
    ast.Await: 'await expression',
    ast.Yield: 'yield expression',
    ast.YieldFrom: 'yield expression',
    ast.Slice: 'slice outside of a subscript',
}
_OPERATORS = {
    ast.Sub: '-', ast.Mult: '*', ast.Div: '/', ast.FloorDiv: '//', ast.Mod: '%',
    ast.Pow: '**', ast.MatMult: '@', ast.LShift: '<<', ast.RShift: '>>',
    ast.BitOr: '|', ast.BitXor: '^', ast.BitAnd: '&',
    ast.UAdd: 'unary +', ast.Not: 'not', ast.Invert: '~',
}
_ALLOWED_FUNCTIONS = ', '.join(f'{name}()' for name in FUNCTIONS)


def _snippet(node: ast.AST) -> str:
    text = ast.unparse(node)
    return repr(text if len(text) <= 60 else text[:57] + '...')


def _dotted_name(node: ast.AST) -> T.Optional[str]:
    """Return ``a.b.c`` for a name or a chain of attributes on a name."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return '.'.join(reversed(parts))


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


class _Checker:
    """Collect grammar violations and referenced names of a parsed expression."""

    def __init__(self):
        self.errors: T.List[str] = []
        self.names: T.Dict[str, None] = {}

    def visit(self, node: ast.AST):
        if isinstance(node, ast.Constant):
            if not (isinstance(node.value, str) or _is_int(node.value)):
                self.errors.append(f'literal {node.value!r} is not allowed (only string and integer literals)')
        elif isinstance(node, (ast.Name, ast.Attribute)):
            name = _dotted_name(node)
            if name is None:
                self.errors.append(f'attribute access is not allowed: {_snippet(node)}')
            elif any(part.startswith('__') for part in name.split('.')):
                self.errors.append(f'dunder name is not allowed: {name!r}')
            else:
                self.names[name] = None
        elif isinstance(node, ast.Subscript):
            self.visit(node.value)
            if isinstance(node.slice, ast.Slice):
                for bound in (node.slice.lower, node.slice.upper, node.slice.step):
                    if bound is not None:
                        self.visit(bound)
            else:
                self.visit(node.slice)
        elif isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Add):
                self.visit(node.left)
                self.visit(node.right)
            else:
                op = _OPERATORS.get(type(node.op), type(node.op).__name__)
                self.errors.append(f"operator '{op}' is not allowed (only '+'): {_snippet(node)}")
        elif isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                self.visit(node.operand)
            else:
                op = _OPERATORS.get(type(node.op), type(node.op).__name__)
                self.errors.append(f"operator '{op}' is not allowed: {_snippet(node)}")
        elif isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id in FUNCTIONS:
                if node.keywords or len(node.args) != 1 or isinstance(node.args[0], ast.Starred):
                    self.errors.append(f'{func.id}() takes exactly one positional argument: {_snippet(node)}')
                else:
                    self.visit(node.args[0])
            elif isinstance(func, ast.Name):
                self.errors.append(f"call to '{func.id}()' is not allowed (allowed functions: {_ALLOWED_FUNCTIONS})")
            else:
                self.errors.append(f'method or attribute call is not allowed: {_snippet(node)} '
                                   f'(allowed functions: {_ALLOWED_FUNCTIONS})')
        else:
            construct = _CONSTRUCTS.get(type(node), f'{type(node).__name__} expression')
            self.errors.append(f'{construct} is not allowed: {_snippet(node)}')


def _parse(expr) -> ast.expr:
    if not isinstance(expr, str):
        raise ExpressionError(f'expression must be a string, got {type(expr).__name__}')
    try:
        return ast.parse(expr.strip(), mode='eval').body
    except SyntaxError as exc:
        raise ExpressionError(f'invalid syntax in {expr!r}: {exc.msg}') from None
    except (ValueError, RecursionError, MemoryError):
        raise ExpressionError(f'expression is too long or too deeply nested: {expr[:60]!r}') from None


def _check(expr) -> T.Tuple[ast.expr, _Checker]:
    body = _parse(expr)
    checker = _Checker()
    try:
        checker.visit(body)
    except RecursionError:
        raise ExpressionError(f'expression is too deeply nested: {expr[:60]!r}') from None
    return body, checker


@functools.lru_cache(maxsize=1024)
def _compile(expr: str) -> T.Tuple[ast.expr, T.Tuple[str, ...]]:
    body, checker = _check(expr)
    if checker.errors:
        raise ExpressionError('; '.join(checker.errors))
    return body, tuple(checker.names)


def check_expression(expr, known_names: T.Optional[T.Iterable[str]] = None) -> T.List[str]:
    """Return a list of problems with ``expr`` (empty if it is valid).

    Grammar violations are always reported. If ``known_names`` is given, names
    that are not in it (context names, probe and part names) are reported too.
    """
    try:
        _, checker = _check(expr)
    except ExpressionError as exc:
        return [str(exc)]
    errors = list(checker.errors)
    if known_names is not None:
        known = set(known_names)
        for name in checker.names:
            if name in known:
                continue
            if '.' in name:
                errors.append(f'unknown name {name!r}: dotted names must refer to a probe part '
                              '(attribute access is not allowed)')
            else:
                errors.append(f'unknown name {name!r}')
    return errors


def expression_references(expr) -> T.Set[str]:
    """Return the names and dotted part names that ``expr`` reads.

    Called function names such as ``rc`` are not included.
    Raises ExpressionError if the expression is invalid.
    """
    if not isinstance(expr, str):
        _parse(expr)
    return set(_compile(expr)[1])


def evaluate_expression(expr, namespace: T.Mapping[str, T.Any]):
    """Evaluate ``expr``, looking up names and dotted part names in ``namespace``."""
    if not isinstance(expr, str):
        _parse(expr)
    body, _ = _compile(expr)
    try:
        return _evaluate(body, namespace)
    except RecursionError:
        raise ExpressionError(f'expression is too deeply nested: {expr[:60]!r}') from None


def _evaluate(node: ast.AST, namespace: T.Mapping[str, T.Any]):
    # Every branch re-checks the node shape, so nothing outside the grammar
    # is evaluated even if _compile's checks are bypassed.
    if isinstance(node, ast.Constant) and (isinstance(node.value, str) or _is_int(node.value)):
        return node.value
    if isinstance(node, (ast.Name, ast.Attribute)):
        name = _dotted_name(node)
        if name is not None and not name.startswith('__') and '.__' not in name:
            if name not in namespace:
                raise ExpressionError(f'unknown name {name!r}')
            return namespace[name]
    elif isinstance(node, ast.Subscript):
        value = _evaluate(node.value, namespace)
        if isinstance(node.slice, ast.Slice):
            bounds = [None if bound is None else _evaluate(bound, namespace)
                      for bound in (node.slice.lower, node.slice.upper, node.slice.step)]
            if not isinstance(value, (str, list, tuple)):
                raise ExpressionError(f'cannot slice {type(value).__name__}: {_snippet(node)}')
            if not all(bound is None or _is_int(bound) for bound in bounds):
                raise ExpressionError(f'slice bounds must be integers: {_snippet(node)}')
            try:
                return value[slice(*bounds)]
            except ValueError as exc:
                raise ExpressionError(f'{exc}: {_snippet(node)}') from None
        key = _evaluate(node.slice, namespace)
        if isinstance(value, dict):
            try:
                return value[key]
            except (KeyError, TypeError):
                raise ExpressionError(f'key {key!r} not found in {_snippet(node.value)}') from None
        if isinstance(value, (str, list, tuple)) and _is_int(key):
            try:
                return value[key]
            except IndexError:
                raise ExpressionError(f'index {key} out of range for {_snippet(node.value)} '
                                      f'(length {len(value)})') from None
        raise ExpressionError(f'cannot index {type(value).__name__} with {type(key).__name__}: {_snippet(node)}')
    elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _evaluate(node.left, namespace)
        right = _evaluate(node.right, namespace)
        if (isinstance(left, str) and isinstance(right, str)) or (_is_int(left) and _is_int(right)):
            return left + right
        raise ExpressionError(f'cannot add {type(left).__name__} and {type(right).__name__}: {_snippet(node)}')
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        operand = _evaluate(node.operand, namespace)
        if _is_int(operand):
            return -operand
        raise ExpressionError(f'cannot negate {type(operand).__name__}: {_snippet(node)}')
    elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in FUNCTIONS
          and len(node.args) == 1 and not node.keywords and not isinstance(node.args[0], ast.Starred)):
        arg = _evaluate(node.args[0], namespace)
        if not isinstance(arg, str):
            raise ExpressionError(f'{node.func.id}() expects a sequence string, got {type(arg).__name__}: '
                                  f'{_snippet(node)}')
        return FUNCTIONS[node.func.id](arg)
    raise ExpressionError(f'unsupported expression: {_snippet(node)}')


def _probe_nodes(probes) -> T.Iterator[T.Tuple[str, dict]]:
    """Yield ``(node_name, config)`` for every probe and part, named as in the DAG."""
    def walk(name, config):
        if not isinstance(config, dict):
            return
        yield name, config
        if 'template' in config and isinstance(config.get('parts'), dict):
            for part_name, part_config in config['parts'].items():
                yield from walk(f'{name}.{part_name}', part_config)

    if isinstance(probes, dict):
        for probe_name, probe_config in probes.items():
            yield from walk(str(probe_name), probe_config)


def check_probe_expressions(probes, check_names: bool = True) -> T.List[str]:
    """Check every ``expr`` in a protocol's ``probes`` mapping.

    Returns one message per problem, naming the probe part. With
    ``check_names``, names must be context names or probe/part names.
    """
    nodes = list(_probe_nodes(probes))
    known = set(CONTEXT_NAMES) | {name for name, _ in nodes} if check_names else None
    errors = []
    for name, config in nodes:
        if 'template' in config or 'expr' not in config:
            continue
        for error in check_expression(config['expr'], known):
            errors.append(f'probe part {name!r}: {error}')
    return errors
