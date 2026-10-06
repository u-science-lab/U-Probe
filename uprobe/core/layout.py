"""Named, contiguous target parts with variable lengths."""
import re
import ast
from string import Formatter

import pandas as pd

from .sampling import normalize_sampling


def validate_layout_references(protocol):
    layout = protocol.get('extracts', {}).get('target_region', {}).get('layout')
    names = set(layout.get('parts', {})) if isinstance(layout, dict) else set()
    if layout is not None and 'target_parts' in protocol.get('probes', {}):
        raise ValueError('target_parts is reserved for layout references')

    def walk(value):
        if isinstance(value, dict):
            if isinstance(value.get('expr'), str):
                for node in ast.walk(ast.parse(value['expr'], mode='eval')):
                    if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and node.value.id == 'target_parts':
                        key = node.slice
                        if isinstance(key, ast.Constant) and (not isinstance(key.value, str) or key.value not in names):
                            raise ValueError(f'Probe expression references missing target part {key.value!r}')
                        if not layout:
                            raise ValueError('target_parts references require extracts.target_region.layout')
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    walk(protocol.get('probes', {}))
    for attr in protocol.get('attributes', {}).values():
        target = attr.get('target', '').replace(':', '.')
        if target.startswith('target_parts.') and target[len('target_parts.'):] not in names:
            raise ValueError(f'Attribute references missing target part {target}')


def compile_layout(layout, length):
    lower, upper, _ = normalize_sampling(length)
    if not isinstance(layout, dict) or not isinstance(layout.get('parts'), dict) or not layout['parts']:
        raise ValueError('layout must define a template and non-empty parts mapping')
    template = layout.get('template')
    if not isinstance(template, str):
        raise ValueError('layout.template must be a string')
    names = []
    try:
        for literal, name, spec, conversion in Formatter().parse(template):
            if literal or name is None or spec or conversion or not re.fullmatch(r'[A-Za-z_]\w*', name):
                raise ValueError('layout.template must contain only named parts, e.g. {left}{gap}{right}')
            names.append(name)
    except ValueError as exc:
        raise ValueError(f'Invalid layout.template: {exc}') from exc
    if len(set(names)) != len(names) or set(names) != set(layout['parts']):
        raise ValueError('layout.template must reference every part exactly once')
    bounds = []
    for name in names:
        part = layout['parts'][name]
        value = part.get('length') if isinstance(part, dict) else None
        if type(value) is int:
            a = b = value
        elif isinstance(value, (list, tuple)) and len(value) == 2:
            a, b = value
        else:
            raise ValueError(f'layout part {name}: length must be an integer or [min, max]')
        if any(type(v) is not int or v < 0 for v in (a, b)) or a > b:
            raise ValueError(f'layout part {name}: lengths must be non-negative integers with min <= max')
        bounds.append((a, b))
    suffix_min = [sum(a for a, _ in bounds[i:]) for i in range(len(bounds) + 1)]
    suffix_max = [sum(b for _, b in bounds[i:]) for i in range(len(bounds) + 1)]
    combinations = {}

    def visit(index, total, sizes):
        if index == len(names):
            combinations.setdefault(total, []).append(tuple(sizes))
            return
        a, b = bounds[index]
        first = max(a, lower - total - suffix_max[index + 1])
        last = min(b, upper - total - suffix_min[index + 1])
        for size in range(first, last + 1):
            visit(index + 1, total + size, [*sizes, size])

    visit(0, 0, [])
    if not combinations:
        raise ValueError('layout part lengths cannot satisfy target_region.length')
    return names, combinations


def expand_layout_candidates(df, compiled):
    """Expand each extracted sequence into its independent, named layouts."""
    names, combinations = compiled
    records = []
    for row in df.to_dict('records'):
        sequence = row['target_region']
        for index, sizes in enumerate(combinations.get(len(sequence), []), 1):
            candidate = row.copy()
            candidate['probe_id'] = f"{row['probe_id']}__layout{index}"
            candidate['layout_id'] = '_'.join(f'{name}:{size}' for name, size in zip(names, sizes))
            offset = 0
            for name, size in zip(names, sizes):
                column = f'target_parts.{name}'
                candidate[column] = sequence[offset:offset + size]
                candidate[column + '.length'] = size
                # Export 1-based inclusive coordinates; empty parts have end = start - 1.
                candidate[column + '.start'] = offset + 1
                candidate[column + '.end'] = offset + size
                offset += size
            records.append(candidate)
    return pd.DataFrame(records)
