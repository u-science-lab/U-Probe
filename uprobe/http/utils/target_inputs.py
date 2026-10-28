"""Normalize target identifiers without changing sequences or gene-name case."""
from copy import deepcopy
import unicodedata

def normalize_target_name(value):
    if not isinstance(value, str):
        raise ValueError('Target names must be strings.')
    def boundary(c):
        return c.isspace() or unicodedata.category(c) == 'Cf'
    start, end = 0, len(value)
    while start < end and boundary(value[start]):
        start += 1
    while end > start and boundary(value[end - 1]):
        end -= 1
    name = value[start:end]
    if not name:
        raise ValueError('Target name cannot be empty or contain only whitespace.')
    if any(c.isspace() or unicodedata.category(c).startswith('C') for c in name):
        raise ValueError(f'Target name {name!r} contains internal whitespace or invisible/control characters. Enter one identifier per target.')
    return name

def normalize_target_inputs(protocol):
    if not isinstance(protocol, dict):
        raise ValueError('Protocol must be a YAML mapping.')
    result = deepcopy(protocol)
    targets = result.get('targets')
    if not isinstance(targets, list) or not targets:
        raise ValueError('Targets must be a non-empty list.')
    seen = set()
    normalized = []
    for entry in targets:
        if isinstance(entry, str):
            name = normalize_target_name(entry)
            cleaned = name
        elif isinstance(entry, dict) and len(entry) == 1:
            raw, sequence = next(iter(entry.items()))
            name = normalize_target_name(raw)
            cleaned = {name: sequence}
        else:
            raise ValueError('Each target must be a name or a single name-to-sequence mapping.')
        if name in seen:
            raise ValueError(f'Duplicate target after whitespace cleanup: {name!r}.')
        seen.add(name)
        normalized.append(cleaned)
    result['targets'] = normalized
    def clean_mapping(mapping, field):
        if not isinstance(mapping, dict):
            raise ValueError(f'{field} must be a mapping.')
        cleaned = {}
        for raw, value in mapping.items():
            name = normalize_target_name(raw)
            if name in cleaned:
                raise ValueError(f'Duplicate target in {field} after whitespace cleanup: {name!r}.')
            cleaned[name] = value
        return cleaned
    if 'encoding' in result:
        result['encoding'] = clean_mapping(result['encoding'], 'encoding')
    post = result.get('post_process') or {}
    for field in ('avoid_otp', 'equal_space'):
        if field in post:
            post[field] = clean_mapping(post[field], field)
    return result
