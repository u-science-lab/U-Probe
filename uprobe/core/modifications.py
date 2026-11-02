"""Chemical annotations kept separate from nucleotide sequences and assembly."""

SUPPORTED_TYPES = frozenset(('phosphate', 'hydroxyl', 'fluorophore', 'disulfide', 'phosphorothioate'))
IUPAC = frozenset('ACGTURYSWKMBDHVN')


def _path_text(path):
    return ''.join(f'[{item}]' if isinstance(item, int) else ('.' if i else '') + str(item)
                   for i, item in enumerate(path))


class ModificationError(ValueError):
    def __init__(self, path, reason):
        self.path = tuple(path)
        self.reason = reason
        super().__init__(f'{_path_text(self.path)}: {reason}')


def normalize_modifications(value, probe_name, sequence_length=None):
    """Return explicit records; internal positions identify 1-based backbone bonds."""
    root = ('probes', probe_name, 'modifications')
    if not isinstance(value, dict):
        raise ModificationError(root, 'must be a mapping of 5prime, 3prime or internal')
    records = []
    occupied = set()
    for position, entries in value.items():
        path = root + (position,)
        if position not in ('5prime', '3prime', 'internal'):
            raise ModificationError(path, 'unknown position; use 5prime, 3prime or internal')
        if position == 'internal':
            if not isinstance(entries, list):
                raise ModificationError(path, 'must be a list of modification objects')
            items = [(item, path + (i,)) for i, item in enumerate(entries)]
        else:
            items = [(entries, path)]
        for item, item_path in items:
            if isinstance(item, str) and position != 'internal':
                item = {'type': item}
            if not isinstance(item, dict):
                raise ModificationError(item_path, 'must be a modification object with type')
            allowed = {'type', 'name'} | ({'after_base'} if position == 'internal' else set())
            for key in item:
                if key not in allowed:
                    raise ModificationError(item_path + (key,), 'unsupported modification field')
            kind = item.get('type')
            if not isinstance(kind, str) or kind not in SUPPORTED_TYPES:
                raise ModificationError(item_path + ('type',), 'unsupported type; use ' + ', '.join(sorted(SUPPORTED_TYPES)))
            if (position == 'internal') != (kind in ('disulfide', 'phosphorothioate')):
                raise ModificationError(item_path + ('type',), f'{kind} is not supported at {position}')
            record = {'probe': probe_name, 'type': kind, 'position': position}
            if kind == 'fluorophore':
                if not isinstance(item.get('name'), str) or not item['name'].strip():
                    raise ModificationError(item_path + ('name',), 'fluorophore requires a non-empty dye name')
                record['name'] = item['name']
            elif 'name' in item:
                raise ModificationError(item_path + ('name',), 'name applies only to fluorophore')
            if position == 'internal':
                after = item.get('after_base')
                if type(after) is not int or after < 1:
                    raise ModificationError(item_path + ('after_base',), 'must be a positive integer (1-based bond)')
                if sequence_length is not None and after >= sequence_length:
                    raise ModificationError(item_path + ('after_base',),
                                            f'must be less than sequence length {sequence_length}')
                record['after_base'] = after
                site = (position, after)
            else:
                if sequence_length is not None and sequence_length < 1:
                    raise ModificationError(item_path, 'cannot modify an empty sequence')
                site = (position,)
            if site in occupied:
                raise ModificationError(item_path, 'duplicate or conflicting modification at the same site')
            occupied.add(site)
            records.append(record)
    return records


def validate_probe_modifications(probes):
    """Validate actual probe annotations and reject declarations on sequence parts."""
    if not isinstance(probes, dict):
        return

    def walk(config, path, top=False):
        if not isinstance(config, dict):
            return
        if 'modifications' in config:
            if not top:
                raise ModificationError(path + ('modifications',), 'modifications are allowed only on actual probes')
            normalize_modifications(config['modifications'], path[1])
        parts = config.get('parts')
        if top and config.get('modifications') and isinstance(parts, dict) and 'modifications' in parts:
            raise ModificationError(path + ('parts', 'modifications'),
                                    'part name conflicts with the modification metadata column')
        if isinstance(parts, dict):
            for name, part in parts.items():
                walk(part, path + ('parts', name))

    for name, config in probes.items():
        walk(config, ('probes', name), top=True)


def validate_dna_sequence(sequence, path):
    """Reject chemistry notation rather than allowing it into DNA calculations."""
    if not isinstance(sequence, str) or set(sequence.upper()) - IUPAC:
        raise ModificationError(path, 'must contain nucleotide/IUPAC bases only; declare chemical modifications in modifications')



def check_modification_yaml_duplicates(text):
    """Reject annotation keys that a normal YAML mapping loader would discard."""
    import yaml
    root = yaml.compose(text, Loader=yaml.SafeLoader)

    def walk(node, path):
        if isinstance(node, yaml.MappingNode):
            seen = set()
            for key_node, value_node in node.value:
                if not isinstance(key_node, yaml.ScalarNode):
                    continue
                key = key_node.value
                child = path + (key,)
                if key in seen and 'modifications' in child:
                    raise ModificationError(child, 'duplicate YAML key')
                seen.add(key)
                walk(value_node, child)
        elif isinstance(node, yaml.SequenceNode):
            for i, item in enumerate(node.value):
                walk(item, path + (i,))

    walk(root, ())


def has_chemical_notation(sequence):
    """Recognize synthesis annotations without changing legacy non-DNA warnings."""
    import re
    return bool(re.search(r'-SS-|/[^/\s]+/|\*|phosphate|phosphorothioate|fluorophore|disulfide|'
                          r'\[[^\]]*(?:phos|fluor|disulf|Cy[0-9]|FAM|HEX|ROX|biotin)[^\]]*\]',
                          sequence, re.IGNORECASE))
