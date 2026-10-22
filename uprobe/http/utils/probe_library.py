"""Built-in web/AI presets: <configured data_dir>/probes/<method>.yaml.

No cache, JSON fallback, metadata wrapper, or CLI-template merge.
"""
import logging

import yaml

from .paths import get_probes_dir

log = logging.getLogger(__name__)


class PresetLoader(yaml.SafeLoader):
    """Catch accidental duplicate keys instead of silently overwriting fields."""


def _mapping(loader, node, deep=False):
    loader.flatten_mapping(node)
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise ValueError(f'Duplicate YAML key {key!r} at line {key_node.start_mark.line + 1}')
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


PresetLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_builtin_probes():
    """Return the existing method -> config API shape, rereading files each call.

    Bad files are logged and skipped independently. Unfinished but well-formed
    methods are preserved; this is not scientific or full protocol validation.
    """
    root = get_probes_dir()
    if not root.is_dir():
        raise FileNotFoundError(f'Built-in probe directory does not exist: {root}')
    result = {}
    for file in sorted(root.glob('*.yaml')):
        try:
            config = yaml.load(file.read_text(encoding='utf-8'), Loader=PresetLoader)
            if not isinstance(config, dict) or not isinstance(config.get('probes'), dict):
                raise ValueError('Expected a configuration with a probes mapping')
            result[file.stem] = config
        except (OSError, ValueError, TypeError, yaml.YAMLError) as exc:
            log.error('Cannot load built-in probe %s: %s', file.name, exc)
    return result


def load_builtin_probe(name):
    """Resolve a selected method from the same fresh catalog used by the web API."""
    return load_builtin_probes()[name]
