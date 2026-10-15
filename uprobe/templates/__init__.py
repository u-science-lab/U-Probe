"""Protocol templates for published probe technologies.

Each ``<name>.yaml`` in this package is a complete U-Probe protocol. Its
top-level ``meta`` block documents the probe architecture, the cited sources
of every parameter, the hybridization conditions of the source protocol and
what the user must replace. The workflow ignores ``meta``.
"""
from importlib import resources

import yaml


def _template_files():
    files = (f for f in resources.files(__name__).iterdir() if f.name.endswith('.yaml'))
    return {f.name[:-len('.yaml')]: f for f in sorted(files, key=lambda f: f.name)}


def template_names():
    return list(_template_files())


def read_template(name):
    """Return the template YAML text unchanged, comments included."""
    files = _template_files()
    if name not in files:
        raise KeyError(f"Unknown template {name!r}. Available: {', '.join(files)}")
    return files[name].read_text(encoding='utf-8')


def load_template(name):
    return yaml.safe_load(read_template(name))


def format_citation(citation):
    text = citation.get('short', '')
    return f"{text}, doi:{citation['doi']}" if citation.get('doi') else text


def list_templates():
    """Return name, technology and primary citation of every template."""
    rows = []
    for name in template_names():
        meta = load_template(name).get('meta') or {}
        citations = meta.get('citations') or [{}]
        rows.append({
            'name': name,
            'technology': meta.get('technology', ''),
            'citation': format_citation(citations[0]),
        })
    return rows
