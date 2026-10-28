"""The web presets in examples/web-presets are reference copies, not library templates."""
from pathlib import Path

import yaml

from uprobe.templates import template_names

PRESETS = Path(__file__).resolve().parent.parent / 'examples' / 'web-presets'


def test_web_presets_are_not_library_templates():
    names = template_names()
    assert not [name for name in names if name.endswith('_web')]
    assert {'mipseq', 'starmap'} <= set(names)


def test_web_preset_copies_parse_and_have_probes():
    files = sorted(PRESETS.glob('*.yaml'))
    assert [f.name for f in files] == ['MiP-seq_web.yaml', 'STARmap_web.yaml']
    for file in files:
        config = yaml.safe_load(file.read_text(encoding='utf-8'))
        assert isinstance(config.get('probes'), dict), file.name
