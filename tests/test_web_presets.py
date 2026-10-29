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
    assert [f.name for f in files] == ['ExSeq_web.yaml', 'MERFISH_web.yaml', 'MiP-seq_web.yaml', 'Oligo-LiveFISH_web.yaml', 'OpenFISH_web.yaml', 'PRISM_web.yaml', 'RAEFISH_web.yaml', 'STARmap_web.yaml', 'TDDN-FISH_web.yaml', 'seqFISH+_web.yaml']
    for file in files:
        config = yaml.safe_load(file.read_text(encoding='utf-8'))
        assert isinstance(config.get('probes'), dict), file.name
