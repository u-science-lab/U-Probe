"""The probe presets and protocol skeleton the built-in agent copies its RNA source from."""
import json
from pathlib import Path

import pytest
import yaml

from uprobe.core.validate import SOURCES, validate_protocol

ROOT = Path(__file__).resolve().parents[1]
PRESETS = json.loads((ROOT / 'examples' / 'probe.template.json').read_text())
SOURCE = 'extracts.target_region.source'


@pytest.mark.parametrize('name', sorted(PRESETS))
def test_example_presets_name_the_region_they_extract(name):
    # source: exon is the CDS of protein-coding genes and draws a validate warning; the example
    # preset is a generic probe, not a whole-mRNA (smFISH-style) design, so it says CDS.
    preset = PRESETS[name]
    source = preset['extracts']['target_region']['source']
    assert source in SOURCES and source != 'exon'
    protocol = {'name': name, 'genome': 'hg38', 'targets': ['GAPDH'],
                'encoding': {'GAPDH': {key: 'ACGTACGTACGTACGTACGT' for key in preset['barcodes']['barcodes']}},
                **{key: preset[key] for key in ('extracts', 'probes', 'attributes', 'post_process') if key in preset}}
    findings = validate_protocol(protocol)
    assert [f for f in findings if f.path == SOURCE] == []
    assert [f for f in findings if f.severity == 'error'] == []


def test_agent_protocol_skeleton_offers_no_exon_default():
    skeleton = yaml.safe_load((ROOT / 'uprobe' / 'http' / 'agent' / 'templates' / 'DEFAULT_PROTOCOL.yaml').read_text())
    assert skeleton['extracts']['target_region']['source'] is None
