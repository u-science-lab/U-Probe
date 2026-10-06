import copy
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from uprobe.core.sampling import normalize_sampling, iter_target_windows
from uprobe.core.gen.fun import generate_target_seqs
from uprobe.core.utils import reverse_complement
from uprobe.core.attributes import add_attributes
from uprobe.core.api import UProbeAPI
from uprobe.core.process import filter_table


@pytest.mark.parametrize('length,step,overlap', [
    (0, 1, None), (True, 1, None), ([46, 40], 1, None),
    ([40], 1, None), ([40, 46, 48], 1, None), ([40, 46.0], 1, None),
    (40, 0, None), (40, -1, None), (40, True, None),
    (40, 1.5, None), (40, None, 40), ([40, 46], None, 39),
])
def test_invalid_sampling(length, step, overlap):
    with pytest.raises(ValueError):
        normalize_sampling(length, step, overlap)


def test_range_boundaries_and_step():
    windows = list(iter_target_windows('ACGTAC', [3, 5], 2))
    assert [(a, b) for a, b, _ in windows] == [(0, 3), (0, 4), (0, 5), (2, 5), (2, 6)]
    assert list(iter_target_windows('AC', [3, 5], 1)) == []
    assert list(iter_target_windows('ACG', [3, 5], 1)) == [(0, 3, 'ACG')]
    assert list(iter_target_windows('ACGTAC', 4, overlap=3)) == list(iter_target_windows('ACGTAC', 4, step=1))
    assert normalize_sampling(40, 1, 20) == (40, 40, 1)


@pytest.fixture
def reference(tmp_path):
    seq = 'ACGT' * 12
    fa = tmp_path / 'reference.fa'
    fa.write_text('>1\n' + seq + '\n')
    gtf = tmp_path / 'reference.gtf'
    gtf.write_text('1\ttest\texon\t1\t48\t.\t+\t.\tgene_name "G"; transcript_id "T";\n'
                   '1\ttest\tUTR\t1\t48\t.\t+\t.\tgene_name "G"; transcript_id "T";\n')
    return seq, str(fa), str(gtf)


@pytest.mark.parametrize('source', ['exon', 'CDS', 'UTR', 'genome', 'custom'])
def test_all_sources_range(reference, source):
    seq, fa, gtf = reference
    if source == 'custom':
        api = object.__new__(UProbeAPI)
        api.protocol = {'extracts': {'target_region': {'source': 'exon', 'length': [40, 46], 'step': 2}}}
        api._parse_targets = lambda: ([], {'G': seq})
        df = api.generate_target_seqs()
    else:
        df = generate_target_seqs(source, ['1:0-48'] if source == 'genome' else ['G'], fa, gtf, [40, 46], step=2)
    expected = list(iter_target_windows(seq, [40, 46], 2))
    assert df.target_region.tolist() == [s for _, _, s in expected]
    assert df.probe_id.is_unique
    if 'start' in df:
        assert df.start.tolist() == [a + 1 for a, _, _ in expected]
        assert df.end.tolist() == [b for _, b, _ in expected]


def test_length_attributes_and_filtering():
    df = pd.DataFrame({'target_region': ['A' * 40, 'C' * 46], 'pad.part1': ['ACG', 'ACGT']})
    config = {'attributes': {'target_length': {'target': 'target_region', 'type': 'length'},
                             'part_length': {'target': 'pad:part1', 'type': 'length'}}}
    df = add_attributes(df, config, {})
    assert df.target_length.tolist() == [40, 46]
    assert df.part_length.tolist() == [3, 4]
    result = filter_table(df, {'target_length': {'condition': 'target_length >= 42 & target_length <= 46'}})
    assert result.target_length.tolist() == [46]
    assert add_attributes(df.iloc[:0].copy(), config, {}).empty


def test_same_start_mapping_does_not_overwrite(monkeypatch, tmp_path):
    from uprobe.core import attributes
    monkeypatch.chdir(tmp_path)
    def align(_, records, *args):
        assert len(records) == 2
        return {key: len(seq) for key, seq in records.items()}
    monkeypatch.setattr(attributes, 'count_n_bowtie2_aligned_genes', align)
    df = pd.DataFrame({'target_region': ['A' * 40, 'A' * 46], 'exon_name': ['E', 'E'], 'start': [1, 1]})
    config = {'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie2'}}}
    result = add_attributes(df, config, {'fasta': str(tmp_path / 'ref.fa'), 'align_index': ['bowtie2']})
    assert result.hits.tolist() == [40, 46]


def test_range_workflow_exports_numeric_length(tmp_path, reference):
    seq, _, _ = reference
    api = object.__new__(UProbeAPI)
    api.protocol = {
        'name': 'range_test', 'genome': 'test', 'targets': [{'G': seq}],
        'encoding': {}, 'extracts': {'target_region': {'source': 'exon', 'length': [40, 46], 'step': 2}},
        'probes': {'pad': {'template': '{part1}', 'parts': {'part1': {'expr': 'rc(target_region)'}}}},
        'attributes': {'target_length': {'target': 'target_region', 'type': 'length'}},
        'post_process': {'filters': {'target_length': {'condition': 'target_length >= 42'}}},
        'summary': {'report_name': ''},
    }
    api.genome = {}
    api.output_dir = tmp_path
    api.build_genome_index = lambda **kwargs: None
    api.validate_targets = lambda **kwargs: True
    api._parse_targets = lambda: ([], {'G': seq})
    result = api.run_workflow(raw_csv=True)
    assert sorted(result.target_length.unique()) == [42, 43, 44, 45, 46]
    assert result.target_length.tolist() == result.target_region.str.len().tolist()
    assert len(list(tmp_path.glob('*.xlsx'))) == 2


@pytest.mark.asyncio
async def test_http_submit_range_and_reject_invalid(monkeypatch):
    import io
    import yaml
    from fastapi import HTTPException, UploadFile
    from uprobe.http.routers import workflow
    from types import SimpleNamespace
    saved = []
    monkeypatch.setattr(workflow, 'update_task_in_db', lambda user, task: saved.append(task))
    async def submit(length, step):
        payload = yaml.safe_dump({'extracts': {'target_region': {'source': 'exon', 'length': length, 'step': step}}}).encode()
        return await workflow.submit_task(UploadFile(filename='test.yaml', file=io.BytesIO(payload)), SimpleNamespace(username='test'))
    assert (await submit([40, 46], 1))['status'] == 'success'
    assert saved[0].parameters.probe_length == [40, 46]
    assert yaml.safe_load(saved[0].yaml_content)['extracts']['target_region']['step'] == 1
    with pytest.raises(HTTPException) as exc:
        await submit([46, 40], 1)
    assert exc.value.status_code == 422
    assert len(saved) == 1


def test_deployed_mip_seq_regression(reference):
    preset_path = Path(os.environ.get('UPROBE_PRESET', '/home/qzhang/uprobe_all/cathe/probe.json'))
    if not preset_path.exists():
        pytest.skip('deployed preset is tested on lc_zq')
    preset = json.loads(preset_path.read_text())['MiP-seq']
    assert preset['extracts']['target_region'] == {'source': 'exon', 'length': 40, 'step': 1}
    seq, fa, gtf = reference
    old = generate_target_seqs('exon', ['G'], fa, gtf, 40, overlap=39)
    new = generate_target_seqs('exon', ['G'], fa, gtf, 40, step=1)
    pd.testing.assert_frame_equal(old, new)
    api = object.__new__(UProbeAPI)
    api.protocol = copy.deepcopy(preset)
    api.protocol['encoding'] = {'G': {'BC1': 'ACGTAC', 'BC2': 'ATCGTA'}}
    probes = api.construct_probes(new)
    assert len(probes) == 9
    for i, target in enumerate(new.target_region):
        assert probes.iloc[i]['mRNA'] == target[:26] + target[27:40]
        assert probes.iloc[i]['pad_probe'] == reverse_complement(target[:13]) + 'ACGTAC' + 'A' + 'ATCGTA' + reverse_complement(target[13:26])
        assert probes.iloc[i]['amp_probe'] == reverse_complement(target[27:40]) + target[26] + reverse_complement('ATCGTA')[:-2]
    combined = pd.concat([new, probes], axis=1)
    # Real fold/Tm/self-match calculations; external alignment is covered separately.
    api.protocol['attributes'].pop('target_mapped_genes')
    result = add_attributes(combined, api.protocol, {})
    assert result[list(api.protocol['attributes'])].notna().all().all()
