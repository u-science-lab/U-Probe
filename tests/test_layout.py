import copy
import io
import itertools

import pandas as pd
import pytest
import yaml

from uprobe.core.api import UProbeAPI
from uprobe.core.layout import compile_layout, expand_layout_candidates, validate_layout_references
from uprobe.core.attributes import add_attributes
from uprobe.core.process import filter_table
from uprobe.core.utils import reverse_complement


def layout(**lengths):
    return {'template': ''.join('{' + name + '}' for name in lengths),
            'parts': {name: {'length': length} for name, length in lengths.items()}}


def test_combinations_match_brute_force():
    names, groups = compile_layout(layout(a=[12, 14], b=[12, 14], c=[12, 14]), [40, 46])
    assert names == ['a', 'b', 'c']
    actual = {sizes for options in groups.values() for sizes in options}
    assert actual == {sizes for sizes in itertools.product(range(12,15), repeat=3) if 40 <= sum(sizes) <= 46}
    assert set(groups) == {40, 41, 42}


@pytest.mark.parametrize('bad,total', [
    ({'template': '{a}{a}', 'parts': {'a': {'length': 20}}}, 40),
    ({'template': '{a}A', 'parts': {'a': {'length': 40}}}, 40),
    ({'template': '{a}', 'parts': {'a': {'length': 40}, 'b': {'length': 0}}}, 40),
    (layout(a=[25,19]), 40), (layout(a=-1), 40), (layout(a=True), 40),
    (layout(a=[0,1.5]), 40), (layout(a=0), 40), (layout(a=20,b=20), [41,46]),
    ({'template': '{a', 'parts': {'a': {'length': 40}}}, 40),
])
def test_invalid_layout(bad,total):
    with pytest.raises(ValueError): compile_layout(bad,total)


def test_order_comes_from_template_and_zero_gap():
    spec = layout(left=19, gap=0, right=21)
    spec['parts'] = dict(reversed(list(spec['parts'].items())))
    compiled = compile_layout(spec, 40)
    df = expand_layout_candidates(pd.DataFrame({'probe_id':['P'], 'target_region':['A'*19+'C'*21]}), compiled)
    assert df.iloc[0]['target_parts.left'] == 'A'*19
    assert df.iloc[0]['target_parts.gap'] == ''
    assert df.iloc[0]['target_parts.right'] == 'C'*21
    for name, start, end in [('left', 1, 19), ('gap', 20, 19), ('right', 20, 40)]:
        row = df.iloc[0]
        assert row[f'target_parts.{name}.start'] == start
        assert row[f'target_parts.{name}.end'] == end
        assert row[f'target_parts.{name}.length'] == end - start + 1
        assert row[f'target_parts.{name}'] == row['target_region'][start - 1:end]


def make_api(spec, seq='ACGT'*12):
    api = object.__new__(UProbeAPI)
    api.protocol = {
        'name':'layout_test', 'genome':'test', 'targets':[{'G':seq}], 'encoding':{},
        'extracts':{'target_region':{'source':'exon','length':[40,46],'step':2,'layout':spec}},
        'probes':{'pad':{'template':'{left}{right}', 'parts':{
            'left':{'expr':"rc(target_parts['left'])"},
            'right':{'expr':"rc(target_parts['right'][-6:])"}}}},
        'attributes':{'left_length':{'target':'target_parts.left','type':'length'},
                      'left_tm':{'target':'pad.left','type':'annealing_temperature'}},
        'post_process':{'filters':{'left_length':{'condition':'left_length >= 21'}}},
        'summary':{'report_name':''},
    }
    api._parse_targets = lambda: ([], {'G':seq})
    return api


def test_dynamic_references_and_filters():
    api = make_api(layout(left=[19,25],gap=[0,2],right=[19,25]))
    targets = api.generate_target_seqs()
    assert targets.probe_id.is_unique
    assert targets.groupby(['start','end']).size().max() > 1
    assert targets['target_parts.gap.length'].min() == 0
    probes = api.construct_probes(targets)
    for i,row in targets.iterrows():
        assert probes.iloc[i]['pad.left'] == reverse_complement(row['target_parts.left'])
        assert probes.iloc[i]['pad.right'] == reverse_complement(row['target_parts.right'][-6:])
    combined = add_attributes(pd.concat([targets,probes],axis=1),api.protocol,{})
    filtered = filter_table(combined, api.protocol['post_process']['filters'])
    assert not filtered.empty and filtered.left_length.min() >= 21
    assert combined.left_tm.notna().all()


def test_layout_full_workflow_and_report(tmp_path):
    api = make_api(layout(left=[19,25],gap=[0,2],right=[19,25]))
    api.genome = {}
    api.output_dir = tmp_path
    api.build_genome_index = lambda **kwargs: None
    api.validate_targets = lambda **kwargs: True
    result = api.run_workflow(raw_csv=True)
    assert not result.empty
    assert len(list(tmp_path.glob('*.xlsx'))) == 2
    api.protocol['summary'] = {'report_name':'rna_report'}
    report = api.generate_report(result)
    text = report['html_reports'][0].read_text()
    assert 'target_parts.left' in text and 'layout_id' in text


def test_missing_reference_fails_before_generation():
    api = make_api(layout(left=20,right=20))
    api.protocol['probes']['pad']['parts']['left']['expr'] = "rc(target_parts['missing'])"
    with pytest.raises(ValueError, match='missing target part'): api.generate_target_seqs()
    with pytest.raises(ValueError): validate_layout_references({'probes':api.protocol['probes']})


def test_attribute_calculation_cached(monkeypatch):
    from uprobe.core import attributes
    calls=[]
    monkeypatch.setattr(attributes,'cal_temp',lambda seq: calls.append(seq) or len(seq))
    df=pd.DataFrame({'s':['ACGT','ACGT','ACG']})
    result=add_attributes(df,{'attributes':{'tm':{'target':'s','type':'annealing_temperature'}}},{})
    assert result.tm.tolist()==[4,4,3] and calls==['ACGT','ACG']


@pytest.mark.parametrize('source', ['exon','CDS','UTR','genome'])
def test_layout_all_genome_sources(tmp_path, source):
    seq='ACGT'*12
    fa=tmp_path/'ref.fa';fa.write_text('>1\n'+seq+'\n')
    gtf=tmp_path/'ref.gtf';gtf.write_text('1\ttest\texon\t1\t48\t.\t-\t.\tgene_name "G"; transcript_id "T";\n'
                                        '1\ttest\tUTR\t1\t48\t.\t-\t.\tgene_name "G"; transcript_id "T";\n')
    api=make_api(layout(left=[19,25],gap=[0,2],right=[19,25]))
    api.genome={'fasta':str(fa),'gtf':str(gtf)}
    api.protocol['extracts']['target_region']['source']=source
    api._parse_targets=lambda: (['1:0-48'] if source=='genome' else ['G'],{})
    result=api.generate_target_seqs()
    assert not result.empty and result.probe_id.is_unique
    for _,row in result.iterrows():
        assert row['target_region']==row['target_parts.left']+row['target_parts.gap']+row['target_parts.right']


def test_alignment_reused_across_layouts(monkeypatch,tmp_path):
    from uprobe.core import attributes
    monkeypatch.chdir(tmp_path)
    def align(_,records,*args):
        assert records=={'0':'ACGT','1':'ACG'}
        return {'0':2,'1':3}
    monkeypatch.setattr(attributes,'count_n_bowtie2_aligned_genes',align)
    result=add_attributes(pd.DataFrame({'s':['ACGT','ACGT','ACG']}),
                          {'attributes':{'hits':{'target':'s','type':'mapped_genes','aligner':'bowtie2'}}},
                          {'fasta':str(tmp_path/'ref.fa'),'align_index':['bowtie2']})
    assert result.hits.tolist()==[2,2,3]


@pytest.mark.asyncio
async def test_http_layout_roundtrip(monkeypatch):
    from fastapi import UploadFile, HTTPException
    from types import SimpleNamespace
    from uprobe.http.routers import workflow
    saved=[]
    monkeypatch.setattr(workflow,'update_task_in_db',lambda user,task:saved.append(task))
    api=make_api(layout(left=[19,25],gap=[0,2],right=[19,25]))
    payload=yaml.safe_dump(api.protocol).encode()
    await workflow.submit_task(UploadFile(filename='test.yaml',file=io.BytesIO(payload)),SimpleNamespace(username='test'))
    assert yaml.safe_load(saved[0].yaml_content)['extracts']['target_region']['layout']==api.protocol['extracts']['target_region']['layout']
    bad=copy.deepcopy(api.protocol)
    bad['extracts']['target_region']['layout']=layout(left=1,right=1)
    with pytest.raises(HTTPException) as exc:
        await workflow.submit_task(UploadFile(filename='bad.yaml',file=io.BytesIO(yaml.safe_dump(bad).encode())),SimpleNamespace(username='test'))
    assert exc.value.status_code==422 and len(saved)==1
