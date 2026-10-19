"""Expression-aware off-target counting and the per-probe off-target lists of mapped_genes."""
import gzip
from pathlib import Path

import pandas as pd
import pytest
import yaml

from uprobe.core import attributes
from uprobe.core.attributes import add_attributes
from uprobe.core.cli import _validate_and_normalize_protocol
from uprobe.core.gene_expression import ExpressionTable, normalize_gene_id, validate_gene_expression
from uprobe.core.offtarget import GeneAnnotation, OffTargetCounter, format_gene_hits, parse_gene_hits
from uprobe.core.report.formatting import ordered_columns
from uprobe.core.report.run_summary import build_run_summary
from uprobe.core.validate import validate_protocol

DATA = Path(__file__).parent / 'data'
# Gene id -> name; ENSG5 and ENSG6 share a name, ENSG7 has none.
GENES = {'ENSG1': 'ACTB', 'ENSG2': 'ACTG1', 'ENSG3': 'POTEF', 'ENSG4': 'ACTBP2', 'ENSG5': 'Y_RNA',
         'ENSG6': 'Y_RNA', 'ENSG7': None}
SETTINGS = {'aligner': 'bowtie', 'strand': 'sense', 'max_mismatches': 2}


@pytest.fixture
def gtf(tmp_path):
    path = tmp_path / 'genome.gtf'
    lines = ['#!genome-build test\n']
    for i, (gene, name) in enumerate(GENES.items()):
        attrs = f'gene_id "{gene}"; transcript_id "T{i}";' + (f' gene_name "{name}";' if name else '')
        lines.append(f'1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "{gene}";\n')
        lines.append(f'1\ttest\texon\t1\t100\t.\t+\t.\t{attrs}\n')
    path.write_text(''.join(lines))
    return path


def write_table(path, text):
    path.write_text(text)
    return str(path)


@pytest.fixture
def table(tmp_path):
    # Versioned ids; ENSG3 and ACTG1 (by name) missing; ENSG4 unexpressed.
    return write_table(tmp_path / 'tpm.tsv', 'gene_id\tgene_name\tTPM\n'
                       'ENSG1.17\tACTB\t2000\nENSG2.9\tACTG1\t800\nENSG4.1\tACTBP2\t0.2\n'
                       'ENSG5.1\tY_RNA\t0\nENSG6.1\tY_RNA\t3\n')


# -- the gene_expression block -----------------------------------------------

@pytest.mark.parametrize('block,message', [
    ([], 'must be a mapping'),
    ({'table': 'x.tsv'}, 'missing min_value'),
    ({'min_value': 1}, 'missing table'),
    ({'table': 'x.tsv', 'min_value': '1'}, 'min_value must be a number'),
    ({'table': 'x.tsv', 'min_value': True}, 'min_value must be a number'),
    ({'table': 'x.tsv', 'min_value': 1, 'min_tpm': 1}, 'unknown field'),
    ({'table': 'x.tsv', 'min_value': 1, 'gene_column': []}, 'gene_column must be'),
    ({'table': 'x.tsv', 'min_value': 1, 'value_column': 3}, 'value_column must be'),
    ({'table': 'x.tsv', 'min_value': 1, 'min_fraction_of_target': '0.05'}, 'min_fraction_of_target must be'),
    ({'table': 'x.tsv', 'min_value': 1, 'min_fraction_of_target': 1.5}, 'min_fraction_of_target must be'),
    ({'table': 'x.tsv', 'min_value': 1, 'min_fraction_of_target': -0.1}, 'min_fraction_of_target must be'),
    ({'table': 'x.tsv', 'min_value': 1, 'min_fraction_of_target': True}, 'min_fraction_of_target must be'),
])
def test_block_is_checked(block, message):
    with pytest.raises(ValueError, match=message):
        validate_gene_expression(block)


def test_block_normalized():
    settings = validate_gene_expression({'table': 'x.tsv', 'gene_column': 'id', 'min_value': 1})
    assert settings == {'table': 'x.tsv', 'gene_columns': ['id'], 'value_column': None, 'min_value': 1.0,
                        'min_fraction_of_target': None}
    settings = validate_gene_expression({'table': 'x.tsv', 'min_value': 1, 'min_fraction_of_target': 0})
    assert settings['min_fraction_of_target'] == 0.0


def test_normalize_gene_id():
    assert normalize_gene_id(' ENSG00000075624.17 ') == 'ENSG00000075624'
    assert normalize_gene_id('ENSMUSG00000029580.15') == 'ENSMUSG00000029580'
    assert normalize_gene_id('ENSG00000182378.14_PAR_Y') == 'ENSG00000182378_PAR_Y'
    assert normalize_gene_id('ACTB') == 'ACTB'


# -- reading tables ------------------------------------------------------------

def test_table_by_id_and_name(table):
    expression = ExpressionTable({'table': table, 'gene_column': ['gene_id', 'gene_name'], 'min_value': 1})
    assert expression.value_column == 'TPM'  # the only numeric column
    assert expression.value('ENSG1') == 2000 and expression.value('ENSG1.3') == 2000
    assert expression.value('ENSG9', 'ACTG1') == 800  # by name
    assert expression.value('ENSG9', 'Y_RNA') == 3  # a name listed twice: the highest value
    assert expression.value('ENSG3', 'POTEF') is None
    assert expression.describe()['table_rows'] == 5


def test_csv_gzip_and_gct(tmp_path):
    path = tmp_path / 'tpm.csv.gz'
    with gzip.open(path, 'wt') as f:
        f.write('symbol,sample_a,note\nACTB,5,x\nGAPDH,0.5,y\n')
    expression = ExpressionTable({'table': str(path), 'min_value': 1})
    assert (expression.gene_columns, expression.value_column) == (['symbol'], 'sample_a')
    assert expression.value('ENSG1', 'ACTB') == 5
    gct = write_table(tmp_path / 'median.gct', '#1.2\n2\t2\nName\tDescription\tLung\tLiver\n'
                      'ENSG00000075624.17\tACTB\t900\t700\nENSG00000111640.15\tGAPDH\t1000\t5\n')
    with pytest.raises(ValueError, match='set value_column: the expression table has 2 numeric columns'):
        ExpressionTable({'table': gct, 'min_value': 1})
    liver = ExpressionTable({'table': gct, 'value_column': 'Liver', 'min_value': 10})
    assert liver.value('ENSG00000111640') == 5
    comments = write_table(tmp_path / 'counts.txt', '# Program:featureCounts\n# more\nGeneid\tTPM\nACTB\t3\n')
    assert ExpressionTable({'table': comments, 'min_value': 1}).value('ENSG1', 'ACTB') == 3


def test_table_errors(tmp_path, table):
    with pytest.raises(FileNotFoundError):
        ExpressionTable({'table': str(tmp_path / 'missing.tsv'), 'min_value': 1})
    with pytest.raises(ValueError, match='gene_column ensembl not in the expression table'):
        ExpressionTable({'table': table, 'gene_column': 'ensembl', 'min_value': 1})
    with pytest.raises(ValueError, match='value_column gene_name of the expression table is not numeric'):
        ExpressionTable({'table': table, 'value_column': 'gene_name', 'min_value': 1})


# -- counting ------------------------------------------------------------------

def test_gene_names_from_gtf(gtf):
    annotation = GeneAnnotation.from_gtf(gtf)
    assert annotation.names == GENES
    assert annotation.label('ENSG2') == 'ACTG1'
    assert annotation.label('ENSG5') == 'Y_RNA(ENSG5)'  # shared name
    assert annotation.label('ENSG7') == 'ENSG7'
    assert annotation.target_keys('ACTB') == {'ENSG1', 'ACTB'}
    assert GeneAnnotation.from_gtf(None).names == {}


def test_format_and_parse():
    annotation = GeneAnnotation(GENES)
    text = format_gene_hits([('ENSG3', 2), ('ENSG2', 0), ('ENSG5', 1)], annotation)
    assert text == 'ACTG1:0, Y_RNA(ENSG5):1, POTEF:2'
    assert parse_gene_hits(text) == ['ACTG1', 'Y_RNA(ENSG5)', 'POTEF']
    assert format_gene_hits([('ENSG3', None), ('ENSG2', None)], annotation) == 'ACTG1, POTEF'
    assert parse_gene_hits('ACTG1, POTEF') == ['ACTG1', 'POTEF']
    assert parse_gene_hits('') == [] and parse_gene_hits(float('nan')) == []


def hits():
    return {'S1': {'ENSG1': 0, 'ENSG2': 0, 'ENSG3': 1, 'ENSG4': 0, 'ENSG5': 2, 'ENSG6': 2},
            'S2': {'ENSG1': 0}, 'S3': {}}


def test_counts_and_lists_without_expression(gtf):
    df = pd.DataFrame({'target': ['ACTB', 'ACTB', 'GAPDH', 'ACTB'], 'seq': ['S1', 'S2', 'S2', 'S3']})
    report = {}
    OffTargetCounter({}, {'gtf': str(gtf)}, report).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    assert df.hits.tolist() == [6, 1, 1, 0]  # every gene hit, the target included, as before
    assert df.hits_off_targets.tolist() == [
        'ACTBP2:0, ACTG1:0, POTEF:1, Y_RNA(ENSG5):2, Y_RNA(ENSG6):2', '', 'ACTB:0', '']
    assert 'hits_unexpressed' not in df
    assert report == {'attributes': {'hits': {'aligner': 'bowtie', 'max_mismatches': 2, 'strand': 'sense',
                                              'off_target_genes': 6}}}


def test_unexpressed_genes_are_not_counted(gtf, table, caplog):
    df = pd.DataFrame({'target': ['ACTB', 'GAPDH'], 'seq': ['S1', 'S1']})
    protocol = {'gene_expression': {'table': table, 'gene_column': 'gene_id', 'min_value': 1}}
    report = {}
    OffTargetCounter(protocol, {'gtf': str(gtf)}, report).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    # ACTB (target) + ACTG1 (800 TPM) + POTEF (not in the table: counted) + Y_RNA ENSG6 (3 TPM)
    assert df.hits.tolist() == [4, 4]
    assert df.hits_off_targets[0] == 'ACTG1:0, POTEF:1, Y_RNA(ENSG6):2'
    assert df.hits_unexpressed[0] == 'ACTBP2:0, Y_RNA(ENSG5):2'
    assert df.hits_off_targets[1] == 'ACTB:0, ACTG1:0, POTEF:1, Y_RNA(ENSG6):2'  # GAPDH's off-target ACTB
    stats = report['attributes']['hits']
    assert (stats['off_target_genes'], stats['below_min_value'], stats['missing_from_table']) == (6, 2, 1)
    assert stats['missing_from_table_examples'] == ['POTEF']
    assert report['gene_expression']['annotation_genes_in_table'] == 5
    assert 'not in the expression table and count as expressed (e.g. POTEF)' in caplog.text


def test_target_counts_even_when_unexpressed(gtf, tmp_path):
    table = write_table(tmp_path / 'tpm.tsv', 'gene\tTPM\nACTB\t0\nACTG1\t0\n')
    df = pd.DataFrame({'target': ['ACTB'], 'seq': ['S1']})
    protocol = {'gene_expression': {'table': table, 'min_value': 1}}
    OffTargetCounter(protocol, {'gtf': str(gtf)}).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    assert df.hits[0] == 5  # ACTB + POTEF, ACTBP2, Y_RNA x2 (not in the table); ACTG1 ignored
    assert df.hits_unexpressed[0] == 'ACTG1:0'


def test_threshold_relative_to_the_target(gtf, table, caplog):
    # GAPDH is neither annotated nor in the table: only min_value applies to its off-targets.
    df = pd.DataFrame({'target': ['ACTB', 'GAPDH'], 'seq': ['S1', 'S1']})
    protocol = {'gene_expression': {'table': table, 'gene_column': 'gene_id', 'min_value': 1,
                                    'min_fraction_of_target': 0.05}}
    report = {}
    OffTargetCounter(protocol, {'gtf': str(gtf)}, report).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    # ACTB (2000 TPM, found by its gene id): threshold max(1, 0.05 x 2000) = 100, so Y_RNA ENSG6 (3) no longer counts.
    assert df.hits.tolist() == [3, 4]
    assert df.hits_off_targets.tolist() == ['ACTG1:0, POTEF:1', 'ACTB:0, ACTG1:0, POTEF:1, Y_RNA(ENSG6):2']
    assert df.hits_unexpressed.tolist() == ['ACTBP2:0, Y_RNA(ENSG5):2, Y_RNA(ENSG6):2', 'ACTBP2:0, Y_RNA(ENSG5):2']
    stats = report['attributes']['hits']
    assert (stats['below_min_value'], stats['below_target_fraction'], stats['missing_from_table']) == (2, 1, 1)
    described = report['gene_expression']
    assert described['min_fraction_of_target'] == 0.05
    assert described['targets'] == {'ACTB': {'value': 2000.0, 'threshold': 100.0}}
    assert described['targets_not_in_table'] == ['GAPDH']
    assert 'target(s) not in the expression table, so only min_value applies' in caplog.text
    final = df.assign(hits_off_targets=df.hits_off_targets)
    summary = build_run_summary({'name': 'x'}, {'final': final, 'off_targets': report})
    warning = next(w for w in summary['warnings'] if w['code'] == 'targets_not_in_expression_table')
    assert warning['targets'] == ['GAPDH']


def test_threshold_fraction_can_exceed_min_value_only(gtf, table):
    df = pd.DataFrame({'target': ['ACTB', 'ACTG1'], 'seq': ['S1', 'S1']})
    protocol = {'gene_expression': {'table': table, 'gene_column': ['gene_id', 'gene_name'], 'min_value': 1,
                                    'min_fraction_of_target': 0.5}}
    report = {}
    OffTargetCounter(protocol, {'gtf': str(gtf)}, report).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    # ACTB: threshold 1000 drops ACTG1 (800); ACTG1: threshold 400 keeps ACTB (2000).
    assert df.hits_off_targets.tolist() == ['POTEF:1', 'ACTB:0, POTEF:1']
    assert report['gene_expression']['targets'] == {'ACTB': {'value': 2000.0, 'threshold': 1000.0},
                                                    'ACTG1': {'value': 800.0, 'threshold': 400.0}}
    assert report['gene_expression']['targets_not_in_table'] == []
    # Without the option every target's threshold is min_value and no target is reported.
    del protocol['gene_expression']['min_fraction_of_target']
    report = {}
    OffTargetCounter(protocol, {'gtf': str(gtf)}, report).annotate(df, 'hits', 'seq', hits(), SETTINGS)
    assert df.hits_off_targets[0] == 'ACTG1:0, POTEF:1, Y_RNA(ENSG6):2'
    assert 'targets' not in report['gene_expression'] and 'below_target_fraction' not in report['attributes']['hits']


def test_table_matching_no_annotated_gene_is_an_error(gtf, tmp_path):
    table = write_table(tmp_path / 'mouse.tsv', 'gene\tTPM\nActb\t5\n')
    df = pd.DataFrame({'target': ['ACTB'], 'seq': ['S1']})
    counter = OffTargetCounter({'gene_expression': {'table': table, 'min_value': 1}}, {'gtf': str(gtf)})
    with pytest.raises(ValueError, match='none of the 7 genes of the genome annotation is in'):
        counter.annotate(df, 'hits', 'seq', hits(), SETTINGS)


def test_add_attributes_lists_off_targets(monkeypatch, gtf, table, tmp_path):
    monkeypatch.chdir(tmp_path)
    by_sequence = {'AAAA': hits()['S1'], 'CCCC': hits()['S2']}
    monkeypatch.setattr(attributes, 'bowtie_gene_hits',
                        lambda _, records, *args: {key: by_sequence[seq] for key, seq in records.items()})
    df = pd.DataFrame({'target': ['ACTB', 'ACTB'], 'target_region': ['AAAA', 'CCCC']})
    protocol = {'gene_expression': {'table': table, 'min_value': 1},
                'attributes': {'target_mapped_genes': {'target': 'target_region', 'type': 'mapped_genes',
                                                       'aligner': 'bowtie', 'strand': 'sense'}}}
    report = {}
    result = add_attributes(df, protocol, {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(gtf),
                                           'align_index': ['bowtie']}, report)
    assert result.target_mapped_genes.tolist() == [4, 1]
    assert result.target_mapped_genes_off_targets.tolist() == ['ACTG1:0, POTEF:1, Y_RNA(ENSG6):2', '']
    assert result.target_mapped_genes_unexpressed.tolist() == ['ACTBP2:0, Y_RNA(ENSG5):2', '']
    assert report['attributes']['target_mapped_genes']['missing_from_table'] == 1


# -- outputs -------------------------------------------------------------------

def test_summary_reports_off_targets_and_missing_genes():
    report = {'gene_expression': {'table': 'tpm.tsv', 'min_value': 1.0},
              'attributes': {'hits': {'aligner': 'bowtie', 'max_mismatches': 2, 'strand': 'sense',
                                      'off_target_genes': 6, 'below_min_value': 2, 'missing_from_table': 2,
                                      'missing_from_table_examples': ['LINC1', 'POTEF'],
                                      'missing_labels': ['LINC1', 'POTEF']}}}
    final = pd.DataFrame({'target': ['ACTB', 'ACTB'], 'hits_off_targets': ['POTEF:1, ACTG1:0', '']})
    summary = build_run_summary({'name': 'x'}, {'final': final, 'off_targets': report})
    entry = summary['off_targets']['attributes']['hits']
    assert entry['final_off_target_genes'] == 2 and entry['final_missing_from_table'] == ['POTEF']
    assert 'missing_labels' not in entry
    assert summary['off_targets']['gene_expression'] == report['gene_expression']
    warning = next(w for w in summary['warnings'] if w['code'] == 'off_target_genes_not_in_expression_table')
    assert warning['attributes'] == {'hits': 2}
    assert build_run_summary({}, {'final': final})['off_targets'] is None


def test_off_target_columns_follow_their_attribute():
    protocol = {'extracts': {'target_region': {}}, 'probes': {'p': {}},
                'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes'},
                               'gc': {'target': 'target_region', 'type': 'gc_content'}}}
    columns = ['probe_id', 'hits_unexpressed', 'hits_off_targets', 'target_region', 'gc', 'hits', 'p']
    assert ordered_columns(columns, protocol) == [
        'probe_id', 'target_region', 'gc', 'hits', 'hits_unexpressed', 'hits_off_targets', 'p']


# -- uprobe run / uprobe validate --------------------------------------------

@pytest.fixture
def rna():
    protocol = yaml.safe_load((DATA / 'RNA_format.yaml').read_text(encoding='utf-8'))
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    return protocol


def test_run_rejects_an_invalid_block(rna):
    rna['gene_expression'] = {'table': 'tpm.tsv'}
    with pytest.raises(ValueError, match='missing min_value'):
        _validate_and_normalize_protocol(rna)


def test_validate_checks_the_block_and_the_table(rna, table, tmp_path):
    rna['gene_expression'] = {'table': table, 'gene_column': 'gene_id', 'value_column': 'TPM', 'min_value': 1,
                              'min_fraction_of_target': 0.05}
    assert validate_protocol(rna) == []

    def finding(path, severity='error'):
        found = [f for f in validate_protocol(rna) if f.path == path and f.severity == severity]
        assert found, path
        return found[0]

    rna['gene_expression']['value_column'] = 'tpm'
    assert 'value_column tpm not in the expression table' in finding('gene_expression.value_column').message
    assert 'TPM' in finding('gene_expression.value_column').suggestion
    rna['gene_expression'].update(value_column='TPM', gene_column='ensembl')
    assert 'gene_column ensembl not in' in finding('gene_expression.gene_column').message
    rna['gene_expression'].update(gene_column='gene_id', table=str(tmp_path / 'missing.tsv'))
    assert 'not found' in finding('gene_expression.table').message
    rna['gene_expression'] = {'table': table}
    assert 'missing min_value' in finding('gene_expression').message
    rna['gene_expression'] = {'table': table, 'min_value': 1, 'min_fraction_of_target': 5}
    assert 'min_fraction_of_target must be a number from 0 to 1' in finding('gene_expression').message
    assert 'min_fraction_of_target: 0.05' in finding('gene_expression').suggestion
    rna['gene_expression'] = {'table': table, 'min_value': 1}
    del rna['attributes']['target_mappedGenes']
    assert 'there is none' in finding('gene_expression', 'warning').message
