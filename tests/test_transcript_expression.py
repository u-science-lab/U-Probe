"""gene_expression.level: transcript: an off-target gene is valued by the transcripts a sequence hits."""
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml

from uprobe.core import attributes
from uprobe.core.attributes import add_attributes, gene_hits_of
from uprobe.core.cli import _validate_and_normalize_protocol
from uprobe.core.gene_expression import ExpressionTable, identifier_level, unversioned, validate_gene_expression
from uprobe.core.offtarget import GeneAnnotation, OffTargetCounter
from uprobe.core.process.off_target_pool import limit_off_target_pool
from uprobe.core.report.formatting import ordered_columns
from uprobe.core.report.run_summary import build_run_summary
from uprobe.core.validate import validate_protocol

DATA = Path(__file__).parent / 'data'
needs_bowtie = pytest.mark.skipif(not (shutil.which('bowtie') and shutil.which('bowtie-build')),
                                  reason='needs Bowtie 1 (bowtie, bowtie-build)')
SETTINGS = {'aligner': 'bowtie', 'strand': 'sense', 'max_mismatches': 2}
# Gene id -> (name, transcripts). Lnc carries a segment of Gapdh on T2A only, the way mouse
# lncRNAs embed Gapdh segments; its expression lies on T2B, which the probes do not hit.
GENES = {'G1': ('Gapdh', ['T1A', 'T1B']), 'G2': ('Lnc', ['T2A', 'T2B']), 'G3': ('Para', ['T3A', 'T3B']),
         'G4': ('Pseudo', ['T4A', 'T4B']), 'G5': ('Other', ['T5A'])}
# Versioned, as salmon lists Ensembl transcripts; T4A and T5A are missing.
TRANSCRIPT_TPM = {'T1A.2': 1500, 'T1B.1': 500, 'T2A.1': 0, 'T2B.3': 50, 'T3A.1': 30, 'T3B.1': 0, 'T4B.1': 5}
# The same sample by gene: the sums over each gene's transcripts.
GENE_TPM = {'G1': 2000, 'G2': 50, 'G3': 30, 'G4': 5}


def hits_of(*records):
    """{transcript.fa record: mismatches} of one sequence; a record given as 'G2_T2A:1' has 1 mismatch."""
    return {record.split(':')[0]: int(record.split(':')[1]) if ':' in record else 0 for record in records}


# The transcripts each window hits, its target's included.
TRANSCRIPT_HITS = {
    'S1': hits_of('G1_T1A', 'G1_T1B', 'G2_T2A:1'),  # Lnc on its unexpressed isoform only
    'S2': hits_of('G1_T1A', 'G2_T2A', 'G2_T2B:2'),  # Lnc on both isoforms: 50 TPM
    'S3': hits_of('G1_T1A', 'G3_T3A:1', 'G3_T3B:1'),  # Para: 30 TPM
    'S4': hits_of('G1_T1A', 'G3_T3B'),  # Para's unexpressed isoform
    'S5': hits_of('G1_T1A', 'G4_T4A:1'),  # a transcript missing from the table
    'S6': hits_of('G1_T1A', 'G4_T4A:1', 'G4_T4B:2'),  # one missing, one at 5 TPM
    'S7': hits_of('G1_T1A', 'G5_T5A'),  # a gene the table does not list at all
}


def gene_hits():
    return {sequence: gene_hits_of(records) for sequence, records in TRANSCRIPT_HITS.items()}


@pytest.fixture
def gtf(tmp_path):
    path = tmp_path / 'genome.gtf'
    path.write_text(''.join(f'1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "{gene}"; transcript_id "{transcript}"; '
                            f'gene_name "{name}";\n'
                            for gene, (name, transcripts) in GENES.items() for transcript in transcripts))
    return path


def write_table(path, header, rows):
    path.write_text(header + ''.join(f'{key}\t{value}\n' for key, value in rows.items()))
    return str(path)


@pytest.fixture
def transcript_table(tmp_path):
    # Like salmon's quant.sf: several numeric columns, so value_column: TPM.
    path = tmp_path / 'quant.sf'
    path.write_text('Name\tLength\tTPM\tNumReads\n'
                    + ''.join(f'{name}\t1000\t{tpm}\t{tpm * 10}\n' for name, tpm in TRANSCRIPT_TPM.items()))
    return str(path)


@pytest.fixture
def gene_table(tmp_path):
    return write_table(tmp_path / 'genes.tsv', 'gene_id\tTPM\n', GENE_TPM)


def by_transcript(table, **options):
    return {'gene_expression': {'table': table, 'level': 'transcript', 'value_column': 'TPM', 'min_value': 1,
                                **options}}


def annotate(protocol, gtf, sequences, target='Gapdh', report=None):
    df = pd.DataFrame({'target': [target] * len(sequences), 'seq': sequences})
    level = protocol.get('gene_expression', {}).get('level')
    counter = OffTargetCounter(protocol, {'gtf': str(gtf)}, {} if report is None else report)
    counter.annotate(df, 'hits', 'seq', gene_hits(), SETTINGS, TRANSCRIPT_HITS if level == 'transcript' else None)
    return df.set_index('seq')


# -- the block and the table -----------------------------------------------------------------------

def test_block_with_level():
    settings = validate_gene_expression({'table': 'q.sf', 'level': 'transcript', 'transcript_column': 'Name',
                                         'min_value': 1})
    assert (settings['level'], settings['transcript_column'], settings['gene_columns']) == ('transcript', 'Name', None)
    assert validate_gene_expression({'table': 'x.tsv', 'min_value': 1})['level'] == 'gene'


@pytest.mark.parametrize('block,message', [
    ({'level': 'isoform'}, 'level must be one of gene, transcript'),
    ({'level': 'transcript', 'gene_column': 'gene_id'}, 'gene_column applies to level: gene'),
    ({'transcript_column': 'Name'}, 'transcript_column applies to level: transcript'),
    ({'level': 'transcript', 'transcript_column': ['Name']}, 'transcript_column must be a column name'),
])
def test_block_with_level_is_checked(block, message):
    with pytest.raises(ValueError, match=message):
        validate_gene_expression({'table': 'x.tsv', 'min_value': 1, **block})


def test_transcripts_are_matched_with_or_without_version(tmp_path):
    table = write_table(tmp_path / 't.tsv', 'transcript_id\tTPM\n', {
        'ENSMUST00000073605.15': 7, 'ENST00000381192.10_PAR_Y': 2, 'NM_008084.3': 4,
        'ENST00000000001.1|ENSG00000000001.1|-|-|A-201|A|1000|protein_coding|': 9,
        'MSTRG.1.1': 1, 'MSTRG.1.2': 3, 'ENSMUST00000000002': 6})
    expression = ExpressionTable(by_transcript(table)['gene_expression'])
    assert expression.transcript_column == 'transcript_id'
    assert expression.transcript_value('ENSMUST00000073605') == 7  # an unversioned annotation
    assert expression.transcript_value('ENSMUST00000073605.14') == 7  # another version
    assert expression.transcript_value('ENSMUST00000000002.3') == 6  # a versioned annotation
    assert expression.transcript_value('ENST00000381192_PAR_Y') == 2
    assert expression.transcript_value('ENST00000381192') is None
    assert expression.transcript_value('NM_008084') == 4
    assert expression.transcript_value('ENST00000000001') == 9  # a GENCODE FASTA header as the name
    # Isoform numbers are not versions: exact ids match, the shared stem does not.
    assert (expression.transcript_value('MSTRG.1.2'), expression.transcript_value('MSTRG.1')) == (3, None)
    assert expression.transcripts_value(['ENSMUST00000073605', 'NM_008084', 'XM_1']) == (11, ['XM_1'])
    assert expression.transcripts_value(['XM_1']) == (None, ['XM_1'])
    assert unversioned('ENSMUST00000073605.15') == 'ENSMUST00000073605'
    assert expression.describe()['level'] == 'transcript' and 'gene_columns' not in expression.describe()


def test_identifier_level():
    assert identifier_level(['ENSMUST00000073605.15', 'ENST00000381192', 'NM_008084.3', 'XR_1']) == 'transcript'
    assert identifier_level(['ENSMUSG00000057666', 'ENSG00000111640.15']) == 'gene'
    assert identifier_level(['Gapdh', 'Actb']) is None and identifier_level([]) is None


def test_gene_annotation_reads_the_transcripts(gtf):
    assert GeneAnnotation.from_gtf(gtf).transcripts == {}  # read only when asked for
    annotation = GeneAnnotation.from_gtf(gtf, transcripts=True)
    assert annotation.names['G2'] == 'Lnc' and annotation.gene_transcripts('G2') == ['T2A', 'T2B']
    assert annotation.transcript_id('G2_T2B') == 'T2B'
    assert annotation.transcript_id('G9_NM_008084.3') == 'NM_008084.3'  # not annotated: after the gene id


def test_gene_ids_with_underscores_keep_their_transcripts(tmp_path):
    gtf = tmp_path / 'gencode.gtf'
    gtf.write_text('X\tt\texon\t1\t9\t.\t+\t.\tgene_id "ENSG1.1_PAR_Y"; transcript_id "ENST1.1_PAR_Y"; '
                   'gene_name "SHOX";\n')
    annotation = GeneAnnotation.from_gtf(gtf, transcripts=True)
    assert annotation.transcript_id('ENSG1.1_PAR_Y_ENST1.1_PAR_Y') == 'ENST1.1_PAR_Y'


# -- counting --------------------------------------------------------------------------------------

def test_expression_on_an_unhit_isoform_does_not_count(gtf, transcript_table, gene_table):
    rows = annotate(by_transcript(transcript_table), gtf, ['S1', 'S2', 'S3', 'S4'])
    # Lnc's 50 TPM lie on T2B: S1 hits only T2A (0 TPM), S2 both. Para's 30 TPM lie on T3A.
    assert rows.hits.to_dict() == {'S1': 1, 'S2': 2, 'S3': 2, 'S4': 1}
    assert rows.hits_off_targets.to_dict() == {'S1': '', 'S2': 'Lnc:0', 'S3': 'Para:1', 'S4': ''}
    assert rows.hits_unexpressed.to_dict() == {'S1': 'Lnc:1', 'S2': '', 'S3': '', 'S4': 'Para:0'}
    assert rows.hits_hit_transcripts.to_dict() == {
        'S1': 'Lnc: T2A=0', 'S2': 'Lnc: T2A=0, T2B=50', 'S3': 'Para: T3A=30, T3B=0', 'S4': 'Para: T3B=0'}
    # By gene, the same sample counts Lnc and Para wherever they are hit.
    rows = annotate({'gene_expression': {'table': gene_table, 'min_value': 1}}, gtf, ['S1', 'S2', 'S3', 'S4'])
    assert rows.hits.to_dict() == {'S1': 2, 'S2': 2, 'S3': 2, 'S4': 2}
    assert 'hits_hit_transcripts' not in rows


def test_the_target_value_is_the_sum_over_all_its_transcripts(gtf, transcript_table, caplog):
    report = {}
    rows = annotate(by_transcript(transcript_table, min_fraction_of_target=0.05), gtf, ['S2', 'S3'], report=report)
    # Gapdh: T1A 1500 + T1B 500, though S3 hits T1A only; threshold max(1, 0.05 x 2000) = 100.
    assert report['gene_expression']['targets'] == {'Gapdh': {'value': 2000.0, 'threshold': 100.0}}
    assert rows.hits.to_dict() == {'S2': 1, 'S3': 1}
    assert rows.hits_unexpressed.to_dict() == {'S2': 'Lnc:0', 'S3': 'Para:1'}
    assert report['attributes']['hits']['below_target_fraction'] == 2
    report = {}
    rows = annotate(by_transcript(transcript_table, min_fraction_of_target=0.01), gtf, ['S2', 'S3'], report=report)
    assert rows.hits.to_dict() == {'S2': 2, 'S3': 2}  # threshold 20: Lnc (50) and Para (30) count
    # A target none of whose transcripts is listed keeps min_value.
    report = {}
    annotate(by_transcript(transcript_table, min_fraction_of_target=0.05), gtf, ['S3'], target='Other',
             report=report)
    assert report['gene_expression']['targets_not_in_table'] == ['Other']
    assert report['gene_expression']['targets'] == {}


@pytest.mark.parametrize('absent', ['count', 'ignore'])
def test_hit_transcripts_missing_from_the_table(gtf, transcript_table, caplog, absent):
    report = {}
    rows = annotate(by_transcript(transcript_table, absent=absent), gtf, ['S5', 'S6', 'S7'], report=report)
    assert rows.hits_hit_transcripts.to_dict() == {'S5': 'Pseudo: T4A=NA', 'S6': 'Pseudo: T4A=NA, T4B=5',
                                                   'S7': 'Other: T5A=NA'}
    if absent == 'count':  # a missing transcript may be expressed: the gene counts
        assert rows.hits.to_dict() == {'S5': 2, 'S6': 2, 'S7': 2}
        assert 'which make them count as expressed (e.g. Other, Pseudo)' in caplog.text
    else:  # missing transcripts add nothing: Pseudo counts on S6 by T4B alone (5 TPM)
        assert rows.hits.to_dict() == {'S5': 1, 'S6': 2, 'S7': 1}
        assert rows.hits_unexpressed.to_dict() == {'S5': 'Pseudo:1', 'S6': '', 'S7': 'Other:0'}
        assert 'which add nothing to their value (gene_expression.absent: ignore)' in caplog.text
    stats = report['attributes']['hits']
    assert (stats['missing_from_table'], stats['missing_transcripts'], stats['absent']) == (2, 2, absent)
    assert stats['missing_from_table_examples'] == ['Other', 'Pseudo']
    summary = build_run_summary({'name': 'x'}, {'final': rows.reset_index(), 'off_targets': report})
    assert summary['off_targets']['gene_expression']['level'] == 'transcript'
    warning = next(w for w in summary['warnings'] if w['code'] == 'off_target_genes_not_in_expression_table')
    assert warning['attributes'] == {'hits': 2} and 'have hit transcripts that are not in' in warning['message']


def test_summary_reports_the_level(gtf, transcript_table, gene_table):
    report = {}
    annotate(by_transcript(transcript_table), gtf, ['S1'], report=report)
    described = report['gene_expression']
    assert (described['level'], described['transcript_column'], described['value_column']) == (
        'transcript', 'Name', 'TPM')
    assert (described['annotation_transcripts'], described['annotation_transcripts_in_table']) == (9, 7)
    assert (described['annotation_genes'], described['annotation_genes_in_table']) == (5, 4)
    report = {}
    annotate({'gene_expression': {'table': gene_table, 'min_value': 1}}, gtf, ['S1'], report=report)
    assert report['gene_expression']['level'] == 'gene' and 'transcript_column' not in report['gene_expression']


def test_tables_that_do_not_fit_the_level(tmp_path):
    gtf = tmp_path / 'ensembl.gtf'
    gtf.write_text('1\tt\texon\t1\t9\t.\t+\t.\tgene_id "ENSMUSG00000057666"; transcript_id "ENSMUST00000073605"; '
                   'gene_name "Gapdh";\n')
    df = pd.DataFrame({'target': ['Gapdh'], 'seq': ['S']})
    genes = write_table(tmp_path / 'genes.tsv', 'gene_id\tTPM\n', {'ENSMUSG00000057666.9': 5, 'ENSMUSG00000000001': 1})
    transcripts = write_table(tmp_path / 'tx.tsv', 'transcript_id\tTPM\n', {'ENSMUST00000073605.15': 5,
                                                                            'ENSMUST00000000001.5': 1})
    with pytest.raises(ValueError, match='none of the 1 transcripts .* look like gene ids: set gene_expression.level: '
                                         'gene'):
        OffTargetCounter(by_transcript(genes), {'gtf': str(gtf)}).annotate(df, 'h', 'seq', {}, SETTINGS, {})
    with pytest.raises(ValueError, match='none of the 1 genes .* look like transcript ids: set '
                                         'gene_expression.level: transcript'):
        OffTargetCounter({'gene_expression': {'table': transcripts, 'min_value': 1}},
                         {'gtf': str(gtf)}).annotate(df, 'h', 'seq', {}, SETTINGS)
    with pytest.raises(ValueError, match="needs the genome's GTF"):
        OffTargetCounter(by_transcript(transcripts), {}).annotate(df, 'h', 'seq', {}, SETTINGS, {})
    with pytest.raises(ValueError, match='needs the transcripts each sequence hits'):
        OffTargetCounter(by_transcript(transcripts), {'gtf': str(gtf)}).annotate(df, 'h', 'seq', {}, SETTINGS)


def test_add_attributes_searches_transcripts_for_a_table_of_transcripts(monkeypatch, gtf, transcript_table, tmp_path):
    monkeypatch.chdir(tmp_path)
    by_sequence = {'AAAA': TRANSCRIPT_HITS['S1'], 'CCCC': TRANSCRIPT_HITS['S2']}
    calls = []

    def search(name):
        def hits(_, records, *args):
            calls.append(name)
            return {key: dict(by_sequence[seq]) for key, seq in records.items()}
        return hits

    for name in ('bowtie_transcript_hits', 'bowtie2_transcript_hits'):
        monkeypatch.setattr(attributes, name, search(name))
    for name in ('bowtie_gene_hits', 'bowtie2_gene_hits'):
        monkeypatch.setattr(attributes, name, lambda *args: pytest.fail('the gene search does not list transcripts'))
    df = pd.DataFrame({'target': ['Gapdh', 'Gapdh'], 'target_region': ['AAAA', 'CCCC']})
    for aligner in ('bowtie', 'bowtie2'):
        protocol = {**by_transcript(transcript_table),
                    'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': aligner}}}
        result = add_attributes(df.copy(), protocol, {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(gtf),
                                                      'align_index': [aligner]})
        assert result.hits.tolist() == [1, 2] and result.hits_off_targets.tolist() == ['', 'Lnc:0']
        assert result.hits_hit_transcripts.tolist() == ['Lnc: T2A=0', 'Lnc: T2A=0, T2B=50']
    assert calls == ['bowtie_transcript_hits', 'bowtie2_transcript_hits']


def test_off_target_columns_follow_their_attribute():
    protocol = {'extracts': {'target_region': {}}, 'attributes': {'hits': {'target': 'target_region'}}}
    columns = ['hits_hit_transcripts', 'hits_unexpressed', 'target_region', 'hits', 'hits_off_targets']
    assert ordered_columns(columns, protocol) == [
        'target_region', 'hits', 'hits_hit_transcripts', 'hits_unexpressed', 'hits_off_targets']


# -- copies, the pool limit and tables of genes work as before -------------------------------------

def test_copies_come_before_transcript_expression(tmp_path, gtf, transcript_table):
    # Lnc copies Gapdh: exempt whatever the transcripts it is hit on.
    target = 'ACGTTGCAAGGCTTACCGATGCATGCCAGTTAGCATCGGATCCATGGTACCAGCTAGGTCATGCAA' * 3
    (tmp_path / 'transcript.fa').write_text(f'>G1_T1A\n{target}\n>G2_T2A\n{target}\n>G3_T3A\nACGT\n')
    df = pd.DataFrame({'target': 'Gapdh', 'seq': ['S2', 'S3'], 'target_start': [1, 50],
                       'target_end': [len(target), 69]})
    df['target_region'] = [target, target[49:69]]
    report = {}
    settings = {**SETTINGS, 'exempt_copies': {'kmer': 30, 'min_shared': 0.5}}
    OffTargetCounter(by_transcript(transcript_table), {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(gtf)},
                     report).annotate(df, 'hits', 'seq', gene_hits(), settings, TRANSCRIPT_HITS)
    rows = df.set_index('seq')
    assert rows.hits_copies['S2'] == 'Lnc:0' and rows.hits_off_targets['S2'] == ''
    assert rows.hits_hit_transcripts['S2'] == ''  # a copy is not valued
    assert rows.hits_off_targets['S3'] == 'Para:1'
    assert report['attributes']['hits']['copy_genes'] == 1


def test_a_table_of_genes_works_as_before(gtf, gene_table):
    counter = OffTargetCounter({'gene_expression': {'table': gene_table, 'min_value': 1}}, {'gtf': str(gtf)})
    df = pd.DataFrame({'target': 'Gapdh', 'seq': ['S1', 'S4']})
    counter.annotate(df, 'hits', 'seq', gene_hits(), SETTINGS)  # no transcript hits needed
    assert counter.annotation.transcripts == {}  # nor the GTF's transcripts
    assert list(df.columns) == ['target', 'seq', 'hits', 'hits_off_targets', 'hits_unexpressed']
    assert df.hits_off_targets.tolist() == ['Lnc:1', 'Para:0']


def test_the_pool_limit_counts_the_genes(gtf, transcript_table):
    df = annotate(by_transcript(transcript_table), gtf, ['S2', 'S3', 'S3', 'S6']).reset_index()
    df['target_start'] = range(len(df))
    kept, report = limit_off_target_pool(df, {'attribute': 'hits', 'max_probes_per_gene': 1})
    # Genes, not transcripts, are limited: the second S3 probe hits Para again.
    assert kept.seq.tolist() == ['S2', 'S3', 'S6'] and report['genes_at_limit'] == {'Para': 1}


# -- uprobe run / uprobe validate ------------------------------------------------------------------

@pytest.fixture
def rna():
    protocol = yaml.safe_load((DATA / 'RNA_format.yaml').read_text(encoding='utf-8'))
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    return protocol


def test_run_rejects_a_gene_column_with_level_transcript(rna):
    rna['gene_expression'] = {'table': 'q.sf', 'level': 'transcript', 'gene_column': 'gene_id', 'min_value': 1}
    with pytest.raises(ValueError, match='gene_column applies to level: gene'):
        _validate_and_normalize_protocol(rna)


def test_validate_checks_the_level(rna, tmp_path):
    table = write_table(tmp_path / 'tx.tsv', 'transcript_id\tgene_id\tTPM\tNumReads\n', {
        'ENSMUST00000073605.15\tENSMUSG00000057666': '5\t10', 'ENSMUST00000000001.5\tENSMUSG00000000001': '1\t2'})

    def findings(severity='error'):
        return {f.path: f for f in validate_protocol(rna) if f.severity == severity}

    rna['gene_expression'] = {'table': table, 'level': 'transcript', 'value_column': 'TPM', 'min_value': 1}
    assert validate_protocol(rna) == []
    rna['gene_expression']['transcript_column'] = 'gene_id'
    found = findings()['gene_expression.level']
    assert 'level is transcript, but column gene_id holds gene ids (e.g. ENSMUSG00000057666)' in found.message
    assert found.suggestion == 'level: gene\ngene_column: gene_id   # in place of transcript_column'
    rna['gene_expression'] = {'table': table, 'value_column': 'TPM', 'min_value': 1}  # level: gene
    found = findings()['gene_expression.level']
    assert 'level is gene, but column transcript_id holds transcript ids' in found.message
    assert found.suggestion == 'level: transcript'
    rna['gene_expression'] = {'table': table, 'level': 'transcript', 'transcript_column': 'Name', 'min_value': 1}
    found = findings()['gene_expression.transcript_column']
    assert 'transcript_column Name not in the expression table' in found.message
    rna['gene_expression'] = {'table': table, 'level': 'transcript', 'min_value': 1}
    assert 'set value_column' in findings()['gene_expression.value_column'].message
    rna['gene_expression'] = {'table': table, 'level': 'isoform', 'min_value': 1}
    found = findings()['gene_expression']
    assert 'level must be one of gene, transcript' in found.message and 'level: gene' in found.suggestion


def test_smfish_template_recommends_a_table_of_transcripts(tmp_path):
    from uprobe.templates import read_template
    text = read_template('smfish_stellaris')
    lines = text[text.index('\n# gene_expression:\n') + 1:].splitlines()
    block = '\n'.join(line[2:] for line in lines[:next(i for i, line in enumerate(lines) if not line.startswith('#'))])
    # Uncommenting the level line makes the example a table of transcripts, valid as it stands.
    example = yaml.safe_load(block.replace('# level: transcript', 'level: transcript'))['gene_expression']
    settings = validate_gene_expression(example)
    assert (settings['level'], settings['value_column'], settings['min_fraction_of_target']) == (
        'transcript', 'TPM', 0.05)
    assert 'by transcript if you have one' in text


# -- with Bowtie 1 ---------------------------------------------------------------------------------

@needs_bowtie
def test_add_attributes_with_bowtie_values_the_hit_isoform(tmp_path, monkeypatch):
    """A toy genome: Lnc's first exon embeds a Gapdh segment; only its other isoform is expressed."""
    from uprobe.core.tools import build_transcripts_index
    import random
    rng = random.Random(7)

    def random_sequence(length):
        return ''.join(rng.choice('ACGT') for _ in range(length))

    gapdh, embedded, lnc_exon2, para = random_sequence(300), random_sequence(100), random_sequence(200), \
        random_sequence(300)
    gapdh = gapdh[:100] + embedded + gapdh[200:]
    chromosome = (random_sequence(50) + gapdh + random_sequence(50) + embedded + random_sequence(50) + lnc_exon2
                  + random_sequence(50) + para + random_sequence(50))
    (tmp_path / 'genome.fa').write_text(f'>1\n{chromosome}\n')
    exons = [('G1', 'Gapdh', 'T1A', 51, 350), ('G2', 'Lnc', 'T2A', 401, 500), ('G2', 'Lnc', 'T2A', 551, 750),
             ('G2', 'Lnc', 'T2B', 551, 750), ('G3', 'Para', 'T3A', 801, 1100)]
    (tmp_path / 'genome.gtf').write_text(''.join(
        f'1\tt\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{tid}"; gene_name "{name}";\n'
        for gene, name, tid, start, end in exons))
    build_transcripts_index(tmp_path / 'genome.gtf', tmp_path / 'genome.fa', tmp_path / 'bowtie_transcript',
                            threads=1, aligner='bowtie')
    table = write_table(tmp_path / 'tx.tsv', 'transcript_id\tTPM\n', {'T1A': 900, 'T2A': 0, 'T2B': 80, 'T3A': 9})
    monkeypatch.chdir(tmp_path)
    windows = [gapdh[i:i + 20] for i in range(0, len(gapdh) - 19, 10)]
    df = pd.DataFrame({'target': 'Gapdh', 'target_region': windows})
    attribute = {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie', 'max_mismatches': 0,
                 'strand': 'sense'}
    genome = {'fasta': str(tmp_path / 'genome.fa'), 'gtf': str(tmp_path / 'genome.gtf'), 'align_index': ['bowtie']}
    inside = [100 <= i * 10 and i * 10 + 20 <= 200 for i in range(len(windows))]
    for level, counted in (('gene', True), ('transcript', False)):
        protocol = {'attributes': {'hits': attribute}, 'gene_expression': {
            'table': table, 'level': level, 'min_value': 1} if level == 'transcript' else {
            'table': write_table(tmp_path / 'g.tsv', 'gene\tTPM\n', {'Gapdh': 900, 'Lnc': 80, 'Para': 9}),
            'min_value': 1}}
        result = add_attributes(df.copy(), protocol, genome)
        # Windows in the embedded segment hit Lnc's T2A (0 TPM); Lnc's 80 TPM lie on T2B.
        assert result.hits_off_targets[inside].eq('Lnc:0' if counted else '').all(), level
        assert (result.hits[[not x for x in inside]] == 1).all()
        if level == 'transcript':
            assert result.hits_unexpressed[inside].eq('Lnc:0').all()
            assert result.hits_hit_transcripts[inside].eq('Lnc: T2A=0').all()


# -- Ensembl 113 mouse with the NIH/3T3 salmon tables (read only; skipped when absent) -------------

# Folders prepared by the probe-design benchmark: the transcriptome (records.tsv, sequences.fa and
# its Bowtie 1 index) and the salmon quantification by gene and by transcript.
REFERENCE = Path(os.environ.get('UPROBE_REFERENCE_DATA', '~/Projects/probe-design-bench/resources/data')).expanduser()


@needs_bowtie
def test_reference_mouse_gapdh_probes_by_gene_and_by_transcript(tmp_path):
    from uprobe.core.process import post_process
    from uprobe.templates import load_template
    transcriptome, salmon = REFERENCE / 'mouse_grcm39_ens113_transcriptome', REFERENCE / 'mouse_nih3t3_ens113_salmon'
    needed = [transcriptome / 'records.tsv', transcriptome / 'sequences.fa', transcriptome / 'bowtie' / 'index.1.ebwt',
              salmon / 'expression.tsv', salmon / 'expression_transcripts.tsv']
    if not all(path.is_file() for path in needed):
        pytest.skip(f'reference data not available: {REFERENCE}')
    records = pd.read_csv(transcriptome / 'records.tsv', sep='\t')
    # transcript.fa records as U-Probe names them from the Ensembl GTF: {gene_id}_{transcript_id}, unversioned.
    record_of = {t: f'{g.split(".")[0]}_{t.split(".")[0]}' for t, g in zip(records.transcript_id, records.gene_id)}
    canonical = records[(records.gene_symbol == 'Gapdh') & (records.ensembl_canonical == 1)].transcript_id.iloc[0]
    sequences = {}
    with open(transcriptome / 'sequences.fa') as f:
        name = None
        for line in f:
            if line.startswith('>'):
                name = line[1:].split(None, 1)[0]
            elif name == canonical:
                sequences[name] = sequences.get(name, '') + line.strip()
    mrna = sequences[canonical]
    # The smFISH template: 20-mers along the reference mRNA, every site within 2 substitutions.
    windows = pd.DataFrame([('Gapdh', mrna[i:i + 20], i + 1, i + 20) for i in range(len(mrna) - 19)],
                           columns=['target', 'target_region', 'target_start', 'target_end'])
    unique = list(dict.fromkeys(windows.target_region))
    found = attributes.bowtie_transcript_hits(str(tmp_path), {str(i): s for i, s in enumerate(unique)}, 'q',
                                              str(transcriptome / 'bowtie' / 'index'), 2, 4, 'sense')
    transcript_hits = {s: {record_of[r]: m for r, m in found[str(i)].items()} for i, s in enumerate(unique)}
    # The hit genes with all their transcripts: the GTF and transcript.fa (for exempt_copies) they need.
    genes = {record.split('_')[0] for hit in transcript_hits.values() for record in hit}
    subset = records[records.gene_id.str.split('.').str[0].isin(genes)]
    (tmp_path / 'genome.gtf').write_text(''.join(
        f'1\tt\texon\t1\t9\t.\t+\t.\tgene_id "{g.split(".")[0]}"; transcript_id "{t.split(".")[0]}"; '
        f'gene_name "{name}";\n' for t, g, name in zip(subset.transcript_id, subset.gene_id, subset.gene_symbol)))
    wanted = {record_of[t] for t in subset.transcript_id}
    with open(transcriptome / 'sequences.fa') as f, open(tmp_path / 'transcript.fa', 'w') as out:
        keep = False
        for line in f:
            if line.startswith('>'):
                keep = record_of[line[1:].split(None, 1)[0]] in wanted
                if keep:
                    out.write(f'>{record_of[line[1:].split(None, 1)[0]]}\n')
            elif keep:
                out.write(line)
    template = load_template('smfish_stellaris')
    settings = {**SETTINGS, 'exempt_copies': template['attributes']['target_mapped_genes']['exempt_copies']}
    recommended = {'value_column': 'tpm', 'min_value': 1, 'min_fraction_of_target': 0.05}
    tables = {'gene': {'table': str(salmon / 'expression.tsv'), 'gene_column': 'gene_id', **recommended},
              'transcript': {'table': str(salmon / 'expression_transcripts.tsv'), 'level': 'transcript',
                             **recommended}}
    results = {}
    for level, block in tables.items():
        df, report = windows.copy(), {}
        counter = OffTargetCounter({'gene_expression': block}, {'fasta': str(tmp_path / 'g.fa'),
                                                                'gtf': str(tmp_path / 'genome.gtf')}, report)
        counter.annotate(df, 'target_mapped_genes', 'target_region',
                         {s: gene_hits_of(r) for s, r in transcript_hits.items()}, settings, transcript_hits)
        final = post_process(df, {'post_process': template['post_process']})
        results[level] = (int((df.target_mapped_genes <= 1).sum()), len(final), df, report)
    (gene_windows, gene_probes, by_gene, _), (tx_windows, tx_probes, by_tx, report) = results.values()
    freed = pd.Series([label for a, b in zip(by_gene.target_mapped_genes_off_targets,
                                             by_tx.target_mapped_genes_off_targets)
                       for label in set(a.split(', ')) - set(b.split(', ')) if label]).str.split(':').str[0]
    print(f'\nGapdh ({len(windows)} windows): gene table {gene_windows} windows without an expressed off-target gene, '
          f'{gene_probes} probes; transcript table {tx_windows} windows, {tx_probes} probes; '
          f'genes no longer counted (windows): {freed.value_counts().head(5).to_dict()}')
    assert report['gene_expression']['level'] == 'transcript'
    assert report['gene_expression']['annotation_transcripts_in_table'] == report['gene_expression'][
        'annotation_transcripts']
    # Both levels come from one quantification (gene TPM = the sum over its transcripts), so the
    # transcript table only ever counts fewer genes, and here frees most of the mRNA.
    assert report['gene_expression']['targets']['Gapdh']['value'] == pytest.approx(
        pd.read_csv(salmon / 'expression.tsv', sep='\t').set_index('gene_symbol').tpm['Gapdh'], abs=0.01)
    assert (by_tx.target_mapped_genes <= by_gene.target_mapped_genes).all()
    # 47, not the 48 requested: the template filters GC to 40-60% and keeps 3 nt between sites
    # (48 with its earlier 2-nt spacing and no GC filter).
    assert tx_windows > 1.5 * gene_windows and tx_probes == 47 > gene_probes
