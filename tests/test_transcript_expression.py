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
