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
