"""uprobe genome add --transcriptome ... --records ...: a transcriptome FASTA and a transcript table as a genome entry."""
import random
import shutil

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from uprobe.core.api import UProbeAPI
from uprobe.core.cli import cli
from uprobe.core.gen.fun import generate_target_seqs, read_gtf
from uprobe.core.genome_paths import index_prefix, index_root
from uprobe.core.tools import build_genome
from uprobe.core.transcriptome import TranscriptomeError, add_transcriptome, read_records
from uprobe.core.validate import validate_protocol

needs_bowtie2 = pytest.mark.skipif(shutil.which('bowtie2-build') is None, reason='bowtie2 is not installed')

# transcript, gene name, gene id, MANE Select, canonical, CDS start, CDS end
ROWS = [
    ('ENST01.3', 'AAA', 'ENSG01', '1', '1', 31, 150),
    ('ENST02.1', 'AAA', 'ENSG01', '0', '0', 31, 90),
    ('ENST03.2', 'BBB', 'ENSG02', 'MANE Select', '1', '', ''),
    ('ENST04.1', 'CC_C', 'ENSG_03', '', '1', '', ''),
]


@pytest.fixture
def inputs(tmp_path):
    rng = random.Random(11)
    seqs = {row[0]: ''.join(rng.choice('ACGT') for _ in range(200 + 20 * i)) for i, row in enumerate(ROWS)}
    seqs['ENST99.1'] = ''.join(rng.choice('ACGT') for _ in range(120))  # not in the table
    fasta = tmp_path / 'tx.fa'
    fasta.write_text(''.join(f'>{name} cdna gene:x\n{seq[:70]}\n{seq[70:]}\n' for name, seq in seqs.items()))
    table = tmp_path / 'tx.tsv'
    table.write_text('transcript_id\tgene_name\tgene_id\tMANE_Select\tcanonical\tcds_start\tcds_end\n' +
                     ''.join('\t'.join(map(str, row)) + '\n' for row in ROWS))
    return fasta, table, seqs


def add(tmp_path, inputs, **options):
    fasta, table, _ = inputs
    return add_transcriptome('tx', fasta, table, out_dir=tmp_path / 'genomes', genomes_yaml=tmp_path / 'genomes.yaml',
                             align_index=['bowtie'], **options)


def test_writes_a_single_exon_annotation_with_tags(tmp_path, inputs):
    result = add(tmp_path, inputs)
    _, _, seqs = inputs
    directory = tmp_path / 'genomes' / 'tx'
    gtf = read_gtf(str(directory / 'annotation.gtf'), extract_fields=['gene_id', 'transcript_id', 'gene_name'])
    exons = gtf[gtf['type'] == 'exon'].set_index('transcript_id')
    assert sorted(exons.index) == ['ENST01.3', 'ENST02.1', 'ENST03.2', 'ENST04.1']  # ENST99.1 is not listed
    assert exons.loc['ENST03.2', 'end'] == len(seqs['ENST03.2']) and set(exons['strand']) == {'+'}
    assert exons.loc['ENST01.3', 'gene_name'] == 'AAA' and exons.loc['ENST04.1', 'gene_name'] == 'CC_C'
    # MANE from a 1 or a value containing "Select"; canonical from 1; CDS lines only where the table gives them.
    info = gtf.groupby('transcript_id')['info'].first()
    assert 'tag "MANE_Select"' in info['ENST01.3'] and 'tag "Ensembl_canonical"' in info['ENST01.3']
    assert 'MANE_Select' in info['ENST03.2'] and 'MANE_Select' not in info['ENST04.1']
    assert 'Ensembl_canonical' not in info['ENST02.1']
    assert gtf[gtf['type'] == 'CDS'].set_index('transcript_id')[['start', 'end']].to_dict('index') == {
        'ENST01.3': {'start': 31, 'end': 150}, 'ENST02.1': {'start': 31, 'end': 90}}
    # transcript.fa carries the <gene_id>_<transcript_id> names the off-target search splits at the first '_'.
    names = [line[1:].strip() for line in (directory / 'transcript.fa').read_text().splitlines() if line[0] == '>']
    assert names == ['ENSG01_ENST01.3', 'ENSG01_ENST02.1', 'ENSG02_ENST03.2', 'ENSG-03_ENST04.1']
    assert result['summary'] == {'fasta_records': 5, 'transcripts': 4, 'genes': 3, 'unlisted': 1, 'table_only': 0,
                                 'gene_ids_changed': 1, 'with_cds': 2}
    sequences = dict(zip(names, [l for l in (directory / 'transcript.fa').read_text().splitlines() if l[0] != '>']))
    assert sequences['ENSG02_ENST03.2'] == seqs['ENST03.2']


def test_entry_keeps_derived_files_out_of_the_fasta_directory(tmp_path, inputs):
    result = add(tmp_path, inputs)
    entry = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())['tx']
    assert entry == result['entry']
    directory = tmp_path / 'genomes' / 'tx'
    assert entry['fasta'] == str(directory / 'transcriptome.fa') and entry['gtf'] == str(directory / 'annotation.gtf')
    assert entry['index_dir'] == str(directory) and entry['align_index'] == ['bowtie'] and entry['jellyfish'] is False
    assert index_root(entry) == directory
    assert index_prefix(entry, 'bowtie', 'transcript') == directory / 'bowtie_transcript' / 'transcriptome'
    # The FASTA the user gave is only linked: nothing is written next to it.
    assert sorted(p.name for p in inputs[0].parent.iterdir() if p.name.startswith('tx.')) == ['tx.fa', 'tx.tsv']
    elsewhere = add(tmp_path, inputs, index_dir=tmp_path / 'fast')
    assert elsewhere['entry']['index_dir'] == str(tmp_path / 'fast') and (tmp_path / 'fast' / 'transcript.fa').is_file()


def test_unlisted_records_and_versions(tmp_path, inputs):
    included = add(tmp_path, inputs, unlisted='include')
    assert included['summary']['transcripts'] == 5 and included['summary']['genes'] == 4
    gtf = read_gtf(str(tmp_path / 'genomes' / 'tx' / 'annotation.gtf'), extract_fields=['gene_name', 'transcript_id'])
    assert gtf.set_index('transcript_id').loc['ENST99.1', 'gene_name'].unique().tolist() == ['ENST99.1']
    # A table without versions still matches versioned FASTA records.
    fasta, table, _ = inputs
    table.write_text(table.read_text().replace('.3\t', '\t').replace('.1\t', '\t').replace('.2\t', '\t'))
    assert add(tmp_path, inputs)['summary']['transcripts'] == 4
    with pytest.raises(TranscriptomeError, match='none of the 5 FASTA records'):
        table.write_text('transcript_id\tgene_name\nOTHER\tX\n')
        add(tmp_path, inputs)


def test_column_detection_and_errors(tmp_path):
    table = tmp_path / 't.tsv'
    table.write_text('tx\tsymbol\ttags\nT1\tG1\tMANE_Select,basic\nT2\tG2\t\n')
    with pytest.raises(TranscriptomeError, match='no transcript id column'):
        read_records(table)
    records = read_records(table, transcript_column='tx')
    assert records['T1']['tags'] == {'MANE_Select', 'basic'} and records['T1']['gene_id'] == 'G1'
    with pytest.raises(TranscriptomeError, match="column 'nope' is not in the table"):
        read_records(table, transcript_column='nope')
    table.write_text('transcript_id\tx\nT1\t1\n')
    with pytest.raises(TranscriptomeError, match='no gene column'):
        read_records(table)


def test_the_entry_serves_gene_targets_and_validates(tmp_path, inputs, monkeypatch):
    add(tmp_path, inputs)
    monkeypatch.chdir(tmp_path)
    entry = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())['tx']
    mane = generate_target_seqs('mRNA', ['AAA', 'BBB'], entry['fasta'], entry['gtf'], min_length=30, step=30,
                                transcripts='mane')
    assert mane.attrs['transcripts']['AAA'] == {'selection': 'mane', 'chosen_by': 'MANE_Select', 'ids': ['ENST01.3']}
    assert set(mane['target']) == {'AAA', 'BBB'}
    seqs = inputs[2]
    first = mane[mane['target'] == 'AAA'].iloc[0]['target_region']
    assert first == seqs['ENST01.3'][:30]
    protocol = {'name': 'p', 'genome': 'tx', 'targets': ['AAA', 'BBB'], 'encoding': {},
                'extracts': {'target_region': {'source': 'mRNA', 'length': 30, 'transcripts': 'mane'}},
                'probes': {'p1': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}},
                'attributes': {'gc': {'target': 'target_region', 'type': 'gc_content'}}}
    (tmp_path / 'p.yaml').write_text(yaml.safe_dump(protocol))
    result = CliRunner().invoke(cli, ['validate', '-p', 'p.yaml', '-g', 'genomes.yaml'])
    assert result.exit_code == 0, result.output
    assert 'OK: 0 error(s)' in result.output


def test_cli_genome_add_and_list(tmp_path, inputs, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fasta, table, _ = inputs
    runner = CliRunner()
    dry = runner.invoke(cli, ['genome', 'add', 'tx', '--transcriptome', str(fasta), '--records', str(table),
                              '--align-index', 'bowtie', '--dry-run'])
    assert dry.exit_code == 0 and 'Dry run' in dry.output and not (tmp_path / 'genomes.yaml').exists()
    done = runner.invoke(cli, ['genome', 'add', 'tx', '--transcriptome', str(fasta), '--records', str(table),
                               '--align-index', 'bowtie', '--existing-index', 'bowtie_transcript=/nowhere/prefix'])
    assert done.exit_code == 0, done.output
    entry = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())['tx']
    assert entry['indexes'] == {'bowtie_transcript': '/nowhere/prefix'}
    bad = runner.invoke(cli, ['genome', 'add', 'tx', '--transcriptome', str(fasta), '--records', str(table),
                              '--existing-index', 'oops'])
    assert bad.exit_code != 0 and 'ALIGNER_KIND=PREFIX' in bad.output
    listing = runner.invoke(cli, ['genome', 'list'])
    assert listing.exit_code == 0 and 'bowtie index: transcript missing' in listing.output


# -- index_dir and existing indexes ---------------------------------------------------------

@needs_bowtie2
def test_index_dir_and_existing_index_pointer(tmp_path, inputs):
    add(tmp_path, inputs, index_dir=tmp_path / 'idx')
    entry = {**yaml.safe_load((tmp_path / 'genomes.yaml').read_text())['tx'], 'align_index': ['bowtie2']}
    build_genome(entry, threads=1)
    prefix = tmp_path / 'idx' / 'bowtie2_transcript' / 'transcriptome'
    assert (tmp_path / 'idx' / 'bowtie2_transcript' / 'transcriptome.1.bt2').is_file()
    assert not (tmp_path / 'genomes' / 'tx' / 'bowtie2_transcript').exists()
    # Point at that index: it is used, nothing is built next to the FASTA, and the validator accepts it.
    pointed = {**entry, 'index_dir': str(tmp_path / 'elsewhere'), 'indexes': {'bowtie2_transcript': str(prefix)}}
    build_genome(pointed, threads=1)
    assert not (tmp_path / 'elsewhere' / 'bowtie2_transcript').exists()
    assert index_prefix(pointed, 'bowtie2', 'transcript') == prefix
    # A pointer to nothing is an error, not a silent build.
    broken = {**pointed, 'indexes': {'bowtie2_transcript': str(tmp_path / 'nope' / 'prefix')}}
    with pytest.raises(ValueError, match='points to'):
        build_genome(broken, threads=1)
    assert not (tmp_path / 'nope').exists()


def test_validate_checks_index_dir_and_pointers(tmp_path, inputs):
    from uprobe.core.validate import validate_protocol
    add(tmp_path, inputs)
    genomes = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())
    protocol = {'name': 'p', 'genome': 'tx', 'targets': ['AAA'], 'encoding': {},
                'extracts': {'target_region': {'source': 'mRNA', 'length': 30}},
                'probes': {'p1': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}},
                'attributes': {'gc': {'target': 'target_region', 'type': 'gc_content'}}}

    def messages(**extra):
        genomes['tx'].update(extra)
        return [(f.path, f.message) for f in validate_protocol(protocol, genomes=genomes) if f.severity == 'error']
    assert messages() == []
    assert 'index_dir must be a directory path' in messages(index_dir=3)[0][1]
    assert 'unknown index' in messages(index_dir=str(tmp_path), indexes={'bowtie_exon': '/x'})[0][1]
    assert 'no complete bowtie index' in messages(indexes={'bowtie_transcript': str(tmp_path / 'none')})[0][1]
    assert 'must map' in messages(indexes=['x'])[0][1]


# -- an existing index must name its records <gene_id>_<transcript_id> --------------------

def fake_index(prefix, aligner='bowtie'):
    from uprobe.core.tools import _BOWTIE_SMALL_SUFFIXES, _BOWTIE2_SMALL_SUFFIXES
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (_BOWTIE_SMALL_SUFFIXES if aligner == 'bowtie' else _BOWTIE2_SMALL_SUFFIXES):
        prefix.with_name(prefix.name + suffix).write_text('index')
    return prefix


PLAIN = ['ENST01.3', 'ENST02.1', 'ENST03.2', 'NM_000003.1']  # a transcript FASTA's own names
NAMED = ['ENSG01_ENST01.3', 'ENSG01_ENST02.1', 'ENSG02_ENST03.2', 'ENSG-03_ENST04.1']  # U-Probe's transcript.fa


def test_existing_index_with_other_record_names_is_refused(tmp_path, inputs, monkeypatch):
    import uprobe.core.offtarget as offtarget
    prefix = fake_index(tmp_path / 'ext' / 'tx')
    monkeypatch.setattr(offtarget, 'index_record_names', lambda aligner, where: PLAIN)
    # genome add: refused before anything is written.
    with pytest.raises(TranscriptomeError, match=r"4 of 4 records are not named <gene_id>_<transcript_id>.*"
                                                 r"ENST01.3 -> gene 'ENST01.3'.*NM_000003.1 -> gene 'NM'"):
        add(tmp_path, inputs, existing_indexes={'bowtie_transcript': str(prefix)})
    assert not (tmp_path / 'genomes.yaml').exists()
    assert not list((tmp_path / 'genomes' / 'tx').glob('*.tmp'))
    # The same index named in genomes.yaml by hand: validate and the run refuse it.
    add(tmp_path, inputs)
    genomes = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())
    genomes['tx']['indexes'] = {'bowtie_transcript': str(prefix)}
    protocol = {'name': 'p', 'genome': 'tx', 'targets': ['AAA'], 'encoding': {},
                'extracts': {'target_region': {'source': 'mRNA', 'length': 30}},
                'probes': {'p1': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}},
                'attributes': {'mg': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie'}}}
    errors = [f.message for f in validate_protocol(protocol, genomes=genomes) if f.severity == 'error']
    assert any('records are not named <gene_id>_<transcript_id>' in message for message in errors), errors
    with pytest.raises(ValueError, match='records are not named'):
        build_genome(genomes['tx'], threads=1, used={'genome': set(), 'transcript': {'bowtie'}})
    # Records named after the annotation are accepted.
    monkeypatch.setattr(offtarget, 'index_record_names', lambda aligner, where: NAMED)
    assert not [f for f in validate_protocol(protocol, genomes=genomes) if f.severity == 'error' and 'records' in f.message]
    build_genome(genomes['tx'], threads=1, used={'genome': set(), 'transcript': {'bowtie'}})
    assert add(tmp_path, inputs, existing_indexes={'bowtie_transcript': str(prefix)})['entry']['indexes']


def test_unchecked_record_names_are_a_warning(tmp_path, inputs, monkeypatch):
    import uprobe.core.offtarget as offtarget
    prefix = fake_index(tmp_path / 'ext' / 'tx')
    monkeypatch.setattr(offtarget, 'index_record_names', lambda aligner, where: None)  # no bowtie-inspect
    add(tmp_path, inputs, existing_indexes={'bowtie_transcript': str(prefix)})
    genomes = yaml.safe_load((tmp_path / 'genomes.yaml').read_text())
    protocol = {'name': 'p', 'genome': 'tx', 'targets': ['AAA'], 'encoding': {},
                'extracts': {'target_region': {'source': 'mRNA', 'length': 30}},
                'probes': {'p1': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}},
                'attributes': {'gc': {'target': 'target_region', 'type': 'gc_content'}}}
    warnings = [f.message for f in validate_protocol(protocol, genomes=genomes) if f.severity == 'warning']
    assert any('bowtie-inspect is not on PATH' in message for message in warnings), warnings


@needs_bowtie2
def test_plain_transcript_index_is_refused_with_bowtie2(tmp_path, inputs):
    import subprocess
    fasta = inputs[0]
    subprocess.run(['bowtie2-build', '-q', str(fasta), str(tmp_path / 'plain')], check=True)
    with pytest.raises(TranscriptomeError, match="ENST01.3 -> gene 'ENST01.3'"):
        add(tmp_path, inputs, existing_indexes={'bowtie2_transcript': str(tmp_path / 'plain')})
    add(tmp_path, inputs)
    built = tmp_path / 'genomes' / 'tx' / 'transcript.fa'
    subprocess.run(['bowtie2-build', '-q', str(built), str(tmp_path / 'named')], check=True)
    entry = add(tmp_path, inputs, existing_indexes={'bowtie2_transcript': str(tmp_path / 'named')})['entry']
    build_genome({**entry, 'align_index': ['bowtie2']}, threads=1)


# -- genome add again: indexes of the previous sequences are removed --------------------------

def test_re_adding_changed_sequences_removes_the_indexes_built_from_them(tmp_path, inputs):
    directory = tmp_path / 'genomes' / 'tx'
    add(tmp_path, inputs)
    transcript_index = fake_index(directory / 'bowtie_transcript' / 'transcriptome')
    kmers = directory / 'jf_genome' / 'transcriptome.k35.jf'
    kmers.parent.mkdir()
    kmers.write_text('index')
    # The same inputs again: transcript.fa is unchanged, its indexes stay.
    assert add(tmp_path, inputs)['removed_indexes'] == []
    assert transcript_index.with_name('transcriptome.1.ebwt').is_file() and kmers.is_file()
    # A changed gene id renames the transcript.fa records: the bowtie index of the old names goes.
    fasta, table, _ = inputs
    table.write_text(table.read_text().replace('ENSG01', 'ENSG09'))
    removed = add(tmp_path, inputs)['removed_indexes']
    assert removed == [str(directory / 'bowtie_transcript')]
    assert not (directory / 'bowtie_transcript').exists() and kmers.is_file()  # the FASTA is the same
    # Another FASTA: the indexes of the old one (jf_genome) go too.
    fake_index(directory / 'bowtie_transcript' / 'transcriptome')
    other = tmp_path / 'other.fa'
    other.write_text(fasta.read_text().replace('>ENST02.1 cdna gene:x\n', '>ENST02.1 cdna gene:x\nA'))
    result = add_transcriptome('tx', other, table, out_dir=tmp_path / 'genomes', genomes_yaml=tmp_path / 'genomes.yaml',
                               align_index=['bowtie'])
    assert sorted(result['removed_indexes']) == [str(directory / 'bowtie_transcript'), str(directory / 'jf_genome')]


@needs_bowtie2
def test_re_added_transcriptome_gets_an_index_of_its_new_records(tmp_path, inputs):
    from uprobe.core.offtarget import index_record_names
    add(tmp_path, inputs)
    entry = {**yaml.safe_load((tmp_path / 'genomes.yaml').read_text())['tx'], 'align_index': ['bowtie2']}
    build_genome(entry, threads=1)
    prefix = index_prefix(entry, 'bowtie2', 'transcript')
    assert 'ENSG01_ENST01.3' in index_record_names('bowtie2', prefix)
    _, table, _ = inputs
    table.write_text(table.read_text().replace('ENSG01', 'ENSG09'))
    add(tmp_path, inputs)
    build_genome(entry, threads=1)
    names = index_record_names('bowtie2', prefix)
    assert 'ENSG09_ENST01.3' in names and not any(name.startswith('ENSG01_') for name in names)
