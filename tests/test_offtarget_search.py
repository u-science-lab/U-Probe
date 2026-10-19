"""Exhaustive off-target search: mapped_genes with aligner: bowtie (Bowtie 1)."""
import random
import shutil
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest
import yaml

from uprobe.core import attributes
from uprobe.core.attributes import add_attributes, bowtie_gene_hits
from uprobe.core.offtarget import BOWTIE2_LOCAL_MIN_LENGTH, mapped_genes_settings
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import validate_protocol

DATA = Path(__file__).parent / 'data'
needs_bowtie = pytest.mark.skipif(not (shutil.which('bowtie') and shutil.which('bowtie-build')),
                                  reason='needs Bowtie 1 (bowtie, bowtie-build)')


def mutate(seq, positions, rng):
    seq = list(seq)
    for i in positions:
        seq[i] = rng.choice([b for b in 'ACGT' if b != seq[i]])
    return ''.join(seq)


def brute_force(transcripts, query, max_mismatches, strand):
    """gene -> fewest mismatches over every window of every transcript, by plain Hamming distance."""
    sites = {'both': [query, reverse_complement(query)], 'sense': [query],
             'antisense': [reverse_complement(query)]}[strand]
    genes = {}
    for name, seq in transcripts.items():
        gene = name.split('_')[0]
        for site in sites:
            for i in range(len(seq) - len(site) + 1):
                mismatches = sum(a != b for a, b in zip(seq[i:i + len(site)], site))
                if mismatches <= max_mismatches and mismatches < genes.get(gene, max_mismatches + 1):
                    genes[gene] = mismatches
    return genes


@pytest.fixture(scope='module')
def toy_transcriptome(tmp_path_factory):
    """Random transcripts with planted sites of 0-3 substitutions on both strands, and its Bowtie 1 index."""
    rng = random.Random(7)
    transcripts = {f'G{g}_T{g}{t}': ''.join(rng.choice('ACGT') for _ in range(400))
                   for g in range(12) for t in range(2)}
    names = sorted(transcripts)
    queries = {}
    for q in range(30):
        source = names[q % len(names)]
        length = 18 + q % 5
        start = rng.randrange(0, 400 - length)
        query = transcripts[source][start:start + length]
        queries[str(q)] = query
        for k in range(4):  # copies with k substitutions, sense or antisense, in other transcripts
            name = rng.choice(names)
            copy = mutate(query, rng.sample(range(length), k), rng)
            if rng.random() < 0.5:
                copy = reverse_complement(copy)
            at = rng.randrange(0, 400 - length)
            transcripts[name] = transcripts[name][:at] + copy + transcripts[name][at + length:]
    root = tmp_path_factory.mktemp('toy')
    (root / 'transcript.fa').write_text(''.join(f'>{n}\n{s}\n' for n, s in transcripts.items()))
    return root, transcripts, queries


@pytest.fixture(scope='module')
def bowtie_index(toy_transcriptome):
    from uprobe.core.tools import build_transcripts_index
    root, _, _ = toy_transcriptome
    return build_transcripts_index(root / 'toy.gtf', root / 'toy.fa', root / 'bowtie_transcript',
                                   threads=1, aligner='bowtie')


@needs_bowtie
@pytest.mark.parametrize('strand', ['both', 'sense', 'antisense'])
@pytest.mark.parametrize('max_mismatches', [0, 2, 3])
def test_bowtie_finds_every_site_brute_force_finds(toy_transcriptome, bowtie_index, tmp_path, strand, max_mismatches):
    _, transcripts, queries = toy_transcriptome
    hits = bowtie_gene_hits(str(tmp_path), queries, 'q', bowtie_index, max_mismatches, 1, strand)
    for key, query in queries.items():
        assert hits[key] == brute_force(transcripts, query, max_mismatches, strand), (key, query)


@needs_bowtie
def test_build_index_creates_only_a_transcript_index_for_bowtie(toy_transcriptome, bowtie_index, monkeypatch):
    from uprobe.core import tools
    root, _, _ = toy_transcriptome
    assert tools._bowtie_index_exists(Path(bowtie_index))
    assert Path(bowtie_index) == root / 'bowtie_transcript' / 'toy'
    build = Mock()
    monkeypatch.setattr(tools, 'build_bowtie_index', build)
    tools.build_genome({'fasta': str(root / 'toy.fa'), 'gtf': str(root / 'toy.gtf'), 'align_index': ['bowtie']})
    build.assert_not_called()  # complete index reused
    assert not (root / 'bowtie_genome').exists()


@needs_bowtie
def test_mapped_genes_attribute_with_bowtie(toy_transcriptome, bowtie_index, tmp_path, monkeypatch):
    root, transcripts, queries = toy_transcriptome
    monkeypatch.chdir(tmp_path)
    df = pd.DataFrame({'target_region': list(queries.values()) + [queries['0']]})
    protocol = {'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                                        'max_mismatches': 2, 'strand': 'sense'}}}
    genome = {'fasta': str(root / 'toy.fa'), 'align_index': ['bowtie']}
    result = add_attributes(df, protocol, genome)
    expected = [len(brute_force(transcripts, q, 2, 'sense')) for q in df.target_region]
    assert result.hits.tolist() == expected


def test_build_genome_indexes_only_the_transcriptome_for_bowtie(monkeypatch, tmp_path):
    from uprobe.core import tools
    def build(source, prefix, threads):
        for suffix in tools._BOWTIE_SMALL_SUFFIXES:
            Path(str(prefix) + suffix).write_text('index')
    build = Mock(side_effect=build)
    monkeypatch.setattr(tools, 'build_bowtie_index', build)
    monkeypatch.setattr(tools, 'build_bowtie2_index', Mock())
    (tmp_path / 'transcript.fa').write_text('>g_t\nACGT\n')
    genome = {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(tmp_path / 'g.gtf'), 'align_index': ['bowtie']}
    tools.build_genome(genome, threads=3)
    build.assert_called_once_with(tmp_path / 'transcript.fa', tmp_path / 'bowtie_transcript' / 'g', 3)
    tools.build_genome(genome)
    build.assert_called_once()  # a complete index is kept
    tools.build_bowtie2_index.assert_not_called()
    assert not (tmp_path / 'bowtie_genome').exists()
    with pytest.raises(ValueError, match='must use bowtie_transcript'):
        tools.build_transcripts_index(tmp_path / 'g.gtf', tmp_path / 'g.fa', tmp_path / 'bowtie2_transcript',
                                      aligner='bowtie')


def test_bowtie2_strand_options(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(attributes._attributes.subp, 'check_call', lambda cmd, shell: calls.append(cmd))
    for strand in ('both', 'sense', 'antisense'):
        attributes.bowtie2_align_se_sen('q.fq', 'idx', str(tmp_path / f'{strand}.sam'), strand=strand)
    assert ['--norc' in c or '--nofw' in c for c in calls] == [False, True, True]
    assert '--norc' in calls[1] and '--nofw' in calls[2]


def test_short_sequences_warn_with_bowtie2(monkeypatch, tmp_path, caplog):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(attributes, 'bowtie2_gene_hits', lambda _, records, *args: {k: {} for k in records})
    warnings = []
    monkeypatch.setattr(attributes.log, 'warning', warnings.append)
    df = pd.DataFrame({'s': ['A' * 20, 'C' * BOWTIE2_LOCAL_MIN_LENGTH]})
    add_attributes(df, {'attributes': {'hits': {'target': 's', 'type': 'mapped_genes', 'aligner': 'bowtie2'}}},
                   {'fasta': str(tmp_path / 'ref.fa'), 'align_index': ['bowtie2']})
    assert len(warnings) == 1 and '1 of 2 sequences are shorter than 22 nt' in warnings[0]


@pytest.mark.parametrize('attr,message', [
    ({'aligner': 'bowtie', 'max_mismatches': 4}, 'max_mismatches must be an integer from 0 to 3'),
    ({'aligner': 'bowtie', 'max_mismatches': 1.5}, 'max_mismatches must be an integer'),
    ({'aligner': 'bowtie2', 'strand': 'plus'}, 'strand must be one of both, sense, antisense'),
    ({'aligner': 'blast'}, 'aligner must be one of bowtie2, bowtie'),
    ({'aligner': 'bowtie', 'exempt_copies': {'kmer': 40}}, r'attributes.hits.exempt_copies.kmer must be'),
])
def test_mapped_genes_settings_reject(attr, message):
    with pytest.raises(ValueError, match=message):
        mapped_genes_settings(attr, 'attributes.hits')


def test_mapped_genes_settings_defaults():
    assert mapped_genes_settings({'aligner': 'bowtie'}) == {'aligner': 'bowtie', 'strand': 'both', 'max_mismatches': 2,
                                                             'exempt_copies': None}
    assert mapped_genes_settings({'aligner': 'bowtie2'})['max_mismatches'] is None
    assert mapped_genes_settings({'aligner': 'bowtie2', 'exempt_copies': True})['exempt_copies'] == {
        'kmer': 30, 'min_shared': 0.5}


# -- uprobe validate -----------------------------------------------------

@pytest.fixture
def rna():
    protocol = yaml.safe_load((DATA / 'RNA_format.yaml').read_text(encoding='utf-8'))
    del protocol['attributes']['target_mappedGenes']['min_mapq']
    return protocol


@pytest.fixture
def genomes(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>1\n' + 'ACGT' * 100 + '\n')
    gtf = tmp_path / 'ref.gtf'
    gtf.write_text('1\ttest\texon\t1\t200\t.\t+\t.\tgene_id "ID1"; gene_name "xxx"; transcript_id "T1";\n')
    return {'test': {'fasta': str(fasta), 'gtf': str(gtf), 'align_index': ['bowtie2'], 'jellyfish': False}}


def paths(findings, severity):
    return {f.path: f for f in findings if f.severity == severity}


def test_validate_accepts_bowtie_and_checks_its_options(rna):
    attr = rna['attributes']['target_mappedGenes']
    attr.update(aligner='bowtie', max_mismatches=2, strand='sense')
    assert validate_protocol(rna) == []
    attr.update(max_mismatches=5, strand='minus')
    errors = paths(validate_protocol(rna), 'error')
    assert 'max_mismatches' in errors['attributes.target_mappedGenes.max_mismatches'].message
    assert errors['attributes.target_mappedGenes.strand'].suggestion.startswith('strand: sense')
    attr.update(aligner='bowtie2', max_mismatches=2, strand='sense')
    warnings = paths(validate_protocol(rna), 'warning')
    assert 'only applies to aligner: bowtie' in warnings['attributes.target_mappedGenes.max_mismatches'].message
    attr['aligner'] = 'bwa'
    assert 'bowtie2 or bowtie' in paths(validate_protocol(rna), 'error')['attributes.target_mappedGenes.aligner'].message


def test_validate_requires_bowtie_in_align_index(rna, genomes):
    rna['attributes']['target_mappedGenes']['aligner'] = 'bowtie'
    finding = paths(validate_protocol(rna, genomes), 'error')['attributes.target_mappedGenes']
    assert 'needs bowtie' in finding.message and 'add bowtie to align_index' in finding.suggestion
    genomes['test']['align_index'] = ['bowtie2', 'bowtie']
    assert 'attributes.target_mappedGenes' not in paths(validate_protocol(rna, genomes), 'error')


def test_validate_warns_about_short_sequences_with_bowtie2(rna):
    rna['extracts']['target_region'].update(length=20, overlap=None, step=10)
    rna['probes'] = {'probe_1': {'template': '{b}', 'parts': {'b': {'expr': 'rc(target_region)'}}}}
    rna['attributes'] = {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie2'},
                         'probe_hits': {'target': 'probe_1.b', 'type': 'mapped_genes', 'aligner': 'bowtie2'}}
    rna['post_process'] = {'filters': {'hits': {'condition': 'hits <= 1'}}}
    rna['summary'] = {'report_name': 'rna_report'}
    warnings = paths(validate_protocol(rna), 'warning')
    assert 'cannot align sequences shorter than 22 nt' in warnings['attributes.hits.aligner'].message
    assert 'probe_1.b (20 nt)' in warnings['attributes.probe_hits.aligner'].message
    assert 'strand: sense' in warnings['attributes.hits.aligner'].suggestion
    assert 'conda install -c bioconda bowtie' in warnings['attributes.hits.aligner'].suggestion
    # 22-30 nt align, but Bowtie2 misses many of their off-target sites.
    rna['extracts']['target_region']['length'] = 22
    warning = paths(validate_protocol(rna), 'warning')['attributes.hits.aligner']
    assert 'misses many off-target sites' in warning.message
    rna['extracts']['target_region']['length'] = 31
    assert 'attributes.hits.aligner' not in paths(validate_protocol(rna), 'warning')
    rna['attributes']['hits']['aligner'] = 'bowtie'
    rna['extracts']['target_region']['length'] = 20
    assert 'attributes.hits.aligner' not in paths(validate_protocol(rna), 'warning')


def test_genome_list_shows_transcript_index_only_for_bowtie(tmp_path):
    from click.testing import CliRunner
    from uprobe.core.cli import cli
    fasta = tmp_path / 'g.fa'
    fasta.write_text('>1\nACGT\n')
    (tmp_path / 'g.gtf').write_text('')
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump({'g': {'fasta': str(fasta), 'gtf': str(tmp_path / 'g.gtf'),
                                             'align_index': ['bowtie2', 'bowtie']}}))
    result = CliRunner().invoke(cli, ['genome', 'list', '-g', str(genomes)])
    assert 'bowtie2 index: genome missing, transcript missing' in result.output
    assert 'bowtie index: transcript missing' in result.output


# -- generated attributes (protocols without `attributes`) ----------------

def installed_tools(monkeypatch, tools):
    """Make shutil.which find exactly `tools`."""
    from uprobe.core import offtarget
    monkeypatch.setattr(offtarget.shutil, 'which', lambda name, *args, **kwargs: f'/opt/bin/{name}'
                        if name in tools else None)


def test_bowtie_unavailable_reasons(tmp_path, monkeypatch):
    from uprobe.core import tools
    from uprobe.core.offtarget import bowtie_unavailable
    genome = {'fasta': str(tmp_path / 'g.fa'), 'align_index': ['bowtie2', 'bowtie']}
    installed_tools(monkeypatch, {'bowtie', 'bowtie-build'})
    assert bowtie_unavailable(genome) == ''
    assert bowtie_unavailable(None) == 'the genome configuration is not known'
    assert bowtie_unavailable({**genome, 'align_index': ['bowtie2']}) == "the genome's align_index has no bowtie"
    installed_tools(monkeypatch, {'bowtie'})
    assert 'index is not built and bowtie-build is not on PATH' in bowtie_unavailable(genome)
    (tmp_path / 'bowtie_transcript').mkdir()
    for suffix in tools._BOWTIE_SMALL_SUFFIXES:
        (tmp_path / 'bowtie_transcript' / f'g{suffix}').write_text('index')
    assert bowtie_unavailable(genome) == ''
    installed_tools(monkeypatch, {'bowtie-build'})
    assert bowtie_unavailable(genome) == 'bowtie (Bowtie 1) is not on PATH'


@pytest.fixture
def rna_defaults(rna):
    """RNA_format.yaml with 20-nt windows and no attributes, post_process or summary: all generated."""
    for key in ('attributes', 'post_process', 'summary'):
        del rna[key]
    rna['extracts']['target_region'] = {'source': 'exon', 'length': 20, 'step': 10}
    rna['probes'] = {'probe_1': {'template': '{b}', 'parts': {'b': {'expr': 'rc(target_region)'}}}}
    return rna


def test_generated_off_target_attribute_uses_bowtie_when_the_genome_can(rna_defaults, genomes, monkeypatch):
    import copy
    from uprobe.core import cli as cli_module
    from uprobe.core.cli import _validate_and_normalize_protocol
    warnings = []
    monkeypatch.setattr(cli_module.log, 'warning', warnings.append)
    installed_tools(monkeypatch, {'bowtie', 'bowtie-build'})
    genome = {**genomes['test'], 'align_index': ['bowtie2', 'bowtie']}
    attr = _validate_and_normalize_protocol(copy.deepcopy(rna_defaults), genome)['attributes']['target_mapped_genes']
    assert attr == {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                    'max_mismatches': 2, 'strand': 'sense'}
    assert warnings == []
    # Without bowtie in align_index (or without a genome entry): Bowtie2, with a warning for 20-nt windows.
    for entry in (genomes['test'], None):
        attr = _validate_and_normalize_protocol(copy.deepcopy(rna_defaults), entry)['attributes']['target_mapped_genes']
        assert attr == {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie2'}
    assert len(warnings) == 2 and "align_index has no bowtie" in warnings[0]
    assert 'the 20-nt windows count no off-target gene' in warnings[0]


def test_validate_reports_the_generated_aligner(rna_defaults, genomes, monkeypatch):
    installed_tools(monkeypatch, {'bowtie', 'bowtie-build'})
    short = 'cannot align sequences shorter than 22 nt'
    finding = paths(validate_protocol(rna_defaults, genomes), 'warning')['attributes']
    assert short in finding.message and 'generated target_mapped_genes then uses aligner: bowtie' in finding.suggestion
    genomes['test']['align_index'] = ['bowtie2', 'bowtie']
    findings = validate_protocol(rna_defaults, genomes)
    assert not any(short in f.message for f in findings)
    assert not paths(findings, 'error')


def test_cli_run_generates_attributes_for_its_genome(tmp_path, monkeypatch):
    from click.testing import CliRunner
    from uprobe.core import cli as cli_module
    seen = []

    def normalize(protocol, genome=None):
        seen.append(genome)
        raise ValueError('stop here')
    monkeypatch.setattr(cli_module, '_validate_and_normalize_protocol', normalize)
    entry = {'fasta': 'g.fa', 'gtf': 'g.gtf', 'align_index': ['bowtie']}
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump({'other': {}, 'g': entry}))
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({'genome': 'g'}))
    result = CliRunner().invoke(cli_module.cli, ['run', '-p', str(protocol), '-g', str(genomes),
                                                 '-o', str(tmp_path / 'out')])
    assert result.exit_code == 1 and seen == [entry]


@needs_bowtie
def test_run_without_attributes_searches_off_targets_with_bowtie(tmp_path):
    # The toy genome (3 genes), with only a Bowtie 1 index: the generated
    # target_mapped_genes must use it, exhaustively and on the sense strand.
    import json
    from click.testing import CliRunner
    from uprobe.core.cli import cli
    from test_run_summary import toy_genome
    genomes = toy_genome(tmp_path, align_index=['bowtie'])
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({
        'name': 'toy', 'genome': 'toy', 'encoding': {}, 'targets': ['ToyGeneA', 'ToyGeneC'],
        'extracts': {'target_region': {'source': 'exon', 'length': 20, 'step': 50}},
        'probes': {'probe': {'template': '{b}', 'parts': {'b': {'expr': 'rc(target_region)'}}}}}))
    result = CliRunner().invoke(cli, ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(tmp_path / 'out'),
                                      '--threads', '2', '--json'])
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout)
    stats = summary['off_targets']['attributes']['target_mapped_genes']
    assert (stats['aligner'], stats['max_mismatches'], stats['strand']) == ('bowtie', 2, 'sense')
    table = pd.read_csv(tmp_path / 'out' / summary['output_files']['tsv'], sep='\t')
    assert len(table) > 0 and (table['target_mapped_genes'] >= 1).all()
