"""mapped_genes with `oligos`: each oligo of a multi-oligo unit searched on its own (Bowtie 1)."""
import random
import shutil

import pandas as pd
import pytest

from uprobe.core.attributes import add_attributes, bowtie_unit_hits
from uprobe.core.offtarget import check_oligo_parts, mapped_genes_settings
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import validate_protocol

needs_bowtie = pytest.mark.skipif(not (shutil.which('bowtie') and shutil.which('bowtie-build')),
                                  reason='needs Bowtie 1 (bowtie, bowtie-build)')

rng = random.Random(11)


def random_seq(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(seq, positions):
    seq = list(seq)
    for i in positions:
        seq[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[seq[i]]
    return ''.join(seq)


# A 52-nt HCR-like window: 25-nt left site, 2-nt gap, 25-nt right site.
LEFT, GAP, RIGHT = random_seq(25), random_seq(2), random_seq(25)
WINDOW = LEFT + GAP + RIGHT
FILLER = 300


def transcript(*pieces):
    """Random filler around the pieces, each piece preceded by 40 random nt."""
    return ''.join(random_seq(40) + piece for piece in pieces) + random_seq(FILLER)


# Off-target genes, named after what they hold.
TRANSCRIPTS = {
    # the target itself
    'TARGET_T1': transcript(WINDOW),
    # each half within 2 substitutions, adjacent: 4 over the window, missed by a whole-window search
    'SPLIT2_T1': transcript(mutate(LEFT, [3, 17]) + GAP + mutate(RIGHT, [5, 20])),
    # the left half only
    'LEFTONLY_T1': transcript(LEFT),
    # both halves within 1 substitution, but 200 nt apart on the transcript
    'APART_T1': transcript(mutate(LEFT, [10]), random_seq(160) + mutate(RIGHT, [10])),
    # both halves on one gene, but on two different transcripts
    'TWOTX_T1': transcript(LEFT), 'TWOTX_T2': transcript(RIGHT),
    # each half with 3 substitutions
    'SPLIT3_T1': transcript(mutate(LEFT, [2, 9, 16]) + GAP + mutate(RIGHT, [2, 9, 16])),
    # the window's reverse complement: the probes have its sequence and cannot bind it
    'ANTI_T1': transcript(reverse_complement(WINDOW)),
    # the whole window with 2 substitutions in the left half (found by both searches)
    'WHOLE2_T1': transcript(mutate(LEFT, [1, 23]) + GAP + RIGHT),
}


@pytest.fixture(scope='module')
def index(tmp_path_factory):
    from uprobe.core.tools import build_transcripts_index
    root = tmp_path_factory.mktemp('oligos')
    (root / 'transcript.fa').write_text(''.join(f'>{n}\n{s}\n' for n, s in TRANSCRIPTS.items()))
    prefix = build_transcripts_index(root / 'toy.gtf', root / 'toy.fa', root / 'bowtie_transcript',
                                     threads=1, aligner='bowtie')
    return root, prefix


def genes(records):
    return {record.split('_')[0]: mismatches for record, mismatches in records.items()}


@needs_bowtie
def test_a_gene_counts_when_every_oligo_hits_one_of_its_transcripts(index, tmp_path):
    _, prefix = index
    unit = (LEFT, RIGHT)
    hits = bowtie_unit_hits(str(tmp_path), [unit], 'u', str(prefix), 2, 1, 'sense')[unit]
    # Worst oligo's mismatches: 2 for SPLIT2 (2 + 2), 2 for WHOLE2 (2 + 0), 1 for APART (1 + 1).
    assert genes(hits) == {'TARGET': 0, 'SPLIT2': 2, 'APART': 1, 'WHOLE2': 2}


@needs_bowtie
def test_max_gap_keeps_only_sites_next_to_each_other(index, tmp_path):
    _, prefix = index
    unit = (LEFT, RIGHT)
    hits = bowtie_unit_hits(str(tmp_path), [unit], 'u', str(prefix), 2, 1, 'sense', max_gap=5)[unit]
    assert genes(hits) == {'TARGET': 0, 'SPLIT2': 2, 'WHOLE2': 2}
    hits = bowtie_unit_hits(str(tmp_path), [unit], 'v', str(prefix), 2, 1, 'sense', max_gap=1)[unit]
    assert genes(hits) == {}  # the sites are 2 nt apart everywhere


@needs_bowtie
def test_strand_both_finds_the_antisense_unit_and_strand_antisense_only_it(index, tmp_path):
    _, prefix = index
    unit = (LEFT, RIGHT)
    both = bowtie_unit_hits(str(tmp_path), [unit], 'b', str(prefix), 2, 1, 'both', max_gap=5)[unit]
    assert set(genes(both)) == {'TARGET', 'SPLIT2', 'WHOLE2', 'ANTI'}
    anti = bowtie_unit_hits(str(tmp_path), [unit], 'a', str(prefix), 2, 1, 'antisense', max_gap=5)[unit]
    assert set(genes(anti)) == {'ANTI'}


@needs_bowtie
def test_short_oligos_and_joined_parts(index, tmp_path):
    # A padlock-like oligo (both arm sites joined, 26 nt) and a 13-nt primer site.
    _, prefix = index
    unit = (WINDOW[:26], WINDOW[39:])
    hits = bowtie_unit_hits(str(tmp_path), [unit], 'p', str(prefix), 2, 1, 'sense', max_gap=20)[unit]
    assert {'TARGET', 'WHOLE2'} <= set(genes(hits))
    assert 'LEFTONLY' not in genes(hits)


def layout_protocol(oligos, **options):
    return {
        'name': 'pairs', 'genome': 'test', 'targets': ['TARGET'], 'encoding': {'TARGET': {}},
        'extracts': {'target_region': {'source': 'exon', 'length': 52, 'step': 1, 'layout': {
            'template': '{left}{gap}{right}',
            'parts': {'left': {'length': 25}, 'gap': {'length': 2}, 'right': {'length': 25}}}}},
        'probes': {'odd': {'template': '{b}', 'parts': {'b': {'expr': "rc(target_parts['left'])"}}},
                   'even': {'template': '{b}', 'parts': {'b': {'expr': "rc(target_parts['right'])"}}}},
        'attributes': {'pair_genes': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                                      'max_mismatches': 2, 'strand': 'sense', 'oligos': oligos, **options},
                       'window_genes': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                                        'max_mismatches': 2, 'strand': 'sense'}},
        'post_process': {'filters': {'pair_genes': {'condition': 'pair_genes <= 1'}}},
    }


@needs_bowtie
def test_mapped_genes_attribute_with_oligos(index, tmp_path, monkeypatch):
    root, _ = index
    monkeypatch.chdir(tmp_path)
    other = random_seq(52)  # no transcript holds it
    df = pd.DataFrame({'target_region': [WINDOW, other, WINDOW]})
    for name, (start, end) in {'left': (0, 25), 'gap': (25, 27), 'right': (27, 52)}.items():
        df[f'target_parts.{name}'] = df.target_region.str[start:end]
    genome = {'fasta': str(root / 'toy.fa'), 'align_index': ['bowtie']}
    report = {}
    result = add_attributes(df, layout_protocol(['left', 'right']), genome, report)
    # The whole-window search misses SPLIT2 (4 substitutions over 52 nt) and APART.
    assert result.window_genes.tolist() == [2, 0, 2]
    assert result.pair_genes.tolist() == [4, 0, 4]
    assert set(result.pair_genes_off_targets[0].split(', ')) == {'TARGET:0', 'APART:1', 'SPLIT2:2', 'WHOLE2:2'}
    assert not any(column.startswith('__') for column in result.columns)
    assert report['attributes']['pair_genes']['oligos'] == [['left'], ['right']]
    result = add_attributes(df.drop(columns=['window_genes', 'pair_genes'], errors='ignore'),
                            layout_protocol(['left', 'right'], max_gap=5), genome)
    assert result.pair_genes.tolist() == [3, 0, 3]


def test_settings():
    settings = mapped_genes_settings({'aligner': 'bowtie', 'oligos': [['arm5', 'arm3'], 'primer'], 'max_gap': 3})
    assert settings['oligos'] == (('arm5', 'arm3'), ('primer',)) and settings['max_gap'] == 3
    assert mapped_genes_settings({'aligner': 'bowtie'})['oligos'] is None


@pytest.mark.parametrize('attr, message', [
    ({'aligner': 'bowtie2', 'oligos': ['left', 'right']}, 'oligos needs aligner: bowtie'),
    ({'aligner': 'bowtie', 'oligos': ['left']}, 'at least two oligos'),
    ({'aligner': 'bowtie', 'oligos': 'left'}, 'at least two oligos'),
    ({'aligner': 'bowtie', 'oligos': ['left', 'left']}, 'more than once'),
    ({'aligner': 'bowtie', 'oligos': ['left', []]}, 'each oligo is a layout part name'),
    ({'aligner': 'bowtie', 'max_gap': 3}, 'max_gap only applies with oligos'),
    ({'aligner': 'bowtie', 'oligos': ['left', 'right'], 'max_gap': -1}, 'max_gap must be a non-negative'),
])
def test_settings_reject(attr, message):
    with pytest.raises(ValueError, match=message):
        mapped_genes_settings(attr, 'attributes.hits')


def test_oligo_parts_must_be_adjacent_layout_parts():
    names = ['arm5_site', 'arm3_site', 'gap', 'primer_site']
    check_oligo_parts((('arm5_site', 'arm3_site'), ('primer_site',)), names)
    with pytest.raises(ValueError, match='adjacent and in layout order'):
        check_oligo_parts((('arm5_site', 'primer_site'), ('arm3_site',)), names)
    with pytest.raises(ValueError, match="'left' is not a part"):
        check_oligo_parts((('left',), ('primer_site',)), names)


def test_validate_checks_oligos():
    def errors(protocol):
        return {f.path: f.message for f in validate_protocol(protocol) if f.severity == 'error'}
    assert errors(layout_protocol(['left', 'right'], max_gap=5)) == {}
    found = errors(layout_protocol(['left', 'middle']))
    assert "'middle' is not a part" in found['attributes.pair_genes.oligos']
    protocol = layout_protocol(['left', 'right'])
    protocol['attributes']['pair_genes']['target'] = 'odd.b'
    assert 'needs target: target_region' in errors(protocol)['attributes.pair_genes.oligos']
    protocol = layout_protocol(['left', 'right'])
    protocol['attributes']['pair_genes']['aligner'] = 'bowtie2'
    assert 'needs aligner: bowtie' in errors(protocol)['attributes.pair_genes.oligos']
    protocol = layout_protocol(['left', 'right'])
    del protocol['attributes']['pair_genes']['oligos']
    protocol['attributes']['pair_genes']['max_gap'] = 4
    assert 'only applies with oligos' in errors(protocol)['attributes.pair_genes.max_gap']


@needs_bowtie
def test_run_with_oligos_on_a_toy_genome(tmp_path):
    # `uprobe run` on the toy genome: layout parts reach the attribute and summary.json lists the oligos.
    import json
    import yaml
    from click.testing import CliRunner
    from uprobe.core.cli import cli
    from test_run_summary import toy_genome
    genomes = toy_genome(tmp_path, align_index=['bowtie'])
    protocol = layout_protocol(['left', 'right'])
    protocol.update(genome='toy', targets=['ToyGeneA', 'ToyGeneC'], encoding={})
    protocol['extracts']['target_region']['step'] = 20
    del protocol['attributes']['window_genes']
    (tmp_path / 'protocol.yaml').write_text(yaml.safe_dump(protocol))
    result = CliRunner().invoke(cli, ['run', '-p', str(tmp_path / 'protocol.yaml'), '-g', str(genomes),
                                      '-o', str(tmp_path / 'out'), '--threads', '2', '--json'])
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout)
    stats = summary['off_targets']['attributes']['pair_genes']
    assert stats['oligos'] == [['left'], ['right']] and stats['max_mismatches'] == 2
    table = pd.read_csv(tmp_path / 'out' / summary['output_files']['tsv'], sep='\t')
    assert len(table) > 0 and (table['pair_genes'] == 1).all()
