"""source: mRNA (the whole mature RNA), extracts.target_region.transcripts and windows across exon
junctions (one selected transcript) on a toy annotation."""
import json
import logging
import random
import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from uprobe.core.cli import cli
from uprobe.core.copies import copy_fractions, target_sequences
from uprobe.core.gen.fun import TranscriptSelectionError, generate_target_seqs, read_gtf, select_transcripts
from uprobe.core.process import post_process
from uprobe.core.utils import reverse_complement
from uprobe.core.validate import render_text, validate_protocol

MANE = 'tag "MANE_Select"; tag "Ensembl_canonical";'
CANONICAL = 'tag "Ensembl_canonical";'
# (gene, transcript, tags, feature, start, end, strand); 1-based inclusive, as in a GTF.
FEATURES = [
    # P, + strand: ENSTP1 (MANE Select) has a 10-nt middle exon; ENSTP2 a long alternative last exon.
    *[('P', 'ENSTP1', MANE, *f, '+') for f in [
        ('exon', 101, 300), ('exon', 401, 410), ('exon', 501, 700),
        ('five_prime_utr', 101, 200), ('CDS', 201, 300), ('CDS', 401, 410), ('CDS', 501, 600),
        ('three_prime_utr', 601, 700)]],
    *[('P', 'ENSTP2', 'tag "basic";', *f, '+') for f in [
        ('exon', 101, 300), ('exon', 801, 1700),
        ('five_prime_utr', 101, 200), ('CDS', 201, 300), ('CDS', 801, 900), ('three_prime_utr', 901, 1700)]],
    # M, - strand, no MANE Select (as outside human): ENSTM3 is its canonical transcript.
    *[('M', 'ENSTM3', CANONICAL, *f, '-') for f in [
        ('exon', 2001, 2200), ('exon', 2401, 2600),
        ('three_prime_utr', 2001, 2100), ('CDS', 2101, 2200), ('CDS', 2401, 2500), ('five_prime_utr', 2501, 2600)]],
    *[('M', 'ENSTM4', '', *f, '-') for f in [('exon', 2001, 2250), ('exon', 2401, 2600), ('CDS', 2151, 2250),
                                             ('CDS', 2401, 2500)]],
    # N, non-coding, GENCODE-style versioned ids.
    *[('N', 'ENSTN5.2', CANONICAL, 'exon', s, e, '+') for s, e in [(3001, 3200), (3301, 3500)]],
    ('N', 'ENSTN6.1', '', 'exon', 3001, 3500, '+'),
    # X, not a target.
    ('X', 'ENSTX7', MANE, 'exon', 3601, 3800, '+'),
    # Q: no tags at all.
    ('Q', 'ENSTQ8', '', 'exon', 3901, 4000, '+'),
]


@pytest.fixture
def annotation(tmp_path):
    rng = random.Random(7)
    genome = ''.join(rng.choice('ACGT') for _ in range(4200))
    (tmp_path / 'g.fa').write_text(f'>1\n{genome}\n')
    lines = [f'1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "G{gene}"; transcript_id "{transcript}"; '
             f'gene_name "{gene}"; {tags}\n' for gene, transcript, tags, feature, start, end, strand in FEATURES]
    (tmp_path / 'g.gtf').write_text(''.join(lines))
    return str(tmp_path / 'g.fa'), str(tmp_path / 'g.gtf'), genome


def spliced(genome, transcript, feature='exon'):
    """The RNA (or CDS) of a transcript, 5' to 3'."""
    parts = sorted((s, e, strand) for _, t, _, f, s, e, strand in FEATURES if t == transcript and f == feature)
    sequence = ''.join(genome[s - 1:e] for s, e, _ in parts)
    return reverse_complement(sequence) if parts[0][2] == '-' else sequence


def genomic(row):
    """1-based genomic interval of an exon window."""
    _, exon_start, exon_end, strand = row['exon_name'].split('_')
    if strand == '+':
        return int(exon_start) + row['start'] - 1, int(exon_start) + row['end'] - 1
    return int(exon_end) - row['end'] + 1, int(exon_end) - row['start'] + 1


def inside(interval, spans):
    return any(start <= interval[0] and interval[1] <= end for start, end in spans)


def ranges(gene, feature):
    return [(s, e) for g, _, _, f, s, e, _ in FEATURES if g == gene and f == feature]


needs_bowtie = pytest.mark.skipif(not (shutil.which('bowtie') and shutil.which('bowtie-build')),
                                  reason='needs Bowtie 1 (bowtie, bowtie-build)')


# -- source: mRNA ------------------------------------------------------------------

def test_mrna_has_utr_windows_and_exon_stays_the_coding_sequence(annotation):
    fasta, gtf, _ = annotation
    windows = {source: generate_target_seqs(source, ['P', 'M', 'N'], fasta, gtf, min_length=20, step=1)
               for source in ('mRNA', 'exon', 'CDS')}
    utrs = ranges('P', 'five_prime_utr') + ranges('P', 'three_prime_utr') + \
        ranges('M', 'five_prime_utr') + ranges('M', 'three_prime_utr')
    for source, df in windows.items():
        coding = df[df.target != 'N']
        located = [genomic(row) for _, row in coding.iterrows()]
        in_utr = sum(inside(interval, utrs) for interval in located)
        if source == 'mRNA':
            assert in_utr > 1000  # 5' and 3' UTR windows of both genes
        else:
            # exon and CDS: only CDS windows, as before mRNA existed.
            assert in_utr == 0
            assert all(inside(interval, ranges('P', 'CDS') + ranges('M', 'CDS')) for interval in located)
    pd.testing.assert_frame_equal(windows['exon'], windows['CDS'])
    # A non-coding gene has no CDS: all three sources tile its exons.
    for source in ('exon', 'CDS'):
        assert windows[source][windows[source].target == 'N'].target_region.tolist() == \
            windows['mRNA'][windows['mRNA'].target == 'N'].target_region.tolist()


def test_all_transcripts_keep_one_window_per_exon(annotation):
    # transcripts: all tiles each exon entry of the merged exon model on its own, as before
    # junction windows existed: P's exons 101-300 (both transcripts), 501-700 and 801-1700.
    fasta, gtf, genome = annotation
    df = generate_target_seqs('mRNA', ['P', 'M'], fasta, gtf, min_length=20, step=1)
    p = df[df.target == 'P']
    assert len(p) == 181 + 181 + 881 and not df.exon_name.str.contains(';').any()
    assert sorted(set(p.exon_name)) == ['1_101_300_+', '1_501_700_+', '1_801_1700_+']
    assert set(p.loc[p.exon_name == '1_101_300_+', 'n_trans']) == {2}
    assert set(p.loc[p.exon_name == '1_101_300_+', 'transcript_names'].map(sorted).map(tuple)) == {('ENSTP1', 'ENSTP2')}
    for _, row in df.iterrows():
        start, end = genomic(row)
        window = genome[start - 1:end]
        assert (window if row['exon_name'].endswith('+') else reverse_complement(window)) == row['target_region']
    # The 10-nt exon gives no window but still counts on the model: exon 501-700 starts at 200 + 10 + 1.
    assert p.loc[p.exon_name == '1_501_700_+', 'target_start'].min() == 211
    # Two selected transcripts are tiled on their merged exon model too.
    both = generate_target_seqs('mRNA', ['P'], fasta, gtf, min_length=20, step=1, transcripts=['ENSTP1', 'ENSTP2'])
    pd.testing.assert_frame_equal(both, p.reset_index(drop=True))


def test_post_process_steps_work_on_mrna_windows(annotation):
    fasta, gtf, _ = annotation
    df = generate_target_seqs('mRNA', ['P', 'M'], fasta, gtf, min_length=20, step=1)
    # Alternative exon entries (101-300 of ENSTP1 and ENSTP2) give windows at the same model positions.
    assert df.duplicated(['target', 'target_start']).any()
    df['hits_off_targets'] = ['OFF:0' if i % 3 == 0 else '' for i in range(len(df))]
    config = {'post_process': {'remove_overlap': {'location_interval': 2},
                               'off_target_pool': {'attribute': 'hits', 'max_probes_per_gene': 2},
                               'equal_space': {'number_desired': 10}}}
    result = post_process(df, config)
    for target, group in result.groupby('target'):
        assert len(group) == 10
        starts = sorted(zip(group.target_start, group.target_end))
        assert all(b[0] > a[1] + 2 for a, b in zip(starts, starts[1:])), target
        assert (group.hits_off_targets != '').sum() <= 2


def test_utr_source_reads_ensembl_utr_features(annotation):
    fasta, gtf, _ = annotation
    df = generate_target_seqs('UTR', ['P'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    assert sorted(set(df.utr_name)) == ['P_101_200', 'P_601_700']
    assert len(df) == 2 * 81


# -- one transcript: windows across exon-exon junctions ------------------------------

def genomic_blocks(row):
    """Genomic intervals (1-based, 5' to 3') of a window, from exon_name and start/end along the named segments."""
    blocks, offset = [], 0
    for name in row['exon_name'].split(';'):
        _, segment_start, segment_end, strand = name.split('_')
        segment_start, segment_end = int(segment_start), int(segment_end)
        first, last = max(row['start'] - offset, 1), min(row['end'] - offset, segment_end - segment_start + 1)
        assert first <= last, f'{name} is named but not covered by window {row["probe_id"]}'
        blocks.append((segment_start + first - 1, segment_start + last - 1) if strand == '+' else
                      (segment_end - last + 1, segment_end - first + 1))
        offset += segment_end - segment_start + 1
    return blocks, strand


def sequence_of(genome, blocks, strand):
    pieces = [genome[start - 1:end] for start, end in blocks]
    return ''.join(pieces) if strand == '+' else ''.join(reverse_complement(piece) for piece in pieces)


@pytest.mark.parametrize('transcripts', ['mane', ['ENSTP1', 'ENSTM3.1', 'ENSTN5']])
def test_one_transcript_is_tiled_across_its_exon_junctions(annotation, transcripts):
    fasta, gtf, genome = annotation
    df = generate_target_seqs('mRNA', ['P', 'M', 'N'], fasta, gtf, min_length=20, step=1, transcripts=transcripts)
    junctions = {'P': 0, 'M': 0, 'N': 0}
    for gene, transcript in (('P', 'ENSTP1'), ('M', 'ENSTM3'), ('N', 'ENSTN5.2')):
        rna = spliced(genome, transcript)
        rows = df[df.target == gene]
        # Every 20-mer of the spliced RNA, in order, at its position on the RNA; ids numbered 5' to 3'.
        assert rows.target_region.tolist() == [rna[i:i + 20] for i in range(len(rna) - 19)]
        assert rows.target_start.tolist() == list(range(1, len(rna) - 18))
        assert (rows.target_end - rows.target_start).eq(19).all()
        assert rows.probe_id.tolist() == [f'{gene}_{n}' for n in range(1, len(rows) + 1)]
        assert set(rows.transcript_names.map(tuple)) == {(transcript,)} and set(rows.n_trans) == {1}
        # exon_name and start/end place each window on the genome, also across junctions.
        for _, row in rows.iterrows():
            blocks, strand = genomic_blocks(row)
            assert sequence_of(genome, blocks, strand) == row['target_region']
            assert row['sub_region'] == f'{row["start"]}_{row["end"]}' and row['end'] - row['start'] == 19
        junctions[gene] = int(rows.exon_name.str.contains(';').sum())
    # 19 windows span each junction; around P's 10-nt exon, 9 of them span both of its junctions.
    assert junctions == {'P': 19 + 19 - 9, 'M': 19, 'N': 19}
    p = df[df.target == 'P'].set_index('target_start')
    assert p.loc[191, 'exon_name'] == '1_101_300_+;1_401_410_+' and (p.loc[191, 'start'], p.loc[191, 'end']) == (191, 210)
    assert p.loc[192, 'exon_name'] == '1_101_300_+;1_401_410_+;1_501_700_+'
    assert p.loc[201, 'exon_name'] == '1_401_410_+;1_501_700_+' and (p.loc[201, 'start'], p.loc[201, 'end']) == (1, 20)
    assert p.loc[211, 'exon_name'] == '1_501_700_+' and (p.loc[211, 'start'], p.loc[211, 'end']) == (1, 20)
    m = df[df.target == 'M'].set_index('target_start')  # - strand: the exon at 2401-2600 comes first
    assert m.loc[1, 'exon_name'] == '1_2401_2600_-' and m.loc[200, 'exon_name'] == '1_2401_2600_-;1_2001_2200_-'


def test_one_transcript_with_length_range_and_step(annotation):
    fasta, gtf, genome = annotation
    df = generate_target_seqs('mRNA', ['P'], fasta, gtf, min_length=[18, 22], step=3, transcripts='mane')
    rna = spliced(genome, 'ENSTP1')
    assert df.target_region.tolist() == [rna[i:i + k] for i in range(0, len(rna) - 17, 3)
                                         for k in range(18, min(22, len(rna) - i) + 1)]
    for _, row in df.iterrows():
        assert rna[row['target_start'] - 1:row['target_end']] == row['target_region']
        assert sequence_of(genome, *genomic_blocks(row)) == row['target_region']


@pytest.mark.parametrize('source', ['CDS', 'exon'])
def test_one_transcript_cds_is_tiled_across_its_junctions(annotation, source):
    # The spliced coding sequence, the 10-nt CDS piece of ENSTP1 included; exon and CDS stay identical.
    fasta, gtf, genome = annotation
    df = generate_target_seqs(source, ['P', 'M', 'N'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    for gene, transcript, feature in (('P', 'ENSTP1', 'CDS'), ('M', 'ENSTM3', 'CDS'), ('N', 'ENSTN5.2', 'exon')):
        sequence = spliced(genome, transcript, feature)
        rows = df[df.target == gene]
        assert rows.target_region.tolist() == [sequence[i:i + 20] for i in range(len(sequence) - 19)]
        assert all(sequence[s - 1:e] == w for s, e, w in zip(rows.target_start, rows.target_end, rows.target_region))
        for _, row in rows.iterrows():
            assert sequence_of(genome, *genomic_blocks(row)) == row['target_region']
    assert len(spliced(genome, 'ENSTP1', 'CDS')) == 210
    assert set(df.loc[df.target == 'P', 'exon_name'].str.split(';').explode()) == {
        '1_201_300_+', '1_401_410_+', '1_501_600_+'}


def test_post_process_steps_work_on_junction_windows(annotation):
    fasta, gtf, genome = annotation
    df = generate_target_seqs('mRNA', ['P', 'M'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    df['hits_off_targets'] = ['OFF:0' if i % 3 == 0 else '' for i in range(len(df))]
    spaced = post_process(df.copy(), {'post_process': {'remove_overlap': {'location_interval': 2}}})
    # Windows 22 nt apart from the 5' end; some of them span junctions.
    for target, group in spaced.groupby('target'):
        assert group.target_start.tolist() == list(range(1, group.target_start.max() + 1, 22)), target
    assert spaced.exon_name.str.contains(';').any()
    config = {'post_process': {'remove_overlap': {'location_interval': 2},
                               'off_target_pool': {'attribute': 'hits', 'max_probes_per_gene': 2},
                               'equal_space': {'number_desired': 10}}}
    result = post_process(df, config)
    for target, group in result.groupby('target'):
        rna = spliced(genome, {'P': 'ENSTP1', 'M': 'ENSTM3'}[target])
        assert len(group) == 10
        assert all(rna[s - 1:e] == w for s, e, w in zip(group.target_start, group.target_end, group.target_region))
        starts = sorted(zip(group.target_start, group.target_end))
        assert all(b[0] > a[1] + 2 for a, b in zip(starts, starts[1:])), target
        assert (group.hits_off_targets != '').sum() <= 2


@needs_bowtie
def test_junction_windows_are_found_in_the_transcriptome(annotation, tmp_path, monkeypatch):
    # mapped_genes searches the spliced transcripts (transcript.fa), so a window across a
    # junction is found on its own transcript like any other window.
    from uprobe.core.attributes import add_attributes
    from uprobe.core.tools import build_transcripts_index
    fasta, gtf, _ = annotation
    build_transcripts_index(Path(gtf), Path(fasta), Path(fasta).parent / 'bowtie_transcript', threads=1,
                            aligner='bowtie')
    monkeypatch.chdir(tmp_path)
    df = generate_target_seqs('mRNA', ['P', 'M'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    protocol = {'attributes': {'hits': {'target': 'target_region', 'type': 'mapped_genes', 'aligner': 'bowtie',
                                        'max_mismatches': 2, 'strand': 'sense', 'exempt_copies': True}}}
    genome = {'fasta': fasta, 'gtf': gtf, 'align_index': ['bowtie']}
    result = add_attributes(df, protocol, genome, {})
    assert result.exon_name.str.contains(';').sum() == 29 + 19
    assert result.hits.eq(1).all() and result.hits_off_targets.eq('').all()


# -- transcripts -------------------------------------------------------------------

def test_mane_canonical_and_fallback(annotation, caplog):
    fasta, gtf, genome = annotation
    caplog.set_level(logging.INFO, logger='uprobe.core.gen.fun')
    df = generate_target_seqs('mRNA', ['P', 'M', 'N'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    assert df.attrs['transcripts'] == {
        'P': {'selection': 'mane', 'chosen_by': 'MANE_Select', 'ids': ['ENSTP1']},
        'M': {'selection': 'mane', 'chosen_by': 'Ensembl_canonical', 'ids': ['ENSTM3']},
        'N': {'selection': 'mane', 'chosen_by': 'Ensembl_canonical', 'ids': ['ENSTN5.2']}}
    # One line for all targets without MANE (every target of a non-human genome), not one per target.
    fallback = [r.getMessage() for r in caplog.records if 'no MANE_Select transcript' in r.getMessage()]
    assert fallback == ['2 of 3 target(s) have no MANE_Select transcript (MANE covers human genes only); '
                        'using their Ensembl_canonical transcript: M ENSTM3, N ENSTN5.2']
    for gene, transcript in (('P', 'ENSTP1'), ('M', 'ENSTM3'), ('N', 'ENSTN5.2')):
        rows = df[df.target == gene]
        assert set(rows.transcript_names.map(tuple)) == {(transcript,)}
        assert set(rows.target_region) <= {spliced(genome, transcript)[i:i + 20] for i in range(2000)}
    canonical = generate_target_seqs('mRNA', ['P', 'M'], fasta, gtf, min_length=20, step=1, transcripts='canonical')
    assert {t: e['chosen_by'] for t, e in canonical.attrs['transcripts'].items()} == \
        {'P': 'Ensembl_canonical', 'M': 'Ensembl_canonical'}
    everything = generate_target_seqs('mRNA', ['P'], fasta, gtf, min_length=20, step=1)
    assert everything.attrs['transcripts'] == {'P': {'selection': 'all', 'chosen_by': 'all',
                                                     'ids': ['ENSTP1', 'ENSTP2']}}
    assert len(everything) > len(df[df.target == 'P'])
    # CDS with one transcript: the windows of its spliced coding sequence, the 10-nt CDS piece included.
    cds = generate_target_seqs('CDS', ['P'], fasta, gtf, min_length=20, step=1, transcripts='mane')
    coding = spliced(genome, 'ENSTP1', 'CDS')
    assert cds.target_end.max() == len(coding) == 210
    assert cds.target_region.tolist() == [coding[i:i + 20] for i in range(191)]


def test_listed_transcripts_match_with_or_without_version(annotation):
    fasta, gtf, _ = annotation
    df = generate_target_seqs('mRNA', ['P', 'N'], fasta, gtf, min_length=20, step=1,
                              transcripts=['ENSTP2.4', 'ENSTN6'])
    assert df.attrs['transcripts'] == {'P': {'selection': 'list', 'chosen_by': 'listed', 'ids': ['ENSTP2']},
                                       'N': {'selection': 'list', 'chosen_by': 'listed', 'ids': ['ENSTN6.1']}}
    assert set(df.transcript_names.map(tuple)) == {('ENSTP2',), ('ENSTN6.1',)}


@pytest.mark.parametrize('targets, setting, problems', [
    (['P'], ['ENSTP9'], ['transcript ENSTP9 is not in the GTF', 'no listed transcript belongs to target P; list '
                         'one of its transcripts (e.g. ENSTP1) or use canonical, mane or all']),
    (['P'], ['ENSTP1', 'ENSTX7'], ['transcript ENSTX7 belongs to X, which is not a target']),
    (['P', 'M'], ['ENSTP1'], ['no listed transcript belongs to target M; list one of its transcripts '
                              '(e.g. ENSTM3) or use canonical, mane or all']),
    (['P', 'Q'], 'canonical', ['no transcript of Q is tagged Ensembl_canonical in the GTF; list its transcript '
                               'ids or use transcripts: all']),
    (['Q'], 'mane', ['no transcript of Q is tagged MANE_Select or Ensembl_canonical in the GTF; list its '
                     'transcript ids or use transcripts: all']),
])
def test_wrong_transcripts_are_errors(annotation, targets, setting, problems):
    fasta, gtf, _ = annotation
    with pytest.raises(TranscriptSelectionError) as error:
        generate_target_seqs('mRNA', targets, fasta, gtf, min_length=20, step=1, transcripts=setting)
    assert error.value.problems == problems


@pytest.mark.parametrize('value', ['cannonical', 'ENSTP1', [], [''], {'P': 'mane'}, 3])
def test_invalid_transcripts_values(annotation, value):
    fasta, gtf, _ = annotation
    with pytest.raises(ValueError, match='transcripts must be all, canonical, mane or a list of transcript ids'):
        select_transcripts(read_gtf(gtf, extract_fields=['gene_name', 'transcript_id']), ['P'], value)


def test_copies_are_judged_against_the_selected_transcript(annotation):
    # A processed copy of ENSTP1, the MANE transcript, is a copy of P's MANE mRNA, but holds
    # less than half of the merged exon model of all P's transcripts (ENSTP2 adds 900 nt).
    fasta, gtf, genome = annotation
    rna = spliced(genome, 'ENSTP1')
    copy = 'ACGT' * 10 + rna + 'TTGCA' * 10
    shares = {}
    for setting in ('all', 'mane'):
        df = generate_target_seqs('mRNA', ['P'], fasta, gtf, min_length=20, step=1, transcripts=setting)
        sequence = target_sequences(df)['P']
        if setting == 'mane':
            # The windows cover the whole spliced mRNA, its 10-nt exon too (junction windows).
            assert sequence == rna.encode()
        shares[setting] = copy_fractions({'P': sequence}, {'P': {'COPY'}}, [('COPY', copy)], 30)[('P', 'COPY')]
    assert shares['mane'] == 1.0 and shares['all'] < 0.5


# -- uprobe validate -----------------------------------------------------------------

def protocol(source='exon', targets=('P',), **region):
    return {'name': 't', 'genome': 'toy', 'targets': list(targets), 'encoding': {},
            'extracts': {'target_region': {'source': source, 'length': 20, 'step': 1, **region}},
            'probes': {'probe': {'template': '{b}', 'parts': {'b': {'expr': 'rc(target_region)'}}}},
            'attributes': {'gc': {'target': 'target_region', 'type': 'gc_content'}},
            'post_process': {'remove_overlap': {'location_interval': 2}, 'equal_space': {'number_desired': 48}},
            'summary': {'report_name': 'rna_report', 'attributes': ['gc']}}


def genomes_of(annotation):
    fasta, gtf, _ = annotation
    return {'toy': {'fasta': fasta, 'gtf': gtf, 'align_index': []}}


SOURCE = 'extracts.target_region.source'
TRANSCRIPTS = 'extracts.target_region.transcripts'


def test_validate_warns_about_source_exon(annotation):
    genomes = genomes_of(annotation)
    findings = validate_protocol(protocol())
    assert [(f.path, f.severity) for f in findings] == [(SOURCE, 'warning')], render_text(findings)
    assert 'only the coding sequence' in findings[0].message and 'source: mRNA' in findings[0].suggestion
    # With the GTF: named for the protein-coding targets only.
    findings = validate_protocol(protocol(targets=('P', 'M', 'N')), genomes, check_targets=True)
    assert [f.message for f in findings if f.path == SOURCE] == [
        "2 target(s) are protein-coding (P, M): source: exon extracts only their coding sequence (CDS "
        "features), without the UTRs"]
    assert validate_protocol(protocol(targets=('N',)), genomes, check_targets=True) == []
    for source in ('mRNA', 'CDS'):
        assert validate_protocol(protocol(source), genomes, check_targets=True) == []


def test_validate_checks_transcripts(annotation):
    genomes = genomes_of(annotation)
    assert validate_protocol(protocol('mRNA', ('P', 'M'), transcripts='mane'), genomes, check_targets=True) == []
    findings = validate_protocol(protocol('mRNA', transcripts='cannonical'))
    assert [(f.path, f.severity) for f in findings] == [(TRANSCRIPTS, 'error')]
    assert findings[0].suggestion.startswith("did you mean 'canonical'?")
    findings = validate_protocol(protocol('mRNA', ('P', 'M'), transcripts=['ENSTP1', 'ENSTX7']), genomes,
                                 check_targets=True)
    assert [(f.path, f.message) for f in findings] == [
        (TRANSCRIPTS, 'transcript ENSTX7 belongs to X, which is not a target'),
        (TRANSCRIPTS, 'no listed transcript belongs to target M; list one of its transcripts (e.g. ENSTM3) or use '
                      'canonical, mane or all')]
    region = protocol('genome', ('1:1-200',), transcripts='mane')
    assert [(f.path, f.severity) for f in validate_protocol(region)] == [(TRANSCRIPTS, 'warning')]


# -- uprobe run: summary.json ----------------------------------------------------------

def run(annotation, tmp_path, **region):
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump(genomes_of(annotation)))
    path = tmp_path / 'protocol.yaml'
    path.write_text(yaml.safe_dump(protocol('mRNA', ('P', 'M'), **region)))
    output = tmp_path / 'out'
    result = CliRunner().invoke(cli, ['run', '-p', str(path), '-g', str(genomes), '-o', str(output), '--json'])
    return result, json.loads((output / 'summary.json').read_text())


def test_summary_reports_the_transcripts_of_each_target(annotation, tmp_path):
    result, summary = run(annotation, tmp_path, transcripts='mane')
    assert result.exit_code == 0, result.output
    assert summary['inputs']['extract_source'] == 'mRNA' and summary['inputs']['transcripts'] == 'mane'
    assert {name: entry['transcripts'] for name, entry in summary['targets'].items()} == {
        'P': {'selection': 'mane', 'chosen_by': 'MANE_Select', 'ids': ['ENSTP1']},
        'M': {'selection': 'mane', 'chosen_by': 'Ensembl_canonical', 'ids': ['ENSTM3']}}
    table = pd.read_csv(tmp_path / 'out' / summary['output_files']['tsv'], sep='\t')
    assert set(table.transcript_names) == {'ENSTP1', 'ENSTM3'}


def test_a_run_parses_the_gtf_once(annotation, tmp_path, monkeypatch):
    # Target validation and extraction each parsed the GTF: about 23 s each, 45 s of a
    # 55-75 s run, on the human Ensembl 113 GTF.
    from uprobe.core import api
    from uprobe.core.gen import fun
    calls = []
    original = fun.read_gtf

    def counting(*args, **kwargs):
        calls.append(args[0] if args else kwargs['gtf'])
        return original(*args, **kwargs)
    monkeypatch.setattr(fun, 'read_gtf', counting)
    monkeypatch.setattr(api, 'read_gtf', counting, raising=False)
    result, summary = run(annotation, tmp_path, transcripts='mane')
    assert result.exit_code == 0, result.output
    assert summary['targets']['M']['final'] > 0
    assert len(calls) == 1


def test_wrong_transcripts_fail_the_run(annotation, tmp_path):
    result, summary = run(annotation, tmp_path, transcripts=['ENSTX7'])
    assert result.exit_code == 1
    assert summary['status'] == 'error' and summary['stopped_at'] == 'extract_target_regions'
    assert 'transcript ENSTX7 belongs to X, which is not a target' in summary['error']


# -- a target whose segments are all shorter than the window -------------------------------------

@pytest.mark.parametrize('source, targets, empty', [
    ('CDS', ('P', 'N'), 'P'),    # P: CDS pieces of 100, 10, 100 and 100 nt
    ('exon', ('P', 'N'), 'P'),   # N (non-coding): exons of 200 and 500 nt
    ('UTR', ('P', 'Q'), 'Q'),    # Q: one 100-nt exon and no UTR; P: an 800-nt 3' UTR
])
def test_a_target_without_long_enough_segments_does_not_stop_the_others(annotation, source, targets, empty,
                                                                         caplog):
    fasta, gtf, _ = annotation
    with caplog.at_level(logging.WARNING):
        windows = generate_target_seqs(source, list(targets), fasta, gtf, 150, step=10)
    other = next(target for target in targets if target != empty)
    assert set(windows['target']) == {other} and len(windows) > 0
    assert f'Gene {empty}' in caplog.text


def test_run_reports_the_target_without_windows(annotation, tmp_path):
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump(genomes_of(annotation)))
    path = tmp_path / 'protocol.yaml'
    path.write_text(yaml.safe_dump(protocol('CDS', ('P', 'N'), length=150, step=10)))
    output = tmp_path / 'out'
    result = CliRunner().invoke(cli, ['run', '-p', str(path), '-g', str(genomes), '-o', str(output), '--json'])
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout)
    assert summary['status'] == 'partial'
    assert summary['targets']['P']['candidate_windows'] == 0 and summary['targets']['N']['final'] > 0
    warning = next(w for w in summary['warnings'] if w['code'] == 'no_candidate_windows')
    assert warning['targets'] == ['P']
