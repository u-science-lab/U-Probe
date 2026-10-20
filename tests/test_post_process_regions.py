import pandas as pd
import pytest

from uprobe.core.gen.fun import segment_offsets
from uprobe.core.process import post_process
from uprobe.core.process.equal_space import equal_space, number_desired
from uprobe.core.process.regions import region_bounds


def test_region_bounds_accepts_every_source_format():
    df = pd.DataFrame({
        'sub_region': ['101-140', '1_40', '5-44', None],
        'start': [None, None, 7, 9],
        'end': [None, None, 46, 48],
    })
    bounds = region_bounds(df)
    assert bounds['start'].tolist() == [101, 1, 7, 9]
    assert bounds['end'].tolist() == [140, 40, 46, 48]


def test_region_bounds_rejects_unparseable_rows():
    with pytest.raises(ValueError):
        region_bounds(pd.DataFrame({'sub_region': ['chr1:1-40']}))


def test_equal_space_orders_exon_sub_regions_by_start():
    # Exon/UTR extraction writes "start_end"; int("2_41") parses as 241, so the old
    # split('-') ordering sorted by the concatenated digits (1_100 -> 1100 after 241).
    df = pd.DataFrame({
        'probe_id': ['a', 'b', 'c'],
        'target': ['G'] * 3,
        'sub_region': ['1_100', '2_41', '50_89'],
    })
    result = equal_space(df, {'number_desired': 2})
    assert result['probe_id'].tolist() == ['a', 'c']


def test_equal_space_on_genome_windows_is_unchanged():
    # Evenly spaced ranks over the windows sorted by start, ends included.
    df = pd.DataFrame({'probe_id': [f'w{i}' for i in range(10)], 'target': ['chr1:0-1000'] * 10,
                       'sub_region': [f'{s}-{s + 39}' for s in (901, 1, 101, 201, 301, 401, 501, 601, 701, 801)]})
    result = equal_space(df, {'chr1:0-1000': {'number_desired': 4}})
    assert result['sub_region'].tolist() == ['1-40', '301-340', '601-640', '901-940']


def test_segment_offsets_follow_the_merged_exon_model():
    # Plus strand: two transcripts share 100-250; the union is 100-300 and 500-600.
    plus = segment_offsets([('1', 100, 300, '+'), ('1', 100, 250, '+'), ('1', 500, 600, '+')])
    assert plus == {('1', 100, 300, '+'): 0, ('1', 100, 250, '+'): 0, ('1', 500, 600, '+'): 201}
    # Minus strand: the 5' end is the highest coordinate.
    minus = segment_offsets([('1', 100, 300, '-'), ('1', 100, 250, '-'), ('1', 500, 600, '-')])
    assert minus == {('1', 500, 600, '-'): 0, ('1', 100, 300, '-'): 101, ('1', 100, 250, '-'): 151}


def exon_windows(starts_per_exon, length=20):
    """Windows numbered within each of three exons, with their merged-model positions."""
    rows = []
    for exon, (offset, starts) in enumerate(starts_per_exon):
        rows += [{'probe_id': f'e{exon}_{s}', 'target': 'G', 'exon_name': f'1_{exon}_x_+', 'start': s,
                  'end': s + length - 1, 'sub_region': f'{s}_{s + length - 1}',
                  'target_start': offset + s, 'target_end': offset + s + length - 1} for s in starts]
    return pd.DataFrame(rows)


def test_equal_space_spreads_exon_windows_along_the_target():
    # Every exon numbers its windows from 1; ordering by start would interleave the exons.
    df = exon_windows([(0, [1, 41, 81]), (100, [1, 41, 81]), (200, [1, 41, 81])]).sample(frac=1, random_state=0)
    result = equal_space(df, {'number_desired': 3})
    assert result['probe_id'].tolist() == ['e0_1', 'e1_41', 'e2_81']
    assert result['target_start'].tolist() == [1, 141, 281]
    # Without the target positions (an older table), the windows are ordered within exons only.
    legacy = equal_space(df.drop(columns=['target_start', 'target_end']), {'number_desired': 3})
    assert legacy['start'].tolist() == [1, 41, 81]


def test_equal_space_keeps_every_probe_of_a_short_target():
    df = exon_windows([(0, [1, 41]), (100, [1])])
    result = equal_space(df, {'G': {'number_desired': 48}})
    assert sorted(result['probe_id']) == sorted(df['probe_id'])
    assert number_desired({'G': {'number_desired': 48}, 'number_desired': 5}, 'G') == 48
    assert number_desired({'number_desired': 5}, 'H') == 5 and number_desired({}, 'G') is None


def test_remove_overlap_orders_genome_windows_numerically():
    # As strings "100-139" < "41-80", so the non-overlapping 41-80 window was dropped.
    df = pd.DataFrame({
        'probe_id': ['a', 'b', 'c'],
        'target': ['chr1:0-200'] * 3,
        'sub_region': ['1-40', '41-80', '100-139'],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a', 'b', 'c']
    assert '_start' not in result.columns


def test_remove_overlap_genome_windows_with_empty_input():
    df = pd.DataFrame(columns=['probe_id', 'target', 'sub_region'])
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result.empty and list(result.columns) == ['probe_id', 'target', 'sub_region']


def test_remove_overlap_keeps_windows_of_each_direct_target():
    # Direct-sequence targets all have transcript_names 'custom'; windows of
    # different targets at the same positions must not remove each other.
    df = pd.DataFrame({
        'probe_id': ['a1', 'a2', 'b1', 'b2'],
        'target': ['A', 'A', 'B', 'B'],
        'transcript_names': ['custom'] * 4,
        'start': [0, 20, 0, 20],
        'end': [39, 59, 39, 59],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a1', 'b1']
    df['start'], df['end'] = [0, 40, 0, 40], [39, 79, 39, 79]
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert result['probe_id'].tolist() == ['a1', 'a2', 'b1', 'b2']


def test_remove_overlap_keeps_windows_of_each_exon():
    # Exon and UTR windows are numbered within their exon or UTR, so windows of
    # different exons of one transcript at the same positions do not overlap.
    # Without target_start (an older table), windows are compared per exon or UTR.
    df = pd.DataFrame({
        'probe_id': ['e1a', 'e1b', 'e1c', 'e2a', 'e2b', 'u1', 'c1'],
        'target': ['G'] * 6 + ['C'],
        'exon_name': ['1_101_200_+'] * 3 + ['1_301_400_+'] * 2 + [None, '.'],
        'utr_name': [None] * 5 + ['1_501_600_+', None],
        'transcript_names': [['T1']] * 5 + [['T1']] + ['custom'],
        'start': [1, 21, 41, 1, 41, 1, 1],
        'end': [40, 60, 80, 40, 80, 40, 40],
    })
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': 0}}})
    assert sorted(result['probe_id']) == ['c1', 'e1a', 'e1c', 'e2a', 'e2b', 'u1']


def dna_windows(tmp_path, targets):
    """Candidate windows of a DNA-mode (source: genome) run, as UProbeAPI extracts them."""
    import random
    from uprobe.core.api import UProbeAPI
    rng = random.Random(5)
    (tmp_path / 'g.fa').write_text('>1\n' + ''.join(rng.choice('ACGT') for _ in range(400)) + '\n')
    api = object.__new__(UProbeAPI)
    api.genome = {'fasta': str(tmp_path / 'g.fa'), 'gtf': None}
    api.protocol = {'targets': targets, 'extracts': {'target_region': {'source': 'genome', 'length': 40, 'step': 20}}}
    return api.generate_target_seqs()


def test_remove_overlap_in_dna_mode_with_regions_and_direct_sequences(tmp_path):
    # Direct sequences add transcript_names and start/end columns that genome
    # windows lack (NaN), so all but the first window of each region were dropped.
    sequence = 'ACGTTGCAAGCTTCGATGCA' * 10
    mixed = dna_windows(tmp_path, ['1:1-200', '1:201-400', {'S': sequence}])
    assert mixed['transcript_names'].isna().any() and mixed['start'].isna().any()
    config = {'post_process': {'remove_overlap': {'location_interval': 0}}}
    result = post_process(mixed, config)
    assert result.groupby('target').size().to_dict() == {'1:1-200': 5, '1:201-400': 5, 'S': 5}
    assert result.loc[result['target'] == '1:1-200', 'sub_region'].tolist() == \
        ['1-40', '41-80', '81-120', '121-160', '161-200']
    # The same windows as runs with only the regions or only the sequence.
    alone = pd.concat([post_process(dna_windows(tmp_path, ['1:1-200', '1:201-400']), config),
                       post_process(dna_windows(tmp_path, [{'S': sequence}]), config)])
    assert sorted(result['probe_id']) == sorted(alone['probe_id'])


def alternative_exon_genome(tmp_path):
    """Gene G (+ strand) and gene H (- strand), each with two transcripts whose exons share sites."""
    import random
    rng = random.Random(3)
    seq = ''.join(rng.choice('ACGT') for _ in range(1000))
    (tmp_path / 'g.fa').write_text(f'>1\n{seq}\n')
    exons = [('G', 'T1', 101, 300, '+'), ('G', 'T1', 501, 600, '+'),
             ('G', 'T2', 101, 250, '+'), ('G', 'T2', 501, 650, '+'),
             ('H', 'T3', 701, 900, '-'), ('H', 'T4', 751, 900, '-')]
    (tmp_path / 'g.gtf').write_text(''.join(
        f'1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; gene_name "{gene}"; transcript_id "{t}";\n'
        for gene, t, start, end, strand in exons))
    return str(tmp_path / 'g.fa'), str(tmp_path / 'g.gtf')


def genomic_interval(row):
    _, exon_start, exon_end, strand = row['exon_name'].split('_')
    if strand == '+':
        return int(exon_start) + row['start'] - 1, int(exon_start) + row['end'] - 1
    return int(exon_end) - row['end'] + 1, int(exon_end) - row['start'] + 1


@pytest.mark.parametrize('location_interval', [0, 2])
def test_remove_overlap_compares_alternative_exon_entries(tmp_path, location_interval):
    # Windows are numbered within each exon entry, and entries of different
    # transcripts (101-300 and 101-250 of G) contain the same sites; they were
    # never compared, so the same window was selected once per entry.
    from uprobe.core.gen.fun import generate_target_seqs
    fasta, gtf = alternative_exon_genome(tmp_path)
    df = generate_target_seqs('exon', ['G', 'H'], fasta, gtf, min_length=40, step=10)
    assert df.duplicated(['target', 'target_start', 'target_end']).any()
    result = post_process(df, {'post_process': {'remove_overlap': {'location_interval': location_interval}}})
    for target, group in result.groupby('target'):
        intervals = sorted(genomic_interval(row) for _, row in group.iterrows())
        assert all(b[0] > a[1] + location_interval for a, b in zip(intervals, intervals[1:])), target
        assert group['target_region'].is_unique
    # The kept windows tile the merged exon model of each gene.
    assert sorted(result.loc[result['target'] == 'G', 'target_start'])[:3] == \
        ([1, 41, 81] if location_interval == 0 else [1, 51, 101])
    # A site in two entries is kept once and lists the transcripts of both.
    for target, transcripts in (('G', 'T1, T2'), ('H', 'T3, T4')):
        first = result[(result['target'] == target) & (result['target_start'] == 1)]
        assert first['transcript_names'].tolist() == [transcripts] and first['n_trans'].tolist() == [2]
    # Model positions 151-200 lie only in T1's exon 101-300.
    only_t1 = result[(result['target'] == 'G') & (result['target_start'] > 150) & (result['target_end'] <= 200)]
    assert len(only_t1) == 1 and only_t1['transcript_names'].tolist() == ['T1']
    assert only_t1['n_trans'].tolist() == [1]
