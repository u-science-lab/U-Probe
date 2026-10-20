"""kmer_count reads the Jellyfish index of its own k-mer length."""
import pandas as pd
import pytest

from uprobe.core import attributes, tools
from uprobe.core.attributes import add_attributes


@pytest.fixture
def jellyfish(monkeypatch, tmp_path):
    """Fake `jellyfish count` (writes the k into the index) and `jellyfish query` (records the index read)."""
    monkeypatch.chdir(tmp_path)
    built, queried = [], []
    def build(fasta, k, out, threads, size):
        built.append(k)
        with open(out, 'w') as handle:
            handle.write(f'000000040{{"key_len":{2 * k},"canonical":true}}')
    def query(outdir, name, records, index, k, threads):
        queried.append((k, str(index)))
        return {key: 0 for key in records}
    monkeypatch.setattr(tools, 'build_jf_index', build)
    monkeypatch.setattr(attributes, 'cal_kmer_count', query)
    return built, queried


def run(tmp_path, *lengths):
    protocol = {'attributes': {f'k{k}': {'target': 's', 'type': 'kmer_count', 'aligner': 'jellyfish', 'kmer_len': k}
                               for k in lengths}}
    genome = {'fasta': str(tmp_path / 'ref.fa'), 'jellyfish': True}
    add_attributes(pd.DataFrame({'s': ['ACGT' * 10]}), protocol, genome, {})


def test_each_kmer_len_gets_its_own_index(jellyfish, tmp_path):
    built, queried = jellyfish
    run(tmp_path, 18, 35)
    index = tmp_path / 'jf_genome'
    assert built == [18, 35]
    assert queried == [(18, str(index / 'ref.k18.jf')), (35, str(index / 'ref.k35.jf'))]
    # A changed kmer_len builds its index; the others are reused.
    run(tmp_path, 35, 25)
    assert built == [18, 35, 25] and queried[-1] == (25, str(index / 'ref.k25.jf'))
    assert sorted(p.name for p in index.iterdir()) == ['ref.k18.jf', 'ref.k25.jf', 'ref.k35.jf']


def test_the_index_of_older_versions_is_used_only_for_its_k(jellyfish, tmp_path):
    built, queried = jellyfish
    (tmp_path / 'jf_genome').mkdir()
    legacy = tmp_path / 'jf_genome' / 'ref.jf'
    legacy.write_text('000000040{"key_len":70,"canonical":true}')  # 35-mers
    run(tmp_path, 35, 18)
    assert queried == [(35, str(legacy)), (18, str(tmp_path / 'jf_genome' / 'ref.k18.jf'))]
    assert built == [18]
    assert tools.jellyfish_kmer_len(legacy) == 35 and tools.jellyfish_kmer_len(tmp_path / 'missing.jf') is None


def test_an_interrupted_build_leaves_no_index(monkeypatch, tmp_path):
    def fail(fasta, k, out, threads, size):
        open(out, 'w').write('partial')
        raise RuntimeError('killed')
    monkeypatch.setattr(tools, 'build_jf_index', fail)
    with pytest.raises(RuntimeError):
        tools.jellyfish_index(tmp_path / 'ref.fa', 35)
    assert list((tmp_path / 'jf_genome').iterdir()) == []
