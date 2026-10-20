"""`uprobe run --threads` is the default of every threaded step."""
import pandas as pd
import yaml
from click.testing import CliRunner

from uprobe.core import api as api_module
from uprobe.core import attributes
from uprobe.core.attributes import add_attributes
from uprobe.core.cli import _generate_default_attributes, cli


def test_attributes_use_the_run_threads_unless_they_set_their_own(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    used = {}
    def search(name):
        def run(outdir, records, task_id, index, *args):
            used[name] = args[-2]  # (..., threads, strand)
            return {key: {} for key in records}
        return run
    monkeypatch.setattr(attributes, 'bowtie2_gene_hits', search('bowtie2'))
    monkeypatch.setattr(attributes, 'bowtie_gene_hits', search('bowtie'))
    monkeypatch.setattr(attributes, 'cal_mapped_sites',
                        lambda outdir, name, records, index, threads: used.setdefault(name, threads) and {})
    protocol = {'attributes': {
        'bt2': {'target': 's', 'type': 'mapped_genes', 'aligner': 'bowtie2'},
        'bt1': {'target': 's', 'type': 'mapped_genes', 'aligner': 'bowtie', 'threads': 2},
        'sites': {'target': 's', 'type': 'mapped_sites', 'aligner': 'bowtie2'},
    }}
    genome = {'fasta': str(tmp_path / 'ref.fa'), 'align_index': ['bowtie', 'bowtie2']}
    frame = pd.DataFrame({'probe_id': ['p1'], 'target': ['t'], 's': ['ACGT' * 6]})
    add_attributes(frame, protocol, genome, {}, 3)
    assert used == {'bowtie2': 3, 'bowtie': 2, 'sites': 3}


def test_jellyfish_index_is_built_with_the_run_threads(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    from uprobe.core import tools
    built = []
    def build(fasta, k, out, threads, size, canonical=True):
        built.append(threads)
        open(out, 'w').write('index')
    monkeypatch.setattr(tools, 'build_jf_index', build)
    monkeypatch.setattr(attributes, 'cal_kmer_count', lambda outdir, name, records, index, k, threads: {
        key: 0 for key in records})
    protocol = {'attributes': {'kmers': {'target': 's', 'type': 'kmer_count', 'aligner': 'jellyfish'}}}
    genome = {'fasta': str(tmp_path / 'ref.fa'), 'jellyfish': True}
    add_attributes(pd.DataFrame({'s': ['ACGT' * 10]}), protocol, genome, {}, 4)
    assert built == [4]


def test_generated_kmer_count_takes_the_run_threads():
    protocol = {'extracts': {'target_region': {'source': 'genome', 'length': 40}},
                'probes': {'p': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}}}
    generated = _generate_default_attributes(protocol)
    assert generated['p_a_kmer_count']['type'] == 'kmer_count' and 'threads' not in generated['p_a_kmer_count']


def test_run_threads_reach_index_building_and_attributes(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(api_module.UProbeAPI, 'build_genome_index',
                        lambda self, threads: seen.setdefault('index', threads))
    def add(df, protocol, genome, report, threads):
        seen['attributes'] = threads
        return df.assign(score=1)
    monkeypatch.setattr(api_module, 'add_attributes', add)
    protocol = tmp_path / 'protocol.yaml'
    protocol.write_text(yaml.safe_dump({
        'name': 't', 'genome': 'g', 'encoding': {}, 'targets': [{'S': 'ACGT' * 20}],
        'extracts': {'target_region': {'source': 'exon', 'length': 40, 'step': 40}},
        'probes': {'p': {'template': '{a}', 'parts': {'a': {'expr': 'rc(target_region)'}}}},
        'attributes': {'score': {'target': 'target_region', 'type': 'gc_content'}},
        'post_process': {'sorts': {'is_ascending': ['score']}}, 'summary': {'report_name': ''}}))
    (tmp_path / 'g.fa').write_text('>1\nACGT\n')
    (tmp_path / 'g.gtf').write_text('1\tt\texon\t1\t4\t.\t+\t.\tgene_name "G"; transcript_id "T";\n')
    genomes = tmp_path / 'genomes.yaml'
    genomes.write_text(yaml.safe_dump({'g': {'fasta': str(tmp_path / 'g.fa'), 'gtf': str(tmp_path / 'g.gtf'),
                                             'align_index': []}}))
    result = CliRunner().invoke(cli, ['run', '-p', str(protocol), '-g', str(genomes), '-o', str(tmp_path / 'out'),
                                      '--threads', '3'])
    assert result.exit_code == 0, result.output
    assert seen == {'index': 3, 'attributes': 3}
