import ast
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

SOURCE = Path(__file__).parents[1] / 'uprobe/core/tools/__init__.py'

class TranscriptIndexTests(unittest.TestCase):
    def setUp(self):
        tree = ast.parse(SOURCE.read_text(encoding='utf-8'))
        tree.body = [n for n in tree.body if not isinstance(n, (ast.Import, ast.ImportFrom))]
        self.ns = dict(Path=Path, get_logger=lambda _: Mock(), __name__='tools')
        exec(compile(tree, str(SOURCE), 'exec'), self.ns)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / 'toy.fa'
        self.fasta.write_text('>chr1\nACGT\n')
        self.transcript = self.root / 'transcript.fa'
        self.transcript.write_text('>gene_transcript\nACGT\n')
        def bowtie(source, prefix, threads):
            for suffix in self.ns['_BOWTIE2_SMALL_SUFFIXES']:
                Path(str(prefix) + suffix).write_text('index')
        def blast(source, prefix, title):
            for suffix in self.ns['_BLAST_SUFFIXES']:
                Path(str(prefix) + suffix).write_text('index')
        self.ns['build_bowtie2_index'] = Mock(side_effect=bowtie)
        self.ns['build_blast_db'] = Mock(side_effect=blast)
        self.ns['extract_trans_seqs'] = Mock()

    def build(self, aligner):
        return self.ns['build_transcripts_index'](self.root/'toy.gtf', self.fasta,
            self.root/f'{aligner}_transcript', threads=2, aligner=aligner)

    def test_blast_ignores_misplaced_bowtie_files(self):
        directory = self.root/'blast_transcript'
        directory.mkdir()
        for suffix in self.ns['_BOWTIE2_SMALL_SUFFIXES']:
            (directory/('toy'+suffix)).write_text('wrong format')
        self.build('blast')
        self.ns['build_blast_db'].assert_called_once()
        self.ns['build_bowtie2_index'].assert_not_called()
        self.assertEqual(self.ns['build_blast_db'].call_args.args[0], self.transcript)
        self.build('blast')
        self.ns['build_blast_db'].assert_called_once()

    def test_bowtie_uses_parent_fasta_and_skips_complete_index(self):
        self.build('bowtie2')
        self.ns['build_bowtie2_index'].assert_called_once_with(
            self.transcript, self.root/'bowtie2_transcript/toy', 2)
        self.build('bowtie2')
        self.ns['build_bowtie2_index'].assert_called_once()

    def test_empty_fasta_is_rejected_before_build(self):
        self.transcript.write_text('')
        with self.assertRaisesRegex(ValueError, 'empty or missing'):
            self.build('blast')
        self.ns['extract_trans_seqs'].assert_called_once()
        self.ns['build_blast_db'].assert_not_called()

    def test_wrong_named_directory_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'must use blast_transcript'):
            self.ns['build_transcripts_index'](self.root/'toy.gtf', self.fasta,
                self.root/'bowtie2_transcript', aligner='blast')

if __name__ == '__main__':
    unittest.main()
