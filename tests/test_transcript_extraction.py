import ast
from pathlib import Path
import os,re,tempfile,unittest
from unittest.mock import Mock
import pandas as pd
from pyfaidx import Fasta

class TranscriptExtractionTests(unittest.TestCase):
    def setUp(self):
        source=Path(__file__).parents[1]/'uprobe/core/gen/fun.py'
        tree=ast.parse(source.read_text())
        tree.body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('read_gtf','extract_trans_seqs')]
        self.ns=dict(pd=pd,os=os,re=re,tempfile=tempfile,Fasta=Fasta,log=Mock(),
            GTF_FIELDS=['chr','source','type','start','end','score','strand','score','info'])
        import typing
        self.ns['t']=typing
        self.ns['reverse_complement']=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
        exec(compile(tree,str(source),'exec'),self.ns)
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.fa=self.root/'genome.fa'; self.gtf=self.root/'genome.gtf'; self.out=self.root/'transcript.fa'
    def run_extract(self,header,rows):
        self.fa.write_text('>'+header+'\nAAAACCCCGGGGTTTT\n')
        self.gtf.write_text('\n'.join(rows)+'\n')
        self.ns['extract_trans_seqs'](self.gtf,self.fa,self.out)
    def exon(self,chrom,start,end,attrs,strand='+'):
        return f'{chrom}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}'
    def test_ncbi_reference_and_geneid(self):
        self.run_extract('NC_000067.7',[self.exon('NC_000067.7',1,4,'transcript_id "XR_1"; db_xref "GeneID:42";')])
        self.assertEqual(self.out.read_text(),'>42_XR_1\nAAAA\n')
    def test_chr_alias_negative_strand_and_gene_name_fallback(self):
        attrs='gene_name "GeneA"; transcript_id "tx1";'
        self.run_extract('chr1',[self.exon('1',1,4,attrs,'-'),self.exon('1',13,16,attrs,'-')])
        self.assertEqual(self.out.read_text(),'>GeneA_tx1\nAAAATTTT\n')
    def test_explicit_gene_id_takes_precedence(self):
        self.run_extract('chr1',[self.exon('chr1',1,4,'gene_id "ENSG1"; transcript_id "tx"; db_xref "GeneID:42";')])
        self.assertTrue(self.out.read_text().startswith('>ENSG1_tx\n'))
    def test_no_matching_reference_does_not_overwrite_file(self):
        self.out.write_text('previous valid output')
        with self.assertRaisesRegex(ValueError,'No usable transcript exons'):
            self.run_extract('chr1',[self.exon('NC_9',1,4,'gene_id "g"; transcript_id "t";')])
        self.assertEqual(self.out.read_text(),'previous valid output')
    def test_missing_gene_identifier_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'Missing gene_id'):
            self.run_extract('NC_1',[self.exon('NC_1',1,4,'transcript_id "t";')])
        self.assertFalse(self.out.exists())

if __name__=='__main__': unittest.main()

