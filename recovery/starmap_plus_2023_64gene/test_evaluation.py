import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent/'common'))
from evaluate_recovery import evaluate_files
from verify_results import verify


def write(path, records):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.report = self.root/'report'
        self.ref_path, self.pair_path, self.generated = [self.root/name for name in ('oligos.csv','pairs.csv','generated.csv.gz')]
        refs, pairs = [], []
        sequences = [('A','A1','AACCGGTT','TTGGCCAA'),('A','A2','ACTGACTG','TCAGTCAG'),('B','B1','AGCTAGCT','CGATCGAT')]
        for gene,pair,pad,pri in sequences:
            pairs.append(dict(gene=gene,pair_id=pair,padlock_full='/5Phos/'+pad,primer_full=pri))
            for role,seq in [('padlock',pad),('primer',pri)]:
                refs.append(dict(gene=gene,pair_id=pair,oligo_id=pair+'_'+role,oligo_type=role,full_sequence='/5Phos/'+seq if role=='padlock' else seq))
        generated = []
        for gene,pair,pad,pri in [('A','new1','AACCGGTT','TTGGCCAA'),('A','new2','ACTGACTG','AAAAAAAA'),
                                  ('B','new3','AGCTAGCT','CCCCCCCC'),('B','new4','GGGGGGGG','CGATCGAT')]:
            for role,seq in [('padlock',pad),('primer',pri)]:
                generated.append(dict(stage='candidate',gene=gene,pair_id=pair,probe_id=pair+'_'+role,oligo_type=role,full_sequence=seq))
        write(self.ref_path,refs)
        write(self.pair_path,pairs)
        write(self.generated,generated)
        evaluate_files(self.ref_path,self.pair_path,self.generated,self.report)

    def tearDown(self):
        self.temp.cleanup()

    def check(self, full=False):
        return verify(self.ref_path,self.pair_path,self.report,self.generated if full else None)

    def mutate(self, name, key, value, index=0):
        path = self.report/name
        with path.open(newline='') as handle:
            records = list(csv.DictReader(handle))
        records[index][key] = value
        write(path,records)

    def test_compact_metrics_and_explicit_scope(self):
        result = self.check()
        all_role = result['recomputed_metrics'][0]
        self.assertAlmostEqual(all_role['SRR_pct'],100*5/6)
        self.assertEqual(all_role['MTR_pct'],87.5)
        self.assertAlmostEqual(all_role['CPSR_pct'],100/3)
        self.assertFalse(result['full_generated_scan_performed'])
        self.assertEqual(result['unverified_fields'],['generated_records'])

    def test_optional_independent_generated_scan(self):
        result = self.check(full=True)
        self.assertTrue(result['full_generated_scan_performed'])
        self.assertEqual(result['generated_records'],{'candidate':8})

    def test_rejects_wrong_full_dna_or_gene_role(self):
        self.mutate('sequence_matches.csv','full_sequence','ACGT')
        with self.assertRaises(ValueError):
            self.check()

    def test_rejects_cross_pair_link(self):
        self.mutate('pair_matches.csv','generated_pair_ids','new2')
        with self.assertRaises(ValueError):
            self.check()

    def test_rejects_wrong_metrics(self):
        self.mutate('summary.csv','MTR_pct','100')
        with self.assertRaises(ValueError):
            self.check()

    def test_rejects_unlinked_probe_ids(self):
        self.mutate('sequence_matches.csv','generated_probe_ids','unlinked_padlock')
        with self.assertRaises(ValueError):
            self.check()

    def test_full_scan_detects_changed_generated_dna(self):
        with gzip.open(self.generated,'rt',newline='') as handle:
            generated = list(csv.DictReader(handle))
        generated[0]['full_sequence']='TTTTTTTT'
        write(self.generated,generated)
        with self.assertRaises(ValueError):
            self.check(full=True)


if __name__ == '__main__':
    unittest.main()
