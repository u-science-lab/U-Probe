import unittest
import csv
import gzip
import tempfile
from pathlib import Path
from sequence_recovery import evaluate, normalize_sequence, iter_csv


def oligo(gene, pair, role, dna, **extra):
    return dict(gene=gene, pair_id=pair, oligo_type=role, full_sequence=dna, **extra)


class RecoveryTests(unittest.TestCase):
    def test_modifications_only_remove_explicit_tokens(self):
        self.assertEqual(normalize_sequence(' /5Phos/ a c\ngt '), ('ACGT', ['5PHOS']))
        for sequence in ('5PhosACGT', 'ACNT', 'AC-GT', 'AC/5Phos', '', 'AC/GT/AA'):
            with self.assertRaises(ValueError):
                normalize_sequence(sequence)

    def test_exact_gene_role_direction_and_dedup(self):
        reference = [oligo('A', 'a', 'padlock', '/5Phos/ACGA'),
                     oligo('A', 'a', 'padlock', 'ACGA'),
                     oligo('A', 'a', 'primer', 'TTCA'),
                     oligo('B', 'b', 'padlock', 'ACGA')]
        generated = [oligo('a', 'new', 'PADLOCK', 'acga', stage='candidate'),
                     oligo('A', 'new', 'padlock', 'ACGA', stage='candidate'),
                     oligo('A', 'new', 'primer', 'TGAA', stage='candidate'),
                     oligo('C', 'new', 'padlock', 'ACGA', stage='candidate')]
        result = evaluate(reference, generated)
        summary = result['summary'][0]
        self.assertEqual(summary['reference_sequences'], 3)
        self.assertEqual(summary['matched_sequences'], 1)
        self.assertEqual(summary['duplicate_reference_records'], 1)
        self.assertAlmostEqual(summary['SRR_pct'], 100/3)
        self.assertEqual(summary['TRR_pct'], 50)
        self.assertEqual(summary['MTR_pct'], 25)
        self.assertEqual(summary['CPSR_pct'], 0)

    def test_pair_requires_same_generated_group_and_stages_separate(self):
        reference = [oligo('A', 'ref', 'padlock', 'ACGA'), oligo('A', 'ref', 'primer', 'TTCA')]
        generated = [oligo('A', 'one', 'padlock', 'ACGA', stage='candidate'),
                     oligo('A', 'two', 'primer', 'TTCA', stage='candidate'),
                     oligo('A', 'complete', 'padlock', 'ACGA', stage='final'),
                     oligo('A', 'complete', 'primer', '/5Phos/TTCA', stage='final')]
        result = evaluate(reference, generated)
        summaries = {r['stage']: r for r in result['summary'] if r['oligo_type'] == 'all'}
        self.assertEqual(summaries['candidate']['SRR_pct'], 100)
        self.assertEqual(summaries['candidate']['CPSR_pct'], 0)
        self.assertEqual(summaries['final']['CPSR_pct'], 100)
        self.assertEqual(result['pair_matches'][0]['all_members_individually_matched'], True)

    def test_exclusion_and_generated_validation(self):
        reference = [oligo('A', 'p', 'padlock', 'ACGA'), oligo('', '', '', '')]
        result = evaluate(reference, [])
        self.assertEqual(result['summary'][0]['SRR_pct'], 0)
        self.assertEqual(result['summary'][0]['excluded_reference_records'], 1)
        with self.assertRaises(ValueError):
            evaluate(reference, [oligo('A', 'p', 'padlock', 'ACGN')])
        with self.assertRaises(ValueError):
            evaluate([oligo('A', 'p', 'padlock', 'ACGA', dna_sequence='ACGT')], [])

    def test_wide_pair_consistency(self):
        reference = [oligo('A', 'p', 'padlock', 'ACGA'), oligo('A', 'p', 'primer', 'TTCA')]
        pair = dict(gene='A', pair_id='p', padlock_full='/5Phos/ACGA', primer_dna='TTCA')
        evaluate(reference, [], [pair])
        with self.assertRaises(ValueError):
            evaluate(reference, [], [dict(pair, primer_dna='AAAA')])

    def test_non_method_specific_role_names(self):
        reference = [oligo('Region1', 'set1', 'encoding', 'ACGA'), oligo('Region1', 'set1', 'readout', 'TTCA')]
        pair = dict(gene='Region1', pair_id='set1', encoding_dna='ACGA', readout_dna='TTCA')
        result = evaluate(reference, reference, [pair], version='example_v1')
        self.assertEqual(result['summary'][0]['CPSR_pct'], 100)
        self.assertEqual(result['pair_matches'][0]['method_version'], 'example_v1')
        self.assertNotIn('padlock_matched', result['pair_matches'][0])

    def test_gzip_stream_and_nonreference_outputs(self):
        reference = [oligo('A', 'p', 'padlock', 'ACGA')]
        generated = [oligo('A', 'p', 'padlock', 'ACGA'), oligo('A', 'new', 'padlock', 'TTAA')]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'generated.csv.gz'
            with gzip.open(path, 'wt', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(generated[0]))
                writer.writeheader()
                writer.writerows(generated)
            result = evaluate(reference, iter_csv(path))
        self.assertEqual(result['summary'][0]['generated_records'], 2)
        self.assertEqual(result['summary'][0]['generated_unique_reference_hits'], 1)
        self.assertEqual(result['summary'][0]['SRR_pct'], 100)
        self.assertEqual(result['summary'][0]['CPSR_pct'], '')
        self.assertEqual(result['per_gene'][0]['CPSR_pct'], '')
        self.assertEqual(result['pair_matches'][0]['cpsr_applicable'], False)


if __name__ == '__main__':
    unittest.main()
