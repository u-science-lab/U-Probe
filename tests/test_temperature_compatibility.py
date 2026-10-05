import ast
from pathlib import Path
import secrets
import unittest
import warnings
import pandas as pd
import primer3

ROOT = Path(__file__).parents[1]

def function(path, name, namespace):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
    exec(compile(tree, str(path), 'exec'), namespace)
    return namespace[name]

class TemperatureCompatibilityTests(unittest.TestCase):
    def test_fisheye_values_and_precision_reach_attribute_table(self):
        namespace = dict(primer3=primer3, pd=pd, secrets=secrets)
        calc = function(ROOT/'uprobe/core/attributes/_attributes.py', 'cal_temp', namespace)
        add = function(ROOT/'uprobe/core/attributes/__init__.py', 'add_attributes', namespace)
        for sequence in ('ACGTACGTACGTA', 'GCGCATATGCGCA', 'ACGTNCGTACGTA'):
            with self.subTest(sequence=sequence), warnings.catch_warnings():
                warnings.simplefilter('ignore')
                expected = primer3.calcTm(sequence, dv_conc=0, formamide_conc=0)
                self.assertEqual(calc(sequence), expected)
                result = add(pd.DataFrame({'target_region':[sequence]}), {'attributes':{'tm':{'target':'target_region','type':'annealing_temperature'}}}, {})
                self.assertEqual(result['tm'].iloc[0], expected)

    def test_known_fisheye_temperature(self):
        calc = function(ROOT/'uprobe/core/attributes/_attributes.py', 'cal_temp', dict(primer3=primer3))
        self.assertAlmostEqual(calc('GTATGACAATGAA'), 27.278686288556173, places=10)

if __name__ == '__main__':
    unittest.main()
