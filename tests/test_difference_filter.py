import pandas as pd
from uprobe.core.attributes import add_attributes
from uprobe.core.process import post_process


def test_absolute_difference_filter_uses_both_computed_values(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    protocol = {
        'attributes': {
            'a_len': {'target': 'a', 'type': 'length'},
            'b_len': {'target': 'b', 'type': 'length'},
            'length_delta': {'type': 'difference', 'of': ['a_len', 'b_len']},
        },
        'post_process': {'filters': {'length_delta': {'condition': 'length_delta <= 2'}}},
    }
    probes = pd.DataFrame({'a': ['AAA', 'AAAAA', 'AA'], 'b': ['AAAAA', 'AAA', 'AAAAA']})
    measured = add_attributes(probes, protocol, {})
    assert measured['length_delta'].tolist() == [2, 2, 3]
    assert post_process(measured, protocol).index.tolist() == [0, 1]
