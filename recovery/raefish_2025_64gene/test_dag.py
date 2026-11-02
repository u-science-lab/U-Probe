#!/usr/bin/env python3
"""Check upstream dependencies with the native DAG without reference matching."""
import argparse
import copy
import json
from pathlib import Path

import yaml
from uprobe.core.gen.probe import DAG


def build(config, context):
    dag = DAG()
    dag.from_config(config)
    dag.run(context)
    assert all(node.done for node in dag.nodes)
    return {node.name: node.result for node in dag.nodes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('design', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.design.read_text())
    target = config['targets'][0]
    context = dict(target_region='A' * 30 + 'C' * 30, target=target,
                   encoding=config['encoding'], probe_index=0,
                   target_parts=dict(upstream='A' * 30, gap='', downstream='C' * 30))
    base = build(config, context)

    changed = copy.deepcopy(context)
    changed['target_parts']['upstream'] = 'G' + 'A' * 29
    result = build(config, changed)
    assert base['padlock'] != result['padlock']
    assert base['encoding_probe'] != result['encoding_probe']
    assert base['priming'] == result['priming']

    changed = copy.deepcopy(context)
    changed['target_parts']['downstream'] = 'T' + 'C' * 29
    result = build(config, changed)
    assert base['padlock'] == result['padlock']
    assert base['encoding_probe'] == result['encoding_probe']
    assert base['priming'] != result['priming']

    changed = copy.deepcopy(context)
    readout = changed['encoding'][target]['R1']
    changed['encoding'][target]['R1'] = ('C' if readout[0] == 'A' else 'A') + readout[1:]
    result = build(config, changed)
    assert base['padlock'] == result['padlock']
    assert base['priming'] == result['priming']
    assert base['encoding_probe'] != result['encoding_probe']

    changed_config = copy.deepcopy(config)
    changed_config['probes']['padlock']['parts']['common5']['expr'] = "'ACTCTCACTTAGGGTTCCAAC'"
    result = build(changed_config, context)
    assert base['padlock'] != result['padlock']
    assert base['priming'] != result['priming']

    evidence = dict(native_dag_dependency_checks=4, passed=4, checks=[
        'padlock binding changes encoding',
        'priming binding remains independently linked',
        'gene readout changes encoding only',
        'padlock RCA-compatible common changes priming',
    ])
    args.output.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(evidence))


if __name__ == '__main__':
    main()
