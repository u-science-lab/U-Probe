"""Stage 1: native genome extraction and DAG construction; no reference oligos."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import logging
import re
from pathlib import Path

import yaml
from uprobe.core.api import UProbeAPI

FIELDS = ['gene', 'oligo_type', 'probe_id', 'full_sequence', 'binding', 'chrom',
          'start', 'end', 'strand', 'native_probe_id', 'sites', 'stage']


def generate(config, fasta, output):
    output.mkdir(parents=True, exist_ok=False)
    cfg = yaml.safe_load(config.read_text())
    unused = output / 'unused.gtf'
    unused.touch()  # Required API path, never read by source: genome.
    genomes = {'ce11': {'fasta': str(fasta.resolve()), 'gtf': str(unused.resolve())}}
    counts = []
    for name in ['uprobe.core.api', 'uprobe.core.gen.fun']:
        logging.getLogger(name).setLevel(logging.WARNING)
    with gzip.open(output / 'candidate.csv.gz', 'wt', newline='', compresslevel=1) as handle:
        writer = csv.DictWriter(handle, FIELDS)
        writer.writeheader()
        for target in cfg['targets']:
            protocol = copy.deepcopy(cfg)
            protocol['targets'] = [target]
            api = UProbeAPI(protocol, genomes, output)
            windows = api.generate_target_seqs()
            total = len(windows)
            valid = windows.target_region.str.fullmatch('[ACGT]{40}')
            windows = windows.loc[valid].reset_index(drop=True)
            built = api.construct_probes(windows)
            assert len(built) == len(windows)
            chrom, extent = target.split(':')
            origin = int(extent.split('-')[0])
            unique = {}
            for window, oligos in zip(windows.to_dict('records'), built.to_dict('records')):
                lo, hi = map(int, window['sub_region'].split('-'))
                start, end = origin + lo - 1, origin + hi - 1
                for role, strand in [('forward', '+'), ('reverse', '-')]:
                    seq = oligos[role]
                    assert re.fullmatch('[ACGT]{120}', seq)
                    native = window['probe_id'] + ':' + role
                    site = [start, end, strand, native]
                    if seq in unique:
                        unique[seq]['sites'].append(site)
                    else:
                        unique[seq] = dict(gene=cfg['encoding'][target]['bin_id'], oligo_type='primary',
                            probe_id=native, full_sequence=seq, binding=seq[60:100], chrom=chrom,
                            start=start, end=end, strand=strand, native_probe_id=native, sites=[site], stage='candidate')
            for row in unique.values():
                row['sites'] = json.dumps(row['sites'], separators=(',', ':'))
                writer.writerow(row)
            counts.append(dict(bin=cfg['encoding'][target]['bin_id'], target=target, windows=total,
                               invalid_windows=int((~valid).sum()), before_dedup=2 * len(windows),
                               candidate=len(unique), duplicates=2 * len(windows) - len(unique)))
            if len(counts) % 10 == 0:
                print('generated bins', len(counts), flush=True)
    manifest = {'calls': ['UProbeAPI.generate_target_seqs', 'UProbeAPI.construct_probes'],
                'source': 'genome', 'reference_oligos_read': False,
                'generation_passes': 1, 'bins': counts,
                'candidate_sha256': hashlib.sha256((output / 'candidate.csv.gz').read_bytes()).hexdigest()}
    (output / 'generation.json').write_text(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generate(args.config, args.fasta, args.output)
