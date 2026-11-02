"""Enumerate native spliced mRNA windows and assemble primary probes with YAML/DAG."""
import argparse
import copy
import csv
import gzip
import json
import re
from pathlib import Path

import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.fun import GTF_RUN_FIELDS, read_gtf


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--codebook', type=Path, required=True)
    p.add_argument('--config', type=Path, default=Path(__file__).with_name('primary_design.yaml'))
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--genes', nargs='+', help='Optional smoke-test subset')
    a = p.parse_args()
    codes = list(csv.DictReader(a.codebook.open()))
    genes = sorted({c['target_id'] for c in codes if not a.genes or c['target_id'] in a.genes})
    cfg = yaml.safe_load(a.config.read_text())
    a.output.mkdir(parents=True, exist_ok=True)
    dest = a.output / 'candidates.csv.gz'
    if dest.exists():
        raise RuntimeError('Use a fresh output directory')
    subset = a.output / 'annotation.gtf'
    with a.gtf.open() as src, subset.open('w') as dst:
        for line in src:
            if line.startswith('#'):
                continue
            m = re.search(r'(?:^|;)\s*gene_name\s+"([^"]+)"', line.split('\t')[-1])
            if m and m.group(1) in genes:
                dst.write(line)
    gtf = read_gtf(str(subset), extract_fields=GTF_RUN_FIELDS)
    genome = {'reference': dict(fasta=str(a.fasta.resolve()), gtf=str(subset.resolve()))}
    details = []
    with gzip.open(dest, 'wt', newline='', compresslevel=1) as handle:
        writer = csv.DictWriter(handle, ['target_id', 'oligo_type', 'probe_id', 'full_sequence', 'stage'])
        writer.writeheader()
        for gene in genes:
            annotation = gtf[gtf.gene_name.eq(gene)]
            txs = sorted(annotation.loc[annotation.type.eq('exon'), 'transcript_id'].dropna().unique())
            if not txs:
                raise RuntimeError('Missing annotation: ' + gene)
            config = copy.deepcopy(cfg)
            config['targets'] = [gene]
            readouts = sorted({c['readout'] for c in codes if c['target_id'] == gene})
            # Instantiate the YAML primary template, changing only its fixed readout.
            config['probes'] = {}
            for i, readout in enumerate(readouts):
                assert re.fullmatch('[ACGT]{30}', readout)
                primary = copy.deepcopy(cfg['probes']['primary'])
                primary['parts']['readout']['template'] = readout
                config['probes'][f'primary_{i}'] = primary
            api = UProbeAPI(config, genome, a.output / 'engine', require_genome=True)
            api._gtf = annotation
            seen, count = set(), 0
            for tx in txs:
                api.protocol['extracts']['target_region']['transcripts'] = [tx]
                windows = api.generate_target_seqs()
                for offset in range(0, len(windows), 300):
                    rows = []
                    for origin, row in enumerate(windows.iloc[offset:offset+300].to_dict('records'), offset):
                        seq = row['target_region'].upper()
                        if seq in seen or set(seq) - set('ACGT'):
                            continue
                        seen.add(seq)
                        row['target_region'] = seq
                        row['recovery_origin'] = f'{gene}:{tx}:window{origin}'
                        rows.append(row)
                    if not rows:
                        continue
                    built = api.construct_probes(pd.DataFrame(rows))
                    for role in config['probes']:
                        for row, seq in zip(rows, built[role]):
                            assert re.fullmatch('[ACGT]{65}', seq)
                            writer.writerow(dict(target_id=gene, oligo_type='primary',
                                probe_id=row['recovery_origin']+':'+role,
                                full_sequence=seq, stage='candidate'))
                            count += 1
            detail = dict(target_id=gene, transcripts=len(txs), unique_windows=len(seen), candidates=count)
            details.append(detail)
            print(json.dumps(detail), flush=True)
    (a.output / 'generation.json').write_text(json.dumps(dict(details=details,
        basic_constraints='ACGT only; deduplicate identical windows within each gene; step 1; all annotated transcripts',
        final_recovery='Not evaluated'), indent=2) + '\n')


if __name__ == '__main__':
    main()
