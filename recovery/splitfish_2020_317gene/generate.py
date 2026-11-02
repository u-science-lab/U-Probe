"""Native split-FISH extraction and DAG assembly. Does not read reference oligos."""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import logging
import multiprocessing as mp
import re
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.fun import read_gtf, GTF_RUN_FIELDS

CONFIG = GTF = GENOMES = OUTPUT = None
FIELDS = ['gene', 'oligo_type', 'probe_id', 'pair_id', 'full_sequence', 'stage']


def run_gene(gene):
    cfg = copy.deepcopy(CONFIG)
    cfg['targets'] = [gene]
    api = UProbeAPI(cfg, GENOMES, OUTPUT/'engine'/gene, require_genome=True)
    api._gtf = GTF[GTF.gene_name.eq(gene)]
    txs = sorted(api._gtf.loc[api._gtf.type.eq('exon'), 'transcript_id'].dropna().unique())
    seen = set()
    count = rejected = 0
    with gzip.open(OUTPUT/'shards'/(gene+'.csv.gz'), 'wt', newline='', compresslevel=1) as f:
        writer = csv.DictWriter(f, FIELDS)
        writer.writeheader()
        for tx in txs:
            api.protocol['extracts']['target_region']['transcripts'] = [tx]
            windows = api.generate_target_seqs()
            for offset in range(0, len(windows), 300):
                rows, origins = [], []
                for row in windows.iloc[offset:offset+300].to_dict('records'):
                    seq = row['target_region'].upper()
                    if seq in seen:
                        continue
                    seen.add(seq)
                    if set(seq)-set('ACGT'):
                        rejected += 1
                        continue
                    row['target_region'] = seq
                    # Both published on-bits at every independent site, without ordinal selection.
                    for bit in cfg['encoding'][gene]['bits']:
                        rows.append(row)
                        origins.append(f"{tx}:{row['probe_id']}:B{bit}")
                if not rows:
                    continue
                built = api.construct_probes(pd.DataFrame(rows))
                assert len(built) == len(origins)
                for pair, item in zip(origins, built.to_dict('records')):
                    # The two roles stay in their native construction row throughout evaluation.
                    for role in ('p1', 'p2'):
                        seq = item[role]
                        assert re.fullmatch('[ACGT]{88}', seq)
                        writer.writerow(dict(gene=gene, oligo_type=role, probe_id=pair+':'+role,
                                             pair_id=pair, full_sequence=seq, stage='candidate'))
                        count += 1
    detail = dict(gene=gene, transcripts=len(txs), unique_windows=len(seen),
                  non_acgt_windows=rejected, candidate_records=count)
    print(json.dumps(detail),flush=True)
    return detail


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    global CONFIG, GTF, GENOMES, OUTPUT
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--fasta',type=Path,required=True)
    p.add_argument('--gtf',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--genes',nargs='+')
    p.add_argument('--workers',type=int,default=8)
    a = p.parse_args()
    for name in ('uprobe.core.api', 'uprobe.core.gen.fun'):
        logging.getLogger(name).setLevel(logging.WARNING)
    CONFIG = yaml.safe_load(a.config.read_text())
    if a.genes:
        assert set(a.genes) <= set(CONFIG['targets'])
        CONFIG['targets'] = a.genes
    OUTPUT = a.output.resolve()
    OUTPUT.mkdir(parents=True,exist_ok=True)
    assert not any(OUTPUT.iterdir()), 'Use a fresh output directory'
    (OUTPUT/'shards').mkdir()
    subset = OUTPUT/'annotation.gtf'
    genes = set(CONFIG['targets'])
    with a.gtf.open() as src, subset.open('w') as dst:
        for line in src:
            if line.startswith('#'):
                continue
            m = re.search(r'gene_name "([^"]+)"',line)
            if m and m[1] in genes:
                dst.write(line)
    GTF = read_gtf(str(subset),extract_fields=GTF_RUN_FIELDS)
    GENOMES = {'reference':dict(fasta=str(a.fasta.resolve()),gtf=str(subset))}
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=mp.get_context('fork')) as pool:
        details = list(pool.map(run_gene,CONFIG['targets']))
    with gzip.open(OUTPUT/'candidates.csv.gz','wb',compresslevel=1) as out:
        for i,gene in enumerate(CONFIG['targets']):
            with gzip.open(OUTPUT/'shards'/(gene+'.csv.gz'),'rb') as src:
                header = src.readline()
                if i == 0:
                    out.write(header)
                shutil.copyfileobj(src,out,4*1024*1024)
    result = dict(genes=details,candidate_records=sum(d['candidate_records'] for d in details),
                  fasta_sha256=sha(a.fasta),gtf_sha256=sha(a.gtf),config_sha256=sha(a.config),
                  final_recovery='Not evaluated',
                  policy='All M4 transcripts; step 1; unique 52-nt ACGT windows within gene; both on-bits per window. Pair IDs namespace unchanged native window IDs by transcript and bit.')
    (OUTPUT/'generation.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
