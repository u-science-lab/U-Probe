#!/usr/bin/env python3
'''Native YAML candidate construction; never reads reference oligos or binding arms.'''
import argparse
import copy
import csv
import gzip
import hashlib
import json
import multiprocessing as mp
import re
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.fun import read_gtf, GTF_RUN_FIELDS

CONFIG = GTF = GENOMES = OUTPUT = None
FIELDS = ['stage', 'gene', 'probe_id', 'oligo_type', 'full_sequence', 'modification_5p']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def run_gene(gene):
    started = time.time()
    cfg = copy.deepcopy(CONFIG)
    cfg['targets'] = [gene]
    api = UProbeAPI(cfg, GENOMES, OUTPUT/'engine', require_genome=True)
    api._gtf = GTF[GTF.gene_name.eq(gene)]
    ids = sorted(api._gtf.loc[api._gtf.type.eq('exon'), 'transcript_id'].dropna().unique())
    if not ids:
        raise RuntimeError('No annotated transcripts for '+gene)
    n = cfg['encoding'][gene]['variant_count']
    if n != len(cfg['encoding'][gene]['variants']) or n < 1:
        raise ValueError('Invalid encoding variant count')
    seen = set()
    sites = rejected = count = 0
    with gzip.open(OUTPUT/'candidates_by_gene'/(gene+'.csv.gz'), 'wt', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for tx in ids:
            api.protocol['extracts']['target_region']['transcripts'] = [tx]
            df = api.generate_target_seqs()
            if df.empty:
                raise RuntimeError('No windows generated for '+tx)
            for offset in range(0, len(df), 300):
                rows = []
                for row in df.iloc[offset:offset+300].to_dict('records'):
                    seq = row['target_region'].upper()
                    if seq in seen:
                        continue
                    seen.add(seq)
                    if set(seq)-set('ACGT'):
                        rejected += 1
                        continue
                    row['target_region'] = seq
                    # Repeat each site for every YAML encoding variant, with no site/codebook association.
                    rows.extend([row]*n)
                    sites += 1
                if not rows:
                    continue
                built = api.construct_probes(pd.DataFrame(rows))
                if len(built) != len(rows):
                    raise RuntimeError('Native construction row mismatch')
                for typ in cfg['probes']:
                    for seq in built[typ]:
                        if not seq or set(seq)-set('ACGT'):
                            raise RuntimeError('Invalid native probe DNA')
                        count += 1
                        writer.writerow(dict(stage='candidate', gene=gene, probe_id=gene+'_'+str(count),
                            oligo_type=typ, full_sequence=seq,
                            modification_5p=cfg['probes'][typ].get('modifications', {}).get('5prime', 'none')))
    detail = dict(gene=gene, transcripts=len(ids), unique_sites=sites, encoding_variants=n,
                  candidate_records=count, non_acgt_windows=rejected, seconds=round(time.time()-started, 3))
    print(json.dumps(detail), flush=True)
    return detail


def main():
    global CONFIG, GTF, GENOMES, OUTPUT
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--genes', nargs='+', help='Optional YAML target subset for one small end-to-end test')
    args = p.parse_args()
    CONFIG = yaml.safe_load(args.config.read_text())
    if args.workers < 1:
        p.error('--workers must be positive')
    if CONFIG['recovery_execution']['encoding_variants'] != 'cartesian':
        p.error('This runner supports the declared Cartesian encoding policy only')
    if args.genes:
        if not set(args.genes) <= set(CONFIG['targets']):
            p.error('--genes must be YAML targets')
        CONFIG['targets'] = args.genes
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', g) for g in CONFIG['targets']):
        p.error('Targets must be gene names, never target sequences')
    OUTPUT = args.output.resolve()
    if OUTPUT.exists() and any(OUTPUT.iterdir()):
        p.error('Output must be empty')
    (OUTPUT/'candidates_by_gene').mkdir(parents=True)
    GENOMES = {CONFIG['genome']: dict(fasta=str(args.fasta.resolve()), gtf=str(args.gtf.resolve()))}
    subset = OUTPUT/'target_annotation.gtf'
    genes = set(CONFIG['targets'])
    with args.gtf.open() as src, subset.open('w') as dst:
        for line in src:
            if line.startswith('#'):
                continue
            match = re.search(r'(?:^|;)\s*gene_name\s+\x22([^\x22]+)\x22', line.split('\t')[-1])
            if match and match.group(1).upper() in genes:
                dst.write(line)
    GTF = read_gtf(str(subset), extract_fields=GTF_RUN_FIELDS)
    GTF['gene_name'] = GTF.gene_name.str.upper()
    if genes-set(GTF.gene_name):
        raise RuntimeError('Missing annotated targets: '+repr(genes-set(GTF.gene_name)))
    started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('fork')) as pool:
        details = list(pool.map(run_gene, CONFIG['targets'], chunksize=1))
    with gzip.open(OUTPUT/'candidates.csv.gz', 'wb', compresslevel=1) as out:
        header = None
        for gene in CONFIG['targets']:
            with gzip.open(OUTPUT/'candidates_by_gene'/(gene+'.csv.gz'), 'rb') as src:
                current = src.readline()
                if header is None:
                    header = current
                    out.write(header)
                elif current != header:
                    raise RuntimeError('Candidate headers differ')
                shutil.copyfileobj(src, out, 4*1024*1024)
    provenance = dict(config_sha256=sha256(args.config), fasta_sha256=sha256(args.fasta),
                      gtf_sha256=sha256(args.gtf), runner_sha256=sha256(__file__))
    result = dict(candidate_records=sum(d['candidate_records'] for d in details),
                  genes=details, final_stage='not_executed', seconds=round(time.time()-started, 3), **provenance)
    (OUTPUT/'run_details.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
