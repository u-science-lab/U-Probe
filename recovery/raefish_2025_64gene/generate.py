#!/usr/bin/env python3
"""Enumerate transcript windows and assemble RAEFISH with the native U-Probe DAG.

Only target identities, fixed functional sequences and gene codebooks are inputs.
Published target-binding sequences and reference oligos are never read.
"""
import argparse
import copy
import csv
import gzip
import hashlib
import json
import logging
import multiprocessing as mp
import re
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.fun import read_gtf, GTF_RUN_FIELDS

METHOD = 'raefish_2025_64gene'
CONFIG = GTF = GENOMES = OUTPUT = TARGETS = None
ROLES = ('padlock', 'priming', 'encoding')
FIELDS = ['method_version', 'target_id', 'probe_set_id', 'oligo_type', 'full_sequence',
          'transcript_id', 'target_start', 'target_end', 'gap', 'orientation', 'target_region',
          'padlock_binding_start', 'priming_binding_start', 'functional_encoding']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as src:
        for block in iter(lambda: src.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run_gene(target_id):
    started = time.time()
    symbol = TARGETS[target_id]['gene_symbol']
    cfg = copy.deepcopy(CONFIG)
    cfg['targets'] = [symbol]
    api = UProbeAPI(cfg, GENOMES, OUTPUT / 'engine', require_genome=True)
    api._gtf = GTF[GTF.gene_id.str.split('.').str[0].eq(target_id)].copy()
    # U-Probe's gene lookup uses gene_name; the stable Ensembl ID selected this subset.
    api._gtf['gene_name'] = symbol
    txs = sorted(api._gtf.loc[api._gtf.type.eq('exon'), 'transcript_id'].dropna().unique())
    contigs = {line.split('\t')[0].replace('chr', '') for line in Path(GENOMES[cfg['genome']]['fasta'] + '.fai').read_text().splitlines()}
    unavailable_txs = sorted(api._gtf.loc[api._gtf.type.eq('exon') & ~api._gtf['chr'].isin(contigs), 'transcript_id'].dropna().unique())
    txs = [tx for tx in txs if tx not in unavailable_txs]
    count = rejected = sites = 0
    seen = set()
    with gzip.open(OUTPUT / 'candidates_by_gene' / (target_id + '.csv.gz'), 'wt', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for tx in txs:
            api.protocol['extracts']['target_region']['transcripts'] = [tx]
            df = api.generate_target_seqs()
            if df.empty:
                exons = api._gtf[api._gtf.type.eq('exon') & api._gtf.transcript_id.eq(tx)]
                if int((exons.end - exons.start + 1).sum()) < 60:
                    continue
                raise RuntimeError('No windows generated: ' + tx)
            for offset in range(0, len(df), 300):
                rows = []
                for row in df.iloc[offset:offset + 300].to_dict('records'):
                    seq = row['target_region'].upper()
                    # Deduplicate complete 60/61-nt windows, never isolated binding arms.
                    if seq in seen:
                        continue
                    seen.add(seq)
                    if set(seq) - set('ACGT'):
                        rejected += 1
                        continue
                    row['target'] = target_id
                    row['target_region'] = seq
                    for part in ('upstream', 'gap', 'downstream'):
                        row['target_parts.' + part] = row['target_parts.' + part].upper()
                    rows.append(row)
                if not rows:
                    continue
                # Both legal adjacent directions are declared before any matching.
                for orientation in ('padlock_upstream', 'padlock_downstream'):
                    api.protocol['encoding'][target_id]['padlock_part'] = 'upstream' if orientation == 'padlock_upstream' else 'downstream'
                    api.protocol['encoding'][target_id]['priming_part'] = 'downstream' if orientation == 'padlock_upstream' else 'upstream'
                    built = api.construct_probes(pd.DataFrame(rows))
                    if len(built) != len(rows):
                        raise RuntimeError('Native DAG row count mismatch')
                    for origin, native in zip(rows, built.to_dict('records')):
                        sites += 1
                        start = int(origin['target_start'])
                        end = int(origin['target_end'])
                        gap = len(origin['target_parts.gap'])
                        pid = f'{target_id}|{tx}|{start}:{end}|{orientation}'
                        common = dict(method_version=METHOD, target_id=target_id, probe_set_id=pid,
                                      transcript_id=tx, target_start=start, target_end=end, gap=gap,
                                      orientation=orientation, target_region=origin['target_region'], functional_encoding=native['encoding_probe.functional'],
                                      padlock_binding_start=start if orientation == 'padlock_upstream' else start + 30 + gap,
                                      priming_binding_start=start + 30 + gap if orientation == 'padlock_upstream' else start)
                        for role, length in zip(ROLES, (80, 107, 150)):
                            sequence = native['encoding_probe' if role == 'encoding' else role]
                            if len(sequence) != length or set(sequence) - set('ACGT'):
                                raise RuntimeError('Invalid native ' + role)
                            writer.writerow(dict(common, oligo_type=role, full_sequence=sequence))
                            count += 1
    detail = dict(target_id=target_id, gene_symbol=symbol, transcripts=len(txs),
                  native_probe_sets=sites, candidate_records=count,
                  non_acgt_windows=rejected, missing_annotation=not txs, unavailable_contig_transcripts=unavailable_txs,
                  seconds=round(time.time() - started, 3))
    print(json.dumps(detail), flush=True)
    return detail


def main():
    global CONFIG, GTF, GENOMES, OUTPUT, TARGETS
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--targets', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--genes', nargs='+', help='Optional stable-ID subset for smoke testing')
    args = p.parse_args()
    CONFIG = yaml.safe_load(args.config.read_text())
    with args.targets.open() as f:
        TARGETS = {r['target_id']: r for r in csv.DictReader(f)}
    selected = args.genes or CONFIG['targets']
    if args.workers < 1 or not set(selected) <= set(CONFIG['targets']):
        p.error('Invalid workers or target subset')
    OUTPUT = args.output.resolve()
    if OUTPUT.exists() and any(OUTPUT.iterdir()):
        p.error('Output must be empty')
    (OUTPUT / 'candidates_by_gene').mkdir(parents=True)
    GENOMES = {CONFIG['genome']: dict(fasta=str(args.fasta.resolve()), gtf=str(args.gtf.resolve()))}
    subset = OUTPUT / 'target_annotation.gtf'
    with args.gtf.open() as src, subset.open('w') as dst:
        for line in src:
            if line.startswith('#'):
                continue
            match = re.search(r'(?:^|;)\s*gene_id\s+"([^".]+)(?:\.[^"]+)?"', line.split('\t')[-1])
            if match and match.group(1) in selected:
                dst.write(line)
    fields = list(dict.fromkeys([*GTF_RUN_FIELDS, 'gene_id']))
    GTF = read_gtf(str(subset), extract_fields=fields)
    logging.getLogger('uprobe.core.api').setLevel(logging.WARNING)
    logging.getLogger('uprobe.core.gen.fun').setLevel(logging.WARNING)
    started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('fork')) as pool:
        details = list(pool.map(run_gene, selected, chunksize=1))
    with gzip.open(OUTPUT / 'candidates.csv.gz', 'wt', compresslevel=1, newline='') as dst:
        for index, target_id in enumerate(selected):
            with gzip.open(OUTPUT / 'candidates_by_gene' / (target_id + '.csv.gz'), 'rt') as src:
                header = next(src)
                if index == 0:
                    dst.write(header)
                for line in src:
                    dst.write(line)
    manifest = dict(method_version=METHOD, native_api=['UProbeAPI.generate_target_seqs', 'UProbeAPI.construct_probes'],
                    config_sha256=sha256(args.config), targets_sha256=sha256(args.targets),
                    fasta_sha256=sha256(args.fasta), gtf_sha256=sha256(args.gtf), runner_sha256=sha256(__file__),
                    candidate_records=sum(r['candidate_records'] for r in details),
                    native_probe_sets=sum(r['native_probe_sets'] for r in details), genes=details,
                    final_selection='not_evaluated', seconds=round(time.time() - started, 3),
                    deduplication='Complete target windows per gene; first sorted transcript origin retained',
                    constraints='ACGT windows only; all transcript windows; gaps 0 and 1; both RNA adjacency orientations')
    (OUTPUT / 'run_details.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
