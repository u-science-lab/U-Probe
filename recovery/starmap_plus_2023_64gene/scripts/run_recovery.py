#!/usr/bin/env python3
"""Portable independent STARmap candidate recovery; no reference oligos are read."""
import argparse, csv, gzip, hashlib, json, multiprocessing as mp, re, time, copy, shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd
import yaml
from uprobe.core.api import UProbeAPI
from uprobe.core.gen.fun import read_gtf, GTF_RUN_FIELDS
from uprobe.core.layout import compile_layout, expand_layout_candidates

OUTPUT = None
FIELDS = ['stage','gene','pair_id','probe_id','oligo_type','full_sequence','modification_5p']
GTF = None
CONFIG = None
GENOMES = None
RUN_FINGERPRINT = None

def file_sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def run_gene(gene):
    started=time.time()
    cfg=copy.deepcopy(CONFIG); cfg['targets']=[gene]
    api=UProbeAPI(cfg,GENOMES,OUTPUT/'engine',require_genome=True)
    api._gtf=GTF[GTF.gene_name.eq(gene)]
    ids=sorted(api._gtf.loc[api._gtf.type.eq('exon'),'transcript_id'].dropna().unique())
    extract_cfg=copy.deepcopy(cfg)
    layout=extract_cfg['extracts']['target_region'].pop('layout')
    extract_cfg['probes']={}
    extraction_api=UProbeAPI(extract_cfg,GENOMES,OUTPUT/'engine',require_genome=True)
    extraction_api._gtf=api._gtf
    compiled=compile_layout(layout,cfg['extracts']['target_region']['length'])
    seen=set();transcript_lengths={};rejected=0
    output=OUTPUT/'candidates_by_gene'/f'{gene}.csv.gz'; output.parent.mkdir(parents=True,exist_ok=True)
    manifest=OUTPUT/'candidate_sites'/f'{gene}.csv.gz'; manifest.parent.mkdir(parents=True,exist_ok=True)
    count=0
    with gzip.open(output,'wt',newline='') as f, gzip.open(manifest,'wt',newline='') as sf:
        writer=csv.DictWriter(f,fieldnames=FIELDS); writer.writeheader()
        sw=csv.writer(sf); sw.writerow(['pair_id','transcript_names','target_start','target_end','layout_id','target_region'])
        for tx in ids:
            extract_cfg['extracts']['target_region']['transcripts']=[tx]
            # Same native API extraction and layout engine, expanded only in bounded batches.
            df=extraction_api.generate_target_seqs()
            if df.empty:continue
            transcript_lengths[tx]=int(df.target_end.max()) if 'target_end' in df else None
            for offset in range(0,len(df),250):
                raw=[]
                for row in df.iloc[offset:offset+250].to_dict('records'):
                    sequence=row['target_region'].upper()
                    # Every layout of this window is expanded below; window dedup is equivalent
                    # to v1's gene/window/layout dedup but avoids storing each layout key.
                    if sequence in seen:continue
                    seen.add(sequence)
                    if set(sequence)-set('ACGT'):rejected+=len(compiled[1][len(sequence)]);continue
                    row['target_region']=sequence;raw.append(row)
                if not raw:continue
                chunk=expand_layout_candidates(pd.DataFrame(raw),compiled)
                built=api.construct_probes(chunk)
                if len(built)!=len(chunk):raise RuntimeError('construction row mismatch')
                for k,(_,row) in enumerate(chunk.iterrows()):
                    pair=f'{gene}_{count+1:09d}'
                    for typ in ('padlock','primer'):
                        seq=built.iloc[k][typ]
                        if not seq or set(seq)-set('ACGT'):raise RuntimeError('invalid generated DNA')
                        writer.writerow(dict(stage='candidate',gene=gene,pair_id=pair,probe_id=pair+'_'+typ,oligo_type=typ,full_sequence=seq,modification_5p='phosphate' if typ=='padlock' else 'none'))
                    sw.writerow([pair,row.get('transcript_names'),row.get('target_start'),row.get('target_end'),row['layout_id'],row['target_region']])
                    count+=1
            del df
    details=dict(gene=gene,transcripts=ids,transcript_max_window_end=transcript_lengths,candidate_pairs=count,non_acgt_layouts_rejected=rejected,seconds=time.time()-started,run_fingerprint=RUN_FINGERPRINT)
    (OUTPUT/'gene_run_details'/f'{gene}.json').write_text(json.dumps(details,indent=2))
    print(json.dumps(details),flush=True)
    return details

def main():
    global GTF,CONFIG,GENOMES,OUTPUT,RUN_FINGERPRINT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--fasta',type=Path,required=True)
    parser.add_argument('--gtf',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.workers<1:parser.error('--workers must be positive')
    CONFIG=yaml.safe_load(args.config.read_text(encoding='utf-8'))
    targets=CONFIG.get('targets',[])
    if not isinstance(targets,list) or not targets or any(not isinstance(g,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',g) for g in targets) or len(targets)!=len(set(targets)):
        parser.error('targets must be unique gene names, not supplied target sequences')
    genome=CONFIG.setdefault('genome','reference')
    GENOMES={genome:{'fasta':str(args.fasta.resolve()),'gtf':str(args.gtf.resolve())}}
    # A resume is valid only for identical configuration, reference bytes and runner.
    provenance={'config_sha256':file_sha256(args.config),'fasta_sha256':file_sha256(args.fasta),'gtf_sha256':file_sha256(args.gtf),'runner_sha256':file_sha256(Path(__file__))}
    RUN_FINGERPRINT=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    OUTPUT=args.output.resolve();results=OUTPUT
    state=results/'run_state.json'
    if results.exists() and any(results.iterdir()):
        if not args.resume:parser.error('output is not empty; choose another directory or use --resume')
        if not state.exists() or json.loads(state.read_text())['run_fingerprint']!=RUN_FINGERPRINT:
            parser.error('resume rejected: configuration, reference or runner does not match output provenance')
    results.mkdir(parents=True,exist_ok=True)
    state.write_text(json.dumps(dict(run_fingerprint=RUN_FINGERPRINT,**provenance),indent=2))
    (results/'gene_run_details').mkdir(exist_ok=True)
    subset=results/'target_annotation.gtf'
    genes=set(CONFIG['targets']); source=args.gtf
    with source.open() as src,subset.open('w') as dst:
        for line in src:
            if line.startswith('#'): dst.write(line); continue
            attrs=line.split('\t')[-1]
            match=re.search(r'(?:^|;)\s*(?:gene_name|gene)\s+"([^"]+)"',attrs)
            if match and match.group(1) in genes: dst.write(line)
    GTF=read_gtf(str(subset),extract_fields=GTF_RUN_FIELDS)
    missing=genes-set(GTF.gene_name)
    if missing: raise RuntimeError('missing genes '+repr(missing))
    started=time.time()
    completed={}
    if args.resume:
        for gene in CONFIG['targets']:
            marker=results/'gene_run_details'/f'{gene}.json'
            if marker.exists() and (results/'candidates_by_gene'/f'{gene}.csv.gz').exists() and (results/'candidate_sites'/f'{gene}.csv.gz').exists():
                detail=json.loads(marker.read_text())
                if detail.get('run_fingerprint')!=RUN_FINGERPRINT:raise RuntimeError('gene completion provenance mismatch: '+gene)
                completed[gene]=detail
    pending=[g for g in CONFIG['targets'] if g not in completed]
    print(json.dumps({'event':'start','workers':args.workers,'resume':args.resume,'completed':list(completed),'pending_genes':len(pending)}),flush=True)
    # A lost worker raises BrokenProcessPool instead of silently losing a Pool.map task.
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('fork')) as pool:
        fresh=list(pool.map(run_gene,pending,chunksize=1))
    completed.update({d['gene']:d for d in fresh})
    details=[completed[g] for g in CONFIG['targets']]
    combined=results/'candidates.csv.gz';temporary=results/'candidates.tmp.csv.gz'
    with gzip.open(temporary,'wb',compresslevel=1) as out:
        header=None
        for gene in CONFIG['targets']:
            with gzip.open(results/'candidates_by_gene'/f'{gene}.csv.gz','rb') as f:
                current=f.readline()
                if header is None:header=current;out.write(header)
                elif current!=header:raise RuntimeError('candidate CSV headers differ')
                shutil.copyfileobj(f,out,4*1024*1024)
    temporary.replace(combined)
    (results/'run_details.json').write_text(json.dumps(dict(start_unix=started,seconds=time.time()-started,workers=args.workers,candidate_pairs=sum(d['candidate_pairs'] for d in details),genes=details,final_stage='not_applicable_no_published_selection_rule',run_fingerprint=RUN_FINGERPRINT),indent=2))

if __name__=='__main__': main()
