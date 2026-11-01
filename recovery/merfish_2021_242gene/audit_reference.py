#!/usr/bin/env python3
"""Audit published transcript identifiers and binding regions against mm10 annotation."""
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def rc(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    targets = list(csv.DictReader((root/'targets.csv').open()))
    oligos = list(csv.DictReader((root/'reference_oligos.csv').open()))
    wanted = {r['target_id'].upper() for r in targets}
    exons = defaultdict(list)
    txgene = {}
    names = {}
    geneids = defaultdict(set)
    with (args.reference_dir/'mm10.gtf').open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] != 'exon':
                continue
            attr = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            gene = attr.get('gene_name', '').upper()
            if gene not in wanted:
                continue
            tx = attr['transcript_id']
            txgene[tx] = gene
            geneids[gene].add(attr['gene_id'])
            names[gene] = attr['gene_name']
            exons[tx].append((fields[0], int(fields[3])-1, int(fields[4]), fields[6]))
    fasta = args.reference_dir/'mm10.fa'
    index = {}
    with Path(str(fasta)+'.fai').open() as handle:
        for line in handle:
            fields = line.split()
            index[fields[0]] = tuple(map(int, fields[1:5]))
    def fetch(handle, chrom, start, end):
        length, offset, bases, width = index[chrom]
        if not (0 <= start < end <= length):
            raise ValueError('Invalid annotated exon coordinates')
        left = offset+(start//bases)*width+start%bases
        right = offset+((end-1)//bases)*width+(end-1)%bases+1
        handle.seek(left)
        return handle.read(right-left).decode().replace('\n', '').replace('\r', '').upper()
    sequences = defaultdict(list)
    with fasta.open('rb') as handle:
        for tx, blocks in exons.items():
            blocks.sort(key=lambda b: b[1])
            seq = ''.join(fetch(handle, chrom, start, end) for chrom, start, end, _ in blocks)
            if blocks[0][3] == '-':
                seq = rc(seq)
            sequences[txgene[tx]].append((tx, seq))
    coverage = []
    for row in targets:
        gene = row['target_id'].upper()
        original = row.get('transcript_id') or row.get('original_transcript_id') or row.get('ensembl_transcript_id') or ''
        same_gene_ids = {tx for tx, _ in sequences[gene]}
        base_hits = sorted(tx for tx in same_gene_ids if tx.split('.')[0] == original.split('.')[0])
        coverage.append(dict(target_id=gene, published_transcript_id=original,
            annotation_gene_name=names.get(gene, ''), annotation_gene_ids=';'.join(sorted(geneids[gene])),
            transcript_count=len(same_gene_ids), exact_transcript_id_present=original in same_gene_ids,
            transcript_id_base_matches=';'.join(base_hits), accession_status='Exact version' if original in same_gene_ids else ('Base ID present; source unversioned' if '.' not in original else 'Version differs') if base_hits else 'Absent'))
    mappings = []
    for row in oligos:
        gene = row['gene'].upper()
        query = rc(row['binding_sequence'].upper())
        hits = [tx for tx, seq in sequences[gene] if query in seq]
        mappings.append(dict(target_id=gene, oligo_id=row['oligo_id'], mapped=bool(hits),
            matching_transcripts=';'.join(hits), reason='' if hits else 'Binding region absent from same-gene annotated mRNA'))
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir/'annotation_coverage.csv', coverage)
    write_csv(args.outdir/'binding_mapping.csv', mappings)
    print(json.dumps(dict(targets=len(targets), annotated_genes=len(names), transcripts=len(exons),
        exact_transcript_versions=sum(r['exact_transcript_id_present'] for r in coverage),
        matching_transcript_bases=sum(bool(r['transcript_id_base_matches']) for r in coverage),
        reference_oligos=len(oligos), mapped_oligos=sum(r['mapped'] for r in mappings))))


if __name__ == '__main__':
    main()
