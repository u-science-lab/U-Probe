#!/usr/bin/env python3
"""Extract the 42-gene targeted ExSeq visual-cortex probe set from Alon et al. 2021.

Source: NIHMS1670746-supplement-Supp_Tables.xlsx (S9 + S10).
Only Python standard library is required.
"""
import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
ASSAY = "Mouse Visual Cortex"
UNIVERSAL = "TCTCGGGAACGCTGAAGACGGC"
NAME_RE = re.compile(
    r"^Mouse_(?P<gene>.+?)_(?P<accession>NM_\d+)_BC_(?P<logical>[0-3]{4})"
    r"_probeBC_(?P<probe_bc>[0-3]{2})$"
)
DNA = re.compile(r"[ACGT]{74}")


def worksheet_rows(xlsx, requested):
    """Stream one worksheet with Excel shared-string support (no openpyxl)."""
    with ZipFile(xlsx) as archive:
        ns = {"m": MAIN[1:-1], "r": REL[1:-1]}
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in si.iter(MAIN + "t"))
                       for si in shared.findall(MAIN + "si")]
        wb = ET.fromstring(archive.read("xl/workbook.xml"))
        rel_xml = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {r.attrib["Id"]: r.attrib["Target"] for r in rel_xml}
        work_path = None
        for sheet in wb.findall("m:sheets/m:sheet", ns):
            if sheet.get("name", "").strip() == requested:
                rel_id = sheet.get(REL + "id")
                target = targets[rel_id]
                work_path = (
                    target.lstrip("/") if target.startswith("/")
                    else ("xl/" + target if not target.startswith("xl/") else target)
                )
                break
        if not work_path:
            raise ValueError(f"Worksheet {requested!r} not found")
        root = ET.fromstring(archive.read(work_path))
        for row in root.findall(MAIN + "sheetData/" + MAIN + "row"):
            values = {"source_row": int(row.get("r"))}
            for cell in row.findall(MAIN + "c"):
                column = re.match(r"[A-Z]+", cell.get("r")).group()
                v = cell.find(MAIN + "v")
                inline = cell.find(MAIN + "is")
                text = (
                    v.text if v is not None
                    else "".join(t.text or "" for t in inline.iter(MAIN + "t"))
                    if inline is not None else ""
                )
                if cell.get("t") == "s" and text:
                    text = strings[int(text)]
                values[column] = text or ""
            yield values


def save_csv(path, rows, cols):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def parse_oligo(row, targets):
    name = (row.get("B") or "").strip()
    match = NAME_RE.fullmatch(name)
    if not match:
        raise ValueError(f"Unrecognized probe name (S10 row {row['source_row']}): {name!r}")
    gene = match["gene"]
    if gene not in targets:
        raise ValueError(f"Probe gene missing in S9: {name}")
    metadata = targets[gene]
    assert match["logical"] == metadata["logical_barcode"], name
    assert metadata["refseq_accession"].split(".")[0] == match["accession"], name

    original = row.get("C", "")
    if not original.startswith("/5Phos/"):
        raise ValueError(f"Missing 5'-phosphate at S10 row {row['source_row']}")
    full = re.sub(r"\s+", "", original[len("/5Phos/"):]).upper()
    if not DNA.fullmatch(full):
        raise ValueError(f"Expected 74 DNA bases at S10 row {row['source_row']}")
    arm_5p, spacer_5p, solid, universal, illumina, spacer_3p, arm_3p = (
        full[:16], full[16:19], full[19:26], full[26:48],
        full[48:55], full[55:58], full[58:74]
    )
    if (spacer_5p, universal, spacer_3p) != ("AAA", UNIVERSAL, "AAA"):
        raise ValueError(f"Unexpected backbone at S10 row {row['source_row']}")
    assert len(arm_5p) == len(arm_3p) == 16
    return {
        "method": "Targeted ExSeq",
        "dataset": "Alon_2021_Mouse_Visual_Cortex_42gene",
        "species": "Mus musculus",
        "gene": gene,
        "refseq_accession": metadata["refseq_accession"],
        "logical_barcode": metadata["logical_barcode"],
        "probe_barcode_id": match["probe_bc"],
        "probe_id": name,
        "oligo_type": "padlock",
        "five_prime_modification": "/5Phos/",
        "full_sequence": full,
        "arm_5p": arm_5p,
        "spacer_5p": spacer_5p,
        "barcode_solid": solid,
        "universal_site": universal,
        "barcode_illumina": illumina,
        "spacer_3p": spacer_3p,
        "arm_3p": arm_3p,
        "source_sheet": "Table S10",
        "source_row": row["source_row"],
        "source_original_sequence": original,
    }



def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, help='Original S9/S10 workbook; re-extract and verify committed inputs')
    parser.add_argument('--reference-dir', required=True, type=Path, help='Directory containing mm10.fa and mm10.gtf')
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    if args.workbook:
        targets = {}
        for row in worksheet_rows(args.workbook, 'Table S9'):
            if row.get('A', '').strip() == ASSAY:
                targets[row['B'].strip().upper()] = dict(refseq_accession=row['C'].strip(), logical_barcode=row['D'].strip(), source_row=row['source_row'])
        oligos = [parse_oligo(row, targets) for row in worksheet_rows(args.workbook, 'Table S10') if row.get('A', '').strip() == ASSAY]
        assert len(targets) == 42 and len(oligos) == 334
        assert len({o['full_sequence'] for o in oligos}) == 334
        target_rows = [dict(target_id=g, annotation_gene_name=g.capitalize(), **v) for g,v in targets.items()]
        ref_rows = [dict(gene=o['gene'],pair_id='',oligo_type='padlock',oligo_id=o['probe_id'],barcode=o['logical_barcode'],full_sequence='/5Phos/'+o['full_sequence'],dna_sequence=o['full_sequence'],binding_sequence=o['arm_3p']+o['arm_5p'],binding_length_nt=32,modification_5p='5phosphate',refseq_accession=o['refseq_accession'],probe_barcode_id=o['probe_barcode_id'],arm_5p=o['arm_5p'],arm_3p=o['arm_3p'],source_sheet='Table S10',source_row=o['source_row']) for o in oligos]
        codes = [dict(target_id=o['gene'],logical_barcode=o['logical_barcode'],probe_barcode_id=o['probe_barcode_id'],solid_barcode_7nt=o['barcode_solid'],illumina_barcode_7nt=o['barcode_illumina']) for o in oligos]
        assert len({(r['target_id'],r['solid_barcode_7nt'],r['illumina_barcode_7nt']) for r in codes}) == len(codes)
        for filename, rows in [('targets.csv',target_rows),('reference_oligos.csv',ref_rows),('codebook.csv',codes)]:
            path=root/filename
            if path.exists():
                current=list(csv.DictReader(path.open()))
                assert current == [{k:str(v) for k,v in row.items()} for row in rows], filename+' differs from workbook'
            else: save_csv(path,rows,list(rows[0]))
    targets=list(csv.DictReader((root/'targets.csv').open()))
    oligos=list(csv.DictReader((root/'reference_oligos.csv').open()))
    wanted={r['target_id'] for r in targets}
    exons=defaultdict(list); txgene={}; geneids=defaultdict(set); names={}
    with (args.reference_dir/'mm10.gtf').open() as handle:
        for line in handle:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if f[2]!='exon': continue
            a=dict(re.findall(r'(\w+) "([^"]*)"',f[8])); gene=a.get('gene_name','').upper()
            if gene not in wanted: continue
            tx=a['transcript_id']; txgene[tx]=gene; geneids[gene].add(a['gene_id']); names[gene]=a['gene_name']
            exons[tx].append((f[0],int(f[3])-1,int(f[4]),f[6]))
    # Read indexed genomic FASTA directly: transcripts are independently reconstructed from GTF exons.
    fasta=args.reference_dir/'mm10.fa'; index={}
    fai=Path(str(fasta)+'.fai')
    if not fai.is_file(): raise ValueError('FASTA index missing; run samtools faidx on the local FASTA')
    for line in fai.open():
        f=line.split(); index[f[0]]=tuple(map(int,f[1:5]))
    def fetch(handle, chrom, start, end):
        if chrom not in index: raise ValueError('GTF chromosome absent from FASTA index: '+chrom)
        length,offset,bases,width=index[chrom]
        if not (0 <= start < end <= length): raise ValueError('GTF exon outside FASTA bounds: '+chrom)
        if bases <= 0 or width < bases: raise ValueError('Invalid FASTA index geometry')
        left=offset+(start//bases)*width+start%bases
        right=offset+((end-1)//bases)*width+(end-1)%bases+1
        handle.seek(left)
        return handle.read(right-left).decode().replace('\n','').replace('\r','').upper()
    sequences=defaultdict(list)
    with fasta.open('rb') as handle:
        for tx,blocks in exons.items():
            blocks.sort(key=lambda b:b[1]); seq=''.join(fetch(handle,c,s,e) for c,s,e,_ in blocks)
            if blocks[0][3]=='-': seq=reverse_complement(seq)
            sequences[txgene[tx]].append((tx,seq))
    assert all(names.get(r['target_id']) == r['annotation_gene_name'] for r in targets), 'Annotation gene-name crosswalk changed'
    coverage=[]; mappings=[]
    for row in targets:
        gene=row['target_id']; coverage.append(dict(target_id=gene,refseq_accession=row['refseq_accession'],annotation_gene_name=names.get(gene,''),annotation_gene_ids=';'.join(sorted(geneids[gene])),transcript_count=len(sequences[gene]),accession_status='RefSeq ID absent from GENCODE GTF; gene-name crosswalk'))
    for row in oligos:
        gene=row['gene']; query=reverse_complement(row['arm_3p']+row['arm_5p'])
        hits=[tx for tx,seq in sequences[gene] if query in seq]
        reason='' if hits else 'Neither 16-nt arm matches any same-gene annotated transcript'
        if not hits and any(reverse_complement(row['arm_3p']) in seq or reverse_complement(row['arm_5p']) in seq for _,seq in sequences[gene]): reason='Individual-arm or isoform coverage without contiguous 32-nt match'
        mappings.append(dict(target_id=gene,oligo_id=row['oligo_id'],mapped=bool(hits),matching_transcripts=';'.join(hits),reason=reason))
    args.outdir.mkdir(parents=True,exist_ok=True)
    save_csv(args.outdir/'annotation_coverage.csv',coverage,list(coverage[0]))
    save_csv(args.outdir/'arm_mapping.csv',mappings,list(mappings[0]))
    print(json.dumps(dict(targets=len(targets),padlocks=len(oligos),annotated_genes=len(geneids),transcripts=len(exons),mapped_padlocks=sum(r['mapped'] for r in mappings),workbook_sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest() if args.workbook else None)))


if __name__ == '__main__':
    import json
    main()
