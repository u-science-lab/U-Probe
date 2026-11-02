import openpyxl,csv,re,pathlib,argparse
ap=argparse.ArgumentParser();ap.add_argument('workbook');ap.add_argument('output_dir');args=ap.parse_args()
p=pathlib.Path(args.output_dir);p.mkdir(parents=True,exist_ok=True)
w=openpyxl.load_workbook(args.workbook,data_only=True)
data=[(int(x[0]),re.sub(r'^PRISM_\d+_','',x[1]).strip().replace('TBX5','Tbx5').replace('MEF2C','Mef2c').replace('Hba-X','Hba-x'),x[2].upper()) for x in list(w['Padlock probe-mouse embryo(64)'].values)[1:] if x[1] and x[2]]
b={int(x[6]):x[8] for x in w['Imaging probe (64 panel)'].values if isinstance(x[6],(int,float)) and x[8]}
def write(path,cols,rows):
 with open(path,'w') as f:
  c=csv.writer(f);c.writerow(cols);c.writerows(rows)
assert len(data)==64 and len({g for n,g,s in data})==64
assert all(len(s)==120 and re.fullmatch('[ACGT]+',s) and s[20:100]==b[n] for n,g,s in data)
write(p/'targets.csv',['gene'],[(g,) for n,g,s in data])
write(p/'codebook.csv',['gene','BC1','BC2','BC3','BC4'],[(g,*[b[n][i:i+20] for i in range(0,80,20)]) for n,g,s in data])
write(p/'reference_oligos.csv',['gene','oligo_type','oligo_id','dna_sequence'],[(g,'padlock',f'PRISM_{n}_{g}',s) for n,g,s in data])
print('INPUT',len(data),'unique',len({s for n,g,s in data}))

write(p/"gene_names.csv",["published_name","gene"],[(x[1],g) for x,(n,g,seq) in zip([x for x in list(w["Padlock probe-mouse embryo(64)"].values)[1:] if x[1] and x[2]],data)])
