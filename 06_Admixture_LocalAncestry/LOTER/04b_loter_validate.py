#!/usr/bin/env python
# Validacion LOTER: held-out de cepas puras para estimar el piso de ruido de
# minor_frac. Referencia = 70% de cada linaje; query = 30% puro + Admix + NoAm.
import allel, numpy as np, pandas as pd, sys, os
import loter.locanc.local_ancestry as lc
np.random.seed(7)
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
meta=pd.read_csv(f"{proj}/results/metadata_final.txt", sep="\t"); meta.columns=["ID","POP"]
sources=["PA","PB1","PB2","PB3","PB4","PB5","PB6"]
outdir=f"{proj}/results/04_loter"; os.makedirs(outdir, exist_ok=True)

cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz",
                   fields=['samples','calldata/GT','variants/CHROM'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT']); chrom=cal['variants/CHROM']
id2pop=dict(zip(meta.ID,meta.POP)); hap0=gt.to_haplotypes()[:,::2]

# split puras: 70% ref / 30% test por linaje
ref_idx={}; test_pure=[]
for p in sources:
    idx=[i for i,s in enumerate(samples) if id2pop.get(s)==p]
    np.random.shuffle(idx); k=max(2,int(round(len(idx)*0.3)))
    test_pure+= idx[:k]; ref_idx[p]=idx[k:]
adm_idx=[i for i,s in enumerate(samples) if id2pop.get(s) in ("Admix","NoAm")]
query_idx=test_pure+adm_idx
qlab=[("PURE:"+id2pop[samples[i]]) if i in set(test_pure) else id2pop[samples[i]] for i in query_idx]

chroms=sorted(set(chrom))
frac=np.zeros((len(query_idx),len(sources))); cnt=np.zeros(len(query_idx))
for c in chroms:
    m=(chrom==c)
    l_H=[hap0[np.ix_(m,ref_idx[p])].T.astype(np.uint8) for p in sources]
    h_adm=hap0[np.ix_(m,query_idx)].T.astype(np.uint8)
    res=np.asarray(lc.loter_local_ancestry(l_H,h_adm,num_threads=4)[0])
    for k in range(res.shape[0]):
        for j in range(len(sources)): frac[k,j]+=np.sum(res[k]==j)
        cnt[k]+=res.shape[1]
    print(f"  {c} ok", file=sys.stderr)
F=frac/cnt[:,None]
df=pd.DataFrame(F,columns=sources); df.insert(0,"ID",[samples[i] for i in query_idx])
df.insert(1,"class",qlab)
df["minor_frac"]=1-df[sources].max(1)
df["group"]=np.where(df["class"].str.startswith("PURE"),"pure(held-out)",
             np.where(df["class"]=="NoAm","NoAm","Admix"))
df.to_csv(f"{outdir}/loter_validation.tsv", sep="\t", index=False)
print("\n=== minor_frac por grupo ===", file=sys.stderr)
print(df.groupby("group").minor_frac.describe()[["count","mean","50%","75%","95%","max".replace("95%","max")]].round(3) if False else
      df.groupby("group")["minor_frac"].describe(percentiles=[.5,.9,.95]).round(3), file=sys.stderr)
