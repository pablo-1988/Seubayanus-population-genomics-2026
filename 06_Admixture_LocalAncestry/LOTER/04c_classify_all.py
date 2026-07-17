#!/usr/bin/env python
# Clasificacion admix/puro para TODAS las cepas via LOTER.
# - Puras: 5-fold CV (cada pura clasificada cuando esta fuera de la referencia).
# - Admix/NoAm: contra todas las referencias puras.
# Umbral: minor_frac >= 0.10 => admixada (validado: puras p95=0.07, max 0.11).
import allel, numpy as np, pandas as pd, sys, os
import loter.locanc.local_ancestry as lc
np.random.seed(7)
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
meta=pd.read_csv(f"{proj}/results/metadata_final.txt", sep="\t"); meta.columns=["ID","POP"]
sources=["PA","PB1","PB2","PB3","PB4","PB5","PB6"]; THR=0.10
outdir=f"{proj}/results/04_loter"

cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz",
                   fields=['samples','calldata/GT','variants/CHROM'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT']); chrom=cal['variants/CHROM']
id2pop=dict(zip(meta.ID,meta.POP)); hap0=gt.to_haplotypes()[:,::2]
chroms=sorted(set(chrom))

def run_loter(ref_idx, query_idx):
    frac=np.zeros((len(query_idx),len(sources))); cnt=np.zeros(len(query_idx))
    for c in chroms:
        m=(chrom==c)
        l_H=[hap0[np.ix_(m,ref_idx[p])].T.astype(np.uint8) for p in sources]
        h_adm=hap0[np.ix_(m,query_idx)].T.astype(np.uint8)
        res=np.asarray(lc.loter_local_ancestry(l_H,h_adm,num_threads=4)[0])
        for k in range(res.shape[0]):
            for j in range(len(sources)): frac[k,j]+=np.sum(res[k]==j)
            cnt[k]+=res.shape[1]
    return frac/cnt[:,None]

rows=[]
# 1) targets (Admix+NoAm) contra todas las puras
allref={p:[i for i,s in enumerate(samples) if id2pop.get(s)==p] for p in sources}
tgt=[i for i,s in enumerate(samples) if id2pop.get(s) in ("Admix","NoAm")]
print("target run...", file=sys.stderr); Ft=run_loter(allref,tgt)
for k,i in enumerate(tgt): rows.append([samples[i],id2pop[samples[i]]]+list(Ft[k]))

# 2) puras via 5-fold
pures={p:[i for i,s in enumerate(samples) if id2pop.get(s)==p] for p in sources}
for p in pures: np.random.shuffle(pures[p])
K=5
for f in range(K):
    ref={}; q=[]
    for p in sources:
        idx=pures[p]; fold=[idx[j] for j in range(len(idx)) if j%K==f]
        ref[p]=[x for x in idx if x not in set(fold)]; q+=fold
    if not q: continue
    print(f"fold {f} ({len(q)} puras)...", file=sys.stderr)
    Fp=run_loter(ref,q)
    for k,i in enumerate(q): rows.append([samples[i],id2pop[samples[i]]]+list(Fp[k]))

df=pd.DataFrame(rows,columns=["ID","label"]+sources)
o=df[sources].values; order=np.argsort(-o,axis=1)
df["top_anc"]=[sources[order[i,0]] for i in range(len(df))]
df["second_anc"]=[sources[order[i,1]] for i in range(len(df))]
df["top_frac"]=o.max(1); df["minor_frac"]=1-df["top_frac"]
df["LOTER_admix"]=df["minor_frac"]>=THR
df=df.sort_values("minor_frac",ascending=False)
df.to_csv(f"{outdir}/loter_classification_all.tsv", sep="\t", index=False)
print("\n=== clasificacion LOTER vs etiqueta ===", file=sys.stderr)
df["label_admix"]=df["label"].isin(["Admix","NoAm"])
print(pd.crosstab(df["label"], df["LOTER_admix"]), file=sys.stderr)
print("\nminor_frac por etiqueta:", file=sys.stderr)
print(df.groupby("label")["minor_frac"].describe(percentiles=[.5,.95]).round(3), file=sys.stderr)
