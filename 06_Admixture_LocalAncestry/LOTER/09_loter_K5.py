#!/usr/bin/env python
# Parte 4 (K5): ancestria local LOTER. Fuentes puras K5 = PA,PB1,PB2,PB3,PB4.
# Query = sub-poblaciones admix (NoAm,SoAm1-4,Hol) del esquema noGen.
# Genomas ~99% homocigotos -> haploides (1 hap por cepa).
import allel, numpy as np, pandas as pd, sys, os
import loter.locanc.local_ancestry as lc
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
vcf=f"{proj}/results/eub_final.vcf.gz"
sets=dict(l.split() for l in open(f"{proj}/results/K5/admix_subpops_noGen/SETS.txt"))
sources=["PA","PB1","PB2","PB3","PB4"]
tgtpops=["NoAm","SoAm1","SoAm2","SoAm3","SoAm4","Hol"]
outdir=f"{proj}/results/K5/loter"; os.makedirs(outdir, exist_ok=True)

print("cargando VCF...", file=sys.stderr)
cal=allel.read_vcf(vcf, fields=['samples','calldata/GT','variants/CHROM','variants/POS'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT'])
chrom=cal['variants/CHROM']; pos=cal['variants/POS']
hap0 = gt.to_haplotypes()[:, ::2]

def pop_mat(idx, mask):
    return hap0[np.ix_(mask, idx)].T.astype(np.uint8)

src_idx={p:[i for i,s in enumerate(samples) if sets.get(s)==p] for p in sources}
tgt_idx=[i for i,s in enumerate(samples) if sets.get(s) in tgtpops]
print("fuentes:",{p:len(src_idx[p]) for p in sources},"| targets:",len(tgt_idx),file=sys.stderr)

chroms=sorted(set(chrom))
anc_per_snp={}; res_frac=np.zeros((len(tgt_idx), len(sources))); counts=np.zeros(len(tgt_idx))
posall=[]
for c in chroms:
    m=(chrom==c)
    l_H=[pop_mat(src_idx[p], m) for p in sources]
    h_adm=pop_mat(tgt_idx, m)
    res=np.asarray(lc.loter_local_ancestry(l_H, h_adm, num_threads=4)[0])  # (n_adm,n_snp) 0..4
    for k in range(res.shape[0]):
        for j in range(len(sources)):
            res_frac[k,j]+=np.sum(res[k]==j)
        counts[k]+=res.shape[1]
    anc_per_snp[c]=res.astype(np.uint8); posall.append(pos[m])
    print(f"  {c}: {res.shape[1]} SNP", file=sys.stderr)

frac=res_frac/counts[:,None]
df=pd.DataFrame(frac, columns=sources)
df.insert(0,"ID",[samples[i] for i in tgt_idx])
df.insert(1,"POP",[sets[samples[i]] for i in tgt_idx])
df["top_anc"]=df[sources].idxmax(1); df["top_frac"]=df[sources].max(1)
df["minor_frac"]=1-df["top_frac"]
df.to_csv(f"{outdir}/loter_ancestry_proportions.tsv", sep="\t", index=False)
np.savez_compressed(f"{outdir}/loter_perSNP.npz",
                    targets=np.array([samples[i] for i in tgt_idx]),
                    sources=np.array(sources),
                    **{c:anc_per_snp[c] for c in chroms},
                    **{f"POS_{c}":posall[i] for i,c in enumerate(chroms)})
print("=== proporciones medias por POP ===", file=sys.stderr)
print(df.groupby("POP")[sources+["minor_frac"]].mean().round(3), file=sys.stderr)
print("OK", file=sys.stderr)
