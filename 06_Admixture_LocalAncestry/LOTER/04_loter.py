#!/usr/bin/env python
# Parte 4: ancestria local con LOTER.
# Fuentes (referencias puras): PA, PB1..PB6. Query: Admix + NoAm.
# Genomas ~99% homocigotos -> se tratan como HAPLOIDES (1 haplotipo por cepa).
import allel, numpy as np, pandas as pd, sys, os
import loter.locanc.local_ancestry as lc
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
vcf=f"{proj}/results/eub_final.vcf.gz"
meta=pd.read_csv(f"{proj}/results/metadata_final.txt", sep="\t"); meta.columns=["ID","POP"]
sources=["PA","PB1","PB2","PB3","PB4","PB5","PB6"]
outdir=f"{proj}/results/04_loter"; os.makedirs(outdir, exist_ok=True)

print("cargando VCF...", file=sys.stderr)
cal=allel.read_vcf(vcf, fields=['samples','calldata/GT','variants/CHROM','variants/POS'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT'])
chrom=cal['variants/CHROM']; pos=cal['variants/POS']
id2pop=dict(zip(meta.ID, meta.POP))
# haploidizar: 1 alelo por cepa (hap 0). Missing (-1) -> se maneja por cromosoma.
hap0 = gt.to_haplotypes()[:, ::2]   # (n_snp, n_sample) alelo del primer haplotipo

def pop_mat(idx, mask):
    # matriz (n_hap, n_snp) para SNPs de un cromosoma
    return hap0[np.ix_(mask, idx)].T.astype(np.uint8)

targets=[s for s in samples if id2pop.get(s) in ("Admix","NoAm")]
src_idx={p:[i for i,s in enumerate(samples) if id2pop.get(s)==p] for p in sources}
tgt_idx=[i for i,s in enumerate(samples) if s in targets]

chroms=sorted(set(chrom))
anc_per_snp={}  # target -> array de ancestria por SNP (global, concatenado)
allpos=[]
res_frac=np.zeros((len(tgt_idx), len(sources)))
counts=np.zeros(len(tgt_idx))
for c in chroms:
    m = (chrom==c)
    l_H=[pop_mat(src_idx[p], m) for p in sources]
    h_adm=pop_mat(tgt_idx, m)
    # LOTER multi-poblacion: loter_local_ancestry devuelve (argmax, counts);
    # argmax = indice de fuente (0..6) por hap x SNP. (loter_smooth es solo 2-vias)
    res = lc.loter_local_ancestry(l_H, h_adm, num_threads=4)
    res = np.asarray(res[0])  # (n_adm, n_snp_c) con valores 0..6
    for k in range(res.shape[0]):
        for j,_ in enumerate(sources):
            res_frac[k,j]+= np.sum(res[k]==j)
        counts[k]+=res.shape[1]
    anc_per_snp[c]=res
    print(f"  {c}: {res.shape[1]} SNP procesados", file=sys.stderr)

# proporciones globales por cepa
frac = res_frac / counts[:,None]
df=pd.DataFrame(frac, columns=sources); df.insert(0,"ID",[samples[i] for i in tgt_idx])
df.insert(1,"POP",[id2pop[samples[i]] for i in tgt_idx])
df["top_anc"]=df[sources].idxmax(1)
df["top_frac"]=df[sources].max(1)
df["minor_frac"]=1-df["top_frac"]
df.to_csv(f"{outdir}/loter_ancestry_proportions.tsv", sep="\t", index=False)
np.savez_compressed(f"{outdir}/loter_perSNP.npz", **{c:anc_per_snp[c] for c in chroms},
                    targets=np.array([samples[i] for i in tgt_idx]),
                    sources=np.array(sources))
print("Proporciones guardadas. Resumen minor_frac:", file=sys.stderr)
print(df[["POP","minor_frac"]].groupby("POP").describe(), file=sys.stderr)
