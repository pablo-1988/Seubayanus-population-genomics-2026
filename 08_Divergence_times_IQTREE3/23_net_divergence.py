#!/usr/bin/env python
# Paso 2 datacion: tiempos de split por divergencia neta de Nei.
#   dxy = Sum_SNPs [p1(1-p2) + p2(1-p1)] / L      (divergencia absoluta por bp)
#   da  = dxy - (pi_x + pi_y)/2                    (remueve el polimorfismo ancestral)
#   T   = da / (2 * rate)   con rate = mu * g = 1.67e-10 * 2920 = 4.876e-7 sust/sitio/ano
# CIs por block-jackknife sobre los 16 cromosomas.
import allel, numpy as np, pandas as pd, sys, os
from itertools import combinations
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
OUT=f"{proj}/results/K5/dating"; os.makedirs(OUT,exist_ok=True)
MU=1.67e-10; G=2920; RATE=MU*G          # 4.876e-7 sust/sitio/ano
L=11863845                              # genoma nuclear (16 cromosomas)
POPS=["PA","PB1","PB2","PB3","PB4","NoAm","Hol"]

pop=dict(l.split() for l in open(f"{OUT}/taxa_pop.tsv"))
print("cargando VCF...",file=sys.stderr)
cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz",
                   fields=['samples','calldata/GT','variants/CHROM'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT'])
CH=cal['variants/CHROM']; chroms=sorted(set(CH))
idx={p:[i for i,s in enumerate(samples) if pop.get(s)==p] for p in POPS}
print({p:len(idx[p]) for p in POPS},file=sys.stderr)

# allele counts por pop y por cromosoma (para el jackknife)
ac_chr={}   # (pop, chrom) -> AlleleCounts
for c in chroms:
    m=(CH==c); g=gt.compress(m,axis=0)
    for p in POPS:
        ac_chr[(p,c)]=g.take(idx[p],axis=1).count_alleles(max_allele=1)

def stats(chrom_subset):
    """pi por pop y dxy/da/T por par, usando solo los cromosomas dados."""
    pi={}
    for p in POPS:
        s=sum(np.nansum(allel.mean_pairwise_difference(ac_chr[(p,c)],fill=0)) for c in chrom_subset)
        pi[p]=s/L
    rows={}
    for a,b in combinations(POPS,2):
        s=sum(np.nansum(allel.mean_pairwise_difference_between(
                ac_chr[(a,c)],ac_chr[(b,c)],fill=0)) for c in chrom_subset)
        dxy=s/L
        da=dxy-(pi[a]+pi[b])/2
        rows[(a,b)]=(dxy,da,da/(2*RATE))
    return pi,rows

pi_all,res_all=stats(chroms)
# --- block jackknife (leave-one-chromosome-out) ---
jk={k:[] for k in res_all}
for c in chroms:
    sub=[x for x in chroms if x!=c]
    _,r=stats(sub)
    for k in jk: jk[k].append(r[k][2])   # T
n=len(chroms)
rows=[]
for (a,b),(dxy,da,T) in res_all.items():
    v=np.array(jk[(a,b)])
    se=np.sqrt((n-1)/n*np.sum((v-v.mean())**2))    # SE jackknife
    rows.append([a,b,dxy,da,T,T-1.96*se,T+1.96*se])
df=pd.DataFrame(rows,columns=["pop1","pop2","dxy","da","T_years","T_lo95","T_hi95"])
df=df.sort_values("T_years")
df.to_csv(f"{OUT}/split_times_da.tsv",sep="\t",index=False)
pd.Series(pi_all).to_csv(f"{OUT}/pi_all_pops.csv",header=["pi"])

print(f"\n=== rate = {RATE:.4g} sust/sitio/ano (mu={MU}, g={G}) ===",file=sys.stderr)
print("=== pi por poblacion ===",file=sys.stderr)
print(pd.Series(pi_all).round(6).to_string(),file=sys.stderr)
print("\n=== tiempos de split (divergencia neta da) ===",file=sys.stderr)
out=df.copy()
out["dxy"]=out["dxy"].round(5); out["da"]=out["da"].round(5)
for c in ["T_years","T_lo95","T_hi95"]: out[c]=out[c].round(0).astype(int)
print(out.to_string(index=False),file=sys.stderr)
