#!/usr/bin/env python
# Parte 4: pintado cromosomico (mosaicos de ancestria) para cepas representativas.
import allel, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
sources=["PA","PB1","PB2","PB3","PB4","PB5","PB6"]
cols=dict(zip(sources, plt.cm.tab10.colors[:7]))
d=np.load(f"{proj}/results/04_loter/loter_perSNP.npz", allow_pickle=True)
targets=list(d["targets"]); chroms=[k for k in d.files if k.startswith("CBS")]
cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz", fields=['variants/CHROM','variants/POS'])
CH=cal['variants/CHROM']; POS=cal['variants/POS']
pos_by={c: POS[CH==c] for c in chroms}
chroms=sorted(chroms)

# cepas representativas (deben estar en targets)
import pandas as pd
cl=pd.read_csv(f"{proj}/results/04_loter/loter_classification_all.tsv",sep="\t").set_index("ID")
pick=[]
for cand in ["ABFM5L","yHKS210","yHDPN421","1200_Ch","yHAB565","CL451.4","yHQL2270","yHAB56"]:
    if cand in targets: pick.append(cand)
pick=pick[:5]

fig,axes=plt.subplots(len(pick),1,figsize=(12,2.1*len(pick)))
if len(pick)==1: axes=[axes]
for ax,strain in zip(axes,pick):
    ti=targets.index(strain)
    for yi,c in enumerate(chroms):
        anc=d[c][ti]; p=pos_by[c]; L=p.max()
        # dibujar segmentos por tramos consecutivos de misma ancestria
        idx=np.where(np.diff(anc)!=0)[0]
        starts=np.r_[0, idx+1]; ends=np.r_[idx, len(anc)-1]
        for s,e in zip(starts,ends):
            ax.barh(yi, p[e]-p[s]+1, left=p[s], height=0.8,
                    color=cols[sources[anc[s]]], edgecolor='none')
    ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms], fontsize=6)
    lab=cl.loc[strain,"label"]; mf=cl.loc[strain,"minor_frac"]
    ax.set_title(f"{strain}  ({lab}, minor_frac={mf:.2f})", fontsize=10, loc="left")
    ax.set_xlabel("posicion (bp)", fontsize=8); ax.invert_yaxis()
fig.legend(handles=[Patch(color=cols[s],label=s) for s in sources],
           ncol=7, loc="lower center", frameon=False, bbox_to_anchor=(0.5,-0.02))
plt.tight_layout(rect=[0,0.03,1,1])
plt.savefig(f"{proj}/results/04_loter/chromosome_painting.png", dpi=150, bbox_inches="tight")
print("Pintado cromosomico:", pick)
