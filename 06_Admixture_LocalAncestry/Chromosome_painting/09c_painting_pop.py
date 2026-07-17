#!/usr/bin/env python
# Pintado cromosomico LOTER (K5) para TODAS las cepas de una sub-poblacion.
# Uso: 09c_painting_pop.py <POP>   (ej. Hol, NoAm)
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
POP=sys.argv[1]
sources=["PA","PB1","PB2","PB3","PB4"]
COL={"PA":"#F4D03F","PB1":"#E74C3C","PB2":"#2E86C1","PB3":"#1ABC9C","PB4":"#34495E"}
d=np.load(f"{proj}/results/K5/loter/loter_perSNP.npz", allow_pickle=True)
targets=list(d["targets"]); chroms=sorted(k for k in d.files if k.startswith("CBS12357"))
pos_by={c: d[f"POS_{c}"] for c in chroms}
prop=pd.read_csv(f"{proj}/results/K5/loter/loter_ancestry_proportions.tsv",sep="\t").set_index("ID")
pick=[s for s in prop[prop.POP==POP].index if s in targets]

fig,axes=plt.subplots(len(pick),1,figsize=(13,2.0*len(pick)))
if len(pick)==1: axes=[axes]
for ax,strain in zip(axes,pick):
    ti=targets.index(strain)
    for yi,c in enumerate(chroms):
        anc=d[c][ti]; p=pos_by[c]
        idx=np.where(np.diff(anc)!=0)[0]
        starts=np.r_[0, idx+1]; ends=np.r_[idx, len(anc)-1]
        for s,e in zip(starts,ends):
            ax.barh(yi, p[e]-p[s]+1, left=p[s], height=0.8,
                    color=COL[sources[anc[s]]], edgecolor='none')
    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels([c.split('_')[1] for c in chroms], fontsize=6)
    r=prop.loc[strain]
    ax.set_title(f"{strain}  ({POP}; " + " ".join(f"{s}={r[s]:.2f}" for s in sources) + ")",
                 fontsize=9, loc="left")
    ax.set_xlabel("posición (bp)", fontsize=8); ax.invert_yaxis()
fig.legend(handles=[Patch(color=COL[s],label=s) for s in sources],
           ncol=5, loc="lower center", frameon=False, bbox_to_anchor=(0.5,-0.005))
plt.suptitle(f"Pintado cromosómico LOTER — sub-población {POP} (K=5, todas las cepas)",
             y=1.0, fontsize=12)
plt.tight_layout(rect=[0,0.02,1,0.99])
out=f"{proj}/results/K5/loter/painting_{POP}_all.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print("guardado:", out, "| cepas:", len(pick))
