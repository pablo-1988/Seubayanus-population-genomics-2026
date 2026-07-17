#!/usr/bin/env python
# Suaviza la ancestria local de LOTER por voto de mayoria en ventanas de W bp
# (reduce el ruido SNP-a-SNP) y pinta cromosomas. NO re-corre LOTER: usa el perSNP.
# Uso: 14_painting_windows.py <POP> <W_kb>   |   14_painting_windows.py DEMO
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
srcs=["PA","PB1","PB2","PB3","PB4"]
COL={"PA":"#F4D03F","PB1":"#E74C3C","PB2":"#2E86C1","PB3":"#1ABC9C","PB4":"#34495E"}
d=np.load(f"{proj}/results/K5/loter/loter_perSNP.npz", allow_pickle=True)
targets=list(d["targets"]); chroms=sorted(k for k in d.files if k.startswith("CBS12357"))
pos_by={c: d[f"POS_{c}"] for c in chroms}
prop=pd.read_csv(f"{proj}/results/K5/loter/loter_ancestry_proportions.tsv",sep="\t").set_index("ID")

def smooth(anc, pos, W):
    """voto de mayoria en ventanas de W bp -> lista de (start,end,ancestry)."""
    if W<=0:  # raw: bloques por tramos consecutivos
        idx=np.where(np.diff(anc)!=0)[0]; st=np.r_[0,idx+1]; en=np.r_[idx,len(anc)-1]
        return [(pos[s],pos[e],anc[s]) for s,e in zip(st,en)]
    out=[]; lo=pos.min(); hi=pos.max()
    edges=np.arange(lo, hi+W, W)
    for a,b in zip(edges[:-1],edges[1:]):
        m=(pos>=a)&(pos<b)
        if not m.any(): continue
        vals=anc[m]; maj=np.bincount(vals,minlength=len(srcs)).argmax()
        out.append((a,b,maj))
    # fusionar ventanas contiguas de misma ancestria
    merged=[]
    for s,e,v in out:
        if merged and merged[-1][2]==v and abs(merged[-1][1]-s)<1e-6:
            merged[-1]=(merged[-1][0],e,v)
        else: merged.append((s,e,v))
    return merged

def paint(pop, Wbp, fname, title):
    pick=[s for s in prop[prop.POP==pop].index if s in targets]
    fig,axes=plt.subplots(len(pick),1,figsize=(13,1.9*len(pick)))
    if len(pick)==1: axes=[axes]
    for ax,strain in zip(axes,pick):
        ti=targets.index(strain)
        for yi,c in enumerate(chroms):
            for s,e,v in smooth(d[c][ti], pos_by[c], Wbp):
                ax.barh(yi, e-s+1, left=s, height=0.8, color=COL[srcs[v]], edgecolor='none')
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms],fontsize=5)
        r=prop.loc[strain]
        ax.set_title(f"{strain} ({pop}; "+" ".join(f"{s}={r[s]:.2f}" for s in srcs)+")",fontsize=8.5,loc="left")
        ax.invert_yaxis()
    axes[-1].set_xlabel("posición (bp)",fontsize=8)
    fig.legend(handles=[Patch(color=COL[s],label=s) for s in srcs],ncol=5,loc="lower center",
               frameon=False,bbox_to_anchor=(0.5,-0.006))
    plt.suptitle(title,y=1.0,fontsize=11)
    plt.tight_layout(rect=[0,0.015,1,0.99])
    plt.savefig(fname,dpi=140,bbox_inches="tight"); plt.close()
    print("guardado:",fname,"| cepas:",len(pick))

if sys.argv[1]=="DEMO":
    # una cepa Hol a 4 resoluciones para ver si la fragmentacion colapsa
    strain="CDFM21L"; ti=targets.index(strain)
    Ws=[0,10000,25000,50000]; labs=["por-SNP (crudo)","ventana 10 kb","ventana 25 kb","ventana 50 kb"]
    fig,axes=plt.subplots(len(Ws),1,figsize=(13,7.5))
    for ax,W,lab in zip(axes,Ws,labs):
        for yi,c in enumerate(chroms):
            for s,e,v in smooth(d[c][ti], pos_by[c], W):
                ax.barh(yi, e-s+1, left=s, height=0.85, color=COL[srcs[v]], edgecolor='none')
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms],fontsize=5)
        ax.set_title(f"{strain} — {lab}",fontsize=10,loc="left"); ax.invert_yaxis()
    fig.legend(handles=[Patch(color=COL[s],label=s) for s in srcs],ncol=5,loc="lower center",
               frameon=False,bbox_to_anchor=(0.5,-0.01))
    plt.suptitle("¿Es artefacto el mosaico fino? Efecto del tamaño de ventana (voto de mayoría)",y=1.0,fontsize=12)
    plt.tight_layout(rect=[0,0.02,1,0.98])
    out=f"{proj}/results/K5/loter/painting_DEMO_ventanas.png"
    plt.savefig(out,dpi=150,bbox_inches="tight"); print("guardado:",out)
else:
    pop=sys.argv[1]; Wkb=float(sys.argv[2]); W=int(Wkb*1000)
    paint(pop, W, f"{proj}/results/K5/loter/painting_{pop}_{int(Wkb)}kb.png",
          f"Pintado cromosómico LOTER — {pop} (voto de mayoría en ventanas de {int(Wkb)} kb)")
