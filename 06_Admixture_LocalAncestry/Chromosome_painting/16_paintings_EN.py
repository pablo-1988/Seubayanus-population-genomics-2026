#!/usr/bin/env python
# Paintings cromosomicos en INGLES / PDF alta calidad (rasterizados para peso razonable).
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.size":10,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
OUT=f"{proj}/results/K5/figures_EN"; import os; os.makedirs(OUT,exist_ok=True)
srcs=["PA","PB1","PB2","PB3","PB4"]
COL={"PA":"#F4D03F","PB1":"#E74C3C","PB2":"#2E86C1","PB3":"#1ABC9C","PB4":"#34495E"}
d=np.load(f"{proj}/results/K5/loter/loter_perSNP.npz",allow_pickle=True)
targets=list(d["targets"]); chroms=sorted(k for k in d.files if k.startswith("CBS12357"))
pos_by={c:d[f"POS_{c}"] for c in chroms}
prop=pd.read_csv(f"{proj}/results/K5/loter/loter_ancestry_proportions.tsv",sep="\t").set_index("ID")

def smooth(anc,pos,W):
    if W<=0:
        idx=np.where(np.diff(anc)!=0)[0]; st=np.r_[0,idx+1]; en=np.r_[idx,len(anc)-1]
        return [(pos[s],pos[e],anc[s]) for s,e in zip(st,en)]
    out=[]; lo=pos.min(); hi=pos.max()
    for a,b in zip(np.arange(lo,hi+W,W)[:-1],np.arange(lo,hi+W,W)[1:]):
        m=(pos>=a)&(pos<b)
        if not m.any(): continue
        out.append((a,b,np.bincount(anc[m],minlength=len(srcs)).argmax()))
    mg=[]
    for s,e,v in out:
        if mg and mg[-1][2]==v and abs(mg[-1][1]-s)<1e-6: mg[-1]=(mg[-1][0],e,v)
        else: mg.append((s,e,v))
    return mg

def paint_pop(pop,W,name,title):
    # Vector (no rasterizado): cada tramo es un rectangulo real -> Illustrator conserva el color.
    # broken_barh agrupa por color -> menos objetos, PDF mas liviano y editable.
    pick=[s for s in prop[prop.POP==pop].index if s in targets]
    fig,axes=plt.subplots(len(pick),1,figsize=(13,1.9*len(pick)))
    if len(pick)==1: axes=[axes]
    for ax,st in zip(axes,pick):
        ti=targets.index(st)
        for yi,c in enumerate(chroms):
            segs={i:[] for i in range(len(srcs))}
            for s,e,v in smooth(d[c][ti],pos_by[c],W):
                segs[v].append((s,e-s+1))
            for v,spans in segs.items():
                if spans:
                    ax.broken_barh(spans,(yi-0.4,0.8),facecolors=COL[srcs[v]],
                                   edgecolors="none",linewidth=0,rasterized=False)
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms],fontsize=5)
        ax.set_ylim(-0.6,len(chroms)-0.4)
        r=prop.loc[st]
        ax.set_title(f"{st} ({pop}; "+" ".join(f"{s}={r[s]:.2f}" for s in srcs)+")",fontsize=8.5,loc="left")
        ax.invert_yaxis()
    axes[-1].set_xlabel("position (bp)",fontsize=8)
    fig.legend(handles=[Patch(facecolor=COL[s],edgecolor="none",label=s) for s in srcs],ncol=5,
               loc="lower center",frameon=False,bbox_to_anchor=(0.5,-0.006))
    plt.suptitle(title,y=1.0,fontsize=11); plt.tight_layout(rect=[0,0.015,1,0.99])
    fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight"); plt.close(); print("->",name)

# representativa (mayor minor_frac por sub-pop)
def paint_repr(name):
    pick=[]
    for p in ["NoAm","SoAm1","SoAm2","SoAm3","SoAm4","Hol"]:
        sub=prop[prop.POP==p].sort_values("minor_frac",ascending=False)
        if len(sub): pick.append(sub.index[0])
    fig,axes=plt.subplots(len(pick),1,figsize=(13,2.2*len(pick)))
    for ax,st in zip(axes,pick):
        ti=targets.index(st)
        for yi,c in enumerate(chroms):
            for s,e,v in smooth(d[c][ti],pos_by[c],0):
                ax.barh(yi,e-s+1,left=s,height=0.8,color=COL[srcs[v]],edgecolor='none',rasterized=True)
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms],fontsize=6)
        r=prop.loc[st]
        ax.set_title(f"{st}  ({r.POP}; minor={r.minor_frac:.2f}; top={r.top_anc} {r.top_frac:.2f})",fontsize=10,loc="left")
        ax.set_xlabel("position (bp)",fontsize=8); ax.invert_yaxis()
    fig.legend(handles=[Patch(color=COL[s],label=s) for s in srcs],ncol=5,loc="lower center",
               frameon=False,bbox_to_anchor=(0.5,-0.01))
    plt.suptitle("Local-ancestry chromosome painting (LOTER, K=5) — representative admixed strains",y=1.0,fontsize=11.5)
    plt.tight_layout(rect=[0,0.02,1,0.99]); fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight",dpi=300)
    plt.close(); print("->",name)

# DEMO ventanas
def demo(name):
    st="CDFM21L"; ti=targets.index(st); Ws=[0,10000,25000,50000]
    labs=["per-SNP (raw)","20 kb window","25 kb window","50 kb window"]; labs[1]="10 kb window"
    fig,axes=plt.subplots(len(Ws),1,figsize=(13,7.5))
    for ax,W,lab in zip(axes,Ws,labs):
        for yi,c in enumerate(chroms):
            for s,e,v in smooth(d[c][ti],pos_by[c],W):
                ax.barh(yi,e-s+1,left=s,height=0.85,color=COL[srcs[v]],edgecolor='none',rasterized=True)
        ax.set_yticks(range(len(chroms))); ax.set_yticklabels([c.split('_')[1] for c in chroms],fontsize=5)
        ax.set_title(f"{st} — {lab}",fontsize=10,loc="left"); ax.invert_yaxis()
    fig.legend(handles=[Patch(color=COL[s],label=s) for s in srcs],ncol=5,loc="lower center",
               frameon=False,bbox_to_anchor=(0.5,-0.01))
    plt.suptitle("Is the fine mosaic an artifact? Effect of window size (majority vote)",y=1.0,fontsize=12)
    plt.tight_layout(rect=[0,0.02,1,0.98]); fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight",dpi=300)
    plt.close(); print("->",name)

paint_repr("07_painting_representative")
demo("08_painting_window_DEMO")
paint_pop("Hol",25000,"09_painting_Hol_25kb","Local-ancestry chromosome painting — Hol (majority vote, 25 kb windows)")
paint_pop("NoAm",25000,"10_painting_NoAm_25kb","Local-ancestry chromosome painting — NoAm (majority vote, 25 kb windows)")
print("Paintings EN listas.")
