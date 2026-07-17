#!/usr/bin/env python
# Figura resumen: SVs (MUM&Co) y su asociacion con TE (anotados por BLAST).
import glob,os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.size":11,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
MUM=f"{proj}/Data/mumandco"; SV=f"{proj}/results/TE_SV"; OUT=f"{proj}/results/K5/figures_EN"

# --- cargar todos los SVs ---
rows=[]
for f in glob.glob(f"{MUM}/CBS_*_output/CBS_*.SVs_all.tsv"):
    st=os.path.basename(f).replace("CBS_","").replace(".SVs_all.tsv","")
    d=pd.read_csv(f,sep="\t"); d["strain"]=st; rows.append(d)
a=pd.concat(rows,ignore_index=True)
a=a[a.ref_chr.str.contains("Chr")&~a.ref_chr.str.contains("mtDNA")]
a["mobile"]=a.SV_type.str.contains("mobile")

# --- null del enriquecimiento ---
null=np.loadtxt(f"{SV}/null_counts.txt")
NB=6422; obs=1644
mu=null.mean(); fold=obs/mu; p=(np.sum(null>=obs)+1)/(len(null)+1)

fig,ax=plt.subplots(1,3,figsize=(15,4.8))

# Panel A: composicion de SV, mobile vs no
order=["insertion_mobile","deletion_mobile","deletion_novel","insertion_novel",
       "duplication","contraction","transloc","inversion"]
vc=a.SV_type.value_counts()
vals=[vc.get(t,0) for t in order]
cols=["#E67E22" if "mobile" in t else "#95A5A6" for t in order]
ax[0].barh(range(len(order))[::-1],vals,color=cols,edgecolor="black",lw=0.5)
ax[0].set_yticks(range(len(order))[::-1]); ax[0].set_yticklabels(order,fontsize=9)
for i,v in zip(range(len(order))[::-1],vals): ax[0].text(v+15,i,str(v),va="center",fontsize=8.5)
ax[0].set_xlabel("SV count (22 strains)")
ax[0].set_title(f"A  Structural variants by type\n{len(a)} SVs; {a.mobile.mean():.0%} TE-mobile (MUM&Co)",
                loc="left",fontweight="bold",fontsize=11)
ax[0].legend(handles=[Patch(fc="#E67E22",ec="black",lw=.5,label="mobile (TE)"),
                      Patch(fc="#95A5A6",ec="black",lw=.5,label="other")],
             frameon=False,fontsize=9,loc="lower right")
ax[0].spines[["top","right"]].set_visible(False); ax[0].set_xlim(0,max(vals)*1.15)

# Panel B: enriquecimiento (null + observado)
ax[1].hist(null,bins=30,color="#BDC3C7",edgecolor="white")
ax[1].axvline(obs,color="#C0392B",lw=2.5)
ax[1].text(obs,ax[1].get_ylim()[1]*0.9,f"  observed\n  {obs} ({100*obs/NB:.0f}%)",
           color="#C0392B",fontsize=10,va="top",fontweight="bold")
ax[1].text(mu,ax[1].get_ylim()[1]*0.55,f"expected\n{mu:.0f} ({100*mu/NB:.1f}%)",
           color="#555",fontsize=9,ha="center")
ax[1].set_xlabel("SV breakpoints overlapping TE (per permutation)")
ax[1].set_ylabel("Permutations")
ax[1].set_title(f"B  Breakpoints enriched in TE\n{fold:.0f}x enrichment; p<{1/len(null):.3f} ({len(null)} perms)",
                loc="left",fontweight="bold",fontsize=11)
ax[1].spines[["top","right"]].set_visible(False)

# Panel C: concordancia mobile vs novel
tags=["mobile","novel"]; nbp={"mobile":4182,"novel":2240}; ov={"mobile":1337,"novel":307}
pcts=[100*ov[t]/nbp[t] for t in tags]
bars=ax[2].bar(tags,pcts,color=["#E67E22","#95A5A6"],edgecolor="black",lw=0.5,width=0.6)
ax[2].axhline(100*mu/NB,color="#C0392B",ls="--",lw=1.3)
ax[2].text(1.4,100*mu/NB+1,f"genome background\n({100*mu/NB:.1f}%)",color="#C0392B",fontsize=8.5,ha="right")
for b,pc,t in zip(bars,pcts,tags): ax[2].text(b.get_x()+b.get_width()/2,pc+0.6,f"{pc:.1f}%",ha="center",fontsize=10)
ax[2].set_ylabel("Breakpoints overlapping TE (%)")
ax[2].set_title("C  MUM&Co 'mobile' flag vs\nindependent TE annotation",loc="left",fontweight="bold",fontsize=11)
ax[2].set_xticks([0,1]); ax[2].set_xticklabels(["mobile\n(MUM&Co)","novel\n(MUM&Co)"])
ax[2].spines[["top","right"]].set_visible(False); ax[2].set_ylim(0,max(pcts)*1.25)

plt.tight_layout()
fig.savefig(f"{OUT}/40_TE_SV_association.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/40_TE_SV_association.png",dpi=300,bbox_inches="tight")

# --- tabla resumen por cepa ---
summ=a.groupby("strain").agg(SVs=("SV_type","size"),mobile=("mobile","sum")).reset_index()
summ["pct_mobile"]=(100*summ["mobile"]/summ["SVs"]).round(1)
summ.to_csv(f"{SV}/SV_summary_per_strain.tsv",sep="\t",index=False)
print("-> 40_TE_SV_association ; SV_summary_per_strain.tsv")
print(f"TE annotation: 152 loci, 0.39% genome | enrichment {fold:.1f}x p<{1/len(null):.3f}")
