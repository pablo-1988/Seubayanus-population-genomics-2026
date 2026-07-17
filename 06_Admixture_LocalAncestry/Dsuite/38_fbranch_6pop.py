#!/usr/bin/env python
# f-branch SOLO 7 poblaciones (PA,PB1,PB2,PB3,PB4,Hol,NoAm; sin SoAm) -> tabla tidy + heatmap anotado.
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import norm
plt.rcParams.update({"pdf.fonttype":42,"font.size":10,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
D=f"{proj}/results/K5/sep_Hol/02_dsuite"; OUT=f"{proj}/results/K5/figures_EN"

lines=[l for l in open(f"{D}/fbranch_stats.txt").read().split("\n")]
zi=[i for i,l in enumerate(lines) if l.startswith("# Z-scores")][0]
def parse(block):
    hdr=block[0].split("\t"); recips=hdr[2:]; rows=[]
    for l in block[1:]:
        if not l.strip(): continue
        c=l.split("\t")
        for r,v in zip(recips,c[2:]): rows.append((c[0],c[1],r,v))
    return pd.DataFrame(rows,columns=["branch","donor_clade","recipient","val"])
fb=parse([l for l in lines[:zi] if l.strip()])
zz=parse([l for l in lines[zi+1:] if l.strip()])
m=fb.merge(zz,on=["branch","donor_clade","recipient"],suffixes=("_fb","_Z"))
m["f_b"]=pd.to_numeric(m["val_fb"],errors="coerce"); m["Z"]=pd.to_numeric(m["val_Z"],errors="coerce")
m=m.dropna(subset=["f_b"]).drop(columns=["val_fb","val_Z"])
m["SE"]=np.where(m["Z"]>0,m["f_b"]/m["Z"],np.nan)
m["p_value"]=np.where(m["Z"]>0,2*(1-norm.cdf(m["Z"].abs())),np.nan)
def stars(z): return "n.s." if (not np.isfinite(z) or z<=0) else ("***" if z>3.29 else "**" if z>2.58 else "*" if z>1.96 else "n.s.")
m["signif"]=m["Z"].apply(stars)
m=m.rename(columns={"recipient":"recipient_pop"})[
    ["branch","donor_clade","recipient_pop","f_b","Z","SE","p_value","signif"]].sort_values("Z",ascending=False)
m.to_csv(f"{D}/fbranch_tidy_6pop.tsv",sep="\t",index=False,float_format="%.5g")
print(f"pares no-nan={len(m)}  f_b>0={(m['f_b']>0).sum()}  Z>3={(m['Z']>3).sum()}")
print(m[m['Z']>3].to_string(index=False,formatters={"f_b":"{:.4f}".format,"Z":"{:.2f}".format,
      "SE":"{:.4f}".format,"p_value":"{:.2e}".format}))

# --- heatmap ---
recips=["PA","PB1","PB2","PB3","PB4","Hol"]
border=["b4","b5","b6","b7","b8","b9","b10","b11"]
blab={}
for _,r in m.drop_duplicates("branch").iterrows():
    d=r["donor_clade"]; blab[r["branch"]]=d if "," not in d else f"{d.split(',')[0]}+ ({r['branch']})"
FB=pd.DataFrame(index=border,columns=recips,dtype=float); ZS=FB.copy()
for _,r in m.iterrows():
    if r["branch"] in border and r["recipient_pop"] in recips:
        FB.loc[r["branch"],r["recipient_pop"]]=r["f_b"]; ZS.loc[r["branch"],r["recipient_pop"]]=r["Z"]
FBv=FB.values.astype(float); ZSv=ZS.values.astype(float)
FBplot=np.where(FBv>0,FBv,np.nan)
cmap=LinearSegmentedColormap.from_list("fb",["#FFF7EC","#FDBB84","#D7301F","#7F0000"]); cmap.set_bad("#F2F2F2")
fig,ax=plt.subplots(figsize=(6.8,5.6))
im=ax.imshow(FBplot,cmap=cmap,vmin=0,vmax=np.nanmax(FBplot),aspect="auto")
ax.set_xticks(range(len(recips))); ax.set_xticklabels(recips)
ax.set_yticks(range(len(border))); ax.set_yticklabels([blab[b] for b in border],fontsize=9)
ax.set_xlabel("Recipient population (gene flow INTO)"); ax.set_ylabel("Donor branch")
for i in range(FBplot.shape[0]):
    for j in range(FBplot.shape[1]):
        z=ZSv[i,j]; v=FBplot[i,j]
        if np.isfinite(v) and np.isfinite(z) and z>1.96:
            st="***" if z>3.29 else "**" if z>2.58 else "*"
            ax.text(j,i,st,ha="center",va="center",fontsize=9,
                    color="white" if v>0.20 else "black",fontweight="bold")
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.set_label("$f_b$ (admixture proportion)")
ax.set_title("f-branch (pure lineages + Hol (no NoAm))\n* Z>1.96, ** Z>2.58, *** Z>3.29 (block jackknife)",
             fontweight="bold",fontsize=11)
plt.tight_layout()
fig.savefig(f"{OUT}/37_fbranch_6pop.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/37_fbranch_6pop.png",dpi=160,bbox_inches="tight")
print("\n-> fbranch_tidy_6pop.tsv , 37_fbranch_6pop")
