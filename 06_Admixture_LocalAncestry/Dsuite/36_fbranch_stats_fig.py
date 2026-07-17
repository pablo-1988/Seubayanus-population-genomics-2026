#!/usr/bin/env python
# (1) Matriz de comparaciones pareadas entre flujos f_b que entran a un receptor
#     (z = (fb1-fb2)/sqrt(SE1^2+SE2^2)); para Hol y NoAm.
# (2) Figura f-branch (heatmap f_b) anotada con significancia (Z-scores).
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import norm
plt.rcParams.update({"pdf.fonttype":42,"font.size":10,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
D=f"{proj}/results/K5/admix_subpops_noGen/02_dsuite"; OUT=f"{proj}/results/K5/figures_EN"
m=pd.read_csv(f"{D}/fbranch_tidy_full.tsv",sep="\t")

# etiqueta de donante legible: si la rama es una hoja usa la pop; si es clado, "clado (bN)"
def donor_label(r):
    d=r["donor_clade"]
    return d if "," not in d else f"{d.split(',')[0]}+.. ({r['branch']})"
m["donor"]=m.apply(donor_label,axis=1)

# ---------- (1) comparaciones pareadas ----------
def pairwise(recipient):
    s=m[(m["recipient_pop"]==recipient)&(m["Z"]>3)].copy().sort_values("f_b",ascending=False)
    s=s.reset_index(drop=True)
    rows=[]
    for i in range(len(s)):
        for j in range(i+1,len(s)):
            a,b=s.loc[i],s.loc[j]
            se=np.sqrt(a["SE"]**2+b["SE"]**2)
            z=(a["f_b"]-b["f_b"])/se
            p=2*(1-norm.cdf(abs(z)))
            rows.append({"recipient":recipient,"donor_A":a["donor"],"f_b_A":a["f_b"],
                         "donor_B":b["donor"],"f_b_B":b["f_b"],"diff":a["f_b"]-b["f_b"],
                         "z_diff":z,"p_diff":p,"signif":"***" if abs(z)>3.29 else "**" if abs(z)>2.58 else "*" if abs(z)>1.96 else "n.s."})
    return pd.DataFrame(rows)

allpw=[]
for rec in ["Hol","NoAm"]:
    pw=pairwise(rec); allpw.append(pw)
    pw.to_csv(f"{D}/fbranch_pairwise_{rec}.tsv",sep="\t",index=False,float_format="%.5g")
    print(f"\n=== Comparaciones pareadas de flujos hacia {rec} (Z_fb>3) ===")
    print(pw.to_string(index=False,formatters={"f_b_A":"{:.3f}".format,"f_b_B":"{:.3f}".format,
          "diff":"{:.3f}".format,"z_diff":"{:.2f}".format,"p_diff":"{:.2e}".format}))
pd.concat(allpw).to_csv(f"{D}/fbranch_pairwise_HolNoAm.tsv",sep="\t",index=False,float_format="%.5g")

# ---------- (2) heatmap f-branch anotado ----------
recips=["PA","NoAm","SoAm1","Hol","PB2","PB3","SoAm4","PB1","SoAm3","PB4","SoAm2"]
branch_order=["b5","b6","b8","b9","b10","b11","b12","b13","b14","b15","b16","b17","b18","b19","b20","b21"]
blab={}
for _,r in m.drop_duplicates("branch").iterrows():
    d=r["donor_clade"]; blab[r["branch"]]=d if "," not in d else f"{d.split(',')[0]}+ ({r['branch']})"
FB=pd.DataFrame(index=branch_order,columns=recips,dtype=float)
ZS=pd.DataFrame(index=branch_order,columns=recips,dtype=float)
for _,r in m.iterrows():
    if r["branch"] in branch_order and r["recipient_pop"] in recips:
        FB.loc[r["branch"],r["recipient_pop"]]=r["f_b"]
        ZS.loc[r["branch"],r["recipient_pop"]]=r["Z"]
FBv=FB.values.astype(float); ZSv=ZS.values.astype(float)
FBplot=np.where(FBv>0,FBv,np.nan)   # 0 = sin senal -> blanco

cmap=LinearSegmentedColormap.from_list("fb",["#FFF7EC","#FDBB84","#D7301F","#7F0000"])
cmap.set_bad("#F2F2F2")
fig,ax=plt.subplots(figsize=(9.5,7.5))
im=ax.imshow(FBplot,cmap=cmap,vmin=0,vmax=np.nanmax(FBplot),aspect="auto")
ax.set_xticks(range(len(recips))); ax.set_xticklabels(recips,rotation=45,ha="right")
ax.set_yticks(range(len(branch_order))); ax.set_yticklabels([blab[b] for b in branch_order],fontsize=8.5)
ax.set_xlabel("Recipient population (gene flow INTO)")
ax.set_ylabel("Donor branch")
# estrellas de significancia
for i in range(FBplot.shape[0]):
    for j in range(FBplot.shape[1]):
        z=ZSv[i,j]; v=FBplot[i,j]
        if np.isfinite(v) and np.isfinite(z) and z>1.96:
            st="***" if z>3.29 else "**" if z>2.58 else "*"
            col="white" if v>0.20 else "black"
            ax.text(j,i,st,ha="center",va="center",fontsize=8,color=col,fontweight="bold")
cb=fig.colorbar(im,ax=ax,fraction=0.045,pad=0.03); cb.set_label("$f_b$ (admixture proportion)")
ax.set_title("f-branch statistic with significance\n(* Z>1.96, ** Z>2.58, *** Z>3.29; block jackknife)",
             fontweight="bold",fontsize=11)
plt.tight_layout()
fig.savefig(f"{OUT}/35_fbranch_significance.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/35_fbranch_significance.png",dpi=160,bbox_inches="tight")
print("\n-> fbranch_pairwise_{Hol,NoAm}.tsv , 35_fbranch_significance")
