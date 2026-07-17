#!/usr/bin/env python
# Paso 3: compara LSD2 (TMRCA del arbol) vs divergencia neta (da) + tabla de reescalado
# por generaciones/ano, y figuras en ingles/PDF.
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"pdf.fonttype":42,"font.size":11})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
D=f"{proj}/results/K5/dating"; OUT=f"{proj}/results/K5/figures_EN"
MU=1.67e-10; G_REF=2920
COL={"PA":"#F4D03F","PB1":"#E74C3C","PB2":"#2E86C1","PB3":"#1ABC9C","PB4":"#34495E",
     "NoAm":"#000000","Hol":"#E91E63"}
PARAPHYLETIC={"PB1"}   # no monofiletica en el arbol ML -> su TMRCA no es un split real

da=pd.read_csv(f"{D}/split_times_da.tsv",sep="\t")
ls=pd.read_csv(f"{D}/tmrca_lsd2.tsv",sep="\t")
crown=pd.read_csv(f"{D}/crown_ages_lsd2.tsv",sep="\t")

m=da.merge(ls,on=["pop1","pop2"])
m["tree_reliable"]=~(m.pop1.isin(PARAPHYLETIC)|m.pop2.isin(PARAPHYLETIC))
m["ratio_tree_over_da"]=m.TMRCA_years/m.T_years
m=m.sort_values("T_years")
m.to_csv(f"{D}/dating_comparison.tsv",sep="\t",index=False)

# --- tabla de reescalado (T escala como 1/g) ---
rows=[]
for _,r in m.iterrows():
    for g in [365,1000,2920]:
        k=G_REF/g
        rows.append([r.pop1,r.pop2,g,r.T_years*k,r.TMRCA_years*k])
resc=pd.DataFrame(rows,columns=["pop1","pop2","gen_per_year","T_split_da_years","TMRCA_tree_years"])
resc.to_csv(f"{D}/rescaling_generations.tsv",sep="\t",index=False)

print("=== comparacion (g=2920) ===")
p=m[["pop1","pop2","da","T_years","T_lo95","T_hi95","TMRCA_years","tree_reliable"]].copy()
for c in ["T_years","T_lo95","T_hi95","TMRCA_years"]: p[c]=p[c].round(0).astype(int)
p["da"]=p["da"].round(5)
print(p.to_string(index=False))
print("\n=== crown age (LSD2) ===")
print(crown.assign(crown_years=crown.crown_years.round(0)).to_string(index=False))

# ================= FIGURA A: comparacion de metodos =================
fig,ax=plt.subplots(figsize=(11,6))
mm=m.sort_values("T_years").reset_index(drop=True)
y=np.arange(len(mm)); h=0.38
ax.barh(y+h/2, mm.T_years, h, color="#2E86C1", label="Net divergence $d_a$ (split time)")
ax.errorbar(mm.T_years, y+h/2, xerr=[mm.T_years-mm.T_lo95, mm.T_hi95-mm.T_years],
            fmt='none', ecolor="black", capsize=2, lw=0.8)
cols=["#95A5A6" if not r else "#E74C3C" for r in mm.tree_reliable]
ax.barh(y-h/2, mm.TMRCA_years, h, color=cols, label="Tree TMRCA (IQ-TREE3 + LSD2)")
ax.set_yticks(y); ax.set_yticklabels([f"{a}–{b}" for a,b in zip(mm.pop1,mm.pop2)],fontsize=9)
ax.set_xlabel(f"Years before present  (μ={MU:.2e}/bp/gen × {G_REF} gen/year)")
ax.set_title("Population divergence times: net divergence vs dated tree",pad=28)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.10),ncol=2,fontsize=9,frameon=False)
ax.set_xlim(0,3700)
ax.text(0.99,0.02,"grey = tree TMRCA unreliable (PB1 is paraphyletic)",transform=ax.transAxes,
        ha="right",va="bottom",fontsize=8,color="#7F8C8D")
ax.invert_yaxis(); ax.grid(axis="x",alpha=0.25)
plt.tight_layout()
fig.savefig(f"{OUT}/26_divergence_times_comparison.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/26_divergence_times_comparison.png",dpi=160,bbox_inches="tight")
print("\n-> 26_divergence_times_comparison")

# ================= FIGURA B: reescalado por generaciones/ano =================
fig,ax=plt.subplots(figsize=(9,5.5))
for g,ls_ in zip([365,1000,2920],["-","--",":"]):
    sub=resc[resc.gen_per_year==g].sort_values("T_split_da_years")
    ax.plot(range(len(sub)),sub.T_split_da_years,ls_,marker="o",ms=4,
            label=f"{g} generations/year")
ax.set_xticks(range(len(mm))); ax.set_xticklabels([f"{a}–{b}" for a,b in zip(mm.pop1,mm.pop2)],
                                                  rotation=90,fontsize=8)
ax.set_ylabel("Split time (years before present)"); ax.set_yscale("log")
ax.set_title("Divergence times scale as 1/g — the generations-per-year assumption dominates")
ax.legend(); ax.grid(alpha=0.25)
plt.tight_layout()
fig.savefig(f"{OUT}/27_dating_generation_sensitivity.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/27_dating_generation_sensitivity.png",dpi=160,bbox_inches="tight")
print("-> 27_dating_generation_sensitivity")
