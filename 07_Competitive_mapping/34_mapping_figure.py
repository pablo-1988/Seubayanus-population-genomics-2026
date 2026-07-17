#!/usr/bin/env python
# Tabla + figura del mapeo competitivo (eubayanus CBS12357 vs cerevisiae S288C).
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype":42,"font.size":11,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"; OUT=f"{proj}/results/K5/figures_EN"
df=pd.read_csv(f"{proj}/results/competitive_mapping/mapping_summary.tsv",sep="\t")

meta={"HYBRID_CBS1483":("CBS1483","Lager hybrid (control)"),
      "Hol_yHRVM108":("yHRVM108","Hol"),"PA_QC18":("QC18","PA"),
      "PB1_CL1005":("CL1005","PB1"),"PB2_CL1101":("CL1101","PB2"),
      "PB3_CL905":("CL905","PB3"),"PB4_CO150":("CO150","PB4"),
      "SoAm1_NR2":("NR2","SoAm1")}
df["strain"]=df["sample"].map(lambda x:meta[x][0])
df["pop"]=df["sample"].map(lambda x:meta[x][1])
order=["HYBRID_CBS1483","PA_QC18","PB1_CL1005","PB2_CL1101","PB3_CL905","PB4_CO150","Hol_yHRVM108","SoAm1_NR2"]
df=df.set_index("sample").loc[order].reset_index()

# --- tabla limpia (CSV) ---
tab=df[["pop","strain","total_reads","pct_mapped",
        "eub_breadth_pct","eub_meandepth","scer_breadth_pct","scer_meandepth",
        "pct_reads_eub","pct_reads_scer"]].copy()
tab.columns=["Population","Strain","Total_reads","Pct_reads_mapped",
             "Eub_genome_covered_pct","Eub_mean_depth","Scer_genome_covered_pct","Scer_mean_depth",
             "Pct_mapped_reads_to_Eub","Pct_mapped_reads_to_Scer"]
tab.to_csv(f"{proj}/results/competitive_mapping/mapping_table.csv",index=False)
print(tab.to_string(index=False))

# --- figura: breadth por genoma (barras agrupadas) ---
fig,ax=plt.subplots(figsize=(11,5.6))
y=np.arange(len(df))[::-1]; h=0.38
lab=[f"{r.pop}\n{r.strain}" for r in df.itertuples()]
ax.barh(y+h/2, df["eub_breadth_pct"], height=h, color="#2E86C1",
        label="S. eubayanus (CBS12357)", edgecolor="white")
ax.barh(y-h/2, df["scer_breadth_pct"], height=h, color="#E67E22",
        label="S. cerevisiae (S288C)", edgecolor="white")
for i,r in zip(y,df.itertuples()):
    ax.text(r.eub_breadth_pct+1, i+h/2, f"{r.eub_breadth_pct:.1f}%  ({r.eub_meandepth:.0f}x)",
            va="center",fontsize=8.5)
    ax.text(r.scer_breadth_pct+1, i-h/2, f"{r.scer_breadth_pct:.1f}%  ({r.scer_meandepth:.1f}x)",
            va="center",fontsize=8.5,
            color="#B9540A" if r.scer_breadth_pct>50 else "grey")
ax.set_yticks(y); ax.set_yticklabels(lab,fontsize=9)
ax.set_xlim(0,118); ax.set_xlabel("Genome covered (breadth, %)  [mean depth in parentheses]")
ax.axhline(len(df)-1.5,color="grey",lw=0.8,ls="--")
ax.set_title("Competitive mapping of Nanopore reads to a concatenated\nS. eubayanus + S. cerevisiae reference",
             fontweight="bold",fontsize=12)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=2,frameon=False)
ax.spines[["top","right"]].set_visible(False)
fig.text(0.5,-0.10,
 "Only the lager hybrid (CBS1483) covers BOTH genomes genome-wide (~96%). All pure and admixed "
 "S. eubayanus strains cover ~99% of the eubayanus genome and <0.6% of the cerevisiae genome "
 "(trace reads pile up on conserved loci only) -> they are not interspecific hybrids.",
 ha="center",fontsize=8.6,wrap=True)
plt.tight_layout(rect=[0,0.02,1,1])
fig.savefig(f"{OUT}/34_competitive_mapping.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/34_competitive_mapping.png",dpi=160,bbox_inches="tight")
print("\n-> 34_competitive_mapping + mapping_table.csv")
