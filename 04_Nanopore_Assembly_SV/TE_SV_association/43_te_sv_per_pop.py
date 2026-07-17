#!/usr/bin/env python
# Enriquecimiento de breakpoints de SV en TE, SEPARADO POR POBLACION.
# Null por permutacion (reubicar cada breakpoint al azar dentro de su cromosoma).
import glob,os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
rng=np.random.default_rng(1)
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.size":11,
                     "font.family":"sans-serif","font.sans-serif":["Arial"]})
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
MUM=f"{proj}/Data/mumandco"; TE=f"{proj}/work/TE"; SV=f"{proj}/results/TE_SV"; OUT=f"{proj}/results/K5/figures_EN"
WIN=250; NPERM=2000

# --- cromosomas y mascara TE (cumsum booleano por cromosoma) ---
chrlen={l.split()[0]:int(l.split()[1]) for l in open(f"{TE}/genome_nuclear.txt")}
cum={c:np.zeros(L+1,dtype=np.int32) for c,L in chrlen.items()}
te=np.zeros(1)
for l in open(f"{TE}/TE_annotation_CBS12357.bed"):
    c,s,e=l.split()[:3]; s,e=int(s),int(e)
    if c in cum: cum[c][s:e]=1
cum={c:np.concatenate([[0],np.cumsum(v)]) for c,v in cum.items()}  # prefix sum len L+2
def in_te(c,s,e):  # overlap si hay algun bp TE en [s,e)
    s=max(0,s); e=min(chrlen[c],e)
    if e<=s: return False
    return (cum[c][e]-cum[c][s])>0

# --- pop por cepa ---
sub={}
for l in open(f"{proj}/results/K5/admix_subpops_noGen/SETS.txt"):
    a=l.rstrip("\n").split("\t")
    if len(a)==2: sub[a[0]]=a[1]
alias={"E12":"C12","S13HH":"S13_HH"}
# override (por el usuario): estas 2 cepas 'xxx' se asignan a su linaje
OVERRIDE={"CL467":"PB1","CL1111":"PB2"}
def popof(st):
    if st in OVERRIDE: return OVERRIDE[st]
    for k in (alias.get(st,st),st,st+".1"):
        if k in sub: return sub[k]
    return "??"

# --- cargar breakpoints por cepa/pop (+ composicion de tipos de SV) ---
recs=[]; alltypes=[]
for f in glob.glob(f"{MUM}/CBS_*_output/CBS_*.SVs_all.tsv"):
    st=os.path.basename(f).replace("CBS_","").replace(".SVs_all.tsv","")
    pop=popof(st); d=pd.read_csv(f,sep="\t")
    d=d[d.ref_chr.str.contains("Chr")&~d.ref_chr.str.contains("mtDNA")]
    alltypes.extend(d.SV_type.tolist())
    for _,r in d.iterrows():
        c=r.ref_chr+"_polished"
        if c not in chrlen: continue
        for p in (r.ref_start,r.ref_stop):
            recs.append((pop,st,c,int(p)))
bp=pd.DataFrame(recs,columns=["pop","strain","chr","pos"])
svtypes=pd.Series(alltypes)
print("breakpoints:",len(bp),"| pops:",bp["pop"].value_counts().to_dict())

def enrich(subdf):
    chrs=subdf["chr"].values; pos=subdf["pos"].values; n=len(subdf)
    obs=sum(in_te(c,p-WIN,p+WIN) for c,p in zip(chrs,pos))
    Ls=np.array([chrlen[c] for c in chrs])
    null=np.empty(NPERM)
    for i in range(NPERM):
        rs=(rng.random(n)*(Ls-2*WIN)).astype(int)  # start aleatorio por cromosoma
        null[i]=sum(in_te(c,s,s+2*WIN) for c,s in zip(chrs,rs))
    mu=null.mean(); fold=obs/mu if mu>0 else np.inf
    p=(np.sum(null>=obs)+1)/(NPERM+1)
    return n,obs,100*obs/n,mu,100*mu/n,fold,p

POPS=["PA","PB1","PB2","PB3","PB4","Hol","NoAm","SoAm1","xxx"]
rows=[]
for pop in POPS:
    s=bp[bp["pop"]==pop]
    if len(s)==0: continue
    nstr=s.strain.nunique()
    n,obs,pobs,mu,pexp,fold,pv=enrich(s)
    rows.append(dict(pop=pop,n_strains=nstr,breakpoints=n,in_TE=obs,pct_TE=round(pobs,1),
                     exp_pct=round(pexp,2),fold=round(fold,1),p=pv))
res=pd.DataFrame(rows)
res["pop"]=res["pop"].replace({"xxx":"Admix"})
res.to_csv(f"{SV}/TE_enrichment_per_pop.tsv",sep="\t",index=False)
print("\n",res.to_string(index=False))

# --- SVs por cepa (para panel C) ---
perstr=(bp.groupby(["pop","strain"]).size()/2).reset_index(name="SVs")
perstr["pop"]=perstr["pop"].replace({"xxx":"Admix"})

# --- figura 2x2 ---
from matplotlib.patches import Patch
CPOP={"PA":"#F4D03F","PB1":"#E74C3C","PB2":"#2E86C1","PB3":"#1ABC9C","PB4":"#34495E",
      "Hol":"#6EC6FF","NoAm":"#7F8C8D","SoAm1":"#9B59B6","Admix":"#D0D3D4"}
fig,ax=plt.subplots(1,3,figsize=(17,5))

# Panel A: composicion de SVs por tipo (mobile vs otro)  [antes fig 40A]
order=["insertion_mobile","deletion_mobile","deletion_novel","insertion_novel",
       "duplication","contraction","transloc","inversion"]
vc=svtypes.value_counts(); vals=[vc.get(t,0) for t in order]
cols_t=["#E67E22" if "mobile" in t else "#95A5A6" for t in order]
ax[0].barh(range(len(order))[::-1],vals,color=cols_t,edgecolor="black",lw=0.5)
ax[0].set_yticks(range(len(order))[::-1]); ax[0].set_yticklabels(order,fontsize=9)
for i,v in zip(range(len(order))[::-1],vals): ax[0].text(v+15,i,str(v),va="center",fontsize=8.5)
mobfrac=svtypes.str.contains("mobile").mean()
ax[0].set_xlabel("SV count (22 strains)")
ax[0].set_title(f"A  Structural variants by type\n{len(svtypes)} SVs; {mobfrac:.0%} TE-mobile (MUM&Co)",
                loc="left",fontweight="bold",fontsize=11)
ax[0].legend(handles=[Patch(fc="#E67E22",ec="black",lw=.5,label="mobile (TE)"),
                      Patch(fc="#95A5A6",ec="black",lw=.5,label="other")],
             frameon=False,fontsize=9,loc="lower right")
ax[0].spines[["top","right"]].set_visible(False); ax[0].set_xlim(0,max(vals)*1.15)

# Panel B: fold enrichment por poblacion
r=res.sort_values("fold")
ax[1].barh(r["pop"],r["fold"],color=[CPOP.get(p,"#888") for p in r["pop"]],edgecolor="black",lw=0.5)
for i,(f,pv) in enumerate(zip(r["fold"],r["p"])):
    sig="***" if pv<0.001 else "**" if pv<0.01 else "*" if pv<0.05 else "n.s."
    ax[1].text(f+0.3,i,f"{f:.0f}x {sig}",va="center",fontsize=9)
ax[1].set_xlabel("Fold enrichment of SV breakpoints in TE")
ax[1].set_title("B  TE-mediated SV enrichment by population",loc="left",fontweight="bold",fontsize=11)
ax[1].axvline(1,color="grey",ls=":",lw=1); ax[1].spines[["top","right"]].set_visible(False)
ax[1].set_xlim(0,r["fold"].max()*1.25)

# Panel C: % breakpoints en TE por poblacion
r2=res.sort_values("pct_TE")
ax[2].barh(r2["pop"],r2["pct_TE"],color=[CPOP.get(p,"#888") for p in r2["pop"]],edgecolor="black",lw=0.5)
ax[2].axvline(res["exp_pct"].mean(),color="#C0392B",ls="--",lw=1.2)
ax[2].text(res["exp_pct"].mean()+0.3,0,f"background\n({res['exp_pct'].mean():.1f}%)",color="#C0392B",fontsize=8)
for i,v in enumerate(r2["pct_TE"]): ax[2].text(v+0.4,i,f"{v:.1f}%",va="center",fontsize=9)
ax[2].set_xlabel("% SV breakpoints overlapping TE")
ax[2].set_title("C  Fraction of breakpoints in TE",loc="left",fontweight="bold",fontsize=11)
ax[2].spines[["top","right"]].set_visible(False); ax[2].set_xlim(0,r2["pct_TE"].max()*1.2)

plt.tight_layout()
fig.savefig(f"{OUT}/41_TE_SV_per_pop.pdf",bbox_inches="tight")
fig.savefig(f"{OUT}/41_TE_SV_per_pop.png",dpi=300,bbox_inches="tight")
print("\n-> 41_TE_SV_per_pop ; TE_enrichment_per_pop.tsv")
