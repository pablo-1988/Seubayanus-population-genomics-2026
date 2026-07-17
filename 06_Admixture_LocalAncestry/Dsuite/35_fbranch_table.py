#!/usr/bin/env python
# Convierte la matriz f-branch (+ Z-scores de Dsuite Fbranch -Z) a tabla tidy (long)
# con f_b, Z, SE (=f_b/Z), p-value (dos colas) y significancia. Permite estadistica
# entre valores (z = (fb1-fb2)/sqrt(SE1^2+SE2^2)).
import numpy as np, pandas as pd
from scipy.stats import norm
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
D=f"{proj}/results/K5/admix_subpops_noGen/02_dsuite"

lines=open(f"{D}/fbranch_stats.txt").read().split("\n")
# separar bloque f_b y bloque Z (marcado por "# Z-scores:")
zi=[i for i,l in enumerate(lines) if l.startswith("# Z-scores")][0]
def parse(block):
    hdr=block[0].split("\t"); recips=hdr[2:]
    rows=[]
    for l in block[1:]:
        if not l.strip(): continue
        c=l.split("\t"); br=c[0]; desc=c[1]; vals=c[2:]
        for r,v in zip(recips,vals):
            rows.append((br,desc,r,v))
    return pd.DataFrame(rows,columns=["branch","branch_descendants","recipient","val"])
fb=parse([l for l in lines[:zi] if l.strip()])
zz=parse([l for l in lines[zi+1:] if l.strip()])
m=fb.merge(zz,on=["branch","branch_descendants","recipient"],suffixes=("_fb","_Z"))
m["f_b"]=pd.to_numeric(m["val_fb"],errors="coerce")
m["Z"]=pd.to_numeric(m["val_Z"],errors="coerce")
m=m.drop(columns=["val_fb","val_Z"]).dropna(subset=["f_b"])

# SE y p-value (dos colas). Donde f_b=0 o Z=0 -> sin senal.
m["SE"]=np.where(m["Z"]>0, m["f_b"]/m["Z"], np.nan)
m["p_value"]=np.where(m["Z"]>0, 2*(1-norm.cdf(m["Z"].abs())), np.nan)
def stars(z):
    if not np.isfinite(z) or z==0: return "n.s."
    return "***" if z>3.29 else "**" if z>2.58 else "*" if z>1.96 else "n.s."
m["signif"]=m["Z"].apply(stars)

# donante legible = poblacion o clado que subtiende la rama
m=m.rename(columns={"branch_descendants":"donor_clade","recipient":"recipient_pop"})
m=m[["branch","donor_clade","recipient_pop","f_b","Z","SE","p_value","signif"]]
m=m.sort_values(["Z"],ascending=False).reset_index(drop=True)

# tabla completa (todos los pares no-nan)
m.to_csv(f"{D}/fbranch_tidy_full.tsv",sep="\t",index=False,float_format="%.5g")
# solo senales significativas (Z>3, ~p<0.0027)
sig=m[m["Z"]>3].copy()
sig.to_csv(f"{D}/fbranch_significant.tsv",sep="\t",index=False,float_format="%.5g")

print(f"Pares totales (no-nan): {len(m)}  |  con f_b>0: {(m['f_b']>0).sum()}  |  significativos (Z>3): {len(sig)}")
print("\n=== f-branch SIGNIFICATIVOS (Z>3) ===")
with pd.option_context("display.width",160,"display.max_rows",60):
    print(sig.to_string(index=False,
          formatters={"f_b":"{:.4f}".format,"Z":"{:.2f}".format,"SE":"{:.4f}".format,"p_value":"{:.2e}".format}))
print("\n-> fbranch_tidy_full.tsv , fbranch_significant.tsv")
