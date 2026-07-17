#!/usr/bin/env python
# IC de las edades de nodo del time tree por block-jackknife (deja fuera 1 cromosoma).
# Topologia FIJA (best_ml.nwk) + tasa FIJA -> la incertidumbre viene de los largos de rama.
# Para cada replica: sub-alineamiento sin el cromosoma c + fconst ajustado -> IQ-TREE -te + LSD2.
import allel, numpy as np, os, sys, subprocess, time
from collections import Counter
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
D=f"{proj}/results/K5/dating"; JK=f"{D}/jackknife"; os.makedirs(JK,exist_ok=True)
IQ=f"{proj}/work/iqtree-3.1.3-macOS-arm/bin/iqtree3"
NUC=[f"CBS12357_Chr{i:02d}_polished" for i in range(1,17)]

# --- columnas del alineamiento por cromosoma ---
cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz",fields=['variants/CHROM','variants/REF'])
CH=cal['variants/CHROM']; REF=cal['variants/REF']

# --- sitios constantes por cromosoma (de la referencia) ---
comp={c:Counter() for c in NUC}
name=None; buf=[]
def flush():
    if name in comp: comp[name].update("".join(buf).upper())
for line in open(f"{proj}/Data/CBS12357.fa"):
    if line.startswith(">"):
        flush(); name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
flush()
refsnp={c:Counter(REF[CH==c]) for c in NUC}
const={c:{b: comp[c][b]-refsnp[c].get(b,0) for b in "ACGT"} for c in NUC}

# --- alineamiento completo ---
names=[]; seqs=[]
for l in open(f"{D}/aln_381.fasta"):
    if l.startswith(">"): names.append(l[1:].strip())
    else: seqs.append(l.strip())
A=np.array([list(s) for s in seqs])       # (381, n_snp)
print("alineamiento:",A.shape,file=sys.stderr)

only=sys.argv[1] if len(sys.argv)>1 else None   # correr solo 1 cromosoma (para cronometrar)
todo=[only] if only else NUC
for c in todo:
    t0=time.time()
    keep=(CH!=c)
    sub=A[:,keep]
    fa=f"{JK}/no_{c}.fasta"
    with open(fa,"w") as fh:
        for n,row in zip(names,sub): fh.write(f">{n}\n{''.join(row)}\n")
    fc={b: sum(const[x][b] for x in NUC if x!=c) for b in "ACGT"}
    fcs=",".join(str(fc[b]) for b in "ACGT")
    # OJO: LSD2 corta el argumento de -w en el primer espacio -> usar rutas RELATIVAS
    # (el proyecto vive en "/Volumes/Extreme SSD/..." que tiene un espacio).
    rel=lambda p: os.path.relpath(p, proj)
    cmd=[IQ,"-s",rel(fa),"-fconst",fcs,"-m","GTR+F+I{0.991}+G4{1.499}",
         "-te",rel(f"{D}/best_ml.nwk"),"-o","CBS7001",
         "--date",rel(f"{D}/dates2.txt"),"--date-options",f"-w {rel(D)}/rate.txt",
         "--date-no-outgroup","-T","12","--prefix",rel(f"{JK}/{c}"),"-redo","-quiet"]
    r=subprocess.run(cmd,capture_output=True,text=True,cwd=proj)
    ok=os.path.exists(f"{JK}/{c}.timetree.nwk") and os.path.getsize(f"{JK}/{c}.timetree.nwk")>0
    print(f"{c}: {'OK' if ok else 'FALLO'}  ({time.time()-t0:.0f}s)",file=sys.stderr)
    if not ok: print(r.stdout[-500:],r.stderr[-500:],file=sys.stderr)
    os.remove(fa)
