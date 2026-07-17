#!/usr/bin/env python
# Paso 0 datacion: VCF -> alineamiento FASTA pseudo-haploide (381 taxa) + conteo de
# sitios constantes (-fconst) desde CBS12357.fa.
# Taxa: puras umbral0.8 (PA,PB1,PB2,PB3,PB4) + NoAm + Hol + outgroup CBS7001.
import allel, numpy as np, os, sys
from collections import Counter
proj="/Volumes/Extreme SSD/Eubayanus-Pop"
OUT=f"{proj}/results/K5/dating"; os.makedirs(OUT,exist_ok=True)

# --- mapa de poblaciones ---
p8=dict(l.split() for l in open(f"{proj}/results/K5/umbral0.8/SETS.txt"))
ng=dict(l.split() for l in open(f"{proj}/results/K5/admix_subpops_noGen/SETS.txt"))
pop={}
for s,p in p8.items():
    if p in ("PA","PB1","PB2","PB3","PB4"): pop[s]=p
for s,p in ng.items():
    if p in ("NoAm","Hol"): pop[s]=p
pop["CBS7001"]="Outgroup"
print("taxa:",len(pop),Counter(pop.values()),file=sys.stderr)

# --- VCF ---
print("cargando VCF...",file=sys.stderr)
cal=allel.read_vcf(f"{proj}/results/eub_final.vcf.gz",
                   fields=['samples','calldata/GT','variants/CHROM','variants/POS',
                           'variants/REF','variants/ALT'])
samples=list(cal['samples']); gt=allel.GenotypeArray(cal['calldata/GT'])
CH=cal['variants/CHROM']; POS=cal['variants/POS']
REF=cal['variants/REF']; ALT=cal['variants/ALT'][:,0]
n_snp=len(POS); print("SNPs:",n_snp,file=sys.stderr)

keep=[i for i,s in enumerate(samples) if s in pop]
print("muestras en el alineamiento:",len(keep),file=sys.stderr)

# --- pseudo-haploide: alelo del primer haplotipo; missing -> N ---
hap=gt.to_haplotypes()[:,::2]          # (n_snp, n_sample) 0/1/-1
hap=hap[:,keep]
names=[samples[i] for i in keep]

seqs=[]
for j,s in enumerate(names):
    a=hap[:,j]
    seq=np.where(a==0, REF, np.where(a==1, ALT, "N"))
    seqs.append((s,"".join(seq)))
    if (j+1)%50==0: print(f"  {j+1}/{len(names)}",file=sys.stderr)

fa=f"{OUT}/aln_381.fasta"
with open(fa,"w") as fh:
    for s,seq in seqs: fh.write(f">{s}\n{seq}\n")
print("FASTA:",fa,file=sys.stderr)

# --- sitios constantes: composicion de bases del genoma nuclear menos los SNPs ---
NUC=set(f"CBS12357_Chr{i:02d}_polished" for i in range(1,17))
comp=Counter(); L=0
name=None; buf=[]
def flush():
    global L
    if name in NUC:
        s="".join(buf).upper(); L+=len(s); comp.update(s)
for line in open(f"{proj}/Data/CBS12357.fa"):
    if line.startswith(">"):
        flush(); name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
flush()
print("largo nuclear:",L,file=sys.stderr)

# restar las bases REF de los sitios variables (que ya estan en el alineamiento)
nuc_mask=np.isin(CH,list(NUC))
refsnp=Counter(REF[nuc_mask])
fconst={b: comp[b]-refsnp.get(b,0) for b in "ACGT"}
tot=sum(fconst.values())+n_snp
print("=== sitios constantes (-fconst) ===",file=sys.stderr)
print("  A,C,G,T =",",".join(str(fconst[b]) for b in "ACGT"),file=sys.stderr)
print(f"  suma constantes {sum(fconst.values())} + SNPs {n_snp} = {tot} (nuclear={L}) "
      f"{'OK' if tot==L else 'MISMATCH (revisar: N/ambiguos en la referencia)'}",file=sys.stderr)

with open(f"{OUT}/fconst.txt","w") as fh:
    fh.write(",".join(str(fconst[b]) for b in "ACGT")+"\n")
# archivo de fechas: todas las puntas contemporaneas = 0
with open(f"{OUT}/dates.txt","w") as fh:
    fh.write(f"{len(names)}\n")
    for s in names: fh.write(f"{s}\t0\n")
# tasa fija: 1.67e-10 /bp/gen x 2920 gen/ano
with open(f"{OUT}/rate.txt","w") as fh: fh.write("4.876e-7\n")
# mapa pop para pasos posteriores
with open(f"{OUT}/taxa_pop.tsv","w") as fh:
    for s in names: fh.write(f"{s}\t{pop[s]}\n")
print("OK",file=sys.stderr)
