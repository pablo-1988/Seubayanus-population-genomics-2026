#!/usr/bin/env bash
# Enriquecimiento de breakpoints de SV (MUM&Co) en regiones TE (anotadas por BLAST).
# Null por permutacion (bedtools shuffle). Concordancia con el flag 'mobile' de MUM&Co.
set -euo pipefail
export PATH="/Users/pablo/.local/share/mamba/envs/anotacion/bin:$PATH"
PROJ="/Volumes/Extreme SSD/Eubayanus-Pop"
TE="$PROJ/work/TE"; MUM="$PROJ/Data/mumandco"; OUT="$PROJ/results/TE_SV"
mkdir -p "$OUT"
GENOME="$TE/genome_nuclear.txt"
TEBED="$TE/TE_annotation_CBS12357.bed"
WIN=250; NPERM=1000

# 1) breakpoints de SV de todas las cepas (ref_start y ref_stop -> 2 breakpoints), +/- WIN
#    columnas MUM&Co: 1 ref_chr 3 ref_start 4 ref_stop 6 SV_type ; agregar _polished
: > "$OUT/sv_breakpoints_all.bed"
for f in "$MUM"/CBS_*_output/CBS_*.SVs_all.tsv; do
  st=$(basename "$f" | sed 's/CBS_//;s/.SVs_all.tsv//')
  awk -v st="$st" -v w=$WIN 'NR>1 && $1 ~ /^CBS12357_Chr/ && $1 !~ /mtDNA/ {
    mob=($6 ~ /mobile/)?"mobile":"novel";
    for(p=3;p<=4;p+=1){ s=$p-w; if(s<0)s=0; print $1"_polished\t"s"\t"$p+w"\t"st"|"$6"|"mob }
  }' "$f" >> "$OUT/sv_breakpoints_all.bed"
done
sort -k1,1 -k2,2n "$OUT/sv_breakpoints_all.bed" -o "$OUT/sv_breakpoints_all.bed"
NB=$(wc -l < "$OUT/sv_breakpoints_all.bed")
echo "breakpoints totales (22 cepas x2): $NB"

# subsets mobile / novel (flag de MUM&Co)
grep -E '[|]mobile$' "$OUT/sv_breakpoints_all.bed" > "$OUT/bp_mobile.bed" || true
grep -E '[|]novel$'  "$OUT/sv_breakpoints_all.bed" > "$OUT/bp_novel.bed"  || true

# 2) OBSERVADO: breakpoints que solapan TE (independiente)
obs=$(bedtools intersect -u -a "$OUT/sv_breakpoints_all.bed" -b "$TEBED" | wc -l | tr -d ' ')
pct=$(awk "BEGIN{printf \"%.1f\",100*$obs/$NB}")
echo "observado (breakpoints en TE): $obs / $NB  (${pct}%)"

# 3) NULL: permutar breakpoints por el genoma, contar solape con TE
echo "corriendo $NPERM permutaciones..."
: > "$OUT/null_counts.txt"
for i in $(seq 1 $NPERM); do
  bedtools shuffle -i "$OUT/sv_breakpoints_all.bed" -g "$GENOME" -chrom \
    | bedtools intersect -u -a - -b "$TEBED" | wc -l
done >> "$OUT/null_counts.txt"

# 4) estadistica (fold + p empirico)
python3 - "$obs" "$NB" "$OUT/null_counts.txt" <<'PY'
import sys,numpy as np
obs=int(sys.argv[1]); nb=int(sys.argv[2])
null=np.loadtxt(sys.argv[3])
mu=null.mean(); sd=null.std()
fold=obs/mu if mu>0 else float('inf')
p=(np.sum(null>=obs)+1)/(len(null)+1)
print(f"\n=== ENRIQUECIMIENTO ===")
print(f"observado: {obs}  ({100*obs/nb:.1f}% de breakpoints)")
print(f"esperado (null): {mu:.1f} +/- {sd:.1f}  ({100*mu/nb:.1f}%)")
print(f"fold enrichment: {fold:.1f}x")
print(f"p empirico (1 cola): {p:.4g}  (n={len(null)} perms)")
PY

# 5) concordancia flag MUM&Co 'mobile' vs anotacion TE independiente
echo ""; echo "=== Concordancia MUM&Co-mobile vs TE-BLAST ==="
for tag in mobile novel; do
  n=$(wc -l < "$OUT/bp_$tag.bed" | tr -d ' ')
  ov=$(bedtools intersect -u -a "$OUT/bp_$tag.bed" -b "$TEBED" | wc -l | tr -d ' ')
  awk -v t=$tag -v n=$n -v o=$ov 'BEGIN{printf "  %-7s breakpoints=%d  en TE=%d (%.1f%%)\n",t,n,o,100*o/n}'
done
echo ""; echo "-> resultados en $OUT/"