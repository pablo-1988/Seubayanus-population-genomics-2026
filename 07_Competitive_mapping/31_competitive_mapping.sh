#!/usr/bin/env bash
# Mapeo competitivo de reads Nanopore contra referencia concatenada
# eubayanus (CBS12357) + cerevisiae (S288C R64). Demuestra que las cepas puras
# mapean solo a eubayanus; una hibrida (W34/70) mapea a AMBOS genomas.
# Uso: bash scripts/31_competitive_mapping.sh
set -euo pipefail
PROJ="/Volumes/Extreme SSD/Eubayanus-Pop"
cd "$PROJ"
REF="Data/concat_eub_scer.fa"
READS_DIR="data/raw-reads-nanopore"
OUT="results/competitive_mapping"
BAMDIR="work/hybrid_mapping"
THREADS=8
mkdir -p "$OUT" "$BAMDIR"

SUMMARY="$OUT/mapping_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_reads\tpct_mapped\treads_eub\treads_scer\tpct_reads_eub\tpct_reads_scer\teub_breadth_pct\teub_meandepth\tscer_breadth_pct\tscer_meandepth" > "$SUMMARY"

shopt -s nullglob
FILES=("$READS_DIR"/*.fastq.gz "$READS_DIR"/*.fq.gz "$READS_DIR"/*.fastq "$READS_DIR"/*.fq)
if [ ${#FILES[@]} -eq 0 ]; then echo "No hay reads en $READS_DIR"; exit 1; fi

for fq in "${FILES[@]}"; do
  s=$(basename "$fq"); s=${s%.gz}; s=${s%.fastq}; s=${s%.fq}
  echo ">> $s"
  bam="$BAMDIR/${s}.bam"
  if [ -s "$bam" ] && [ -s "${bam}.bai" ]; then
    echo "   (bam existe, salto mapeo)"
  else
    minimap2 -t "$THREADS" -ax map-ont "$REF" "$fq" 2>"$BAMDIR/${s}.mm2.log" \
      | samtools sort -@ "$THREADS" -o "$bam" -
    samtools index "$bam"
  fi

  # idxstats: reads mapeados (primarios) por contig -> agrupar por prefijo
  samtools idxstats "$bam" > "$OUT/${s}.idxstats.tsv"
  reads_eub=$(awk '$1 ~ /^CBS12357_/ {s+=$3} END{print s+0}' "$OUT/${s}.idxstats.tsv")
  reads_scer=$(awk '$1 ~ /^Scer_/ {s+=$3} END{print s+0}' "$OUT/${s}.idxstats.tsv")

  # flagstat: total y % mapeado
  samtools flagstat "$bam" > "$OUT/${s}.flagstat.txt"
  total=$(awk 'NR==1{print $1}' "$OUT/${s}.flagstat.txt")
  mapped=$(grep -m1 "primary mapped" "$OUT/${s}.flagstat.txt" | awk '{print $1}')
  [ -z "$mapped" ] && mapped=$(grep -m1 " mapped (" "$OUT/${s}.flagstat.txt" | awk '{print $1}')

  # coverage por contig (samtools coverage: $5=covbases $7=meandepth $3=endpos=longitud)
  # breadth = sum(covbases)/sum(len); meandepth ponderado por longitud. Solo cromosomas (sin mtDNA).
  samtools coverage "$bam" > "$OUT/${s}.coverage.tsv"
  eub_breadth=$(awk 'NR>1 && $1 ~ /^CBS12357_Chr/ {cb+=$5; ln+=$3} END{if(ln>0) printf "%.4f", 100*cb/ln; else print 0}' "$OUT/${s}.coverage.tsv")
  scer_breadth=$(awk 'NR>1 && $1 ~ /^Scer_Chr/ {cb+=$5; ln+=$3} END{if(ln>0) printf "%.4f", 100*cb/ln; else print 0}' "$OUT/${s}.coverage.tsv")
  eub_depth=$(awk 'NR>1 && $1 ~ /^CBS12357_Chr/ {dsum+=$7*$3; ln+=$3} END{if(ln>0) printf "%.2f", dsum/ln; else print 0}' "$OUT/${s}.coverage.tsv")
  scer_depth=$(awk 'NR>1 && $1 ~ /^Scer_Chr/ {dsum+=$7*$3; ln+=$3} END{if(ln>0) printf "%.2f", dsum/ln; else print 0}' "$OUT/${s}.coverage.tsv")

  tot_mapped_genomes=$((reads_eub + reads_scer))
  pct_mapped=$(awk -v m="$mapped" -v t="$total" 'BEGIN{if(t>0) printf "%.2f", 100*m/t; else print 0}')
  pe=$(awk -v e="$reads_eub" -v t="$tot_mapped_genomes" 'BEGIN{if(t>0) printf "%.2f", 100*e/t; else print 0}')
  ps=$(awk -v c="$reads_scer" -v t="$tot_mapped_genomes" 'BEGIN{if(t>0) printf "%.2f", 100*c/t; else print 0}')

  echo -e "${s}\t${total}\t${mapped}\t${pct_mapped}\t${reads_eub}\t${reads_scer}\t${pe}\t${ps}\t${eub_breadth}\t${eub_depth}\t${scer_breadth}\t${scer_depth}" >> "$SUMMARY"
done

echo "=== RESUMEN ==="
column -t -s $'\t' "$SUMMARY"
