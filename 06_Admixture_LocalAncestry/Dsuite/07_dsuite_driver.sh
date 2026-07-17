#!/usr/bin/env bash
# Corre Dsuite (Dtrios + Fbranch + heatmap) para una carpeta (K,umbral) dada.
# Uso: 07_dsuite_driver.sh <carpeta>
set -euo pipefail
source "/Volumes/Extreme SSD/Eubayanus-Pop/scripts/env.sh"
D="$1"; cd "$PROJ"
mkdir -p "$D/02_dsuite"
echo "[tree] $D"
"$RSCRIPT" scripts/07_tree.R "$D" > "$D/02_dsuite/tree.log" 2>&1
cp "$D/lineage_tree.nwk" "$D/02_dsuite/"
echo "[Dtrios] $D"
"$DSUITE" Dtrios "$VCF" "$D/SETS.txt" -t "$D/02_dsuite/lineage_tree.nwk" \
    -o "$D/02_dsuite/dsuite" 2> "$D/02_dsuite/dtrios.log"
echo "[Fbranch] $D"
"$DSUITE" Fbranch "$D/02_dsuite/lineage_tree.nwk" "$D/02_dsuite/dsuite_tree.txt" \
    > "$D/02_dsuite/fbranch.txt" 2> "$D/02_dsuite/fbranch.log"
( cd "$D/02_dsuite" && "$PYENV" "$PROJ/work/Dsuite/utils/dtools.py" fbranch.txt \
    lineage_tree.nwk --outgroup Outgroup >/dev/null 2>&1 ) || true
echo "[done] $D"
