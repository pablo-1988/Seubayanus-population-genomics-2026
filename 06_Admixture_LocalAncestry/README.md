# 06 — Admixture, introgression and local ancestry

Scripts added for the 2026 revision. They test and quantify gene flow among
*S. eubayanus* lineages (PA, PB1–PB4) and the admixed sub-populations
(Holarctic/Hol, NoAm, SoAm1–4), using four complementary approaches.

## Dsuite/ — D-statistics and f-branch
Formal tests of introgression with Dsuite (ABBA-BABA / Patterson's D, f4-ratio,
f-branch). `07_dsuite_driver.sh` runs `Dsuite Dtrios` + `Dsuite Fbranch`.
The `35–39_fbranch_*` scripts parse the f-branch matrix, add block-jackknife
Z-scores/SE/p-values (`Fbranch -Z`), build tidy tables and the annotated
heatmaps, and re-run the analysis on population subsets (7-pop without SoAm;
6-pop with Hol; 6-pop with NoAm).

## LOTER/ — local ancestry inference
Parameter-free local-ancestry inference with LOTER, using the five pure lineages
as reference panels and the admixed strains as targets. Includes held-out
validation of the minor-ancestry threshold (`04b`) and per-strain classification
(`04c`); `09_loter_K5.py` is the K=5 run.

## Network_SplitsTree/ — NeighborNet split networks
Reticulate phylogenetic networks (NeighborNet). Split systems are computed in
SplitsTree (headless) and laid out/plotted in R with phangorn/ape. Scripts cover
the full cohort, pure + Hol, and pure + Hol + NoAm + SoAm1.

## Chromosome_painting/ — ancestry mosaics
Chromosome paintings of the per-SNP LOTER assignments, including the
window-smoothing test (per-SNP → 10/25/50 kb majority vote) that shows the mosaic
is not an inference artefact.

> Note: paths inside the scripts point to the local project layout and should be
> adjusted to your environment. Tool versions are listed in the manuscript methods.
