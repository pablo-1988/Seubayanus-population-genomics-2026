# 08 — Divergence-time estimation (IQ-TREE 3 + LSD2)

Time-calibrated phylogeny and independent net-divergence dating of the
*S. eubayanus* lineages.

- `22_build_alignment.py` — pseudo-haploid SNP alignment + constant-site counts
  (`-fconst`) from the reference (critical for correct branch lengths).
- `23_net_divergence.py` — Nei's net divergence (d_a) split times with a
  leave-one-chromosome-out block jackknife.
- `24_extract_node_ages.R` — TMRCA / crown ages from the LSD2 time tree.
- `25_dating_compare.py` — LSD2 vs net-divergence comparison.
- `26_population_timetree.R`, `28_timetree_panelC.R` — population time-tree figures.
- `27_dating_jackknife.py` — jackknife over chromosomes.

ML tree with IQ-TREE 3 (GTR+F+I+G4), dated with LSD2 using a fixed clock
(mutation rate 1.67e-10 /bp/gen x generations-per-year). Node ages are reported
as a range over the generations-per-year assumption; see the manuscript methods.
