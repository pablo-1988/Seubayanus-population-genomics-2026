# TE–SV association

Tests whether structural-variant (SV) breakpoints from MUM&Co are enriched in
transposable-element (TE) regions of CBS12357.

- TE annotation of CBS12357 by BLAST of the *S. cerevisiae* Ty/LTR library
  (built elsewhere) — mostly solo LTRs (~0.4% of the genome).
- `41_te_sv_enrichment.sh` — SV breakpoints vs TE, permutation null
  (`bedtools shuffle`); concordance with MUM&Co's `mobile` flag.
- `42_te_sv_figure.py`, `43_te_sv_per_pop.py` — figures (overall and per lineage).

SV breakpoints are ~23x enriched in TE regions (p<0.001); the enrichment is
present in all lineages (mechanism is species-wide), while the Holarctic lineage
carries a higher total SV burden.
