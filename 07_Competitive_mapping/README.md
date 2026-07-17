# 07 — Competitive mapping (interspecific-hybridity test)

Assembly-free test that the strains are *not* interspecific hybrids. Raw Oxford
Nanopore reads are mapped competitively (`minimap2 -ax map-ont`) against a
concatenated two-species reference: *S. eubayanus* CBS12357 + *S. cerevisiae*
S288C (R64 / GCF_000146045.2), with a distinct contig-name prefix per species.

- `31_competitive_mapping.sh` — mapping + per-genome breadth of coverage and mean
  depth (nuclear chromosomes) from `samtools coverage`.
- `34_mapping_figure.py` — summary table and figure.

Breadth (fraction of each genome covered), not the proportion of reads, is the
decisive statistic: only a true hybrid (lager control CBS1483) covers both genomes
genome-wide, whereas pure/admixed *S. eubayanus* strains cover ~99% of the
*S. eubayanus* genome and <1% of the *S. cerevisiae* genome.
