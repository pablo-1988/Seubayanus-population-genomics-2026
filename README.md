# Seubayanus-population-genomics-2026

Scripts and analysis code for the paper:

> **"Host-shift adaptation shapes genome architecture in *Saccharomyces eubayanus*"**
> Pablo Villarreal-Díaz, Francisco Cubillos et al. *(Manuscript under review, 2026)*

This repository contains all scripts used for genomic, population-genetic, and phenotypic analyses of *Saccharomyces eubayanus* strains described in the paper. Scripts are organized by analysis section, numbered in the logical order in which the analyses are run.

Sections **01–05** cover the core pipeline (variant calling, population structure, population-genetics statistics, Nanopore assembly / structural variants, and phenotyping). Sections **06–08** were added for the 2026 revision and provide direct, independent evidence of gene flow and genome dynamics: admixture and local-ancestry inference, an assembly-free interspecific-hybridity test, and time-calibrated divergence estimation.

---

## Repository Structure

```
.
├── 01_GATK_variant_calling/          # short-read variant calling (GATK4) + base pop-structure
├── 02_Phylogenomics_PopStructure/    # PCA, ADMIXTURE CV, TreeMix
├── 03_PopGenome_statistics/          # π, Tajima's D, FST, Dxy
├── 04_Nanopore_Assembly_SV/          # ONT assembly (LRSDAY), MUM&Co SVs
│   └── TE_SV_association/             #   + TE–structural-variant enrichment (revision)
├── 05_Phenotyping_Maltose/           # growth phenotyping in maltose
├── 06_Admixture_LocalAncestry/       # Dsuite, LOTER, NeighborNet, chromosome painting (revision)
│   ├── Dsuite/
│   ├── LOTER/
│   ├── Network_SplitsTree/
│   └── Chromosome_painting/
├── 07_Competitive_mapping/           # two-species competitive read mapping (revision)
└── 08_Divergence_times_IQTREE3/      # IQ-TREE 3 + LSD2 dating; net divergence (revision)
```

**Workflow order.** `01` produces the SNP call set that feeds `02` (structure), `03` (diversity/differentiation) and, in the revision, `06` (admixture/local ancestry) and `08` (dating). `04` produces the long-read assemblies and structural variants used by its `TE_SV_association/` sub-analysis; the raw long reads are also used by `07`. `05` links the genomic lineages to maltose-growth phenotypes.

---

## Section 1 — GATK Variant Calling (`01_GATK_variant_calling/`)

**Analysis:** Illumina short-read sequencing data from ~120 *S. eubayanus* strains were aligned to the CBS12357 reference genome and processed through the GATK4 Best Practices pipeline to generate a high-quality VCF file. The pipeline also includes downstream population structure analyses.

| Script | Description |
|--------|-------------|
| `GATK_pipeline.sh` | Complete pipeline: GATK4 variant calling (HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs), VCF filtering, phylogenetic tree (IQTree, GTR+ASC), population structure (STRUCTURE, ADMIXTURE, fineSTRUCTURE, fastStructure), gene flow (TreeMix), and gene orthology (OrthoFinder) |

**Key tools:** GATK4, Picard, vcftools, IQTree2, PLINK, ADMIXTURE, STRUCTURE, fineSTRUCTURE4, fastStructure (Docker), TreeMix, OrthoFinder, gffread

---

## Section 2 — Phylogenomics & Population Structure (`02_Phylogenomics_PopStructure/`)

**Analysis:** Principal component analysis (PCA) of genome-wide SNPs, cross-validation error plot to select the optimal number of ancestral populations (K), and population graph with migration events (TreeMix).

| Script | Description |
|--------|-------------|
| `PCA_snps.R` | PCA visualization using PLINK eigenvectors; generates plots for all strains and pure strains (admixed excluded) |
| `CV_error_plot.R` | Cross-validation error plot from ADMIXTURE output to identify optimal K |
| `pipeline_treemix.sh` | TreeMix pipeline: VCF filtering, LD pruning, format conversion, multi-migration runs (m=0–8), and consensus tree with PHYLIP Consense |

**R packages:** tidyverse, ggrepel, ggplot2
**Key tools:** PLINK, ADMIXTURE, TreeMix, PHYLIP, vcftools

---

## Section 3 — Population Genetics Statistics (`03_PopGenome_statistics/`)

**Analysis:** Genome-wide population genetics statistics including nucleotide diversity (π), between-population diversity (Dxy), Tajima's D, and pairwise FST across all *S. eubayanus* populations.

| Script | Description |
|--------|-------------|
| `PopGenome.R` | Calculates π, Tajima's D, pairwise FST, and Dxy using the PopGenome R package on whole-genome VCF data; populations: PA, NoAm, PB1–PB6, Hol-admixed |
| `Fst.R` | Visualizes pairwise FST matrix as heatmaps (diverging RdBu and sequential Blues color scales) |

**R packages:** PopGenome, tidyverse, dplyr, pheatmap, RColorBrewer
**Input format:** VCF files (one per chromosome or whole genome, bgzipped + tabix-indexed for PopGenome)

---

## Section 4 — Oxford Nanopore Assembly & Structural Variants (`04_Nanopore_Assembly_SV/`)

**Analysis:** A subset of strains were sequenced by Oxford Nanopore Technology (ONT). *De novo* genome assemblies and annotations were produced using the LRSDAY pipeline. Structural variants (SVs) between assembled strains and the reference (CBS12357) were detected with MUM&Co.

### 4a — LRSDAY Genome Assembly Pipeline

The LRSDAY (Long Read Sequencing Data Analysis Workflow) pipeline was followed for genome assembly, polishing, scaffolding, and annotation. Example strain: **UCD646** (Holarctic-admixed).

| Script | Description |
|--------|-------------|
| `LRSDAY.00.Long_Reads_Preprocessing.sh` | Adapter trimming (Porechop) and quality/length filtering (Filtlong) |
| `LRSDAY.00.Nanopore_Reads_Basecalling_and_Demultiplexing.sh` | Basecalling and demultiplexing from FAST5 using Guppy |
| `LRSDAY.00.Prepare_Sc_Reference_Genome.sh` | Prepare reference genome (S288C) for downstream steps |
| `LRSDAY.00.Retrieve_Sample_Illumina_Reads.sh` | Download Illumina short reads from SRA (fastq-dump) |
| `LRSDAY.00.Summary_Report_for_Long_Reads.sh` | QC report (read length, N50, quality) using NanoPlot |
| `LRSDAY.01.Long-read-based_Genome_Assembly.sh` | *De novo* assembly with Flye (or Canu, wtdbg2) |
| `LRSDAY.02.Long-read-based_Assembly_Polishing.sh` | Long-read polishing with Racon + Medaka |
| `LRSDAY.03.Short-read-based_Assembly_Polishing.sh` | Short-read polishing with Pilon (Illumina) |
| `LRSDAY.04.Reference-guided_Assembly_Scaffolding.sh` | Scaffolding against reference genome using Ragout |
| `LRSDAY.05.Centromere_Identity_Profiling.sh` | Centromere identification with Exonerate |
| `LRSDAY.06.Mitochondrial_Genome_Assembly_Improvement.sh` | Mitochondrial genome extraction and circularization |
| `LRSDAY.07.Supervised_Final_Assembly.1.sh` | Generate contig modification template for manual curation |
| `LRSDAY.07.Supervised_Final_Assembly.2.sh` | Apply modifications and produce final assembly |
| `LRSDAY.08.Centromere_Annotation.sh` | Centromere re-annotation in final assembly |
| `LRSDAY.09.Nuclear_Gene_Annotation.sh` | Nuclear gene annotation with MAKER |
| `LRSDAY.11.TE_Annotation.sh` | Transposable element annotation with RepeatMasker |
| `LRSDAY.12.Core_X_Element_Annotation.sh` | Core X element identification (nhmmer, HMM-based) |
| `LRSDAY.13.Y_Prime_Element_Annotation.sh` | Y-prime element identification (BLAT) |
| `LRSDAY.14.Gene_Orthology_Identification.sh` | Gene ortholog assignment using ProteinOrtho |
| `LRSDAY.15.Annotation_Integration.sh` | Integrate all annotations into a final GFF3 file |

**LRSDAY documentation:** https://github.com/yjx1217/LRSDAY
**Key tools:** Guppy, Porechop, Filtlong, NanoPlot, Flye, Racon, Medaka, Pilon, Ragout, Exonerate, MAKER, RepeatMasker, HMMER, BLAT, ProteinOrtho, BEDtools, MUMmer4

### 4b — MUM&Co Structural Variant Analysis

| Script | Description |
|--------|-------------|
| `mumandco_v3.8.sh` | Detects all classes of SVs between a query and reference genome using MUMmer4 whole-genome alignment |
| `MUMandCo_nonredundant_population_dataset.sh` | Creates a non-redundant population-level SV dataset by validating calls across multiple assemblies per strain |
| `SV_analysis.R` | Visualizes SV count distributions per strain (dot plots with mean ± SD, heatmaps, by population group) |
| `sv_count.R` | Heatmap of total SV counts per strain ordered by population group |
| `sv_eubayanus.R` | Comprehensive SV analysis: population-specific SVs, chromosomal ideograms (RIdeogram), gene overlap, and GO enrichment |
| `distribución_de_variantes.R` | Identifies SVs unique to each population group from MUM&Co TSV files; histogram of genomic positions |
| `genetic_distance_vs_SV.R` | Correlation between genetic distance from CBS12357 and total SV count; Spearman/Pearson tests, linear models |

**R packages:** tidyverse, ggplot2, reshape2, pheatmap, viridis, RIdeogram, GenomicRanges, rtracklayer, Biostrings, clusterProfiler, org.Sc.sgd.db, readxl, dplyr
**Key tools:** MUMmer4, MUM&Co v3.8

### 4c — TE–Structural-Variant Association (`04_Nanopore_Assembly_SV/TE_SV_association/`) *(2026 revision)*

**Analysis:** Tests whether MUM&Co structural-variant breakpoints are enriched in transposable-element (TE) regions of CBS12357. TEs were annotated by BLAST of the *S. cerevisiae* Ty/LTR library against CBS12357 (predominantly solo LTRs, ~0.4% of the genome). Enrichment of SV breakpoints in TE is tested against a permutation null, and cross-checked against MUM&Co's own `mobile` flag, overall and per lineage.

| Script | Description |
|--------|-------------|
| `41_te_sv_enrichment.sh` | SV breakpoints vs TE regions; permutation null (`bedtools shuffle`); concordance with the MUM&Co `mobile` flag |
| `42_te_sv_figure.py` | Overall figure: SV composition, breakpoint enrichment, mobile-vs-novel concordance |
| `43_te_sv_per_pop.py` | Per-lineage enrichment and fraction of breakpoints in TE |

**Result:** SV breakpoints are ~23× enriched in TE regions (p < 0.001); the enrichment is present in every lineage (a species-wide mechanism), while the Holarctic lineage carries a higher total SV burden.

**Key tools:** BLAST+, bedtools; Python (numpy, pandas, scipy, matplotlib)

---

## Section 5 — Phenotyping in Maltose (`05_Phenotyping_Maltose/`)

**Analysis:** Growth phenotyping of *S. eubayanus* strains in glucose and maltose media at different sugar concentrations. Lag phase duration and maximum optical density (ODmax) were measured across biological triplicates.

| Script | Description |
|--------|-------------|
| `fenotipos_maltose.R` | Generates annotated heatmaps of lag phase and ODmax across conditions (Glucose2, Maltose2, Maltose5, Maltose10, Maltose20); z-score normalized, hierarchical clustering, color-coded by population |

**R packages:** readxl, dplyr, tidyr, pheatmap, RColorBrewer

---

## Section 6 — Admixture & Local Ancestry (`06_Admixture_LocalAncestry/`) *(2026 revision)*

**Analysis:** Direct, independent evidence of gene flow among the pure lineages (PA, PB1–PB4) and the admixed sub-populations (Holarctic/Hol, NoAm, SoAm1–4), using four complementary approaches. This section replaces the reliance on model-based clustering with genome-wide, method-diverse evidence of introgression and mosaic ancestry.

### 6a — Dsuite: D-statistics and f-branch (`Dsuite/`)

| Script | Description |
|--------|-------------|
| `07_dsuite_driver.sh` | Runs `Dsuite Dtrios` (Patterson's D, f4-ratio, BBAA/ABBA/BABA) and `Dsuite Fbranch` over the population tree |
| `35_fbranch_table.py` | Parses the f-branch matrix; adds block-jackknife Z-scores, SE and p-values (`Fbranch -Z`); tidy tables |
| `36_fbranch_stats_fig.py` | Annotated f-branch heatmap + pairwise comparisons between gene-flow estimates |
| `37_fbranch_7pop.py` | f-branch restricted to pure lineages + Hol + NoAm (no SoAm) |
| `38_fbranch_6pop.py` | f-branch for pure lineages + Hol (no NoAm) |
| `39_fbranch_6pop_NoAm.py` | f-branch for pure lineages + NoAm (no Hol) |

### 6b — LOTER: local ancestry (`LOTER/`)

| Script | Description |
|--------|-------------|
| `04_loter.py` | Parameter-free local-ancestry inference; pure lineages as reference panels, admixed strains as targets |
| `04b_loter_validate.py` | Held-out validation of the minor-ancestry threshold |
| `04c_classify_all.py` | Per-strain admixture classification from LOTER proportions |
| `09_loter_K5.py` | LOTER run under the K = 5 lineage scheme |

### 6c — NeighborNet split networks (`Network_SplitsTree/`)

| Script | Description |
|--------|-------------|
| `05_network.R`, `10_network_K5.R`, `10b_network_all.R` | NeighborNet networks (full cohort / K = 5 scheme) |
| `10c_network_from_splits.R` | Layout/plot from SplitsTree-computed split systems (phangorn/ape) |
| `18_network_pure_Hol.R`, `19_network_pure_Hol_NoAm_SoAm1.R` | Networks for selected taxon subsets |

### 6d — Chromosome painting (`Chromosome_painting/`)

| Script | Description |
|--------|-------------|
| `04d_painting.py`, `09b_painting_K5.py`, `09c_painting_pop.py` | Chromosome paintings of per-SNP LOTER assignments |
| `14_painting_windows.py` | Window-smoothing test (per-SNP → 10/25/50 kb majority vote): the mosaic is not an inference artefact |
| `16_paintings_EN.py` | Publication-quality paintings (vector PDF) for Hol and NoAm |

**Key tools:** Dsuite v0.5, LOTER, ADMIXTOOLS 2 (f2/f3/qpAdm), SplitsTree App v6.5.1; R (phangorn, ape); Python (numpy, pandas, scipy, matplotlib)

---

## Section 7 — Competitive Mapping / Interspecific-Hybridity Test (`07_Competitive_mapping/`) *(2026 revision)*

**Analysis:** An assembly-free test that the strains are *not* interspecific hybrids. Raw Oxford Nanopore reads are mapped competitively against a concatenated two-species reference (*S. eubayanus* CBS12357 + *S. cerevisiae* S288C, R64 / GCF_000146045.2). Breadth of coverage per genome — not the proportion of reads — is the decisive statistic: only a true hybrid (lager control CBS1483) covers both genomes genome-wide, whereas pure/admixed *S. eubayanus* strains cover ~99% of the *S. eubayanus* genome and < 1% of the *S. cerevisiae* genome.

| Script | Description |
|--------|-------------|
| `31_competitive_mapping.sh` | `minimap2 -ax map-ont` to the two-species reference; per-genome breadth and mean depth from `samtools coverage` |
| `34_mapping_figure.py` | Summary table and per-genome coverage figure |

**Key tools:** minimap2, samtools; Python (pandas, matplotlib)

---

## Section 8 — Divergence-Time Estimation (`08_Divergence_times_IQTREE3/`) *(2026 revision)*

**Analysis:** Time-calibrated phylogeny of the lineages with IQ-TREE 3 + LSD2, plus an independent estimate of pairwise split times from Nei's net divergence. Because gene flow violates a strict bifurcating tree, the two estimators are reported together. Branch lengths use explicit constant sites (`-fconst`); dating uses a fixed clock (mutation rate 1.67e-10 /bp/gen × generations-per-year), and node ages are reported as a range over the generations-per-year assumption.

| Script | Description |
|--------|-------------|
| `22_build_alignment.py` | Pseudo-haploid SNP alignment + constant-site (A/C/G/T) counts for `-fconst` |
| `23_net_divergence.py` | Nei's net divergence (d_a) split times; leave-one-chromosome-out block jackknife |
| `24_extract_node_ages.R` | TMRCA / crown ages from the LSD2 time tree |
| `25_dating_compare.py` | LSD2 vs net-divergence comparison |
| `26_population_timetree.R`, `28_timetree_panelC.R` | Population time-tree figures |
| `27_dating_jackknife.py` | Chromosome block jackknife for the dating |

**Key tools:** IQ-TREE 3, LSD2; Python (numpy, pandas, matplotlib); R (ape)

---

## Dependencies Summary

### R packages
```r
install.packages(c("tidyverse", "ggplot2", "ggrepel", "reshape2", "dplyr", "tidyr",
                   "pheatmap", "viridis", "RColorBrewer", "MetBrewer", "hrbrthemes",
                   "ggthemes", "readxl", "readr", "tibble", "devtools", "remotes"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "Biostrings", "IRanges",
                       "org.Sc.sgd.db", "clusterProfiler", "enrichplot", "DOSE"))

# GitHub packages
remotes::install_github("pievos101/PopGenome")
devtools::install_github("zwdzwd/RIdeogram")
```

### Bioinformatics tools
| Tool | Version | Usage |
|------|---------|-------|
| GATK | ≥ 4.1.9 | Variant calling |
| Picard | any | Read group manipulation |
| vcftools | ≥ 0.1.16 | VCF filtering |
| IQTree | ≥ 2.0 | Phylogenetic tree |
| PLINK | ≥ 1.9 | LD pruning |
| ADMIXTURE | ≥ 1.3 | Population structure |
| STRUCTURE | ≥ 2.3 | Bayesian clustering |
| fineSTRUCTURE | ≥ 4.1 | Fine-scale structure |
| fastStructure | 1.0 | Structure (variational Bayes) |
| TreeMix | ≥ 1.13 | Population graph |
| OrthoFinder | ≥ 2.5 | Ortholog identification |
| MUMmer | ≥ 4.0 | Genome alignment |
| MUM&Co | 3.8 | SV detection |
| Guppy | any | ONT basecalling |
| Flye | ≥ 2.8 | Genome assembly |
| Medaka | ≥ 1.4 | Long-read polishing |
| Pilon | ≥ 1.23 | Short-read polishing |
| Ragout | ≥ 2.2 | Scaffolding |
| MAKER | ≥ 2.31 | Gene annotation |
| RepeatMasker | any | TE annotation |
| HMMER | ≥ 3.3 | HMM-based search |
| NanoPlot | any | ONT QC |
| Dsuite | 0.5 (r58) | D-statistics, f4-ratio, f-branch |
| ADMIXTOOLS 2 | 2.0.10 | f2, f3, qpAdm |
| LOTER | — | Local ancestry inference |
| SplitsTree App | 6.5.1 | NeighborNet split networks |
| minimap2 | 2.30 | Competitive two-species read mapping |
| samtools | ≥ 1.22 | BAM processing, per-genome coverage |
| BLAST+ | ≥ 2.13 | TE / MAL homology annotation |
| bedtools | ≥ 2.30 | Interval overlap, permutation null |
| IQ-TREE | 3.1.3 | Maximum-likelihood phylogeny |
| LSD2 | 2.4.4 | Least-squares molecular dating |

---

## Citation

If you use these scripts, please cite:

> Villarreal-Díaz P, Cubillos F et al. (2026). *Host-shift adaptation shapes genome architecture in Saccharomyces eubayanus*. *(Under review)*

---

## Contact

Pablo Villarreal-Díaz — GitHub: [@pablo-1988](https://github.com/pablo-1988)
