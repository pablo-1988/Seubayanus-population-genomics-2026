# Seubayanus-population-genomics-2026

Scripts and analysis code for the paper:

> **"Host-shift adaptation shapes genome architecture in *Saccharomyces eubayanus*"**
> Pablo Villarreal-Díaz, Francisco Cubillos et al. *(Manuscript under review, 2026)*

This repository contains all scripts used for genomic, population-genetic, and phenotypic analyses of *Saccharomyces eubayanus* strains described in the paper. Scripts are organized by analysis section.

---

## Repository Structure

```
.
├── 01_GATK_variant_calling/
├── 02_Phylogenomics_PopStructure/
├── 03_PopGenome_statistics/
├── 04_Nanopore_Assembly_SV/
└── 05_Phenotyping_Maltose/
```

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

---

## Section 5 — Phenotyping in Maltose (`05_Phenotyping_Maltose/`)

**Analysis:** Growth phenotyping of *S. eubayanus* strains in glucose and maltose media at different sugar concentrations. Lag phase duration and maximum optical density (ODmax) were measured across biological triplicates.

| Script | Description |
|--------|-------------|
| `fenotipos_maltose.R` | Generates annotated heatmaps of lag phase and ODmax across conditions (Glucose2, Maltose2, Maltose5, Maltose10, Maltose20); z-score normalized, hierarchical clustering, color-coded by population |

**R packages:** readxl, dplyr, tidyr, pheatmap, RColorBrewer

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

---

## Citation

If you use these scripts, please cite:

> Villarreal-Díaz P, Cubillos F et al. (2026). *Host-shift adaptation shapes genome architecture in Saccharomyces eubayanus*. *(Under review)*

---

## Contact

Pablo Villarreal-Díaz — GitHub: [@pablo-1988](https://github.com/pablo-1988)
