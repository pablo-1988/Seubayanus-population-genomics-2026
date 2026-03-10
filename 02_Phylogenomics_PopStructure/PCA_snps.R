# =============================================================================
# Principal Component Analysis (PCA) of SNP Data — S. eubayanus
# =============================================================================
# Description:
#   Performs PCA visualization of genome-wide SNP variation across S. eubayanus
#   strains using eigenvectors and eigenvalues computed by PLINK (--pca).
#   Two analyses are produced:
#     1. All strains (including admixed individuals)
#     2. Pure strains only (admixed individuals excluded)
#   Each plot shows PC1 vs PC2, with points color-coded by population assignment.
#
# Input:
#   - cichlids.eigenvec : PLINK PCA eigenvectors file (samples × PCs)
#   - cichlids.eigenval : PLINK PCA eigenvalues file (variance explained per PC)
#   - data_vcf471_popK7.txt : tab-separated metadata with population assignments
#
# Output:
#   - PCA scatter plots (PC1 vs PC2) for all strains and pure strains
#   - data_pca_puras.csv : exported PCA coordinates for pure strains
#
# Usage:
#   Set your working directory to the folder containing the input files before
#   running this script. Input files were generated using PLINK v1.9 --pca flag.
#
# Dependencies:
#   install.packages(c("tidyverse", "ggrepel"))
# =============================================================================

library(tidyverse)
library(ggrepel)


# ===========================================================================
# SECTION 1: PCA of ALL strains
# ===========================================================================

# Read PLINK PCA eigenvectors (each row = one sample, columns = PC values)
pca <- read_table("./cichlids.eigenvec", col_names = FALSE)

# Read eigenvalues (proportion of variance explained by each PC)
eigenval <- scan("./cichlids.eigenval")

# Remove the redundant first column (PLINK outputs sample ID twice)
pca <- pca[,-1]

# Set column names: first column = individual ID, rest = PC1, PC2, ...
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Load population metadata (sample names + population assignments)
data_pop_pca <- read.table("../../data_vcf471_popK7.txt", header = T, sep = "\t")
as.data.frame(data_pop_pca) -> pop_pca

# Merge PCA coordinates with population metadata
pca <- as.tibble(data.frame(pca, pop_pca$pop))

# Calculate percentage of variance explained by each PC
pve <- data.frame(PC = 2:21, pve = eigenval/sum(eigenval)*100)

# Bar chart of variance explained per PC
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
   ylab("Percentage variance explained (%)") + theme_light()

# Cumulative variance explained
cumsum(pve$pve)

# PCA scatter plot — all strains, color-coded by population
ggplot()+ geom_point(data=pca, aes(PC1, PC2, fill =data_pop_ordenado$pop), size= 2, colour="black", pch = 21) +
  coord_equal() + theme_light()+xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_fill_manual(values= c("#FBFCFC","#95A5A6", "#F4D03F","#E74C3C", "#2E86C1","#1ABC9C","#34495E","#ABEBC6")) +
  theme_gray() + labs(fill = "Populations")


# ===========================================================================
# SECTION 2: PCA of PURE strains only (admixed individuals excluded)
# ===========================================================================

# Read eigenvectors for pure strains (same file, subset applied afterwards)
pca_puras <- read_table("./cichlids.eigenvec", col_names = FALSE)
eigenval_puras <- scan("./cichlids.eigenval")

# Remove redundant first column
pca_puras <- pca_puras[,-1]

# Set column names
names(pca_puras)[1] <- "ind"
names(pca_puras)[2:ncol(pca_puras)] <- paste0("PC", 1:(ncol(pca)-1))

# Export PCA coordinates for pure strains to CSV
write.csv(pca_puras, "data_pca_puras.csv", row.names = FALSE)

# Subset metadata to exclude admixed samples (Pop != "Admix")
data_pop_ordenado_PCA_puras <- subset(data_pop_ordenado_PCA, Pop != "Admix")
as.data.frame(data_pop_ordenado_PCA_puras) -> pop_pca_puras

# Merge PCA with filtered metadata (skip first two rows if needed)
pca_puras <- as.tibble(data.frame(pca_puras[-(1:2),], pop_pca_puras$Pop))

# Calculate percentage of variance explained for pure-strain PCA
pve_puras <- data.frame(PC = 2:21, pve = eigenval_puras/sum(eigenval)*100)

# Bar chart of variance explained
ggplot(pve_puras, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained (%)") + theme_light()

# Cumulative variance explained
cumsum(pve_puras$pve)

# PCA scatter plot — pure strains only, color-coded by population
ggplot()+ geom_point(data=pca_puras, aes(PC1, PC2, fill =data_pop_ordenado_PCA_puras$Pop), size= 2, colour="black", pch = 21) +
  coord_equal() + theme_light()+xlab(paste0("PC1 (", signif(pve_puras$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve_puras$pve[2], 3), "%)")) +
  scale_fill_manual(values= c( "#34495E","#95A5A6", "#F4D03F", "#2ECC71","#E74C3C", "#2E86C1", "#1ABC9C","#9B59B6" )) +
  theme_gray() + labs(fill = "Populations")
