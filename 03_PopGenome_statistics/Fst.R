# =============================================================================
# Pairwise FST Heatmap — S. eubayanus Population Differentiation
# =============================================================================
# Description:
#   Visualizes the pairwise FST (fixation index) matrix between S. eubayanus
#   populations as a heatmap. FST measures genetic differentiation between
#   populations; values range from 0 (no differentiation) to 1 (complete).
#   Two color schemes are provided: diverging (RdBu) and sequential (Blues).
#
# Input:
#   - fst.txt : tab-delimited pairwise FST matrix with population names as
#               row and column headers (e.g., produced by vcftools --weir-fst-pop)
#
# Output:
#   - Two heatmap figures (displayed interactively or saved with ggsave/pdf)
#
# Dependencies:
#   install.packages(c("tidyverse","reshape2","MetBrewer","dplyr",
#                      "RColorBrewer","hrbrthemes","viridis","pheatmap","ggthemes"))
# =============================================================================

library(tidyverse)
library(reshape2)
library(MetBrewer)
library(dplyr)
library(RColorBrewer)
library(hrbrthemes)
library(viridis)
library(pheatmap)
library(ggthemes)

# Define color palettes
color1 <- rev(RColorBrewer::brewer.pal(10, "PiYG"))   # Pink-Yellow-Green (not used below)
color2 <- rev(RColorBrewer::brewer.pal(5, "GnBu"))    # Green-Blue (not used below)
color3 <- rev(RColorBrewer::brewer.pal(5, "RdBu"))    # Diverging Red-Blue (high vs low FST)

# Load pairwise FST matrix (rows = populations, cols = populations)
fst_eubayanus <- read.table("fst.txt", sep = "\t", header = T, row.names = 1)

# --- Heatmap 1: Diverging Red-Blue color scale ---
# Red = high FST (more differentiated), Blue = low FST (less differentiated)
color3 <- rev(RColorBrewer::brewer.pal(5, "RdBu"))

pheatmap(fst_eubayanus, cluster_rows = F, cluster_cols = F, cellwidth = 20,
         cellheight = 15,
         fontsize = 10, color = color3)

# --- Heatmap 2: Sequential Blues color scale ---
# Darker blue = higher FST (greater differentiation)
color4 <- RColorBrewer::brewer.pal(5, "Blues")

pheatmap(fst_eubayanus,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 20,
         cellheight = 15,
         fontsize = 10,
         color = color4)
