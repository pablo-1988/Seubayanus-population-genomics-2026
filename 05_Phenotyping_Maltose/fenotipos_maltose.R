# =============================================================================
# Maltose Phenotyping Heatmaps — S. eubayanus Growth Kinetics
# =============================================================================
# Description:
#   Analyzes and visualizes growth phenotypes of S. eubayanus strains in
#   glucose and maltose media at different sugar concentrations. Two growth
#   parameters are analyzed:
#     1. Lag phase duration (hours) — time before exponential growth begins
#     2. Maximum optical density (ODmax) — proxy for maximum biomass yield
#   Data are averaged across biological triplicates, z-score normalized per
#   condition (column-wise), and displayed as annotated heatmaps ordered by
#   population group (PA, PB1–PB6, Hol-admix).
#
# Input:
#   - Lag.xlsx   : Excel file with columns: Strain, Pop, Condition, Lag
#   - odmax.xlsx : Excel file with columns: Strain, Pop, Condition, ODmax
#   Conditions tested: Glucose2, Maltose2, Maltose5, Maltose10, Maltose20
#   (numbers indicate % sugar concentration)
#
# Output:
#   - Heatmap 1: Lag phase (z-score normalized) across sugar conditions
#   - Heatmap 2: ODmax (z-score normalized) across sugar conditions
#   Both heatmaps are annotated with population colors on the left side
#
# Dependencies:
#   install.packages(c("readxl", "dplyr", "tidyr", "pheatmap", "RColorBrewer"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# =============================================================================
# PART 1: Lag Phase Heatmap
# =============================================================================

# Load lag phase data (Strain, Pop, Condition, Lag)
df <- read_excel("Lag.xlsx")

# Validate that required columns are present
stopifnot(all(c("Strain","Pop","Condition","Lag") %in% names(df)))

# Average across biological triplicates per strain × condition combination
lag_summary <- df %>%
  group_by(Strain, Pop, Condition) %>%
  summarise(Lag = mean(Lag, na.rm = TRUE), .groups = "drop")

# Define the desired order of sugar conditions on x-axis
cond_levels <- c("Glucose2","Maltose2","Maltose5","Maltose10","Maltose20")

# Create a unique row identifier combining strain and population
lag_summary2 <- lag_summary %>%
  mutate(
    Condition = factor(Condition, levels = cond_levels),
    RowID = paste(Strain, Pop, sep = " | ")
  ) %>%
  arrange(Pop, Strain, Condition)  # Sort by population then strain (alphabetical)

# Pivot to wide format: rows = strains, columns = conditions
mat <- lag_summary2 %>%
  select(RowID, Condition, Lag) %>%
  pivot_wider(
    names_from  = Condition,
    values_from = Lag,
    values_fn   = mean     # Average any remaining duplicates
  ) %>%
  as.data.frame()

rownames(mat) <- mat$RowID
mat$RowID <- NULL

# Ensure all condition columns are present and in the correct order
for (cl in cond_levels) if (!cl %in% colnames(mat)) mat[[cl]] <- NA_real_
mat <- mat[, cond_levels, drop = FALSE]

# Build row annotation: maps each row (strain) to its population group
ann_row <- lag_summary2 %>%
  distinct(RowID, Pop, Strain) %>%
  arrange(Pop, Strain) %>%
  select(RowID, Pop) %>%
  as.data.frame()
rownames(ann_row) <- ann_row$RowID
ann_row$RowID <- NULL

# Reorder the matrix rows to match the annotation order
mat <- mat[rownames(ann_row), , drop = FALSE]

# Row labels: show only the strain name (without the population suffix)
row_labels <- lag_summary2 %>%
  distinct(RowID, Strain) %>%
  .[rownames(mat), "Strain", drop = TRUE]

# Define population color palette
pop_levels <- unique(ann_row$Pop)
named_cols <- c(
  "PA"="#f6e58d","PB1"="#e17055","PB2"="#0984e3","PB3"="#00b894",
  "PB4"="#74b9ff","PB5"="#6c5ce7","PB6"="#fd79a8","Hol-admix"="#636e72","REF"="#2d3436"
)
pop_cols <- setNames(
  named_cols[intersect(names(named_cols), pop_levels)],
  intersect(names(named_cols), pop_levels)
)
# Assign fallback colors to any populations not in the predefined palette
if(length(setdiff(pop_levels, names(pop_cols)))>0){
  extra <- setdiff(pop_levels, names(pop_cols))
  pop_cols <- c(pop_cols, setNames(grDevices::rainbow(length(extra)), extra))
}
ann_colors <- list(Pop = pop_cols)

# Define heatmap color gradient (YlGnBu: yellow = low lag, blue = high lag)
pal <- brewer.pal(9, "YlGnBu")
hm_cols <- colorRampPalette(pal)(101)

# Plot lag phase heatmap with z-score normalization per condition (column)
pheatmap(
  mat,
  cellwidth = 15,
  cellheight = 20,
  cluster_rows = TRUE,      # Cluster strains by lag phase pattern
  cluster_cols = FALSE,     # Keep conditions in the defined order
  border_color = NA,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  color = hm_cols,
  main = "Lag phase by condition (z-score per column)",
  fontsize = 10,
  angle_col = 45,
  na_col = "grey90",
  scale = "row"             # z-score normalization applied per row (strain)
)

# =============================================================================
# PART 2: Maximum Optical Density (ODmax) Heatmap
# =============================================================================

# Load ODmax data (Strain, Pop, Condition, ODmax)
df_odmax <- read_excel("odmax.xlsx")

# Validate required columns
stopifnot(all(c("Strain","Pop","Condition","ODmax") %in% names(df_odmax)))

# Average across biological triplicates
lag_summary_od <- df_odmax %>%
  group_by(Strain, Pop, Condition) %>%
  summarise(ODmax = mean(ODmax, na.rm = TRUE), .groups = "drop")

# Define condition order (same as lag phase analysis)
cond_levels <- c("Glucose2","Maltose2","Maltose5","Maltose10","Maltose20")

# Create unique row identifier and sort by population then strain
lag_summary2_od <- lag_summary_od %>%
  mutate(
    Condition = factor(Condition, levels = cond_levels),
    RowID = paste(Strain, Pop, sep = " | ")
  ) %>%
  arrange(Pop, Strain, Condition)

# Pivot to wide format
mat_od <- lag_summary2_od %>%
  select(RowID, Condition, ODmax) %>%
  pivot_wider(
    names_from  = Condition,
    values_from = ODmax,
    values_fn   = mean
  ) %>%
  as.data.frame()

rownames(mat_od) <- mat_od$RowID
mat_od$RowID <- NULL

# Ensure all condition columns are present and ordered correctly
for (cl in cond_levels) if (!cl %in% colnames(mat_od)) mat_od[[cl]] <- NA_real_
mat_od <- mat_od[, cond_levels, drop = FALSE]

# Build row annotation for ODmax heatmap
ann_row_od <- lag_summary2_od %>%
  distinct(RowID, Pop, Strain) %>%
  arrange(Pop, Strain) %>%
  select(RowID, Pop) %>%
  as.data.frame()
rownames(ann_row_od) <- ann_row_od$RowID
ann_row_od$RowID <- NULL

# Reorder ODmax matrix rows
mat_od <- mat_od[rownames(ann_row_od), , drop = FALSE]

# Row labels (strain name only)
row_labels <- lag_summary2_od %>%
  distinct(RowID, Strain) %>%
  .[rownames(mat_od), "Strain", drop = TRUE]

# Reuse the same population color palette
pop_levels <- unique(ann_row_od$Pop)
named_cols <- c(
  "PA"="#f6e58d","PB1"="#e17055","PB2"="#0984e3","PB3"="#00b894",
  "PB4"="#74b9ff","PB5"="#6c5ce7","PB6"="#fd79a8","Hol-admix"="#636e72","REF"="#2d3436"
)
pop_cols <- setNames(
  named_cols[intersect(names(named_cols), pop_levels)],
  intersect(names(named_cols), pop_levels)
)
if(length(setdiff(pop_levels, names(pop_cols)))>0){
  extra <- setdiff(pop_levels, names(pop_cols))
  pop_cols <- c(pop_cols, setNames(grDevices::rainbow(length(extra)), extra))
}
ann_colors <- list(Pop = pop_cols)

# Reuse the same color gradient
pal <- brewer.pal(9, "YlGnBu")
hm_cols <- colorRampPalette(pal)(101)

# Plot ODmax heatmap with z-score normalization per condition (row)
pheatmap(
  mat_od,
  cellwidth = 15,
  cellheight = 20,
  cluster_rows = TRUE,      # Cluster strains by ODmax pattern
  cluster_cols = FALSE,     # Keep conditions in the defined order
  border_color = NA,
  annotation_row = ann_row_od,
  annotation_colors = ann_colors,
  color = hm_cols,
  main = "Maximum Optical Density (ODmax) by condition",
  fontsize = 10,
  angle_col = 45,
  na_col = "grey90",
  cutree_rows = 4,          # Cut dendrogram into 4 groups
  scale = "row"             # z-score normalization applied per row (strain)
)
