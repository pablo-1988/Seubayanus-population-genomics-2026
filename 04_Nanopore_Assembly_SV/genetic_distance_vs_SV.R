# =============================================================================
# Genetic Distance vs. Structural Variant Count — S. eubayanus
# =============================================================================
# Description:
#   Examines the relationship between genetic distance from the reference genome
#   (CBS12357) and the total number of structural variants (SVs) detected in
#   each S. eubayanus strain. Includes scatter plots with per-population
#   regression lines, Spearman and Pearson correlation tests, and linear models
#   with and without population as a covariate. PB1 is set as the reference
#   population for ANCOVA-style comparisons.
#
# Input:
#   - gen_dis-sv.xlsx   : Excel file with columns: strain, pop, gen_dist, total_sv
#   - gen_dis-sv_2.xlsx : Second dataset with corrected/updated statistics
#
# Output:
#   - Scatter plots with regression lines per population
#   - Correlation test results (Spearman and Pearson) printed to console
#   - Linear model summaries printed to console
#
# Dependencies:
#   install.packages(c("readxl", "ggplot2", "dplyr"))
# =============================================================================

library(readxl)
library(ggplot2)
library(dplyr)

# Load dataset 1: genetic distance vs. total SV count per strain
datos <- read_excel("gen_dis-sv.xlsx")

# Inspect structure and first rows
str(datos)
head(datos)

# Scatter plot: genetic distance (x) vs. total SVs (y), colored by population
# with per-population linear regression lines
p2 <- ggplot(datos, aes(x = gen_dist, y = total_sv, color = pop)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    x = "Genetic distance from reference (CBS12357)",
    y = "Total number of structural variants",
    color = "Population"
  )

p2

# --- Correlation Tests ---

# Spearman rank correlation (non-parametric; robust to outliers)
cor_spear <- cor.test(datos$gen_dist, datos$total_sv, method = "spearman")
cor_spear

# Pearson correlation (parametric; assumes normality)
cor_pear <- cor.test(datos$gen_dist, datos$total_sv, method = "pearson")
cor_pear

# --- Linear Models ---

# Simple linear model: total SVs ~ genetic distance (no population effect)
modelo <- lm(total_sv ~ gen_dist, data = datos)
summary(modelo)

# Linear model including population as an additive covariate
modelo_pop <- lm(total_sv ~ gen_dist + pop, data = datos)
summary(modelo_pop)

# --- Set PB1 as the Reference Population ---
# Redefine population factor so PB1 is the baseline for comparisons

datos$pop <- as.factor(datos$pop)

# Set PB1 as the reference level
datos$pop <- relevel(datos$pop, ref = "PB1")

# Verify: PB1 should appear first in levels
levels(datos$pop)

# Re-fit linear model with PB1 as reference population
modelo_pop_PB1 <- lm(total_sv ~ gen_dist + pop, data = datos)
summary(modelo_pop_PB1)


# ===========================================================================
# DATASET 2: Updated/corrected statistics
# ===========================================================================

datos2 <- read_excel("gen_dis-sv_2.xlsx")

# Inspect dataset 2
str(datos2)
head(datos2)

# Scatter plot for dataset 2 with per-population regression lines
p3 <- ggplot(datos2, aes(x = gen_dist, y = total_sv, color = pop)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    x = "Genetic distance from reference (CBS12357)",
    y = "Total number of structural variants",
    color = "Population"
  )

p3
