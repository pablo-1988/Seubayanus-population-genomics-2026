# =============================================================================
# Structural Variant Count Heatmap — S. eubayanus
# =============================================================================
# Description:
#   Visualizes the total number of structural variants (SVs) detected in each
#   S. eubayanus strain using heatmap-style tiles. SVs were detected by MUM&Co
#   using CBS12357 as the reference genome. Two heatmap layouts are produced:
#     1. 2D heatmap (strain × population) — shows which population each strain
#        belongs to and its SV count
#     2. Linear heatmap (1D) — strains ordered by population group and SV count,
#        with a single row of tiles
#   Populations: PA (Patagonia), PB1-PB5 (Patagonian subpopulations), Hol-ADMX
#   (Holarctic-admixed)
#
# Input:
#   - Hardcoded data frame with strain names, population assignments, and total
#     SV counts (derived from MUM&Co summary outputs)
#
# Output:
#   - 2D heatmap and linear heatmap (displayed interactively)
#
# Dependencies:
#   install.packages(c("ggplot2", "viridis", "dplyr"))
# =============================================================================

library(ggplot2)
library(viridis)
library(dplyr)

# Strain-level SV count data
# pop = population assignment; id = strain name; total_sv = number of SVs
df <- data.frame(
  pop = c("PA", "PA", "PB1", "PB2", "PB2", "PB3", "PB2", "PB1", "PB1",
          "PB3", "PB2", "PB1", "PB3", "PB5", "PB4", "PB4", "PA", "PB3", "Hol-ADMX",
          "Hol-ADMX", "Hol-ADMX"),
  id = c("AA3", "AA4", "CL1005", "CL1101", "CL1111", "CL216", "CL248",
         "CL450", "CL467", "CL607", "CL715", "CL815", "CL905", "CO150", "E12",
         "NR2", "QC18", "S13HH", "UCD646", "UCD650", "yHRVM108"),
  total_sv = c(136, 131, 112, 132, 123, 135, 121, 111, 105, 147, 159, 119, 119,
               95, 124, 124, 123, 132, 226, 229, 213)
)

# Order strains by ascending SV count for the x-axis
df$id <- factor(df$id, levels = df$id[order(df$total_sv)])

# --- Heatmap 1: 2D (strain on x-axis, population on y-axis) ---
ggplot(df, aes(x = id, y = pop, fill = total_sv)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Total SVs", option = "plasma", direction = -1) +
  labs(title = "Total SV count heatmap per strain",
       x = "Strain",
       y = "Population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- Heatmap 2: Linear (1D) — ordered by population and SV count ---
# Reorder strains by population group, then by SV count within each group
pop_order <- c("Hol-ADMX", "PA", "PB1", "PB2", "PB3", "PB4", "PB5")
df$pop <- factor(df$pop, levels = pop_order)

df <- df %>%
  arrange(pop, total_sv) %>%
  mutate(id = factor(id, levels = id))

# Linear heatmap (single row) ordered by population then SV count
ggplot(df, aes(x = id, y = 1, fill = total_sv)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Total SVs", option = "plasma", direction = -1) +
  labs(title = "Linear SV count heatmap (ordered by population and count)",
       x = "Strain",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
