# =============================================================================
# Cross-Validation Error Plot for ADMIXTURE / fastStructure Analysis
# =============================================================================
# Description:
#   Visualizes the cross-validation (CV) error across different values of K
#   (number of ancestral populations) from ADMIXTURE or fastStructure runs.
#   The optimal K is identified as the value with the minimum CV error.
#
# Input:
#   - CV error values obtained from ADMIXTURE output (cv_values vector)
#     These are typically parsed from *.log files using:
#     grep "CV error" *.log | awk '{print $3,$4}' | sed 's/(K=//;s/)//'
#
# Output:
#   - Line plot of CV error vs. K with the optimal K highlighted
#
# Dependencies:
#   install.packages(c("ggplot2", "tidyverse", "dplyr"))
# =============================================================================

library(ggplot2)
library(tidyverse)
library(dplyr)

# CV error values for K = 2 to 10 (obtained from ADMIXTURE log files)
cv_values <- c(0.6636590, 0.4956412, 0.4227702, 0.3892754, 0.3378730, 0.3259265, 0.2799016, 0.3225133, 0.3429151)
K_values <- 2:10  # Range of K values tested

# Plot CV error as a connected line chart
plot(K_values, cv_values, type = "b", col = "#95A5A6", pch = 19, lty = 1,
     xlab = "Number of K clusters", ylab = "Cross-Validation Error",
     main = "Cross-Validation Error by K")

# Identify and highlight the K with the minimum CV error
min_cv_index <- which.min(cv_values)
points(K_values[min_cv_index], cv_values[min_cv_index], col = "#34495E", pch = 19, cex = 1.5)

# Add text label indicating the minimum CV value
text(K_values[min_cv_index], cv_values[min_cv_index], labels = paste("Min CV:", cv_values[min_cv_index]), pos = 3, col = "#34495E")
