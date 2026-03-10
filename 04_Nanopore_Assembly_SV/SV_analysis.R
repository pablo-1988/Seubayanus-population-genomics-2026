# =============================================================================
# Structural Variant Summary Visualization — S. eubayanus
# =============================================================================
# Description:
#   Visualizes the distribution of structural variant (SV) counts across
#   S. eubayanus strains using multiple plot types:
#     1. Tile heatmaps (strain × SV type) using a continuous color scale
#     2. Jitter plots showing individual SV values per strain
#     3. Dot plots with mean ± SD by strain, color-coded by population group
#        (three display options: error bars, yellow mean, population colors)
#     4. Heatmap (pheatmap) of total SV matrix across all strains
#     5. Deletion-specific dot plot with population color coding
#
#   Populations: PA (Patagonia), PB1–PB5 (Patagonian subgroups),
#   Admixed (Holarctic-admixed: UCD646, UCD650, yHRVM108, CBS12357)
#
# Input:
#   - SV_eubayanus.txt  : tab-separated SV count matrix (strains × SV types)
#                         Generated from MUM&Co output by summarizing per-strain SVs
#   - SV_eubayanus2.txt : secondary SV count matrix (for pheatmap section)
#   - del.txt           : deletion-specific count matrix (strains × populations)
#
# Output:
#   - Multiple interactive plots (or save with pdf()/png() wrappers)
#
# Dependencies:
#   install.packages(c("tidyverse","reshape2","MetBrewer","dplyr","RColorBrewer",
#                      "hrbrthemes","viridis","pheatmap","ggthemes"))
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


# Load SV count matrix (rows = strains, columns = SV types)
total_sv <- read.table("SV_eubayanus.txt", sep = "\t", header = T)

# Reshape to long format for ggplot (variable = SV type, value = count)
melted_sv <- melt(total_sv)
head(melted_sv)


ggplot(data = melted_sv, aes(x=value, y=variable, fill=value)) + geom_point() +
  geom_tile() + scale_fill_viridis(discrete=FALSE)  + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ theme_base()



ggplot(data = melted_sv, aes(x=value, y=X, fill=X)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  geom_jitter(size=2,position=position_jitter(0), alpha= 0.3) + 
  theme_base()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

####################################################################################################
library(ggplot2)
library(dplyr)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108","CBS12357", "AA3", "AA4","QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Set strain order as factor levels for correct x-axis ordering
melted_sv$X <- factor(melted_sv$X, levels = individual_order)

# Calculate mean and standard deviation of SV counts per strain (across SV types)
summary_stats <- melted_sv %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Assign population group to each strain
summary_stats$group <- NA  # Initialize group column

# Map each strain to its population group
summary_stats$group[summary_stats$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Admixed"
summary_stats$group[summary_stats$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats$group[summary_stats$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats$group[summary_stats$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats$group[summary_stats$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats$group[summary_stats$X %in% c("E12", "NR2")] <- "PB4"
summary_stats$group[summary_stats$X %in% c("CO150")] <- "PB5"


# --- Option 1: Dot plot with horizontal error bars (mean ± SD) ---
ggplot(data = melted_sv, aes(y = X, x = value)) +
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) +
  # Promedio (puntos negros) y barras de error horizontales
  geom_point(data = summary_stats, aes(y = X, x = mean_value), inherit.aes = FALSE, color = "black", size = 2) +
  geom_errorbarh(data = summary_stats, aes(y = X, xmin = mean_value - sd_value, xmax = mean_value + sd_value), 
                 inherit.aes = FALSE, height = 0, color = "black") +  # Ajuste de barra horizontal sin altura
  # Personalizar eje x
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  # Tema y ajustes finales
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) +
  labs(
    y = "Individuos",
    x = "Valor",
    title = "Promedios y puntos individuales"
  )


# --- Option 2: Dot plot with yellow mean points (no error bars) ---
ggplot(data = melted_sv, aes(y = X, x = value)) +
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) +
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats, aes(y = X, x = mean_value), 
             inherit.aes = FALSE, color = "black", size = 3, shape = 21, fill = "yellow", stroke = 0.5) +  # Color y borde
  # Personalizar eje x
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  # Tema y ajustes finales
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) +
  labs(
    y = "Individuos",
    x = "Valor",
    title = "Promedios y puntos individuales"
  )


# --- Option 3: Dot plot with mean points color-coded by population group ---

ggplot(data = melted_sv, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c( "Admixed"= "white", "PA" = "#F4D03F", 
                                "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                                "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Personalizar eje x
  scale_x_continuous(breaks = c(0, 100, 200, 300)) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Total"
  )



# =============================================================================
# Total SV count heatmap using pheatmap
# =============================================================================
# Reads the full SV matrix and displays it as a pheatmap with no clustering
color1 <- rev(RColorBrewer::brewer.pal(10, "PiYG"))
color2 <- RColorBrewer::brewer.pal(9, "GnBu")
color3 <- RColorBrewer::brewer.pal(10, "RdBu")


sv <- read.table("SV_eubayanus2.txt", sep = "\t", header = T, row.names = 1)

pheatmap(sv, 
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 20, 
         cellheight = 15, 
         main = "Total Structural Variation", 
         fontsize = 10, 
         color = color2, 
         border_color = "black")






# =============================================================================
# Deletion-specific analysis — dot plot with population color coding
# =============================================================================

del_sv <- read.table("del.txt", sep = "\t", header = T)

del_long <- melt(del_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(del_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108","CBS12357", "AA3", "AA4","QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
del_long$X <- factor(del_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_del <- del_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_del$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_del$group[summary_stats_del$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_del$group[summary_stats_del$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_del$group[summary_stats_del$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_del$group[summary_stats_del$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_del$group[summary_stats_del$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_del$group[summary_stats_del$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_del$group[summary_stats_del$X %in% c("CO150")] <- "PB5"



ggplot(data = del_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_del, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Personalizar eje x
  scale_x_continuous(breaks = c(0, 50, 100, 150)) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Deletion"
  )






####################################################################################################
#######Insertions####################################################################################

ins_sv <- read.table("ins.txt", sep = "\t", header = T)

ins_long <- melt(ins_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(ins_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108","CBS12357", "AA3", "AA4","QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
ins_long$X <- factor(ins_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_ins <- ins_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_ins$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_ins$group[summary_stats$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_ins$group[summary_stats$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_ins$group[summary_stats$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_ins$group[summary_stats$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_ins$group[summary_stats$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_ins$group[summary_stats$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_ins$group[summary_stats$X %in% c("CO150")] <- "PB5"



ggplot(data = ins_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_ins, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Personalizar eje x
  scale_x_continuous(breaks = c(0, 50,100, 200, 300)) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Insertions"
  )


####################################################################################################
#######Duplications####################################################################################

dup_sv <- read.table("dup.txt", sep = "\t", header = T)

dup_long <- melt(dup_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(dup_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108","CBS12357", "AA3", "AA4","QC18" , "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
dup_long$X <- factor(dup_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_dup <- dup_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_dup$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_dup$group[summary_stats$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_dup$group[summary_stats$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_dup$group[summary_stats$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_dup$group[summary_stats$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_dup$group[summary_stats$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_dup$group[summary_stats$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_dup$group[summary_stats$X %in% c("CO150")] <- "PB5"



ggplot(data = dup_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_dup, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Duplications"
  ) + scale_x_continuous(breaks = c(0, 10, 20))



####################################################################################################
#######Contractions####################################################################################

contr_sv <- read.table("contr.txt", sep = "\t", header = T)

contr_long <- melt(contr_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(contr_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108", "CBS12357", "AA3", "AA4", "QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
contr_long$X <- factor(contr_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_contr <- contr_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_contr$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_contr$group[summary_stats_contr$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_contr$group[summary_stats_contr$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_contr$group[summary_stats_contr$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_contr$group[summary_stats_contr$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_contr$group[summary_stats_contr$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_contr$group[summary_stats_contr$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_contr$group[summary_stats_contr$X %in% c("CO150")] <- "PB5"




ggplot(data = contr_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_contr, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Contractions"
  )

####################################################################################################
#######Inversions####################################################################################

inv_sv <- read.table("inv.txt", sep = "\t", header = T)

inv_long <- melt(inv_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(inv_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108", "CBS12357", "AA3", "AA4", "QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1","E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
inv_long$X <- factor(inv_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_inv <- inv_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_inv$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_inv$group[summary_stats$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_inv$group[summary_stats$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_inv$group[summary_stats$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_inv$group[summary_stats$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_inv$group[summary_stats$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_inv$group[summary_stats$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_inv$group[summary_stats$X %in% c("CO150")] <- "PB5"

ggplot(data = inv_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_inv, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Inversions"
  ) + scale_x_continuous(breaks = c(0, 2, 4, 6))

####################################################################################################
#######Translocations####################################################################################

trans_sv <- read.table("trans.txt", sep = "\t", header = T)

trans_long <- melt(trans_sv, id.vars = "X", variable.name = "variable", value.name = "value")
head(trans_long)

# Define strain display order (Holarctic-admixed first, then PA, PB1–PB5)
individual_order <- c(
  "UCD646", "UCD650", "yHRVM108","CBS12357", "AA3", "AA4", "QC18", "CL450.1", 
  "CL1005.1", "CL815.1", "CL467.1", "CL715.1", "CL1111.1", "CL248.1", "CL1101.1", 
  "CL607.1", "CL905.1", "S13-HH", "CL216.1", "E12", "NR2", "CO150"
)

# Asegúrate de que 'X' sea un factor y asignar el orden deseado
trans_long$X <- factor(trans_long$X, levels = individual_order)

# Calcular el promedio y desviación estándar por individuo
summary_stats_trans <- trans_long %>%
  group_by(X) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

# Crear la columna de grupo con los valores correspondientes
summary_stats_trans$group <- NA  # Inicializar columna

# Asignar el grupo a cada individuo según su nombre
summary_stats_trans$group[summary_stats$X %in% c("UCD646", "UCD650", "yHRVM108", "CBS12357")] <- "Holartic"
summary_stats_trans$group[summary_stats$X %in% c("AA3", "AA4", "QC18")] <- "PA"
summary_stats_trans$group[summary_stats$X %in% c("CL450.1", "CL1005.1", "CL815.1", "CL467.1")] <- "PB1"
summary_stats_trans$group[summary_stats$X %in% c("CL715.1", "CL1111.1", "CL248.1", "CL1101.1")] <- "PB2"
summary_stats_trans$group[summary_stats$X %in% c("CL607.1", "CL905.1", "S13-HH", "CL216.1")] <- "PB3"
summary_stats_trans$group[summary_stats$X %in% c("E12", "NR2")] <- "PB4"
summary_stats_trans$group[summary_stats$X %in% c("CO150")] <- "PB5"



ggplot(data = trans_long, aes(y = X, x = value)) + 
  # Puntos individuales alineados en una línea horizontal
  geom_point(color = "gray", size = 2, alpha = 0.7) + 
  # Promedio (puntos coloreados con borde negro)
  geom_point(data = summary_stats_trans, aes(y = X, x = mean_value, fill = group), 
             inherit.aes = FALSE, size = 3, shape = 21, stroke = 0.5) +  # Color y borde
  # Colores por grupo
  scale_fill_manual(values = c("Holartic" = "white", "PA" = "#F4D03F", 
                               "PB1" = "#E74C3C", "PB2" = "#2E86C1", 
                               "PB3" = "#1ABC9C", "PB4" = "#a9cce3", "PB5"= "#34495E")) + 
  # Tema y ajustes finales
  theme_bw() + 
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    panel.grid.major = element_line(color = "gray90"), # Líneas de la cuadrícula más suaves
    panel.grid.minor = element_blank()  # Eliminar cuadrícula menor
  ) + 
  labs(
    y = "",
    x = "varian counts",
    title = "Translocations"
  ) + scale_x_continuous(breaks = c(0, 3, 6,9))
