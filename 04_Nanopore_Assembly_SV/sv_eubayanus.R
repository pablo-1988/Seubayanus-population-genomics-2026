# =============================================================================
# Chromosomal Ideogram and Gene Overlap Analysis of Structural Variants
# S. eubayanus — Holarctic Population Focus
# =============================================================================
# Description:
#   Comprehensive analysis of structural variants (SVs) detected by MUM&Co in
#   S. eubayanus strains using CBS12357 as the reference genome. This script:
#     1. Loads MUM&Co TSV output files (*.SVs_all.tsv) across all strains
#     2. Maps strains to population groups (Holartic, PA, PB1–PB5)
#     3. Identifies population-specific SVs (found in only one group)
#     4. Builds chromosomal ideograms using RIdeogram with gene density and
#        SV position overlays (focused on Holarctic-admixed strains)
#     5. Finds genes overlapping SV coordinates using GenomicRanges
#     6. Performs GO enrichment analysis on affected genes (org.Sc.sgd.db)
#
# Input:
#   - *.SVs_all.tsv files in the working directory (MUM&Co output per strain)
#     Key columns: ref_chr, ref_start, ref_stop, SV_type, size
#   - CBS12357_polished_20170509.fa : reference genome FASTA (for karyotype)
#   - CBS12357.final.gff3.gz : gene annotation of reference genome
#
# Output:
#   - chromosomes.pdf : chromosomal ideogram with SV positions and gene density
#   - genes_afectados_por_SVs.csv : table of genes overlapping Holarctic SVs
#   - genes_GO_annotation.csv : GO terms for affected genes
#   - df_plot.txt : SV marker coordinates for ideogram (reusable)
#
# Dependencies:
#   install.packages(c("tidyverse","readr","dplyr","devtools","tibble","ggplot2"))
#   devtools::install_github("zwdzwd/RIdeogram")
#   BiocManager::install(c("Biostrings","GenomicRanges","rtracklayer",
#                          "org.Sc.sgd.db","clusterProfiler","enrichplot","DOSE"))
#
# Example strain: UCD646 (Holarctic-admixed, sequenced by Oxford Nanopore)
# =============================================================================

# Install devtools if not already available
install.packages("devtools")

library(tidyverse)
library(readr)
library(dplyr)
library(devtools)
library(RIdeogram)
library(Biostrings)
library(tibble)
library(GenomicRanges)
library(rtracklayer)



# Verify working directory and list available TSV files
getwd()
list.files(pattern = "\\.tsv$")
tsv_files <- list.files(pattern = "\\.tsv$", full.names = TRUE)
print(tsv_files)

df_all <- purrr::map_dfr(tsv_files, function(f) {
  read_tsv(f, col_types = cols()) %>%
    mutate(
      strain = basename(f) %>%
        tools::file_path_sans_ext() %>%
        # quita la parte ".SVs_all" si ese es el sufijo
        sub("\\.SVs_all$", "", .)
    )
})

# Verify data was loaded correctly
glimpse(df_all)
unique(df_all$strain)  # Check that all expected strain names are present

# Define strain-to-population mapping (CBS_* naming from MUM&Co v3.8 output)

strain2group <- tribble(
  ~strain,        ~group,
  "CBS_UCD646",   "Holartic",
  "CBS_UCD650",   "Holartic",
  "CBS_yHRVM108", "Holartic",
  "CBS_AA3",      "PA",
  "CBS_AA4",      "PA",
  "CBS_CL450",    "PB1",
  "CBS_CL1005",   "PB1",
  "CBS_CL815",    "PB1",
  "CBS_CL467",    "PB1",
  "CBS_CL715",    "PB2",
  "CBS_CL1111",   "PB2",
  "CBS_CL248",    "PB2",
  "CBS_CL1101",   "PB2",
  "CBS_CL607",    "PB3",
  "CBS_CL905",    "PB3",
  "CBS_S13HH",    "PB3",
  "CBS_CL216",    "PB3",
  "CBS_E12",      "PB4",
  "CBS_NR2",      "PB4",
  "CBS_CO150",    "PB5"
)

# Join population group to the variant table
df_all <- df_all %>%
  left_join(strain2group, by = "strain")

# Create a unique SV identifier: chromosome_start_stop_type
df_all <- df_all %>%
  mutate(
    sv_id = paste(ref_chr, ref_start, ref_stop, SV_type, sep = "_")
  )

# Recode SV types into broader categories: deletion, insertion, other
df_all <- df_all %>%
  mutate(
    SV_type = case_when(
      SV_type %in% c("deletion_mobile", "deletion_novel")   ~ "deletion",
      SV_type %in% c("insertion_mobile", "insertion_novel") ~ "insertion",
      TRUE                                                   ~ SV_type
    )
  )

# Count how many population groups each SV appears in
sv_group_counts <- df_all %>%
  distinct(sv_id, group) %>%   # ensure each sv_id–group pair counted once
  count(sv_id, name = "n_groups")

# Retain SVs found in only one population group (population-specific)
unique_sv_ids <- sv_group_counts %>% 
  filter(n_groups == 1) %>% 
  pull(sv_id)

df_all <- df_all %>%
  mutate(
    SV_type = case_when(
      SV_type %in% c("deletion_mobile", "deletion_novel")   ~ "deletion",
      SV_type %in% c("insertion_mobile", "insertion_novel") ~ "insertion",
      TRUE                                                   ~ SV_type
    )
  )

df_unique <- df_all %>%
  filter(
    sv_id %in% unique_sv_ids,
    !is.na(group),
    size < 10000
  ) %>%
  distinct(sv_id, group, ref_chr, ref_start, ref_stop, SV_type, size)

# Count population-specific SVs by group and type
df_counts <- df_unique %>% 
  count(group, SV_type, name = "n") %>% 
  arrange(group, desc(n))

print(df_counts)

# Group SV types into broader categories for visualization (already recoded above)
df_unique <- df_unique %>%
  mutate(
    SV_group = case_when(
      SV_type %in% c("deletion_mobile", "deletion_novel")   ~ "deletion",
      SV_type %in% c("insertion_mobile", "insertion_novel") ~ "insertion",
      TRUE                                                   ~ SV_type
    )
  )

# Summary: count population-specific SVs by group and category
df_counts_grouped <- df_unique %>%
  count(group, SV_group, name = "n") %>%
  arrange(group, desc(n))

print(df_counts_grouped)
view(df_counts_grouped)

# =============================================================================
# Holarctic-specific SV analysis and ideogram visualization
# =============================================================================

# Filter to Holarctic-admixed population only (UCD646, UCD650, yHRVM108)
df_holartic <- df_unique %>%
  filter(group == "Holartic")

# Count Holarctic-specific SVs by type
df_holartic %>%
  count(SV_group) %>%
  arrange(desc(n))

# 3. (Opcional) Hacer un histograma de posiciones solo para Holartic
df_holartic %>%
  mutate(midpoint = (ref_start + ref_stop) / 2) %>%
  ggplot(aes(x = midpoint, fill = SV_group)) +
  geom_histogram(bins = 100, colour = "black", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Posición (bp)",
    y     = "Número de SVs únicas",
    fill  = "Tipo SV (agrupado)",
    title = "SVs únicas en el grupo Holartic"
  ) +
  theme_bw()

# ------------------------------------------
# Flujo de trabajo en R para superposición de SVs con genes
# ------------------------------------------
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# 2. Importar archivo GFF3 (anotación del genoma)
gff_file <- "./CBS12357.final.gff3.gz"  # Reemplaza con tu ruta
gff <- import(gff_file)

# Filtrar solo elementos tipo "gene"
genes <- gff[gff$type == "gene"]

# ------------------------------------------
# Filtro SV +- 5000
# ------------------------------------------
df_plot %>% dplyr::select(Chr, Start, End) %>% mutate(df_svs, end_1= Start + 5000, start_2 = End -5000)  -> df_svs

# Convertir a GRanges
svs_gr_2 <- GRanges(
  seqnames = df_svs$Chr,
  ranges = IRanges(start = df_svs$Start, end = df_svs$end_1)
)
# Convertir a GRanges
svs_gr_3 <- GRanges(
  seqnames = df_svs$Chr,
  ranges = IRanges(start = df_svs$start_2, end = df_svs$End))

# 4. Encontrar superposiciones con genes
overlaps_2 <- findOverlaps(svs_gr_2, genes)
genes_affected_2 <- genes[subjectHits(overlaps_2)]
df_genes_affected_2 <- as.data.frame(genes_affected_2)

# 4. Encontrar superposiciones con genes
overlaps_3 <- findOverlaps(svs_gr_3, genes)
genes_affected_3 <- genes[subjectHits(overlaps_3)]
df_genes_affected_3 <- as.data.frame(genes_affected_3)

rbind(df_genes_affected_2, df_genes_affected_3) %>% as_tibble() %>% distinct(Name, .keep_all = T) -> carlos2_df_5000

genes_afected_5000 <- carlos2_df_5000 %>%
  select(seqnames, start, end, width, strand, source, Name)


write_csv(genes_afected_5000, "genes_afected_5000.csv")

# ------------------------------------------
# Filtro SV +- 10000
# ------------------------------------------
df_plot %>% dplyr::select(Chr, Start, End) %>% mutate(df_svs, end_1= Start + 5000, start_2 = End -5000)  -> df_svs

# Convertir a GRanges
svs_gr_2 <- GRanges(
  seqnames = df_svs$Chr,
  ranges = IRanges(start = df_svs$Start, end = df_svs$end_1)
)
# Convertir a GRanges
svs_gr_3 <- GRanges(
  seqnames = df_svs$Chr,
  ranges = IRanges(start = df_svs$start_2, end = df_svs$End))

# 4. Encontrar superposiciones con genes
overlaps_2 <- findOverlaps(svs_gr_2, genes)
genes_affected_2 <- genes[subjectHits(overlaps_2)]
df_genes_affected_2 <- as.data.frame(genes_affected_2)

# 4. Encontrar superposiciones con genes
overlaps_3 <- findOverlaps(svs_gr_3, genes)
genes_affected_3 <- genes[subjectHits(overlaps_3)]
df_genes_affected_3 <- as.data.frame(genes_affected_3)

rbind(df_genes_affected_2, df_genes_affected_3) %>% as_tibble() %>% distinct(Name, .keep_all = T) -> carlos2_df_5000

genes_afected_10000 <- carlos2_df_10000 %>%
  select(seqnames, start, end, width, strand, source, Name)


write_csv(genes_afected_10000, "genes_afected_10000.csv")




##___ SV desde el centromero solo holarticas___##

library(tidyverse)

# Tabla de centrómeros
centros <- tribble(
  ~Chromosome, ~centromere,
  "Chr01", 42487, "Chr02", 480957, "Chr03", 265450, "Chr04", 252556,
  "Chr05", 308165, "Chr06", 98338, "Chr07", 51772, "Chr08", 86433,
  "Chr09", 48405, "Chr10", 654402, "Chr11", 597425, "Chr12", 314739,
  "Chr13", 736102, "Chr14", 129798, "Chr15", 374237, "Chr16", 591038
)

# Preparar datos de SVs
svs <- df_unique %>%
  mutate(
    chrom = str_extract(ref_chr, "Chr\\d+"),
    center = ref_start + (ref_stop - ref_start) / 2
  ) %>%
  left_join(centros, by = c("chrom" = "Chromosome")) %>%
  mutate(
    dist_from_centromere = abs(center - centromere)
  )

# Histograma global
ggplot(svs, aes(x = dist_from_centromere)) +
  geom_histogram(binwidth = 10000, fill = "darkgreen", color = "white") +
  theme_minimal() +
  labs(title = "Distribución global de SVs desde el centrómero",
       x = "Distancia absoluta desde el centrómero (bp)",
       y = "Número de variantes")




ggplot(svs, aes(x = dist_from_centromere)) +
  geom_histogram(binwidth = 10000, fill = "steelblue", color = "white") +
  facet_wrap(~ SV_group, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Distribución de SVs desde el centrómero por tipo de variante",
    x = "Distancia absoluta desde el centrómero (bp)",
    y = "Número de variantes"
  )

ggplot(svs, aes(x = dist_from_centromere, fill = SV_group)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Densidad de SVs desde el centrómero por tipo de variante",
    x = "Distancia absoluta desde el centrómero (bp)",
    y = "Densidad"
  )


#_________________________________Para todas las variantes en todas las poblaciones.

svs_all <- df_all %>%
  mutate(
    chrom = str_extract(ref_chr, "Chr\\d+"),
    center = ref_start + (ref_stop - ref_start) / 2
  ) %>%
  left_join(centros, by = c("chrom" = "Chromosome")) %>%
  mutate(
    dist_from_centromere = abs(center - centromere)
  )

# Histograma global
ggplot(svs_all, aes(x = dist_from_centromere)) +
  geom_histogram(binwidth = 10000, fill = "darkgreen", color = "white") +
  theme_minimal() +
  labs(title = "Distribución global de SVs desde el centrómero",
       x = "Distancia absoluta desde el centrómero (bp)",
       y = "Número de variantes")


#___________________________________________________________________________________________________
#Corregir distancia al centrommetro considerando tamaño de cada uno de los chro.


# Coordenadas de centrómeros (posición central del centrómero en cada cromosoma)
centromeres <- tribble(
  ~chrom, ~centromere,
  "Chr01", 42487, "Chr02", 480957, "Chr03", 265450, "Chr04", 252556,
  "Chr05", 308165, "Chr06", 98338, "Chr07", 51772, "Chr08", 86433,
  "Chr09", 48405, "Chr10", 654402, "Chr11", 597425, "Chr12", 314739,
  "Chr13", 736102, "Chr14", 129798, "Chr15", 374237, "Chr16", 591038
)

# Longitudes totales de cada cromosoma
chrom_lengths <- tribble(
  ~chrom, ~length,
  "Chr01", 230218, "Chr02", 813184, "Chr03", 316620, "Chr04", 1531933,
  "Chr05", 576874, "Chr06", 270161, "Chr07", 1090940, "Chr08", 562643,
  "Chr09", 439888, "Chr10", 745751, "Chr11", 666816, "Chr12", 1078177,
  "Chr13", 924431, "Chr14", 784333, "Chr15", 1091291, "Chr16", 948066
)

# Unir en una sola tabla
genome_info <- left_join(centromeres, chrom_lengths, by = "chrom")





svs <- df_all %>%
  mutate(
    chrom = str_extract(ref_chr, "Chr\\d+"),
    center = ref_start + (ref_stop - ref_start) / 2
  ) %>%
  left_join(genome_info, by = "chrom") %>%
  mutate(
    arm = if_else(center < centromere, "left", "right"),
    arm_start = if_else(arm == "left", 0, centromere),
    arm_end = if_else(arm == "left", centromere, length),
    arm_length = arm_end - arm_start,
    rel_position = (center - arm_start) / arm_length
  ) %>%
  filter(rel_position >= 0 & rel_position <= 1)

ggplot(svs, aes(x = rel_position, fill = SV_type)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Densidad de SVs desde el centrómero",
    x = "Distancia relativa (0 = centrómero, 1 = telómero del brazo)",
    y = "Densidad"
  )

ggplot(svs, aes(x = rel_position)) +
  geom_histogram(binwidth = 0.01, fill = "darkgrey") +
  theme_bw() +
  labs(
    x = "",
    y = "SV count"
  )


library(ggplot2)

ggplot(svs, aes(x = rel_position)) +
  # 1) Densidad con color de relleno y semitransparencia
  geom_density(fill = "#2C3E50", alpha = 0.6, size = 0.8) +
  # 2) Rug plot para marcar cada SV en el eje x
  geom_rug(sides = "b", length = unit(0.02, "npc"), alpha = 0.3) +
  labs( x        = "Posición relativa en el genoma",
    y        = "Densidad",
  ) +
  # 5) Escalas personalizadas
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # 6) Tema limpio con ajustes finos
  theme_bw(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle   = element_text(size = 12, hjust = 0.5),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(color = "gray30"),
    panel.grid.minor = element_blank()
  )




ggplot(svs, aes(x = rel_position)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue") +
  facet_wrap(~ SV_type, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "",
    y = "SV count"
  )


ggplot(svs, aes(x = rel_position, fill = SV_type)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Densidad de variantes estructurales desde el centrómero",
    x = "Distancia relativa (0 = centrómero, 1 = telómero del brazo)",
    y = "Densidad"
  )


svs %>%
  mutate(region = case_when(
    rel_position < 0.25 ~ "pericentromérica",
    rel_position > 0.75 ~ "subtelomérica",
    TRUE ~ "intermedia"
  )) %>%
  count(SV_group, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0)



svs %>%
  mutate(region = case_when(
    rel_position < 0.25 ~ "pericentromérica",
    rel_position > 0.75 ~ "subtelomérica",
  )) %>%
  count(SV_group, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0)




#____________________________________________________________________________________________
# % ins and  del


library(dplyr)
library(stringr)
library(tidyr)

# 1) detectar columna que contiene 'sv' y 'type' (ignore case)
sv_col <- grep("sv[_ ]?type|svtype", names(svs), value = TRUE, ignore.case = TRUE)
if(length(sv_col) == 0) stop("No encontré ninguna columna con 'sv' y 'type' en el nombre. Usa names(svs) para revisar.")
sv_col <- sv_col[1]

# 2) crear columna sv_type de forma robusta y normalizar etiquetas
svs2 <- svs %>%
  mutate(sv_type = as.character(.data[[sv_col]])) %>%       # crea columna nueva a partir del nombre detectado
  mutate(
    sv_type_l = tolower(sv_type),
    sv_type_l = str_replace_all(sv_type_l, "\\s+", " "),
    sv_type_l = str_replace_all(sv_type_l, "[,;|/\\\\]+", ";")  # normaliza separadores a ;
  )

# 3) resumen por fila (cada fila = 1 evento)
summary_byrow <- svs2 %>%
  mutate(sv_type_clean = case_when(
    str_detect(sv_type_l, "ins|insert") ~ "insertion",
    str_detect(sv_type_l, "del|delet")   ~ "deletion",
    TRUE                                 ~ "other"
  )) %>%
  count(sv_type_clean) %>%
  mutate(pct_of_total = round(n / sum(n) * 100, 2))

print(summary_byrow)

# 4) si algunas celdas contienen múltiples tipos (ej. "ins;del"), contar eventos individuales
events <- svs2 %>%
  mutate(sv_type_split = str_split(sv_type_l, ";")) %>%
  unnest(sv_type_split) %>%
  mutate(sv_type_split = str_trim(sv_type_split))

events_summary <- events %>%
  mutate(sv_type_clean = case_when(
    str_detect(sv_type_split, "ins|insert") ~ "insertion",
    str_detect(sv_type_split, "del|delet")   ~ "deletion",
    TRUE                                     ~ "other"
  )) %>%
  count(sv_type_clean) %>%
  mutate(
    pct_of_total_events = round(n / sum(n) * 100, 2)
  )

print(events_summary)

# 5) (opcional) porcentaje relativo solo entre insertion y deletion (entre sí)
insdel <- events_summary %>% filter(sv_type_clean %in% c("insertion", "deletion"))
insdel <- insdel %>% mutate(pct_within_insdel = round(n / sum(n) * 100, 2))
print(insdel)



#____________________________________________________________________________________________
# % ins and  del





