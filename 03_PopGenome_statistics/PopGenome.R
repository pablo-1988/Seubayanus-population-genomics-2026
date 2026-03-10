# =============================================================================
# Population Genetics Statistics — S. eubayanus
# Using the PopGenome R Package
# =============================================================================
# Description:
#   Calculates a comprehensive set of population genetics statistics across
#   S. eubayanus populations using whole-genome SNP data in VCF format.
#   Analyses include:
#     - Nucleotide diversity (pi) per population
#     - Tajima's D neutrality test per population
#     - Pairwise FST (fixation index) between populations
#     - Between-population nucleotide diversity (Dxy)
#   Results are computed genome-wide and/or in sliding windows.
#
# Input:
#   - ./whole/    : directory containing VCF files (one per chromosome or whole genome)
#                   VCFs must be bgzipped and tabix-indexed for PopGenome
#   - popstacks.txt : tab-separated file with columns: Ind, Pop
#                     Defines the population assignment for each sample
#                     (admixed samples excluded from analysis)
#
# Output:
#   - Nucleotide diversity, Tajima's D, FST values (printed to console or
#     exported to data frames for downstream visualization in R)
#
# Dependencies:
#   install.packages("remotes")
#   remotes::install_github("pievos101/PopGenome")
#   install.packages(c("tidyverse", "dplyr"))
#
# Note:
#   PopGenome may require older versions of R (≤ 4.2) or specific Bioconductor
#   packages. If installation fails, try the modified version:
#   remotes::install_github("jrxFive/PopGenome@modified-version")
#
# Populations analyzed:
#   PA, NoAm (North American), PB1, PB2, PB3, PB4, PB5, PB6, Hol-ADM (admixed)
#   Genome size: ~12.5 Mb | SNPs analyzed: ~574,685
# =============================================================================

# Installation options (try in order if previous fails)
install.packages("devtools")
devtools::install_github("jrxFive/PopGenome@modified-version")
remotes::install_github("jrxFive/PopGenome@modified-version")
install.packages(c("bit", "bit64", "ff"))
install.packages("PopGenome_2.7.5.tar.gz", repos = NULL, type = "source")
install.packages("~/Downloads/PopGenome_2.7.5.tar.gz", repos = NULL, type = "source")
install.packages("remotes")
remotes::install_github("pievos101/PopGenome")
library(PopGenome)
library(tidyverse)
library(dplyr)

# Load population metadata (sample IDs + population assignments)
pop_data <- read.table("./popstacks.txt", sep = "\t", header = T)
pop_data

# Exclude admixed samples from population genetics analyses
pop_data <- subset(pop_data, Pop != "Admix")
pop_data

# Load VCF data into PopGenome GENOME.class object
# The 'whole' directory should contain VCF files (one per chromosome or a single whole-genome VCF)
GENOME.class <- readData("./whole", format = "VCF")

# Extract sample lists per population from metadata
pop_data %>% filter(Pop == "PA") %>% pull(Ind) %>% as.character() -> PA
PA ##pop1

pop_data %>% filter(Pop == "NoAm") %>% pull(Ind) %>% as.character() -> NoAm
NoAm ##pop2

pop_data %>% filter(Pop == "PB1") %>% pull(Ind) %>% as.character()-> PB1
PB1 ##pop3

pop_data %>% filter(Pop == "PB2") %>% pull(Ind) %>% as.character() -> PB2
PB2 ##pop4

pop_data %>% filter(Pop == "PB3") %>% pull(Ind) %>% as.character()-> PB3
PB3 ##pop5

pop_data %>% filter(Pop == "PB4") %>% pull(Ind) %>% as.character() -> PB4
PB4 ##pop6

pop_data %>% filter(Pop == "PB5") %>% pull(Ind) %>% as.character()-> PB5
PB5 ##pop7

pop_data %>% filter(Pop == "PB6") %>% pull(Ind) %>% as.character()-> PB6
PB6 ##pop8

pop_data %>% filter(Pop == "Hol-ADM") %>% pull(Ind) %>% as.character()-> HOLadmix
HOLadmix ##pop9



# Assign populations to the GENOME.class object (diploid organisms)
GENOME.class <-set.populations( GENOME.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix), diploid = T)

# Verify population assignments
GENOME.class@populations

# Set genome parameters for downstream calculations
genome_size <- 12500000  # Total genome size in bp (~12.5 Mb for S. eubayanus)
GENOME.class@n.sites2    # Check number of analyzed sites
snp_number <- 574685     # Total number of SNPs in the dataset

# Pairwise 
GENOME.class<- F_ST.stats(GENOME.class, mode="nucleotide")
pairwise.FST <- t(GENOME.class@nuc.F_ST.pairwise)
head(pairwise.FST)

# Diversity within 
GENOME.class <- diversity.stats(GENOME.class)
nucdiv <- GENOME.class@nuc.diversity.within/genome_size
head(nucdiv)

#Neutrality statistics
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
get.neutrality(GENOME.class)[[1]]
GENOME.class@Tajima.D

#______________________________________________________________________________
# SLIDE WINDOW

WHOLE = concatenate.regions(GENOME.class)
slide <- sliding.window.transform(WHOLE,10000, 10000, type=2)
length(slide@region.names)


# Parametros de la ventana
slide <- sliding.window.transform(GENOME.class,10000, 10000, type=1)
length(slide@region.names) # Confirmar numero de ventanas

# Construir data frame con los datos de la ventana
genome.pos <- sapply(slide@region.names, function(x){
  split <- strsplit(x," ")[[1]][c(1,3)] 
  val <- mean(as.numeric(split))
  return(val)
})
genome.pos <- as_tibble(data.frame(genome.pos))

# Pi estadistico
slide <- diversity.stats(slide, pi = TRUE)
# Fst estadistico
slide <- F_ST.stats(slide, mode = "nucleotide")

#Extraer datos
# Diversidad nucleotidica
nd <- slide@nuc.diversity.within/1000
pops <- c("PA","NoAm", "PB1", "PB2", "PB3", "PB4", "PB5", "PB6", "HOLadmix")
colnames(nd) <- paste0(pops, "_pi")
# Fst
fst <- t(slide@nuc.F_ST.pairwise)
x <- colnames(fst)
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)
x <- sub("pop4", pops[4], x)
x <- sub("pop5", pops[5], x)
x <- sub("pop6", pops[6], x)
x <- sub("pop7", pops[7], x)
x <- sub("pop8", pops[8], x)
x <- sub("pop9", pops[9], x)
x <- sub("/", "_", x)
x
# Dxy pairwise diversidad nucleotidica abosluta
dxy <- get.diversity(slide, between = T) [[2]]/500
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

# Unir todos ls data set
slide_data <- as_tibble(data.frame(nd, fst, dxy))

# Visulizacion

# Promedio de cada estadistico
slide_data %>% select(contains("pi")) %>% summarise_all(mean) -> pi_promedio
slide_data %>% select(contains("fst")) %>% summarise_all(mean) -> fst_promedio
slide_data %>% select(contains("dxy")) %>% summarise_all(mean) -> dxy_promedio

# distribucion entre poblaciones pi
pi_g <- slide_data %>% select(contains("pi")) %>% gather(key = "Pop", 
                                                         value = "pi")
pi_acrossGenome <- ggplot(pi_g, aes(Pop, pi, fill = Pop)) + geom_boxplot() + 
  theme_gray() + xlab(NULL) + ggtitle("Diversidad nucleotidica (Pi) Whole Genome")+
  scale_fill_manual(values= c( "#95A5A6", "#F4D03F","#E74C3C", "#2E86C1", "#1ABC9C","#9B59B6","#34495E", "#ABEBC6", "#D6EAF8"))+
  theme_gray() + labs(fill = "Populations")
pi_acrossGenome





# distribucion entre poblaciones fst
fst_g <- slide_data %>% select(contains("fst")) %>% gather(key = "Pop", 
                                                           value = "fst")
fst_acrossGenome <- ggplot(fst_g, aes(Pop, fst)) + geom_boxplot() + 
  theme_light() + xlab(NULL) + ggtitle("Divergencia entre poblaciones (Fst) Whole Genome")
fst_acrossGenome
# distribucion entre poblaciones dxy
dxy_g <- slide_data %>% select(contains("dxy")) %>% gather(key = "Pop", 
                                                           value = "dxy")
dxy_acrossGenome <- ggplot(dxy_g, aes(Pop, dxy)) + geom_boxplot() + 
  theme_light() + xlab(NULL) + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, size = 10))+ 
  ggtitle("Diversidad absoluta (Dxy) Whole Genome")
dxy_acrossGenome

# Distribucion entre poblaciones (un grafico por cada comparacion) para Fst
Fst_dist <- ggplot(slide_data, aes(y = "", x = PA1_PA2_fst)) + 
  geom_line(colour = "red")
Fst_dist <- Fst_dist + xlab("Position (Mpb)") + ylab(expression(italic(F)[ST]))
Fst_dist + theme_light()




# Seleccionar data de inter?ss
PB2HOLadmix <- slide_data %>% select(PB2_HOLadmix_dxy)
mutate(PB2HOLadmix, ventana = 1:41) -> PB2HOLadmix

ggplot(PB2HOLadmix, aes(x = ventana, y= PB2_HOLadmix_dxy)) + geom_point()

##########################################################################################################

GENOME.class <-set.populations(GENOME.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6), diploid = T)

# split the data in 1kb consecutive windows
slide <- sliding.window.transform(GENOME.class,500,500, type=2, whole.data = TRUE)
# total number of windows
length(slide@region.names)

# Statistics
slide <- diversity.stats(slide)
nucdiv <- slide@nuc.diversity.within
# the values have to be normalized by the number of nucleotides in each window
nucdiv <- nucdiv/28
head(nucdiv)


library(ggplot2)
library(ggthemes)

# Supongamos que nucdiv es una matriz donde las filas representan las ventanas y las columnas las poblaciones
# Por ejemplo, si tienes 8 poblaciones, las columnas serán de 1 a 8

# Convertir nucdiv a un dataframe
nucdiv_df <- data.frame(nucdiv)

# Añadir una columna para el número de ventana
nucdiv_df$Window <- 1:nrow(nucdiv_df)

# Derretir el dataframe para que sea compatible con ggplot2
library(reshape2)
nucdiv_melted <- melt(nucdiv_df, id.vars = "Window", variable.name = "Population", value.name = "Nucleotide_Diversity")

# Graficar la diversidad nucleotídica por ventana para cada población

ggplot(nucdiv_melted, aes(x = Window, y = Nucleotide_Diversity, color = Population)) +
  geom_point() + scale_color_manual(values= c("#F4D03F","#95A5A6","#E74C3C", "#2E86C1", "#1ABC9C","#9B59B6","#34495E" ,"#ABEBC6"))+
  labs(x = "Window", y = "Nucleotide Diversity", color = "Population") +
  theme_base()


##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
VCF_split_into_scaffolds("whole/snp_data.recode.vcf", "vcf")


GENOME.class <- readData("./vcf/chr2/", format = "VCF")
GENOME.class@region.names

pop_data %>% filter(Ind == "UCD646") %>% pull(Ind) %>% as.character() -> UCD646
pop_data %>% filter(Ind == "UCD650") %>% pull(Ind) %>% as.character() -> UCD650
pop_data %>% filter(Ind == "yHRVM108") %>% pull(Ind) %>% as.character() -> yHRVM108
pop_data %>% filter(Ind == "yHRVM107") %>% pull(Ind) %>% as.character() -> yHRVM107
pop_data %>% filter(Ind == "yHCT81") %>% pull(Ind) %>% as.character() -> yHCT81
pop_data %>% filter(Ind == "yHAB47") %>% pull(Ind) %>% as.character() -> yHAB47
pop_data %>% filter(Ind == "CDFM21L") %>% pull(Ind) %>% as.character() -> CDFM21L
pop_data %>% filter(Ind == "ABFM5L") %>% pull(Ind) %>% as.character() -> ABFM5L

GENOME.class <-set.populations( GENOME.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw<-sliding.window.transform(GENOME.class,100,100,type=2, whole.data = T)

length(sw@region.names) 

sw <-diversity.stats.between(sw,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("Pop Hol Admix")
   

##UCD646
sw_ucd <-diversity.stats.between(sw,new.populations = list(UCD646,PB2,PB3,PB5))

as_tibble(sw_ucd@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_ucd@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("UCD645")
  
##UCD659
sw_ucd650<-diversity.stats.between(sw,new.populations = list(UCD650,PB2,PB3,PB5))

as_tibble(sw_ucd650@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_ucd650@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("UCD650")

##CDFM21L
sw_CDFM21L<-diversity.stats.between(sw,new.populations = list(CDFM21L,PB2,PB3,PB5))

as_tibble(sw_CDFM21L@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_CDFM21L@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("CDFM21L")


##ABFM5L
sw_ABFM5L<-diversity.stats.between(sw,new.populations = list(ABFM5L,PB2,PB3,PB5))

as_tibble(sw_ABFM5L@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_ABFM5L@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("ABFM5L")

#yHRVM107
sw_yHRVM107<-diversity.stats.between(sw,new.populations = list(yHRVM107,PB2,PB3,PB5))

as_tibble(sw_yHRVM107@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_yHRVM107@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("yHRVM107")


#yHRVM108
sw_yHRVM108<-diversity.stats.between(sw,new.populations = list(yHRVM108,PB2,PB3,PB5))

as_tibble(sw_yHRVM108@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_yHRVM108@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("yHRVM108")


#yHCT81
sw_yHCT81<-diversity.stats.between(sw,new.populations = list(yHCT81,PB2,PB3,PB5))

as_tibble(sw_yHCT81@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_yHCT81@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("yHCT81")



#yHAB47
sw_yHAB47<-diversity.stats.between(sw,new.populations = list(yHAB47,PB2,PB3,PB5))

as_tibble(sw_yHAB47@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/10, `pop1/pop3` = `pop1/pop3`/10, `pop1/pop4`=`pop1/pop4`/10, x = 1:nrow(sw_yHAB47@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"),labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + ggtitle("yHAB47")


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

#CHR1

CHR1.class <- readData("./vcf/chr1/", format = "VCF")
CHR1.class@region.names

CHR1.class <-set.populations(CHR1.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr1<-sliding.window.transform(CHR1.class,100,100,type=2, whole.data = T)
length(sw_chr1@region.names)

sw_chr1 <-diversity.stats.between(sw_chr1,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr1@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr1@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR1") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR2

CHR2.class <- readData("./vcf/chr2/", format = "VCF")
CHR2.class@region.names

CHR2.class <-set.populations(CHR2.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr2<-sliding.window.transform(CHR2.class,100,100,type=2, whole.data = T)
length(sw_chr2@region.names) 

sw_chr2 <-diversity.stats.between(sw_chr2,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr2@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr2@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR2") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))



#CHR3

CHR3.class <- readData("./vcf/chr3/", format = "VCF")
CHR3.class@region.names

CHR3.class <-set.populations(CHR3.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr3<-sliding.window.transform(CHR3.class,100,100,type=2, whole.data = T)
length(sw_chr3@region.names) 

sw_chr3 <-diversity.stats.between(sw_chr3,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr3@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr3@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR3") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))



#CHR4

CHR4.class <- readData("./vcf/chr4/", format = "VCF")
CHR4.class@region.names

CHR4.class <-set.populations(CHR4.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr4<-sliding.window.transform(CHR4.class,100,100,type=2, whole.data = T)
length(sw_chr4@region.names) 

sw_chr4 <-diversity.stats.between(sw_chr4,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr4@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr4@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR4") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR5

CHR5.class <- readData("./vcf/chr5/", format = "VCF")
CHR5.class@region.names

CHR5.class <-set.populations(CHR5.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr5<-sliding.window.transform(CHR5.class,100,100,type=2, whole.data = T)
length(sw_chr5@region.names) 

sw_chr5 <-diversity.stats.between(sw_chr5,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr5@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr5@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR5") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR6

CHR6.class <- readData("./vcf/chr6/", format = "VCF")
CHR6.class@region.names

CHR6.class <-set.populations(CHR6.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr6<-sliding.window.transform(CHR6.class,100,100,type=2, whole.data = T)
length(sw_chr6@region.names) 

sw_chr6 <-diversity.stats.between(sw_chr6,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr6@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr6@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR6") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR7

CHR7.class <- readData("./vcf/chr7/", format = "VCF")
CHR7.class@region.names

CHR7.class <-set.populations(CHR7.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr7<-sliding.window.transform(CHR7.class,100,100,type=2, whole.data = T)
length(sw_chr7@region.names) 

sw_chr7 <-diversity.stats.between(sw_chr7,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr7@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr7@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR7") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR8

CHR8.class <- readData("./vcf/chr8/", format = "VCF")
CHR8.class@region.names

CHR8.class <-set.populations(CHR8.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr8<-sliding.window.transform(CHR8.class,100,100,type=2, whole.data = T)
length(sw_chr8@region.names) 

sw_chr8 <-diversity.stats.between(sw_chr4,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr8@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr8@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR8") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR9

CHR9.class <- readData("./vcf/chr9/", format = "VCF")
CHR9.class@region.names

CHR9.class <-set.populations(CHR9.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr9<-sliding.window.transform(CHR9.class,100,100,type=2, whole.data = T)
length(sw_chr9@region.names) 

sw_chr9 <-diversity.stats.between(sw_chr9,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr9@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr9@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR9") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR10

CHR10.class <- readData("./vcf/chr10/", format = "VCF")
CHR10.class@region.names

CHR10.class <-set.populations(CHR4.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr10<-sliding.window.transform(CHR10.class,100,100,type=2, whole.data = T)
length(sw_chr10@region.names) 

sw_chr10 <-diversity.stats.between(sw_chr10,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr10@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr10@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR10") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR11

CHR11.class <- readData("./vcf/chr11/", format = "VCF")
CHR11.class@region.names

CHR11.class <-set.populations(CHR11.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr11<-sliding.window.transform(CHR11.class,100,100,type=2, whole.data = T)
length(sw_chr11@region.names) 

sw_chr11 <-diversity.stats.between(sw_chr11,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr11@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr11@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR11") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR12

CHR12.class <- readData("./vcf/chr12/", format = "VCF")
CHR12.class@region.names

CHR12.class <-set.populations(CHR4.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr12<-sliding.window.transform(CHR4.class,100,100,type=2, whole.data = T)
length(sw_chr12@region.names) 

sw_chr12 <-diversity.stats.between(sw_chr12,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr12@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr12@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR12") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR13

CHR13.class <- readData("./vcf/chr13/", format = "VCF")
CHR13.class@region.names

CHR13.class <-set.populations(CHR13.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr13<-sliding.window.transform(CHR13.class,100,100,type=2, whole.data = T)
length(sw_chr13@region.names) 

sw_chr13 <-diversity.stats.between(sw_chr13,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr13@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr13@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR13") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))



#CHR14

CHR14.class <- readData("./vcf/chr14/", format = "VCF")
CHR14.class@region.names

CHR14.class <-set.populations(CHR14.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr14<-sliding.window.transform(CHR14.class,100,100,type=2, whole.data = T)
length(sw_chr14@region.names) 

sw_chr14 <-diversity.stats.between(sw_chr14,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr14@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr14@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR14") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR15

CHR15.class <- readData("./vcf/chr15/", format = "VCF")
CHR15.class@region.names

CHR15.class <-set.populations(CHR4.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr15<-sliding.window.transform(CHR15.class,100,100,type=2, whole.data = T)
length(sw_chr15@region.names) 

sw_chr15 <-diversity.stats.between(sw_chr15,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr15@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr15@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR15") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#CHR16

CHR16.class <- readData("./vcf/chr16/", format = "VCF")
CHR16.class@region.names

CHR16.class <-set.populations(CHR16.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix, UCD646, UCD650, yHAB47, yHCT81, yHRVM107, yHRVM108, CDFM21L, ABFM5L))

sw_chr16<-sliding.window.transform(CHR16.class,100,100,type=2, whole.data = T)
length(sw_chr16@region.names) 

sw_chr16 <-diversity.stats.between(sw_chr16,new.populations = list(HOLadmix,PB2,PB3,PB5))


as_tibble(sw_chr16@nuc.diversity.between) %>% mutate(`pop1/pop2` = `pop1/pop2`/100, `pop1/pop3` = `pop1/pop3`/100, `pop1/pop4`=`pop1/pop4`/100, x = 1:nrow(sw_chr16@nuc.diversity.between)) %>%
  select(-4) %>% pivot_longer(cols = 1:3, names_to = "comp", values_to = "div") %>% ggplot(aes(x=x,y=div, color = comp)) + geom_point(size= 1) +
  theme_bw() + scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population") + 
  ggtitle("CHR16") + scale_y_continuous(breaks = seq(0, 0.05, by= 0.01))


#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################


VCF_split_into_scaffolds("./whole2/snp_data.recode.vcf", "vcf2")

GENOME <- readData("./vcf2/", format = "VCF")
GENOME@region.names

sw<-sliding.window.transform(GENOME,10,10,type=1,whole.data = T)

ADM_VIL<-c("CL210.1")

PB3_VIL<-c("CL204.3","CL211.3","CL212.2","CL215.1","CL216.1","CL218.1","CL220.3","CL221.1")

PB1_VOS<-c("CL605.1","CL606.1","CL619.1","CL621.1")


sw <-diversity.stats.between(sw,new.populations = list(ADM_VIL,PB3_VIL,PB1_VOS))


as_tibble(sw@nuc.diversity.between) %>% mutate(`pop1/pop2`=`pop1/pop2`/10, `pop1/pop3`=`pop1/pop3`/10,x=1:nrow(sw@nuc.diversity.between)) %>% 
  select(-3) %>%pivot_longer(cols = 1:2,names_to = "comp",values_to = "div") %>% 
  ggplot(aes(x=x,y=div,color=comp))+geom_point(size=1)+ theme_minimal()


HOLad= c("ABFM5L"  , "CDFM21L"  ,  "UCD646"  , "UCD650" ,  "yHAB47"  ,"yHCT81" ,  "yHRVM107" ,"yHRVM108")

HOL_IRLANDA = c("UCD646")

HOL_USA = c("yHRVM108")

HOL_China = c("CDFM21L")

PB2= c("CL1101.1","CL1108.1","CL1109.1","CL248.1","CL701.1","CL702.1","CL703.2","CL704.2","CL705.1","CL706.2","CL710.1","CL711.2","CL715.1","CL910.1","CL915.1","CO104","CO105","CO108","CO109","CO110","CO114","CO122","CO123","CO126","CO16","CO17","CO18","CO19","CO20","CO21","CO22","CO23","CO24","CO25","CO26","CO27","CO28","CO40","CO41","CO42","CO43","CO44","CO45","CO46","CO47","CO49","CO50","CO51","CO52","CO53","CO54","CO55","CO56","CO57","CO59","CO60","CO75","CO76","CO80","CO89","CO95","yHAB124","yHAB125","yHAB126","yHAB66","yHAB67","yHAB74","yHCT100","yHCT107","yHCT117","yHCT119","yHCT70","yHCT71","yHCT88","yHCT89","yHCT91","yHCT97","yHQL2248","yHQL2250","yHQL2265","yHQL726","yHQL732","yHQL733","yHQL736","yHQL738","yHQL743","yHQL744","yHQL747","yHQL758")

PB3= c("CL1104.1","CL1110.1","CL204.3","CL211.3","CL212.2","CL215.1","CL216.1","CL218.1","CL220.3","CL221.1","CL601.1","CL604.1","CL607.1","CL608.1","CL610.1","CL611.1","CL620.1","CL902.1","CL903.1","CL904.1","CL905.1","CL906.1","CL907.1","CL909.1","CL916.1","CO100","CO101","CO102","CO103","CO106","CO107","CO111","CO112","CO113","CO115","CO116","CO118","CO119","CO120","CO121","CO124","CO125","CO127","CO128","CO129","CO130","CO132","CO135","CO136","CO139","CO140","CO141","CO142","CO144","CO145","CO147","CO151","CO152","CO153","CO155","CO156","CO157","CO61","CO62","CO66","CO69","CO72","CO79","CO84","CO85","CO86","CO92","CO97","CO98","CO99","S13_HH","yHCT105","yHCT115","yHQL2253","yHQL2273")

PB5= c("CO146","CO149","CO150","yHAB54","yHAB55","yHAB57","yHAB59","yHAB62","yHAB65","yHAB68","yHAB75","yHAB77","yHCT108","yHCT120","yHCT61","yHCT87","yHCT93","yHCT94", "yHCT95")

sw <-diversity.stats.between(sw,new.populations = list(HOLad,PB2,PB3,PB5))


as_tibble(sw@nuc.diversity.between) %>%
  mutate(`pop1/pop2`=`pop1/pop2`/10, `pop1/pop3`=`pop1/pop3`/10,`pop1/pop4`=`pop1/pop4`/10, x=1:nrow(sw@nuc.diversity.between)) %>%
  select(-c(4,5,6)) %>%pivot_longer(cols = 1:3,names_to = "comp",values_to = "div") %>%
  ggplot(aes(x=x,y=div,color=comp))+geom_point(size=1)+ theme_bw() + ggtitle("Población completa") +
  scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population")
  

sw_irlanda <-diversity.stats.between(sw,new.populations = list(HOL_IRLANDA,PB2,PB3,PB5))


as_tibble(sw_irlanda@nuc.diversity.between) %>%
  mutate(`pop1/pop2`=`pop1/pop2`/10, `pop1/pop3`=`pop1/pop3`/10,`pop1/pop4`=`pop1/pop4`/10, x=1:nrow(sw_irlanda@nuc.diversity.between)) %>%
  select(-c(4,5,6)) %>%pivot_longer(cols = 1:3,names_to = "comp",values_to = "div") %>%
  ggplot(aes(x=x,y=div,color=comp))+geom_point(size=1)+ theme_bw() + ggtitle("UCD646") +
  scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population")

sw_irlanda@
data.frame(sw_irlanda@n.sites, 1:length(sw_irlanda@n.sites)) %>% filter(sw_irlanda.n.sites < 0)


sw_USA <-diversity.stats.between(sw,new.populations = list(HOL_USA,PB2,PB3,PB5))


as_tibble(sw_USA@nuc.diversity.between) %>%
  mutate(`pop1/pop2`=`pop1/pop2`/10, `pop1/pop3`=`pop1/pop3`/10,`pop1/pop4`=`pop1/pop4`/10, x=1:nrow(sw_USA@nuc.diversity.between)) %>%
  select(-c(4,5,6)) %>%pivot_longer(cols = 1:3,names_to = "comp",values_to = "div") %>%
  ggplot(aes(x=x,y=div,color=comp))+geom_point(size=1)+ theme_bw() + ggtitle("yHRVM108") +
  scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population")


sw_China <-diversity.stats.between(sw,new.populations = list(HOL_China,PB2,PB3,PB5))


as_tibble(sw_China@nuc.diversity.between) %>%
  mutate(`pop1/pop2`=`pop1/pop2`/10, `pop1/pop3`=`pop1/pop3`/10,`pop1/pop4`=`pop1/pop4`/10, x=1:nrow(sw_China@nuc.diversity.between)) %>%
  select(-c(4,5,6)) %>%pivot_longer(cols = 1:3,names_to = "comp",values_to = "div") %>%
  ggplot(aes(x=x,y=div,color=comp))+geom_point(size=1)+ theme_bw() + ggtitle("CDFM21L") +
  scale_color_manual(values= c("#2E86C1", "#1ABC9C","#34495E"), labels = c("PB2", "PB3", "PB5")) + labs(x = "Window", y = "Dxy", color = "Population")



##############################################
##############################################


# distribucion entre poblaciones pi
pi_g <- slide_data %>% select(contains("pi")) %>% gather(key = "Pop", 
                                                         value = "pi")
pi_acrossGenome <- ggplot(pi_g, aes(Pop, pi, fill = Pop)) + geom_boxplot() + 
  theme_gray() + xlab(NULL) + ggtitle("Diversidad nucleotidica (Pi) Whole Genome")+
  scale_fill_manual(values= c( "#95A5A6", "#F4D03F","#E74C3C", "#2E86C1", "#1ABC9C","#9B59B6","#34495E", "#ABEBC6", "#D6EAF8"))+
  theme_gray() + labs(fill = "Populations")
pi_acrossGenome



# Librerías necesarias
library(ggplot2)
library(tidyr)
library(ggthemes)

# Datos en formato matriz
tajima_d <- data.frame(
  Population = paste0("pop ", 1:9),
  Value = c(-0.1866564, -2.487714, 1.133235, -1.044703, -0.3563752, 0.3463621, 0.6988107, 0.5777154, -1.511066)
)

nucdiv <- data.frame(
  Population = paste0("pop ", 1:9),
  Value = c(9.224372e-05, 1.060571e-06, 0.0002890043, 9.016082e-05, 3.018856e-05, 0.0001111895, 0.0001676285, 1.144444e-05, 0.00014912)
)

# Gráfico de Tajima's D
tajima_plot <- ggplot(tajima_d, aes(x = Population, y = Value)) +
  geom_bar(stat = "identity", fill = "#66C2A5", color = "black") +
  theme_minimal(base_family = "Arial") +
  labs(x = "Population", y = "Tajima's D", title = "Tajima's D across Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico de diversidad nucleotídica
nucdiv_plot <- ggplot(nucdiv, aes(x = Population, y = Value)) +
  geom_bar(stat = "identity", fill = "#8DA0CB", color = "black") +
  theme_minimal(base_family = "Arial") +
  labs(x = "Population", y = "Nucleotide Diversity", title = "Nucleotide Diversity across Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mostrar gráficos
print(tajima_plot)
print(nucdiv_plot)

# Filtrar para eliminar pop 9
tajima_d <- tajima_d %>% 
  filter(Population != "pop 9")

nucdiv <- nucdiv %>% 
  filter(Population != "pop 9")



#####objeto de pop con los individuos de la pop
pop_data %>% filter(Pop == "PA") %>% pull(Ind) %>% as.character() -> PA
PA ##pop1

pop_data %>% filter(Pop == "NoAm") %>% pull(Ind) %>% as.character() -> NoAm
NoAm ##pop2

pop_data %>% filter(Pop == "PB1") %>% pull(Ind) %>% as.character()-> PB1
PB1 ##pop3

pop_data %>% filter(Pop == "PB2") %>% pull(Ind) %>% as.character() -> PB2
PB2 ##pop4

pop_data %>% filter(Pop == "PB3") %>% pull(Ind) %>% as.character()-> PB3
PB3 ##pop5

pop_data %>% filter(Pop == "PB4") %>% pull(Ind) %>% as.character() -> PB4
PB4 ##pop6

pop_data %>% filter(Pop == "PB5") %>% pull(Ind) %>% as.character()-> PB5
PB5 ##pop7

pop_data %>% filter(Pop == "PB6") %>% pull(Ind) %>% as.character()-> PB6

# Reordenar las poblaciones (pop 6 en la posición de pop 4)
new_order <- c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 8","pop 7", "pop 6")

tajima_d$Population <- factor(tajima_d$Population, levels = new_order)
nucdiv$Population <- factor(nucdiv$Population, levels = new_order)

# Gráfico de Tajima's D
tajima_plot <- ggplot(tajima_d, aes(x = Population, y = Value)) +
  geom_bar(stat = "identity", fill = "#66C2A5", color = "black") +
  theme_minimal(base_family = "Arial") +
  labs(x = "Population", y = "Tajima's D", title = "Tajima's D across Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Gráfico de diversidad nucleotídica
nucdiv_plot <- ggplot(nucdiv, aes(x = Population, y = Value)) +
  geom_bar(stat = "identity", fill = "#8DA0CB", color = "black") +
  theme_minimal(base_family = "Arial") +
  labs(x = "Population", y = "Nucleotide Diversity", title = "Nucleotide Diversity across Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mostrar gráficos
print(tajima_plot)
print(nucdiv_plot)

# Definir colores personalizados
custom_colors <- c("#F4D03F", "#95A5A6", "#E74C3C", "#2E86C1", "#1ABC9C", "#ABEBC6", "#34495E", "#9B59B6")

# Gráfico de diversidad nucleotídica con puntos y colores personalizados
nucdiv_plot <- ggplot(nucdiv, aes(x = Population, y = Value, color = Population)) +
  geom_point(size = 4) + # Tamaño de los puntos
  scale_color_manual(values = custom_colors) + # Aplicar colores personalizados
  theme_minimal(base_family = "Arial") +
  labs(
    x = "Population", 
    y = "Nucleotide Diversity", 
    title = "Nucleotide Diversity across Populations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none" # Ocultar la leyenda
  ) + theme_base()

# Mostrar gráfico
print(nucdiv_plot)

# Gráfico de Tajima's D como barras con colores personalizados
tajima_barplot <- ggplot(tajima_d, aes(x = Population, y = Value, fill = Population)) +
  geom_bar(stat = "identity", color = "black", width = 0.5) + # Barras con borde negro
  scale_fill_manual(values = custom_colors) + # Aplicar colores personalizados
  theme_minimal(base_family = "Arial") +
  labs(
    x = "Population", 
    y = "Tajima's D", 
    title = "Tajima's D across Populations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Texto del eje X inclinado y más pequeño
    legend.position = "none" # Ocultar la leyenda
  ) + theme_base()

# Mostrar gráfico
print(tajima_barplot)



#______________________________

# Crear el dataframe de diversidad nucleotídica (nucdiv)
nuc_dis <- data.frame(
  Population = c("PA", "PB1", "PB2", "PB3", "PB4", "PB5", "PB6"),
  Value = c(9.224372e-05, 2.890043e-04, 9.016082e-05, 
            3.018856e-05, 1.111895e-04, 1.676285e-04, 1.144444e-05)
)

# Crear un vector con las distancias geográficas (en km) para cada población
distancias <- c(365, 967, 507, 293, 140, 50, 60)  # Ejemplo de distancias, reemplazar con tus valores reales

# Añadir las distancias al dataframe
nuc_dis$Distancia_km <- distancias

# Ver el dataframe actualizado
print(nuc_dis)


# Crear el gráfico de dispersión
ggplot(nuc_dis, aes(x = Distancia_km, y = Value)) +
  geom_point(aes(color = Population), size = 4) +  # Puntos con colores por población
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Línea de regresión
  labs(title = "Correlación entre la distancia geográfica y la diversidad nucleotídica",
       x = "Distancia geográfica (km)",
       y = "Diversidad nucleotídica (π)") +
  theme_minimal()  # Tema minimalista

library(viridis)

ggplot(nuc_dis, aes(x = Distancia_km, y = Value)) +
  geom_point(aes(color = Population), size = 4) +  # Puntos con colores por población
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Línea de regresión
  geom_density_2d_filled(aes(fill = after_stat(level)), color = "white") +  # Densidad 2D con relleno y mapeo de color por nivel
  scale_fill_viridis() +  # Cambiar la escala de color para la densidad (usando scale_fill_viridis)
  labs(title = "Correlación entre la distancia geográfica y la diversidad nucleotídica",
       x = "Distancia geográfica (km)",
       y = "Diversidad nucleotídica (π)") +
  theme_minimal()  # Tema minimalista





############################################################################################
############### TIEMPOS DE DIVERGENCIA #####################################################
############################################################################################


pop_data <- read.table("./popstacks.txt", sep = "\t", header = T)
pop_data

pop_data <- subset(pop_data, Pop != "Admix")
pop_data

# Abrir VCF en genoma completo o por cromosoma.

GENOME.class <- readData("./whole2", format = "VCF")

GENOME.class@n.biallelic.sites
GENOME.class@region.names
GENOME.class@n.sites2

# Pedir resumen de datos 
get.sum.data(GENOME.class)


#####objeto de pop con los individuos de la pop
pop_data %>% filter(Pop == "PA") %>% pull(Ind) %>% as.character() -> PA
PA ##pop1

pop_data %>% filter(Pop == "NoAm") %>% pull(Ind) %>% as.character() -> NoAm
NoAm ##pop2

pop_data %>% filter(Pop == "PB1") %>% pull(Ind) %>% as.character()-> PB1
PB1 ##pop3

pop_data %>% filter(Pop == "PB2") %>% pull(Ind) %>% as.character() -> PB2
PB2 ##pop4

pop_data %>% filter(Pop == "PB3") %>% pull(Ind) %>% as.character()-> PB3
PB3 ##pop5

pop_data %>% filter(Pop == "PB4") %>% pull(Ind) %>% as.character() -> PB4
PB4 ##pop6

pop_data %>% filter(Pop == "PB5") %>% pull(Ind) %>% as.character()-> PB5
PB5 ##pop7

pop_data %>% filter(Pop == "PB6") %>% pull(Ind) %>% as.character()-> PB6
PB6 ##pop8

pop_data %>% filter(Pop == "Hol-ADM") %>% pull(Ind) %>% as.character()-> HOLadmix
HOLadmix ##pop9


pop_data %>% filter(Pop == c("Hol-ADM","PB6", "PB5", "PB4", "PB3","PB2","PB1","NoAm")) %>% pull(Ind) %>% as.character()-> POP1
pop_data %>% 

GENOME.class <-set.populations( GENOME.class,list(PA, NoAm, PB1, PB2, PB3, PB4, PB5, PB6, HOLadmix))
#Check poblaciones diploides
GENOME.class@populations
genome_size <- 12500000 # Añadir tamaño del genoma/cromosoma
GENOME.class@n.sites2
snp_number <- 1330729 # Añadir n° snps

# Neutrality stats 
eubayanus_genome <- neutrality.stats(GENOME.class)
eubayanus_genome@theta_Watterson/snp_number


pop_data %>% filter(Pop %in% c("PB2","Hol-ADM", "PB4", "PB5", "PB1", "PB3", "PB6", "NoAm", "PA" )) %>% 
  pull(Ind) %>% 
  as.character() -> POP0



GENOME.class <-set.populations( GENOME.class,list(POP8, POP7,PB6, POP3,POP4, POP5, POP2, POP1,POP0), diploid = T)
#Check poblaciones diploides
GENOME.class@populations
genome_size <- 12500000 # Añadir tamaño del genoma/cromosoma
GENOME.class@n.sites2
snp_number <- 1330729 # Añadir n° snps

# Neutrality stats 
eubayanus_genome <- neutrality.stats(GENOME.class)
eubayanus_genome@theta_Watterson/snp_number
