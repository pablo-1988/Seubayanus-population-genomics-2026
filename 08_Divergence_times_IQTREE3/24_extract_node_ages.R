#!/usr/bin/env Rscript
# Extrae edades de nodos (TMRCA) del time tree de LSD2 para cada par de poblaciones
# y para el crown de cada poblacion. Edades en anos antes del presente.
suppressMessages({library(ape)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- file.path(proj,"results/K5/dating")
tr <- read.tree(file.path(D,"dated3.timetree.nwk"))
tp <- read.table(file.path(D,"taxa_pop.tsv"), sep="\t", col.names=c("ID","POP"))
p <- setNames(tp$POP, tp$ID)
pop <- p[tr$tip.label]
cat("hojas en el time tree:", Ntip(tr), " (outgroup excluido)\n")
print(table(pop, useNA="ifany"))

# El .nwk de LSD2 queda en unidades de SUSTITUCION (ultrametrico). Se pasa a anos
# dividiendo por la tasa fijada (mu * g). Chequeo: H/RATE debe dar el tMRCA de LSD2.
RATE <- 4.876e-7
dep <- node.depth.edgelength(tr)
H <- max(dep[1:Ntip(tr)])              # root-to-tip en subs/sitio
cat("root-to-tip (subs/sitio):", signif(H,6), "\n")
cat("tMRCA eubayanus (anos):", round(H/RATE,1), " [LSD2 reporto 2534.8]\n\n")
age <- function(node) (H - dep[node]) / RATE   # anos antes del presente

POPS <- c("PA","PB1","PB2","PB3","PB4","NoAm","Hol")
tips <- lapply(POPS, function(q) tr$tip.label[which(pop==q)]); names(tips) <- POPS

# crown de cada poblacion
crown <- sapply(POPS, function(q){
  t <- tips[[q]]; if(length(t)<2) return(NA)
  age(getMRCA(tr, t))
})
cat("=== crown age por poblacion (anos) ===\n"); print(round(crown,0))

# TMRCA por par
rows <- list()
for(i in seq_along(POPS)) for(j in seq_along(POPS)) if(i<j){
  a<-POPS[i]; b<-POPS[j]
  t <- c(tips[[a]], tips[[b]])
  rows[[length(rows)+1]] <- data.frame(pop1=a, pop2=b, TMRCA_years=age(getMRCA(tr,t)))
}
res <- do.call(rbind, rows)
res <- res[order(res$TMRCA_years),]
write.table(res, file.path(D,"tmrca_lsd2.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(POP=POPS, crown_years=crown), file.path(D,"crown_ages_lsd2.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
cat("\n=== TMRCA por par (LSD2, anos) ===\n")
print(data.frame(res, row.names=NULL), digits=5)
