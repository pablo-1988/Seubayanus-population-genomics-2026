#!/usr/bin/env Rscript
# Parte 5: red de reticulacion NeighborNet con TODAS las cepas eubayanus (incl. admix).
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- as.matrix(fread(file.path(proj,"results/05_network/eub_distance.csv"), header=TRUE), rownames=1)
meta <- fread(file.path(proj,"results/metadata_final.txt")); setnames(meta,c("ID","POP"))
# subsample: todas Admix+NoAm + hasta 12 por linaje puro (legibilidad + evita overflow)
set.seed(7); allids <- rownames(D)
popall <- meta$POP[match(allids, meta$ID)]
keep <- c()
for(p in c("PA","PB1","PB2","PB3","PB4","PB5","PB6")){
  ids <- allids[popall==p]; keep <- c(keep, if(length(ids)>5) sample(ids,5) else ids)
}
adm <- allids[popall=="Admix"]; keep <- c(keep, if(length(adm)>25) sample(adm,25) else adm)
keep <- c(keep, allids[popall=="NoAm"])   # todas las Holarticas
D <- D[keep, keep]
pop <- meta$POP[match(rownames(D), meta$ID)]
cat("taxones en la red:", nrow(D), "\n")

cols <- c(PA="#1b9e77", PB1="#d95f02", PB2="#7570b3", PB3="#e7298a",
          PB4="#66a61e", PB5="#e6ab02", PB6="#a6761d", NoAm="#000000", Admix="#999999")
tipcol <- cols[pop]

nnet <- neighborNet(as.dist(D))
ntax <- nrow(D)
png(file.path(proj,"results/05_network/neighbornet.png"), width=2600, height=2200, res=220)
plot(nnet, "2D", tip.color=tipcol, cex=0.55, col.edge="grey55", font=1)
legend("topleft", legend=names(cols), text.col=cols, pch=19, col=cols, bty="n", cex=1.1)
title(sprintf("NeighborNet — S. eubayanus (subset n=%d; incl. todas las Holarticas)", ntax))
dev.off()

# version resaltando SOLO admix vs puras
tc2 <- ifelse(pop %in% c("Admix"), "#d7191c", ifelse(pop=="NoAm","#2b83ba","grey70"))
png(file.path(proj,"results/05_network/neighbornet_admix.png"), width=2600, height=2200, res=220)
plot(nnet, "2D", tip.color=tc2, cex=0.55, col.edge="grey70", font=1)
legend("topleft", legend=c("Admix (rojo)","NoAm/Holartico (azul)","Puras (gris)"),
       text.col=c("#d7191c","#2b83ba","grey40"), bty="n", cex=1.1)
title(sprintf("Reticulaciones y cepas admixadas (NeighborNet, n=%d)", ntax))
dev.off()
cat("Redes generadas.\n")
