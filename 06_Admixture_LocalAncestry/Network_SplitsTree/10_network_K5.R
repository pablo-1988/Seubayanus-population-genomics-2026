#!/usr/bin/env Rscript
# Parte 5 (K5): NeighborNet coloreado por esquema K=5 (puros + sub-admix).
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- as.matrix(fread(file.path(proj,"results/05_network/eub_distance.csv"), header=TRUE), rownames=1)
sets <- fread(file.path(proj,"results/K5/admix_subpops_noGen/SETS.txt"), header=FALSE)
setnames(sets, c("ID","POP")); s2p <- setNames(sets$POP, sets$ID)
outdir <- file.path(proj,"results/K5/network"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

set.seed(7); allids <- rownames(D); popall <- s2p[allids]
pures <- c("PA","PB1","PB2","PB3","PB4")
subadm <- c("NoAm","SoAm1","SoAm2","SoAm3","SoAm4","Hol")
keep <- c()
for(p in pures){ ids <- allids[which(popall==p)]; keep <- c(keep, if(length(ids)>6) sample(ids,6) else ids) }
# foco Holartico: todas NoAm+Hol; SoAm submuestreadas a 10 para legibilidad/velocidad
for(p in c("NoAm","Hol")){ keep <- c(keep, allids[which(popall==p)]) }
for(p in c("SoAm1","SoAm2","SoAm3","SoAm4")){ ids <- allids[which(popall==p)]
  keep <- c(keep, if(length(ids)>10) sample(ids,10) else ids) }
D <- D[keep, keep]; pop <- s2p[rownames(D)]
cat("taxones:", nrow(D), "\n"); print(table(pop))

cols <- c(PA="#F4D03F", PB1="#E74C3C", PB2="#2E86C1", PB3="#1ABC9C", PB4="#34495E",
          NoAm="#000000", SoAm1="#8E44AD", SoAm2="#E67E22", SoAm3="#16A085",
          SoAm4="#C0392B", Hol="#2980B9")
tipcol <- cols[pop]
nnet <- neighborNet(as.dist(D)); ntax <- nrow(D)
png(file.path(outdir,"neighbornet_K5.png"), width=2800, height=2300, res=230)
plot(nnet, "2D", tip.color=tipcol, cex=0.5, col.edge="grey60", font=1)
legend("topleft", legend=names(cols), text.col=cols, pch=19, col=cols, bty="n", cex=1.0)
title(sprintf("NeighborNet — S. eubayanus K=5 (n=%d; sub-admix completas)", ntax))
dev.off()
cat("Red K5 generada.\n")
