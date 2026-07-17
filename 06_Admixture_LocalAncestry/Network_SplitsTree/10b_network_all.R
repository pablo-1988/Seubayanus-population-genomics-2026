#!/usr/bin/env Rscript
# NeighborNet con TODAS las cepas (471), lineas finas para ver reticulaciones.
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- as.matrix(fread(file.path(proj,"results/05_network/eub_distance.csv"), header=TRUE), rownames=1)
sets <- fread(file.path(proj,"results/K5/admix_subpops_noGen/SETS.txt"), header=FALSE)
setnames(sets, c("ID","POP")); s2p <- setNames(sets$POP, sets$ID)
outdir <- file.path(proj,"results/K5/network"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
pop <- s2p[rownames(D)]
cat("taxones:", nrow(D), "\n"); print(table(pop, useNA="ifany"))

cols <- c(PA="#F4D03F", PB1="#E74C3C", PB2="#2E86C1", PB3="#1ABC9C", PB4="#34495E",
          NoAm="#000000", SoAm1="#8E44AD", SoAm2="#E67E22", SoAm3="#16A085",
          SoAm4="#C0392B", Hol="#2980B9", xxx="#BDC3C7")
tipcol <- cols[pop]; tipcol[is.na(tipcol)] <- "#BDC3C7"

cat("construyendo NeighborNet (471 taxones, puede tardar)...\n")
nnet <- neighborNet(as.dist(D)); ntax <- nrow(D)
png(file.path(outdir,"neighbornet_all.png"), width=4000, height=4000, res=300)
plot(nnet, "2D", tip.color=tipcol, cex=0.28, edge.width=0.25, col.edge="grey45", font=1)
legend("topleft", legend=names(cols), text.col=cols, pch=19, col=cols, bty="n", cex=1.1)
title(sprintf("NeighborNet — S. eubayanus todas las cepas (n=%d)", ntax))
dev.off()
saveRDS(nnet, file.path(outdir,"nnet_all.rds"))
cat("Red completa generada.\n")
