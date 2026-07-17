#!/usr/bin/env Rscript
# Dibuja la NeighborNet de las 471 cepas usando los splits YA calculados por
# SplitsTree6 (evita el nnls/overflow de phangorn). Lineas finas + tips por K5.
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
nex  <- file.path(proj,"results/K5/network/eub_all_splits.nex")
D <- as.matrix(fread(file.path(proj,"results/05_network/eub_distance.csv"), header=TRUE), rownames=1)
taxa <- rownames(D)                      # orden 1..471 = orden del NEXUS de distancias
L <- readLines(nex)

# --- CYCLE ---
cyc_line <- grep("^CYCLE", L, value=TRUE)[1]
cycle <- as.integer(strsplit(sub(";","",sub("^CYCLE ","",cyc_line))," +")[[1]])
cycle <- cycle[!is.na(cycle)]

# --- SPLITS matrix ---
i0 <- grep("BEGIN SPLITS;", L); i1 <- grep("MATRIX", L); i1 <- i1[i1>i0][1]
iend <- grep("^;", L); iend <- iend[iend>i1][1]
rows <- L[(i1+1):(iend-1)]
wt <- numeric(0); sp <- list()
for(r in rows){
  parts <- strsplit(r, "\t")[[1]]
  if(length(parts) < 3) next
  wt <- c(wt, as.numeric(trimws(parts[2])))
  tx <- as.integer(strsplit(trimws(gsub(",","",parts[3]))," +")[[1]])
  sp[[length(sp)+1]] <- tx[!is.na(tx)]
}
cat("splits:", length(sp), " taxa:", length(taxa), " cycle:", length(cycle), "\n")

# --- objeto splits de phangorn ---
attr(sp,"weights") <- wt
attr(sp,"labels")  <- taxa
attr(sp,"cycle")   <- cycle
class(sp) <- "splits"

net <- as.networx(sp)                    # layout circular, sin nnls

# --- colores K5 ---
sets <- fread(file.path(proj,"results/K5/admix_subpops_noGen/SETS.txt"), header=FALSE)
setnames(sets,c("ID","POP")); s2p <- setNames(sets$POP, sets$ID)
cols <- c(PA="#F4D03F", PB1="#E74C3C", PB2="#2E86C1", PB3="#1ABC9C", PB4="#34495E",
          NoAm="#000000", SoAm1="#8E44AD", SoAm2="#E67E22", SoAm3="#16A085",
          SoAm4="#C0392B", Hol="#2980B9", xxx="#BDC3C7")
pop <- s2p[net$tip.label]; tipcol <- cols[pop]; tipcol[is.na(tipcol)] <- "#BDC3C7"

outdir <- file.path(proj,"results/K5/network")
png(file.path(outdir,"neighbornet_all.png"), width=4200, height=4200, res=320)
plot(net, "2D", tip.color=tipcol, cex=0.30, edge.width=0.45, col.edge="grey40",
     font=1, show.tip.label=TRUE)
legend("topleft", legend=names(cols), text.col=cols, pch=19, col=cols, bty="n", cex=1.2)
title(sprintf("NeighborNet — S. eubayanus todas las cepas (n=%d; SplitsTree6, fit 99.3%%)", length(taxa)))
dev.off()

# version sin etiquetas para ver la topologia de reticulaciones
png(file.path(outdir,"neighbornet_all_nolabels.png"), width=4200, height=4200, res=320)
plot(net, "2D", show.tip.label=FALSE, edge.width=0.45, col.edge="grey40")
legend("topleft", legend=names(cols), text.col=cols, pch=19, col=cols, bty="n", cex=1.2)
title("NeighborNet — reticulaciones (sin etiquetas)")
dev.off()
cat("Redes generadas en", outdir, "\n")
