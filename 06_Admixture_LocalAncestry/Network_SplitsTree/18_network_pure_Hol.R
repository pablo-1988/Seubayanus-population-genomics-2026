#!/usr/bin/env Rscript
# NeighborNet: todas las puras (PA/PB1-4) + grupo Holartico (Hol), n=299.
# Splits calculados por SplitsTree6; phangorn solo hace layout+plot.
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
OUT  <- file.path(proj,"results/K5/figures_EN"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
nex  <- file.path(proj,"results/K5/network/pure_Hol_splits.nex")
D <- as.matrix(fread(file.path(proj,"results/05_network/eub_distance.csv"),header=TRUE),rownames=1)
sets <- fread(file.path(proj,"results/K5/admix_subpops_noGen/SETS.txt"),header=FALSE)
setnames(sets,c("ID","POP")); s2p <- setNames(sets$POP,sets$ID)
L <- readLines(nex)
# taxa en el orden del NEXUS de entrada (bloque Taxa)
tx0 <- grep("TAXLABELS",L); tx1 <- grep("^;",L); tx1 <- tx1[tx1>tx0][1]
taxa <- sub("^\\[[0-9]+\\] ","",trimws(L[(tx0+1):(tx1-1)]))
taxa <- gsub("^'|'$","",taxa)   # SplitsTree envuelve las etiquetas en comillas

cyc <- as.integer(strsplit(sub(";","",sub("^CYCLE ","",grep("^CYCLE",L,value=TRUE)[1]))," +")[[1]]); cyc<-cyc[!is.na(cyc)]
i0<-grep("BEGIN SPLITS;",L); i1<-grep("MATRIX",L); i1<-i1[i1>i0][1]; iend<-grep("^;",L); iend<-iend[iend>i1][1]
wt<-numeric(0); sp<-list()
for(r in L[(i1+1):(iend-1)]){ p<-strsplit(r,"\t")[[1]]; if(length(p)<3) next
  wt<-c(wt,as.numeric(trimws(p[2]))); txi<-as.integer(strsplit(trimws(gsub(",","",p[3]))," +")[[1]]); sp[[length(sp)+1]]<-txi[!is.na(txi)] }
attr(sp,"weights")<-wt; attr(sp,"labels")<-taxa; attr(sp,"cycle")<-cyc; class(sp)<-"splits"
net <- as.networx(sp)

cols <- c(PA="#F4D03F",PB1="#E74C3C",PB2="#2E86C1",PB3="#1ABC9C",PB4="#34495E",Hol="#E91E63")
pop <- s2p[net$tip.label]; tc <- cols[pop]; tc[is.na(tc)]<-"grey70"
# Hol resaltado con puntos mas grandes
cexv <- ifelse(pop=="Hol", 1.3, 0.5); cexv[is.na(pop)] <- 0.5
cat("taxa:",length(taxa),"| ",paste(names(table(pop)),table(pop),collapse=" "),"\n")

pdf(file.path(OUT,"21_neighbornet_pure_Hol_labels.pdf"),width=13,height=13)
plot(net,"2D",tip.color=tc,cex=0.42,edge.width=0.5,col.edge="grey45",font=1)
legend("topleft",legend=names(cols),text.col=cols,pch=19,col=cols,bty="n",cex=1.2)
title(sprintf("NeighborNet — pure lineages + Holarctic group (n=%d; SplitsTree6, fit 99.8%%)",length(taxa)))
dev.off(); cat("-> 21_neighbornet_pure_Hol_labels\n")

# version con puntos coloreados: reemplazar cada etiqueta por un disco y colorearlo
net2 <- net; net2$tip.label <- rep("●", length(net$tip.label))
pdf(file.path(OUT,"22_neighbornet_pure_Hol_dots.pdf"),width=13,height=13)
plot(net2,"2D",tip.color=tc,cex=0.75,edge.width=0.5,col.edge="grey55",font=1)
legend("topleft",legend=names(cols),text.col=cols,pch=19,col=cols,bty="n",cex=1.2)
title(sprintf("NeighborNet - pure lineages + Holarctic group (n=%d; tips as dots)",length(taxa)))
dev.off(); cat("-> 22_neighbornet_pure_Hol_dots\n")
