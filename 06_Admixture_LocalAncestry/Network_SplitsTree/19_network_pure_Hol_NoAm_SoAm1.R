#!/usr/bin/env Rscript
# NeighborNet: puras (PA/PB1-4) + Hol + NoAm + SoAm1, n=330.
# Splits de SplitsTree6; phangorn solo layout+plot.
suppressMessages({library(phangorn); library(data.table)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
OUT  <- file.path(proj,"results/K5/figures_EN"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
nex  <- file.path(proj,"results/K5/network/pure_Hol_NoAm_SoAm1_splits.nex")
sets <- fread(file.path(proj,"results/K5/admix_subpops_noGen/SETS.txt"),header=FALSE)
setnames(sets,c("ID","POP")); s2p <- setNames(sets$POP,sets$ID)
L <- readLines(nex)
tx0<-grep("TAXLABELS",L); tx1<-grep("^;",L); tx1<-tx1[tx1>tx0][1]
taxa<-sub("^\\[[0-9]+\\] ","",trimws(L[(tx0+1):(tx1-1)])); taxa<-gsub("^'|'$","",taxa)
cyc<-as.integer(strsplit(sub(";","",sub("^CYCLE ","",grep("^CYCLE",L,value=TRUE)[1]))," +")[[1]]); cyc<-cyc[!is.na(cyc)]
i0<-grep("BEGIN SPLITS;",L); i1<-grep("MATRIX",L); i1<-i1[i1>i0][1]; iend<-grep("^;",L); iend<-iend[iend>i1][1]
wt<-numeric(0); sp<-list()
for(r in L[(i1+1):(iend-1)]){ p<-strsplit(r,"\t")[[1]]; if(length(p)<3) next
  wt<-c(wt,as.numeric(trimws(p[2]))); txi<-as.integer(strsplit(trimws(gsub(",","",p[3]))," +")[[1]]); sp[[length(sp)+1]]<-txi[!is.na(txi)] }
attr(sp,"weights")<-wt; attr(sp,"labels")<-taxa; attr(sp,"cycle")<-cyc; class(sp)<-"splits"
net<-as.networx(sp)

cols<-c(PA="#F4D03F",PB1="#E74C3C",PB2="#2E86C1",PB3="#1ABC9C",PB4="#34495E",
        Hol="#E91E63",NoAm="#000000",SoAm1="#8E44AD")
pop<-s2p[net$tip.label]; tc<-cols[pop]; tc[is.na(tc)]<-"grey70"
cat("taxa:",length(taxa),"| ",paste(names(table(pop)),table(pop),collapse=" "),"\n")

pdf(file.path(OUT,"23_neighbornet_pure_Hol_NoAm_SoAm1_labels.pdf"),width=13,height=13)
plot(net,"2D",tip.color=tc,cex=0.4,edge.width=0.5,col.edge="grey45",font=1)
legend("topleft",legend=names(cols),text.col=cols,pch=19,col=cols,bty="n",cex=1.2)
title(sprintf("NeighborNet - pure lineages + Hol + NoAm + SoAm1 (n=%d; SplitsTree6, fit 99.5%%)",length(taxa)))
dev.off(); cat("-> 23_labels\n")

net2<-net; net2$tip.label<-rep("●",length(net$tip.label))
pdf(file.path(OUT,"24_neighbornet_pure_Hol_NoAm_SoAm1_dots.pdf"),width=13,height=13)
plot(net2,"2D",tip.color=tc,cex=0.75,edge.width=0.5,col.edge="grey55",font=1)
legend("topleft",legend=names(cols),text.col=cols,pch=19,col=cols,bty="n",cex=1.2)
title(sprintf("NeighborNet - pure lineages + Hol + NoAm + SoAm1 (n=%d; tips as dots)",length(taxa)))
dev.off(); cat("-> 24_dots\n")
