#!/usr/bin/env Rscript
# Figura: arbol de POBLACIONES con el tiempo de divergencia en cada nodo.
# Colapsa el time tree de LSD2 (380 cepas) a 7 poblaciones.
# Nota: 3 cepas PB1 caen fuera del clado principal de PB1 (parafilia menor) -> se
# excluyen del colapso y se documentan en el pie de figura.
suppressMessages({library(ape)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- file.path(proj,"results/K5/dating"); OUT <- file.path(proj,"results/K5/figures_EN")
RATE <- 4.876e-7
tr <- read.tree(file.path(D,"dated3.timetree.nwk"))
tp <- read.table(file.path(D,"taxa_pop.tsv"),sep="\t",col.names=c("ID","POP"))
p <- setNames(tp$POP, tp$ID); pop <- p[tr$tip.label]

# --- identificar el clado principal de PB1 (excluye las 3 cepas divergentes) ---
pb1 <- tr$tip.label[which(pop=="PB1")]
mr <- getMRCA(tr, pb1)
ch <- tr$edge[tr$edge[,1]==mr, 2]
side <- lapply(ch, function(k) if(k<=Ntip(tr)) tr$tip.label[k] else extract.clade(tr,k)$tip.label)
has3 <- sapply(side, function(s) any(pop[s]=="PB3"))
pb1_main <- intersect(side[[which(has3)[1]]], pb1)      # 59 cepas, hermanas de PB3
pb1_out  <- intersect(side[[which(!has3)[1]]], pb1)     # 3 cepas divergentes
cat("PB1 clado principal:",length(pb1_main)," | cepas fuera:",length(pb1_out),
    "(",paste(pb1_out,collapse=", "),")\n")

reps <- c(PA=tr$tip.label[which(pop=="PA")][1],
          PB1=pb1_main[1],
          PB2=tr$tip.label[which(pop=="PB2")][1],
          PB3=tr$tip.label[which(pop=="PB3")][1],
          PB4=tr$tip.label[which(pop=="PB4")][1],
          NoAm=tr$tip.label[which(pop=="NoAm")][1],
          Hol=tr$tip.label[which(pop=="Hol")][1])
pt <- keep.tip(tr, as.character(reps))
pt$tip.label <- setNames(names(reps), as.character(reps))[pt$tip.label]
# ramas de longitud 0 -> politomia real (radiacion rapida); evita nodos duplicados
pt <- di2multi(pt, tol = 1e-9)
cat("Topologia:", write.tree(pt), "\n")

dep <- node.depth.edgelength(pt); H <- max(dep[1:Ntip(pt)])
ages <- (H - dep) / RATE
nodeages <- ages[(Ntip(pt)+1):length(ages)]
cat("Edades de nodos (anos):", paste(round(nodeages),collapse=", "), "\n")

nn <- c(PA=66, PB1=length(pb1_main), PB2=95, PB3=96, PB4=32, NoAm=21, Hol=8)
crown <- read.table(file.path(D,"crown_ages_lsd2.tsv"),header=TRUE,sep="\t")
cr <- setNames(crown$crown_years, crown$POP)
cols <- c(PA="#F4D03F",PB1="#E74C3C",PB2="#2E86C1",PB3="#1ABC9C",
          PB4="#34495E",NoAm="#000000",Hol="#E91E63")

pdf(file.path(OUT,"28_population_timetree.pdf"), width=10.5, height=6)
par(mar=c(5,2,4,11), xpd=NA)
pt2 <- pt; pt2$edge.length <- pt$edge.length / RATE      # ramas en anos
plot(pt2, edge.width=3, label.offset=60, cex=1.15, tip.color=cols[pt2$tip.label],
     x.lim=c(-120, 3300))
lastPP <- get("last_plot.phylo", envir=.PlotPhyloEnv)
tl <- pt2$tip.label
for(i in seq_along(tl)){
  q <- tl[i]
  text(lastPP$xx[i]+430, lastPP$yy[i],
       sprintf("n=%d   crown %d yr", nn[q], round(cr[q])),
       adj=0, cex=0.78, col="grey35")
}
nodelabels(paste0(round(nodeages)," yr"), node=(Ntip(pt2)+1):(Ntip(pt2)+pt2$Nnode),
           frame="none", adj=c(1.1,-0.6), cex=0.88, col="#B03A2E", font=2)
axisPhylo(backward=TRUE, cex.axis=0.9)
mtext("Years before present", side=1, line=2.8, cex=1.05)
title(expression(paste(italic("S. eubayanus"),
      " population time tree (IQ-TREE3 + LSD2;  ",
      mu, " = 1.67e-10/bp/gen x 2920 gen/year)")), cex.main=1.05, line=1.5)
mtext(paste0("Node labels = divergence time. 3 PB1 strains fall outside the main PB1 clade ",
             "(minor paraphyly) and are excluded from the collapse."),
      side=1, line=4.1, cex=0.72, col="grey35")
dev.off()
cat("-> 28_population_timetree.pdf\n")
