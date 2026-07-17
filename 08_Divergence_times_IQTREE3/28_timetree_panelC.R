#!/usr/bin/env Rscript
# Figura estilo "panel C": arbol poblacional datado con eje en Kya, sombreado de epocas
# (Pleistoceno tardio / Holoceno) y RANGO de tiempo en cada nodo.
# El rango sigue la convencion de la figura previa del paper: extremos por generaciones/ano
#   lento  g=365  (tiempos mas antiguos)   <->   rapido g=2920 (tiempos mas recientes)
# Uso: 28_timetree_panelC.R [escala]   escala = "slow" (default, eje a g=365) o "fast"
suppressMessages({library(ape)})
proj <- "/Volumes/Extreme SSD/Eubayanus-Pop"
D <- file.path(proj,"results/K5/dating"); OUT <- file.path(proj,"results/K5/figures_EN")
RATE <- 4.876e-7          # g=2920
G_FAST <- 2920; G_SLOW <- 365
scale_mode <- ifelse(length(commandArgs(TRUE))>0, commandArgs(TRUE)[1], "slow")
HOLOCENE <- 11.7          # Kya, inicio del Holoceno

tr <- read.tree(file.path(D,"dated3.timetree.nwk"))
tp <- read.table(file.path(D,"taxa_pop.tsv"),sep="\t",col.names=c("ID","POP"))
p <- setNames(tp$POP, tp$ID); pop <- p[tr$tip.label]

# clado principal de PB1 (3 cepas quedan fuera)
pb1 <- tr$tip.label[which(pop=="PB1")]; mr <- getMRCA(tr,pb1)
ch <- tr$edge[tr$edge[,1]==mr,2]
side <- lapply(ch, function(k) if(k<=Ntip(tr)) tr$tip.label[k] else extract.clade(tr,k)$tip.label)
has3 <- sapply(side, function(s) any(pop[s]=="PB3"))
pb1_main <- intersect(side[[which(has3)[1]]], pb1)

reps <- c(PA=tr$tip.label[which(pop=="PA")][1], PB1=pb1_main[1],
          PB2=tr$tip.label[which(pop=="PB2")][1], PB3=tr$tip.label[which(pop=="PB3")][1],
          PB4=tr$tip.label[which(pop=="PB4")][1], NoAm=tr$tip.label[which(pop=="NoAm")][1],
          Hol=tr$tip.label[which(pop=="Hol")][1])
pt <- keep.tip(tr, as.character(reps))
pt$tip.label <- setNames(names(reps), as.character(reps))[pt$tip.label]
# NO se colapsa a politomia: la topologia ML esta RESUELTA (ramas internas 0.00045-0.00047 subs
# en el arbol de sustituciones). LSD2 aplasta algunas a longitud 0 al forzar ultrametricidad,
# asi que varios nodos datan al MISMO tiempo (radiacion rapida) pero el ORDEN si esta resuelto:
# PA se separa primero, luego NoAm, luego el resto.

# edades en anos a g=2920 y a g=365
dep <- node.depth.edgelength(pt); H <- max(dep[1:Ntip(pt)])
age_fast <- (H - dep)/RATE                      # g=2920
age_slow <- age_fast * (G_FAST/G_SLOW)          # g=365  (x8)
nid <- (Ntip(pt)+1):(Ntip(pt)+pt$Nnode)

# escala del eje (Kya)
KY_TRUE <- if(scale_mode=="slow") age_slow/1000 else age_fast/1000   # edades REALES en Kya
H_KY <- max(KY_TRUE[(Ntip(pt)+1)])                                    # ~altura (raiz)
H_KY <- KY_TRUE[Ntip(pt)+1]

# --- SOLO PARA EL DIBUJO: separar los nodos que datan al mismo tiempo ---
# LSD2 les da longitud 0 (radiacion no resoluble en tiempo). Se les impone una separacion
# minima para que no se superpongan. Las EDADES REPORTADAS en las etiquetas siguen siendo
# las reales; la linea punteada marca el tiempo verdadero de esos nodos.
MINSEP <- H_KY * 0.035
disp <- KY_TRUE
ord <- reorder(pt, "cladewise")$edge          # preorden: padres antes que hijos
for(i in seq_len(nrow(ord))){
  par_i <- ord[i,1]; ch_i <- ord[i,2]
  if(ch_i > Ntip(pt)){                        # solo ramas internas
    if(disp[par_i] - disp[ch_i] < MINSEP) disp[ch_i] <- disp[par_i] - MINSEP
  }
}
disp[1:Ntip(pt)] <- 0                          # puntas en el presente
pt2 <- pt
pt2$edge.length <- disp[pt$edge[,1]] - disp[pt$edge[,2]]
XMAX <- max(node.depth.edgelength(pt2))

cols <- c(PA="#F4D03F",PB1="#E74C3C",PB2="#2E86C1",PB3="#1ABC9C",
          PB4="#34495E",NoAm="#000000",Hol="#E91E63")

NT <- Ntip(pt2)
XL <- c(-XMAX*0.30, XMAX*1.30)      # espacio a la izquierda para la etiqueta de la raiz
YL <- c(0.4, NT + 1.3)              # espacio arriba para los rotulos de epoca
pdf(file.path(OUT,"29_population_timetree_panelC.pdf"), width=11.5, height=6.4)
par(mar=c(5,1,3,9), xpd=NA)
plot(pt2, edge.width=2.2, label.offset=XMAX*0.02, cex=1.05,
     tip.color=cols[pt2$tip.label], x.lim=XL, y.lim=YL, plot=FALSE)
# --- sombreado de epocas (x=0 es la raiz, x=XMAX el presente) ---
x_hol <- XMAX - HOLOCENE            # posicion del limite del Holoceno
ytop <- NT + 0.55
if(x_hol > 0){
  rect(0, YL[1], x_hol, ytop, col="#D6EAF0", border=NA)         # Pleistoceno tardio
  rect(x_hol, YL[1], XMAX, ytop, col="#DCDCDC", border=NA)      # Holoceno
  text(x_hol/2, NT+1.0, "Late Pleistocene", cex=0.9, col="grey25")
  text(x_hol+(XMAX-x_hol)/2, NT+1.0, "Holocene", cex=0.9, col="grey25")
} else {
  rect(0, YL[1], XMAX, ytop, col="#DCDCDC", border=NA)
  text(XMAX/2, NT+1.0, "Holocene", cex=0.9, col="grey25")
}
par(new=TRUE)
plot(pt2, edge.width=2.2, label.offset=XMAX*0.02, cex=1.05,
     tip.color=cols[pt2$tip.label], x.lim=XL, y.lim=YL)
# --- linea punteada en el tiempo REAL de los nodos simultaneos (antes de expandir) ---
key <- round(age_fast[nid])
dup <- key[duplicated(key)][1]
if(!is.na(dup)){
  x_true <- XMAX - KY_TRUE[nid[which(key==dup)[1]]]   # posicion del tiempo verdadero
  segments(x_true, YL[1]+0.2, x_true, NT+0.35, lty=3, col="grey45", lwd=1.2)
}
# --- nodos separados para legibilidad; las etiquetas llevan la edad REAL ---
nodelabels(node=nid, pch=19, cex=1.3, col="black")
lab <- sprintf("%.2f - %.2f KYA", age_slow[nid]/1000, age_fast[nid]/1000)
pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
xx <- pp$xx; yy <- pp$yy
# Etiquetas: por defecto arriba-izquierda del nodo. Si dos nodos quedan muy juntos en
# vertical, la del indice par se pone abajo-derecha, usando la separacion horizontal.
for(k in seq_along(nid)){
  n <- nid[k]
  crowded <- any(abs(yy[nid[-k]] - yy[n]) < 0.9)
  if(crowded && k %% 2 == 0){
    text(xx[n]+XMAX*0.012, yy[n]-0.42, lab[k], adj=c(0,0.5), cex=0.75, col="black")
  } else {
    text(xx[n]-XMAX*0.012, yy[n]+0.42, lab[k], adj=c(1,0.5), cex=0.75, col="black")
  }
}
axisPhylo(backward=TRUE, cex.axis=0.9)
mtext("Time (Kya)", side=1, line=2.7, cex=1.05)
title(expression(paste(italic("S. eubayanus"), " population divergence times (IQ-TREE3 + LSD2)")),
      cex.main=1.1, line=1.2)
mtext("Node ranges: slow clock (365 gen/yr) - fast clock (2920 gen/yr). Dotted line = true age of the 3 simultaneous nodes (expanded here for legibility).",
      side=1, line=3.9, cex=0.75, col="grey30")
dev.off()
cat("-> 29_population_timetree_panelC.pdf  (escala:",scale_mode,")\n")
cat("Nodos (slow - fast KYA):\n")
for(i in seq_along(nid)) cat("  ", lab[i], "\n")
