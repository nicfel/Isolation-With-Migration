library(graphics)
library(ape)
library("gridExtra")
library(ggplot2)
library(phytools)

# clear workspace
rm(list = ls())
tr <- read.nexus("mcc/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv_.species.trees")
tr <- ladderize(tr, right = F)


require(grDevices)


pdf("../../text/figures/GreatApes_SpeciesTree.pdf",width=6,height=5,paper='special') 

## set up the plot region:
op <- par(bg = "white")
plot(c(0, 10), c(-0.0003, 0.013), type = "n", xlab = "", ylab = "", axes = F)

scaler = 1500


nh = nodeHeights(tr)
absnh = sort(abs(unique(nh[,1]-max(nh))))


node_locations = data.frame(x = c(1,2,3,4,5,6, 5.5, 4.75, 3.875, 2.9375, 1.9688),
                            ystart = c(0,0,0,0,0,0, absnh[1],absnh[2],absnh[3],absnh[4],absnh[5]),
                            yend = c(absnh[5],absnh[4],absnh[3],absnh[2],absnh[1],absnh[1], absnh[2],absnh[3],absnh[4],absnh[5],0.01),
                            Ne = c(2.1714062911164328E-4, 1.817155688876022E-4, 9.169949697418392E-5, 9.334923676685362E-5, 3.745617301163301E-4,1.4657067038628102E-4,
                                   1.4657067038628102E-4,1.4657067038628102E-4,1.4657067038628102E-4,1.4657067038628102E-4, 1.4657067038628102E-4))
rect(node_locations$x-node_locations$Ne*scaler, node_locations$ystart, node_locations$x+node_locations$Ne*scaler, node_locations$yend,
     lwd=1)
j <- seq(0,4)
segments(node_locations$x[1+j]-node_locations$Ne[1+j]*scaler, absnh[5-j],
         node_locations$x[10-j]+node_locations$Ne[10-j]*scaler, absnh[5-j],
         lwd=1)

k <- seq(3,9)

node_locations$yend[9] = 0.0015

node_locations$x = node_locations$x+3
node_locations$ystart = node_locations$ystart*5+0.0045
node_locations$yend = node_locations$yend*5+0.0045

rect(node_locations$x[k]-node_locations$Ne[k]*scaler, node_locations$ystart[k], node_locations$x[k]+node_locations$Ne[k]*scaler, node_locations$yend[k],
     lwd=2)
j <- seq(2,4)

segments(node_locations$x[1+j]-node_locations$Ne[1+j]*scaler, absnh[5-j]*5+0.0045,
         node_locations$x[10-j]+node_locations$Ne[10-j]*scaler, absnh[5-j]*5+0.0045,
         lwd=2)
axis(2, pos=0.25, at=c(0,0.005,0.01))
axis(2, pos=5.25, at=c(0.0045,0.0085,0.0125), labels=c("0.0", "0.0008", "0.0016"))

text(x=1,y=-0.0003,"Gorilla", cex=0.6)
text(x=2,y=-0.0003,"Human", cex=0.6)
text(x=3,y=-0.0003,"Bonobo", cex=0.6)
text(x=5,y=-0.0003,"Chimpanzee Species", cex=0.6)

text(x=6,y=0.0042,"Bonobo", cex=0.6)
text(x=7,y=0.0042,"Western", cex=0.6)
text(x=8,y=0.0042,"Central", cex=0.6)
text(x=9,y=0.0042,"Eastern", cex=0.6)


par(op)
dev.off()
