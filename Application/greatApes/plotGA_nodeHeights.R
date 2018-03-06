library(graphics)
library(ape)
library("gridExtra")
library(ggplot2)
library(phytools)

# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)





tr <- read.nexus("mcc/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv_mcc.species.trees")
tr <- ladderize(tr, right = F)
plot(tr,show.tip.label=T, direction="downwards",edge.width=2,root.edge = T)



trees <- list.files(path="./mcc/greatapes", pattern="\\.trees$", full.names = TRUE)
for (i in seq(1,length(trees)-1)){
  print(i)
  g_tree = read.nexus(trees[i])
  for (j in seq(1, length(g_tree$tip.label))){
    new.type = data.frame(label = g_tree$tip.label[j], type=strsplit(g_tree$tip.label[j], split="_")[[1]][1])
    if (j == 1){
      type = new.type
    }else{
      type = rbind(type, new.type)
    }
  }
  type <- type[order(type$type),] 
  nh = nodeHeights(g_tree)

  # get the nodes that are coalescent events between types and their heights
  c <- 1
  for (a in seq(1, length(type$label))){
    if (a <length(type$label)){
      for (b in seq(a+1, length(type$label))){
        if (type$type[a]!=type$type[b]){
          n = findMRCA(g_tree, tips=c(type$label[a], type$label[b]), type=c("node","height"))
          new.i_node = data.frame(nr = n,
                                  t1 = type$type[a], 
                                  t2 = type$type[b], 
                                  height=abs(nodeheight(g_tree,n)-max(nh)),
                                  nodecolor="gene")
          if (c==1){
            i_node = new.i_node
            c = 2
          }else{
            i_node = rbind(i_node, new.i_node)
          }
        }
      }
    }
  }
  new.ancestral_nodes = unique(i_node)
  if (i == 1){
    ancestral_nodes = new.ancestral_nodes
  }else{
    ancestral_nodes = rbind(ancestral_nodes, new.ancestral_nodes)
  }
}

# get the mean speciation times
nh = nodeHeights(tr)
for (a in seq(1, length(tr$tip.label))){
  if (a < length(tr$tip.label)){
    for (b in seq(a+1, length(tr$tip.label))){
      if (tr$tip.label[a]!=tr$tip.label[b]){
        n = findMRCA(tr, tips=c(tr$tip.label[a], tr$tip.label[b]), type=c("node","height"))
        new.i_node = data.frame(nr = n,
                                t1 = tr$tip.label[a], 
                                t2 = tr$tip.label[b], 
                                height=abs(nodeheight(tr,n)-max(nh)),
                                nodecolor="species")
        if (c==1){
          i_node = new.i_node
          c = 2
        }else{
          i_node = rbind(i_node, new.i_node)
        }
      }
    }
  }
}
new.ancestral_nodes = unique(i_node)
ancestral_nodes = rbind(ancestral_nodes, new.ancestral_nodes)





species <- unique(type$type)
c = 1
# get all combinations of species
for (a in seq(1, length(species))){
  if (a < length(species)){
    for (b in seq(a+1, length(species))){
      h <- ancestral_nodes$height[which(ancestral_nodes$t1==species[a] & ancestral_nodes$t2==species[b])]
      for (i in seq(1,length(h))){
        new.h = data.frame(x = h[i], t1 = species[a], t2 = species[b])
        if (c==1){
          heights = new.h
          c = 2
        }else{
          heights = rbind(heights, new.h)
        }
      }
    }
  }
}

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
black <- rgb(red=0, green=0,blue=0)

plot_heights <- ggplot(data=ancestral_nodes) + 
  geom_point(aes(t1,height, color=nodecolor)) +
  facet_grid(. ~ t2, scales="free_x", space = "free") + 
  ylab(" substitutions per site")+
  xlab("")+
  scale_colour_manual("",values = c("species" = col1, "gene" = black)) + guides(colour = guide_legend(override.aes = list(size=2)))
plot(plot_heights)
ggsave(plot=plot_heights,"../../text/figures/greatapes/GreatApes_NodeHeights.pdf",width=10, height=8)



