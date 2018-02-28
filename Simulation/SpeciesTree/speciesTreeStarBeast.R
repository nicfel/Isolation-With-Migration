######################################################
######################################################
# Analyses inferred topology and node heights of
# species trees
######################################################
######################################################
library(ggplot2)
# needed to get the node heights
library(phytools)
# needed to read the trees
library(ape)

library(grid)
library(gridExtra)
library(coda)


# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all files with 1 gene
log <- list.files(path="./StarBeastLog", pattern="*.species.trees", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)

# Read In the file with the mutation rates
mut <- read.table("Sequences/relative_rates.txt", header=TRUE, sep="\t")



# Read In Data ---------------------------------
mig_rate <- data.frame()
counter <-1
for (i in seq(1,length(log),1)){
  print(i)
  
  # Make the filenames for all possible migration rates
  trees <- log[i]
  
  tmp <- strsplit(trees,'_')
  tmp2 <- strsplit(tmp[[1]][3], "[.]")
  tmp3 <- strsplit(tmp2[[1]][1], "[S]")
  rep_nr <- tmp2[[1]][1]
  
  
  t <- read.table(gsub("species.trees","log",log[i]), header=TRUE, sep="\t")
  
  # read the true species tree
  true_species_tree <- read.tree(file="","(((A:5,B:5):5,C:10):10,D:20);")
  # get node heights
  true_node_height <- abs(unique(nodeHeights((true_species_tree))[,1])-20)
  
  # read in the estimated species trees
  est_species_tree <- read.nexus(trees)
  
  true_topology <- 0
  
  dvar <- 1
  for (j in 100:length(est_species_tree)){
    # compare the tree topologies to the truth
    is_equal <- all.equal.phylo(true_species_tree,est_species_tree[[j]],
                                use.edge.length = FALSE, use.tip.label = TRUE)
    
    nh <- nodeHeights((est_species_tree[[j]]))
    
    enh <- abs(unique(nh[,1])-max(nh))
    
    if(is_equal){
      true_topology <- true_topology+1

      if (dvar==1){
        dvar<-2

        cond_true_top <- data.frame(n1=enh[1],n2=enh[2],n3=enh[3])
      }else{
        new.cond_true_top <- data.frame(n1=enh[1],n2=enh[2],n3=enh[3])
        cond_true_top <- rbind(cond_true_top, new.cond_true_top)
      }
    }
  }
  
  ess <- effectiveSize(t)
  if (min(ess[2:3])<100){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess[2:3]),trees))
  }else{
    if (counter==1){
      counter <- 2
      tree_infos <- data.frame(true_topology=true_topology, 
                               n1=median(cond_true_top$n1), 
                               n2=median(cond_true_top$n2),
                               n3=median(cond_true_top$n3),
                               fname=trees)
      mutation <- data.frame(true=mut[rep_nr,]$gene1, est=mean(t$mutationRate.gene1))
      for (gene_nr in seq(2,50)){
        gene_name <- paste("gene",gene_nr,sep="")
        gene_name_log <- paste("mutationRate.gene",gene_nr,sep="")
        new.mutation <- data.frame(true=mut[rep_nr,gene_name], est=mean(t[,gene_name_log])) 
        mutation <- rbind(mutation, new.mutation)
      }
    }else{
      new.tree_infos <- data.frame(true_topology=true_topology, 
                               n1=median(cond_true_top$n1), 
                               n2=median(cond_true_top$n2),
                               n3=median(cond_true_top$n3),
                               fname=trees)
      tree_infos <- rbind(tree_infos, new.tree_infos)
      for (gene_nr in seq(1,50)){
        gene_name <- paste("gene",gene_nr,sep="")
        gene_name_log <- paste("mutationRate.gene",gene_nr,sep="")
        new.mutation <- data.frame(true=mut[rep_nr,gene_name], est=mean(t[,gene_name_log])) 
        mutation <- rbind(mutation, new.mutation)
      }
      
    }
  }
}

write.table(tree_infos, file='StarBeast_param.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(mutation, file='StarBeast_mutation.tsv', quote=FALSE, sep='\t', col.names = NA)







