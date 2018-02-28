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
library("OutbreakTools")


# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all files with 1 gene
log <- list.files(path="./log", pattern="*.trees", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)

# Read In the file with the mutation rates
mut <- read.table("Sequences/relative_rates.txt", header=TRUE, sep="\t")

# include the function reading the node labels as well
source("read.nexus.nfm.R")


# Read In Data ---------------------------------
mig_rate <- data.frame()
counter <-1
nenr <- 1
for (i in seq(1,length(log),1)){
  print(i)
  
  # Make the filenames for all possible migration rates
  trees <- log[i]
  
  tmp <- strsplit(trees,'_')
  tmp2 <- strsplit(tmp[[1]][3], "[.]")
  tmp3 <- strsplit(tmp2[[1]][1], "[S]")
  rep_nr <- tmp2[[1]][1]
  
  # read log file
  t <- read.table(gsub("species.trees","log",log[i]), header=TRUE, sep="\t")

  # read the true species tree
  true_species_tree <- read.tree(file="","(((A:0.0025,B:0.0025):0.0025,C:0.005):0.005,D:0.01);")
  # get node heights
  true_node_height <- abs(unique(nodeHeights((true_species_tree))[,1])-0.002)
  
  # read in the estimated species trees
  est_species_tree_with_nodes <- read.annotated.nexus(trees)
  
  true_topology <- 0
  
  
  bpp_name <- gsub("log","BPP", trees)
  bpp_name <- gsub("IM_50Genes","BPP", bpp_name)
  bpp_name <- gsub(".species.trees",".txt", bpp_name)
  
  
  BPP_Input <- strsplit(paste(readLines(bpp_name), collapse=" "), split="#")
  Ne_string = strsplit(BPP_Input[[1]][2], split=" ")
  Ne_true = c(as.numeric(Ne_string[[1]][4]), as.numeric(Ne_string[[1]][5]), as.numeric(Ne_string[[1]][6]), as.numeric(Ne_string[[1]][7]),
              as.numeric(Ne_string[[1]][10]), as.numeric(Ne_string[[1]][9]), as.numeric(Ne_string[[1]][8]))
  
  m1 = strsplit(BPP_Input[[1]][3], split=" ")
  m2 = strsplit(BPP_Input[[1]][4], split=" ")
  m3 = strsplit(BPP_Input[[1]][5], split=" ")
  m4 = strsplit(BPP_Input[[1]][6], split=" ")
  m5 = strsplit(BPP_Input[[1]][7], split=" ")
  m6 = strsplit(BPP_Input[[1]][8], split=" ")
  m7 = strsplit(BPP_Input[[1]][9], split=" ")
  Ne_color = c("sampled","sampled","sampled","sampled","first ancestral","second ancestral","root")
  
  
  m_true = c(as.numeric(m2[[1]][2]), as.numeric(m3[[1]][2]), as.numeric(m4[[1]][2]),
             as.numeric(m1[[1]][3]), as.numeric(m3[[1]][3]), as.numeric(m4[[1]][3]),
             as.numeric(m1[[1]][4]), as.numeric(m2[[1]][4]), as.numeric(m4[[1]][4]),
             as.numeric(m1[[1]][5]), as.numeric(m2[[1]][5]), as.numeric(m3[[1]][5]),
             as.numeric(m7[[1]][4]), as.numeric(m7[[1]][5]), as.numeric(m3[[1]][8]), as.numeric(m4[[1]][8]),
             as.numeric(m6[[1]][5]), as.numeric(m4[[1]][7]))
  m_color = c("sampled","sampled","sampled","sampled","sampled","sampled","sampled","sampled","sampled","sampled","sampled","sampled",
              "first ancestral","first ancestral","first ancestral","first ancestral",
              "second ancestral","second ancestral")
  

  dvar <- 1
  for (j in 100:length(est_species_tree_with_nodes)){
    # compare the tree topologies to the truth
    is_equal <- all.equal.phylo(true_species_tree,est_species_tree_with_nodes[[j]],
                                use.edge.length = FALSE, use.tip.label = TRUE)
    
    if(is_equal){
      true_topology <- true_topology + 1
      
      annotated_tree <- est_species_tree_with_nodes[[j]]
      for (k in seq(1, length(annotated_tree$annotations))){
        if (annotated_tree$annotations[[k]]$species==0){
          A_nr = k
        }
        if (annotated_tree$annotations[[k]]$species==1){
          B_nr = k
        }
        if (annotated_tree$annotations[[k]]$species==2){
          C_nr = k
        }
        if (annotated_tree$annotations[[k]]$species==3){
          D_nr = k
        }
        if (annotated_tree$annotations[[k]]$species>3 && length(annotated_tree$annotations[[k]]$to)==2){
          AB_nr = k
        }
        if (annotated_tree$annotations[[k]]$species>3 && length(annotated_tree$annotations[[k]]$to)==1){
          ABC_nr = k
        }
      }
        
      # get the node name of every branch
      nodes = data.frame(A=A_nr,
                         B=B_nr,
                         C=C_nr,
                         D=D_nr,
                         AB=AB_nr,
                         ABC=ABC_nr)
      
      # get the correct Ne for every branch
      species_nr <- data.frame(A=annotated_tree$annotations[[nodes$A]]$species,
                               B=annotated_tree$annotations[[nodes$B]]$species,
                               C=annotated_tree$annotations[[nodes$C]]$species,
                               D=annotated_tree$annotations[[nodes$D]]$species,
                               AB=annotated_tree$annotations[[nodes$AB]]$species,
                               ABC=annotated_tree$annotations[[nodes$ABC]]$species,
                               root=annotated_tree$root.annotation$species)
      
      # get the correct Ne for every branch
      new.Ne.cond_true_top <- data.frame(A=annotated_tree$annotations[[nodes$A]]$Ne,
                                         B=annotated_tree$annotations[[nodes$B]]$Ne,
                                         C=annotated_tree$annotations[[nodes$C]]$Ne,
                                         D=annotated_tree$annotations[[nodes$D]]$Ne,
                                         AB=annotated_tree$annotations[[nodes$AB]]$Ne,
                                         ABC=annotated_tree$annotations[[nodes$ABC]]$Ne,
                                         root=annotated_tree$root.annotation$Ne)
      
      nh = nodeHeights(annotated_tree)
      new.height.cond_true_top <- data.frame(AB=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("A", "B"), type=c("node","height")))-max(nh)),
                                             ABC=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("A", "C"), type=c("node","height")))-max(nh)),
                                             root=max(nh))
      
      # node names to get height differences
      node_child_heights <- data.frame(A=0,
                                       B=0,
                                       C=0,
                                       D=0,
                                       AB=new.height.cond_true_top$AB,
                                       ABC=new.height.cond_true_top$ABC)
      
      # node names to get height differences
      node_names_for_heights <- data.frame(A="A",
                                           B="B",
                                           C="C",
                                           D="D",
                                           AB="A",
                                           ABC="A")
      
      # get all the routes of migration
      for (a in seq(1, length(annotated_tree$annotations))){
        for (b in seq(1, length(annotated_tree$annotations[[a]]$to))){
          
          to_ind = which(species_nr==annotated_tree$annotations[[a]]$to[[b]])
          
          # get the distance between the different nodes
          name1 <- names(nodes[which(nodes==a)])
          name2 <- names(species_nr[to_ind])
          
          max_childs <- max(node_child_heights[,name1], node_child_heights[,name2])
          if (name1=="D" || name2=="D"){
            par_height = max(nh)
          }else{
            par_height <- abs(nodeheight(annotated_tree, findMRCA(annotated_tree, tips=c(as.character(node_names_for_heights[,name1]), as.character(node_names_for_heights[,name2])), type=c("node","height")))-max(nh))
          }
          if (max_childs>par_height){
            fsdafsad
          }
          rel_rate <- annotated_tree$annotations[[a]]$rates[[b]]*(par_height-max_childs)
          if (rel_rate==0){
            fsd
          }
          
          new.migration = data.frame(from=name1, to=name2, rate=rel_rate)
          if (a == 1 && b == 1){
            new.migration.cond_true_top = new.migration
          }else{
            new.migration.cond_true_top = rbind(new.migration.cond_true_top, new.migration)
          }
        }
      }
      if (j == length(est_species_tree_with_nodes)){
        possible_routes = data.frame(from=new.migration.cond_true_top$from, to=new.migration.cond_true_top$to)
      }
      
      if (dvar==1){
        dvar<-2
        Ne.cond_true_top <- new.Ne.cond_true_top
        migration.cond_true_top <- new.migration.cond_true_top
        height.cond_true_top <- new.height.cond_true_top
      }else{
        Ne.cond_true_top <- rbind(Ne.cond_true_top, new.Ne.cond_true_top)
        migration.cond_true_top <- rbind(migration.cond_true_top, new.migration.cond_true_top)
        height.cond_true_top <- rbind(height.cond_true_top, new.height.cond_true_top)
      }
    }
  }
  

  all_names = c("A","B","C","D","AB","ABC","root")
  ess <- effectiveSize(t)
  if (min(ess[2:3])<50){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
  }else{
    new.tree_infos <- data.frame(true_topology=true_topology, 
                                 n1=median(height.cond_true_top$AB), 
                                 n2=median(height.cond_true_top$ABC),
                                 n3=median(height.cond_true_top$root),
                                 fname=trees)
    for (i in seq(1,length(Ne_true))){
      name <- all_names[i]
      rates.new <- data.frame(est=median(Ne.cond_true_top[, name]), true=Ne_true[i]/8, col=Ne_color[i], rate="Ne")
      if(nenr==1){
        rates = rates.new
        nenr <- 2
      }else{
        rates <- rbind(rates, rates.new)
      }
    }
    for (i in seq(1,length(m_true))){
      name <- paste("m",i,sep="")
      rates.new <- data.frame(est=median(t[, name]), true=m_true[i], col=m_color[i], rate="migration")
      rates <- rbind(rates, rates.new)
    }
    if (counter==1){
      counter <- 2
      tree_infos <- new.tree_infos
      mutation <- data.frame(true=mut[rep_nr,]$gene1, est=mean(t$clockRate.gene1))
      for (gene_nr in seq(2,50)){
        gene_name <- paste("gene",gene_nr,sep="")
        gene_name_log <- paste("clockRate.gene",gene_nr,sep="")
        new.mutation <- data.frame(true=mut[rep_nr,gene_name], est=mean(t[,gene_name_log])) 
        mutation <- rbind(mutation, new.mutation)
      }
    }else{
      tree_infos <- rbind(tree_infos, new.tree_infos)
      for (gene_nr in seq(1,50)){
        gene_name <- paste("gene",gene_nr,sep="")
        gene_name_log <- paste("clockRate.gene",gene_nr,sep="")
        new.mutation <- data.frame(true=mut[rep_nr,gene_name], est=mean(t[,gene_name_log])) 
        mutation <- rbind(mutation, new.mutation)
      }
    }
  }
}


write.table(tree_infos, file='AIM_nodeHeights.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(rates, file='AIM_rates.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(mutation, file='AIM_mutation.tsv', quote=FALSE, sep='\t', col.names = NA)

