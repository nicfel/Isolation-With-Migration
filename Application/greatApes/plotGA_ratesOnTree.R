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
library(coda)
library(grid)
library(gridExtra)
library("OutbreakTools")

# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# get the names of all files with 1 gene
log <- list.files(path="./mcc", pattern="*.trees", full.names = TRUE)


# include the function reading the node labels as well
source("read.nexus.nfm.R")


# Read In Data ---------------------------------
mig_rate <- data.frame()
counter <- 1

i = 1

# Make the filenames for all possible migration rates
trees <- "./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.species.trees"

tmp <- strsplit(trees,'_')
tmp2 <- strsplit(tmp[[1]][3], "[.]")
tmp3 <- strsplit(tmp2[[1]][1], "[S]")
rep_nr <- tmp2[[1]][1]

# read log file
#t <- read.table(gsub("species.trees","log",log[i]), header=TRUE, sep="\t")

# read the true species tree
true_species_tree <- read.tree(file="","(((((Pts:2.0616386643443712E-4,Ptt:2.0616386643443712E-4):5.620471937279749E-5,Ptv:2.623685858072346E-4):7.145726484084223E-4,Ppa:9.76941234215657E-4):0.005169704744388228,Hsa:0.006146645978603885):0.001888109949209401,Ggg:0.008034755927813286);")
# get node heights
true_node_height <- abs(unique(nodeHeights((true_species_tree))[,1])-0.002)

est_species_tree_with_nodes <- read.annotated.nexus(trees)
true_topology <- 0

tree_has_true_topology = c()


dvar <- 1
for (j in 1:length(est_species_tree_with_nodes)){
  print(j)
  # compare the tree topologies to the truth
  is_equal <- all.equal.phylo(true_species_tree,est_species_tree_with_nodes[[j]],
                              use.edge.length = FALSE, use.tip.label = TRUE)
  
  if(is_equal){
    tree_has_true_topology[j] = 1
  }else{
    tree_has_true_topology[j] = 0
  }
  
  
  if(is_equal){
    annotated_tree <- est_species_tree_with_nodes[[j]]
    
    
    for (k in seq(1, length(annotated_tree$annotations))){
      if (annotated_tree$annotations[[k]]$species==0){
        ggg_nr = k
      }
      if (annotated_tree$annotations[[k]]$species==1){
        hsa_nr = k
      }
      if (annotated_tree$annotations[[k]]$species==2){
        ppa_nr = k
      }
      if (annotated_tree$annotations[[k]]$species==3){
        pts_nr = k
      }
      if (annotated_tree$annotations[[k]]$species==4){
        ptt_nr = k
      }
      if (annotated_tree$annotations[[k]]$species==5){
        ptv_nr = k
      }
      if (annotated_tree$annotations[[k]]$species>5 && length(annotated_tree$annotations[[k]]$to)==4){
        ptts_nr = k
      }
      if (annotated_tree$annotations[[k]]$species>5 && length(annotated_tree$annotations[[k]]$to)==3){
        chimp_nr = k
      }
      if (annotated_tree$annotations[[k]]$species>5 && length(annotated_tree$annotations[[k]]$to)==2){
        cbo_nr = k
      }
      if (annotated_tree$annotations[[k]]$species>5 && length(annotated_tree$annotations[[k]]$to)==1){
        cbh_nr = k
      }
      
      
    }
    
    # get the node name of every branch
    nodes = data.frame(Ptts=ptts_nr,
                      chimp=chimp_nr,
                      cbo=cbo_nr,
                      cbh=cbh_nr,
                      Ggg=ggg_nr,
                      Hsa=hsa_nr,
                      Ppa=ppa_nr,
                      Ptv=ptv_nr,
                      Pts=pts_nr,
                      Ptt=ptt_nr)
    

    # get the correct Ne for every branch
    species_nr <- data.frame(Ggg=annotated_tree$annotations[[nodes$Ggg]]$species,
                                   Hsa=annotated_tree$annotations[[nodes$Hsa]]$species,
                                   Ppa=annotated_tree$annotations[[nodes$Ppa]]$species,
                                   Ptv=annotated_tree$annotations[[nodes$Ptv]]$species,
                                   Pts=annotated_tree$annotations[[nodes$Pts]]$species,
                                   Ptt=annotated_tree$annotations[[nodes$Ptt]]$species,
                                   Ptts=annotated_tree$annotations[[nodes$Ptts]]$species,
                                   chimp=annotated_tree$annotations[[nodes$chimp]]$species,
                                   cbo=annotated_tree$annotations[[nodes$cbo]]$species,
                                   cbh=annotated_tree$annotations[[nodes$cbh]]$species,
                                   root=annotated_tree$root.annotation$species)
    
    
    
    
    
    # get the correct Ne for every branch
    new.Ne.cond_true_top <- data.frame(Ggg=annotated_tree$annotations[[nodes$Ggg]]$Ne,
                                   Hsa=annotated_tree$annotations[[nodes$Hsa]]$Ne,
                                   Ppa=annotated_tree$annotations[[nodes$Ppa]]$Ne,
                                   Ptv=annotated_tree$annotations[[nodes$Ptv]]$Ne,
                                   Pts=annotated_tree$annotations[[nodes$Pts]]$Ne,
                                   Ptt=annotated_tree$annotations[[nodes$Ptt]]$Ne,
                                   Ptts=annotated_tree$annotations[[nodes$Ptts]]$Ne,
                                   chimp=annotated_tree$annotations[[nodes$chimp]]$Ne,
                                   cbo=annotated_tree$annotations[[nodes$cbo]]$Ne,
                                   cbh=annotated_tree$annotations[[nodes$cbh]]$Ne,
                                   root=annotated_tree$root.annotation$Ne)
    
    nh = nodeHeights(annotated_tree)
    new.height.cond_true_top <- data.frame(Ptts=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("Pts", "Ptt"), type=c("node","height")))-max(nh)),
                                          chimp=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("Pts", "Ptv"), type=c("node","height")))-max(nh)),
                                          cbo=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("Pts", "Ppa"), type=c("node","height")))-max(nh)),
                                          cbh=abs(nodeheight(annotated_tree,findMRCA(annotated_tree, tips=c("Pts", "Hsa"), type=c("node","height")))-max(nh)),
                                          root=max(nh))
    
    # node names to get height differences
    node_child_heights <- data.frame(Ggg=0,
                                         Hsa=0,
                                         Ppa=0,
                                         Ptv=0,
                                         Pts=0,
                                         Ptt=0,
                                         Ptts=new.height.cond_true_top$Ptts,
                                         chimp=new.height.cond_true_top$chimp,
                                         cbo=new.height.cond_true_top$cbo,
                                         cbh=new.height.cond_true_top$cbh)
    
    
    # node names to get height differences
    node_names_for_heights <- data.frame(Ggg="Ggg",
                                         Hsa="Hsa",
                                         Ppa="Ppa",
                                         Ptv="Ptv",
                                         Pts="Pts",
                                         Ptt="Ptt",
                                         Ptts="Pts",
                                         chimp="Pts",
                                         cbo="Pts",
                                         cbh="Pts")
    
    
    # get all the routes of migration
    for (a in seq(1, length(annotated_tree$annotations))){
      for (b in seq(1, length(annotated_tree$annotations[[a]]$to))){
        
        to_ind = which(species_nr==annotated_tree$annotations[[a]]$to[[b]])
        
        # get the distance between the different nodes
        name1 <- names(nodes[which(nodes==a)])
        name2 <- names(species_nr[to_ind])
        
        max_childs <- max(node_child_heights[,name1], node_child_heights[,name2])
        if (name1=="Ggg" || name2=="Ggg"){
          par_height = max(nh)
        }else{
          par_height <- abs(nodeheight(annotated_tree, findMRCA(annotated_tree, tips=c(as.character(node_names_for_heights[,name1]), as.character(node_names_for_heights[,name2])), type=c("node","height")))-max(nh))
        }
        if (max_childs>par_height){
          fsdafsad
        }
        rel_rate <- annotated_tree$annotations[[a]]$rates[[b]]*(par_height-max_childs)

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


# get the x position for each node
n.x <- data.frame(Ggg=node.height(true_species_tree)[which(true_species_tree$tip.label=="Ggg")],
                  Hsa=node.height(true_species_tree)[which(true_species_tree$tip.label=="Hsa")],
                  Ppa=node.height(true_species_tree)[which(true_species_tree$tip.label=="Ppa")],
                  Ptv=node.height(true_species_tree)[which(true_species_tree$tip.label=="Ptv")],
                  Pts=node.height(true_species_tree)[which(true_species_tree$tip.label=="Pts")],
                  Ptt=node.height(true_species_tree)[which(true_species_tree$tip.label=="Ptt")],
                  Ptts=node.height(true_species_tree)[findMRCA(annotated_tree, tips=c("Pts", "Ptt"), type=c("node","height"))],
                  chimp=node.height(true_species_tree)[findMRCA(annotated_tree, tips=c("Pts", "Ptv"), type=c("node","height"))],
                  cbo=node.height(true_species_tree)[findMRCA(annotated_tree, tips=c("Pts", "Ppa"), type=c("node","height"))],
                  cbh=node.height(true_species_tree)[findMRCA(annotated_tree, tips=c("Pts", "Hsa"), type=c("node","height"))],
                  root=node.height(true_species_tree)[findMRCA(annotated_tree, tips=c("Pts", "Ggg"), type=c("node","height"))])
n.x <- abs(n.x-7)
                  
n.y.start <- data.frame(Ggg=0,
                  Hsa=0,
                  Ppa=0,
                  Ptv=0,
                  Pts=0,
                  Ptt=0,
                  Ptts=median(height.cond_true_top$Ptts),
                  chimp=median(height.cond_true_top$chimp),
                  cbo=median(height.cond_true_top$cbo),
                  cbh=median(height.cond_true_top$cbh),
                  root=median(height.cond_true_top$root))

n.y.end <- data.frame(Ggg=n.y.start$root,
                        Hsa=n.y.start$cbh,
                        Ppa=n.y.start$cbo,
                        Ptv=n.y.start$chimp,
                        Pts=n.y.start$Ptts,
                        Ptt=n.y.start$Ptts,
                        Ptts=n.y.start$chimp,
                        chimp=n.y.start$cbo,
                        cbo=n.y.start$cbh,
                        cbh=n.y.start$root,
                        root=n.y.start$root*1.25)
# get the sum of mean Ne's
scaler = n.y.end$root/30;


species_colors = c("Ggg", "Hsa", "Ppa", "Ptv", "Pts", "Ptt", "Ptts", "chimp", "cbo", "cbh", "root")



nam = names(n.y.end)
for (i in seq(1,length(nam))){
  new.branch <- data.frame(off= median(Ne.cond_true_top[,nam[i]])/2,
                          x1=(n.x[,nam[i]]*scaler),
                           x2=(n.x[,nam[i]]*scaler),
                           ystart=n.y.start[,nam[i]],
                           yend=n.y.end[,nam[i]])
  new.rect.data <- data.frame(x1=(n.x[,nam[i]]*scaler - median(Ne.cond_true_top[,nam[i]])/2),
                              x2=(n.x[,nam[i]]*scaler + median(Ne.cond_true_top[,nam[i]])/2),
                              ystart=n.y.start[,nam[i]],
                              yend=n.y.end[,nam[i]], 
                              color = nam[i])
  new.rect.data.lower <- data.frame(x1=(n.x[,nam[i]]*scaler - quantile(Ne.cond_true_top[,nam[i]], 0.025)/2),
                              x2=(n.x[,nam[i]]*scaler + quantile(Ne.cond_true_top[,nam[i]], 0.025)/2),
                              ystart=n.y.start[,nam[i]],
                              yend=n.y.end[,nam[i]], 
                              color =  nam[i])
  new.rect.data.upper <- data.frame(x1=(n.x[,nam[i]]*scaler - quantile(Ne.cond_true_top[,nam[i]], 0.975)/2),
                                    x2=(n.x[,nam[i]]*scaler + quantile(Ne.cond_true_top[,nam[i]], 0.975)/2),
                                    ystart=n.y.start[,nam[i]],
                                    yend=n.y.end[,nam[i]], 
                                    color =  nam[i])
  if (i==1){
    branch <- new.branch
    rect.data <- new.rect.data
    rect.data.lower <- new.rect.data.lower
    rect.data.upper <- new.rect.data.upper
  }else{
    branch <- rbind(branch, new.branch)
    rect.data <- rbind(rect.data, new.rect.data)
    rect.data.lower <- rbind(rect.data.lower, new.rect.data.lower)
    rect.data.upper <- rbind(rect.data.upper, new.rect.data.upper)
  }
}

# add the nodes to the branches
uni <- sort(unique(branch$yend))
for (i in seq(1,length(uni)-1)){
  vals = which(branch$yend==uni[i])
  if (branch[vals[1],"x1"] < branch[vals[2],"x1"]){
    new.edge <- data.frame(x1=branch[vals[1],"x1"]-branch[vals[1],"off"],
                           x2=branch[vals[2],"x1"]+branch[vals[2],"off"],
                           ystart=uni[i],
                           yend=uni[i])
  }else{
    new.edge <- data.frame(x1=branch[vals[2],"x1"]-branch[vals[2],"off"],
                           x2=branch[vals[1],"x1"]+branch[vals[1],"off"],
                           ystart=uni[i],
                           yend=uni[i])
  }
  if (i==1){
    edge <- new.edge
  }else{
    edge <- rbind(edge, new.edge)
  }
}

dir_exist = FALSE
# build the migration rate arrows
for (i in seq(1, length(possible_routes$from))){
  name1 = as.character(possible_routes[i,1])
  name2 = as.character(possible_routes[i,2])
  start = max(n.y.start[,name1],n.y.start[,name2])
  end = min(n.y.end[,name1],n.y.end[,name2])
  y = (start+end)/2
  
  if (name1=="Ptts"){
    names1 = list(c("Pts", "Ptt", "Ptts"), c("Ptts"))
  }else if(name1=="chimp"){
    names1 = list(c("Pts", "Ptt", "Ptv", "Ptts", "chimp"),c("chimp"))
  }else if(name1=="cbo"){
    names1 = list(c("Pts", "Ptt", "Ptv", "Ppa", "Ptts", "chimp", "cbo"), c("cbo"))
  }else if(name1=="cbh"){
    names1 =  list(c("Pts", "Ptt", "Ptv", "Ppa", "Hsa", "Ptts", "chimp", "cbo", "cbh"), c("cbh"))
  }else{
    names1 = list(c(name1), c(name1))
  }
  
  if (name2=="Ptts"){
    names2 = list(c("Pts", "Ptt", "Ptts"), c("Ptts"))
  }else if(name2=="chimp"){
    names2 = list(c("Pts", "Ptt", "Ptv", "Ptts", "chimp"),c("chimp"))
  }else if(name2=="cbo"){
    names2 = list(c("Pts", "Ptt", "Ptv", "Ppa", "Ptts", "chimp", "cbo"), c("cbo"))
  }else if(name2=="cbh"){
    names2 =  list(c("Pts", "Ptt", "Ptv", "Ppa", "Hsa", "Ptts", "chimp", "cbo", "cbh"), c("cbh"))
  }else{
    names2 = list(c(name2),c(name2))
  }
  
  for (l in seq(1, length(names1))){
    for (a in seq(1, length(names1[[l]]))){
      for (b in seq(1, length(names2[[l]]))){
        # consider the direction
        new.rates <- migration.cond_true_top[intersect(which(migration.cond_true_top$from==names1[[l]][a]),(which(migration.cond_true_top$to==names2[[l]][b]))), "rate"]
        # don't consider the direction
        rate.dir1 <- migration.cond_true_top[intersect(which(migration.cond_true_top$from==names1[[l]][a]),(which(migration.cond_true_top$to==names2[[l]][b]))), "rate"]
        rate.dir2 <- migration.cond_true_top[intersect(which(migration.cond_true_top$from==names2[[l]][b]),(which(migration.cond_true_top$to==names1[[l]][a]))), "rate"]
        
        new.rate.dir = rate.dir1 + rate.dir2
        
        if (a==1 && b==1){
          if (l==1){
            old.rates = list()
            old.rates.dir = list()
          }
          old.rates[[l]] = new.rates
          old.rates.dir[[l]] = new.rate.dir
        }else{
          old.rates[[l]] = old.rates[[l]] + new.rates
          old.rates.dir[[l]] = old.rates.dir[[l]] + new.rate.dir
        }
      }
    }
    if (l==1){
      rate = list()
      rate.dir = list()
    }
    rate[[l]] =  sum(old.rates[[l]]>0)/length(old.rates[[l]])
    rate.dir[[l]] = sum(old.rates.dir[[l]]>0)/length(old.rates.dir[[l]]) #/ (0.95^(2*length(names1)*length(names2)))
    
  }

  if (n.x[,name1]>n.x[,name2]){
    new.migration_rate <- data.frame(x1=n.x[,name1]*scaler - median(Ne.cond_true_top[,name1])/2, x2=n.x[,name2]*scaler + median(Ne.cond_true_top[,name2])/2, y1=y, y2=y, 
                                     rate.clade=rate[[1]] , rate = rate[[2]],  
                                     prior=(0.95^(length(names1)*length(names2))), from=name1, to=name2, color = species_colors[species_nr[,name1]+1])
    new.migration_rate.dir <- data.frame(x1=n.x[,name1]*scaler - median(Ne.cond_true_top[,name1])/2, x2=n.x[,name2]*scaler + median(Ne.cond_true_top[,name2])/2, y1=y, y2=y, 
                                         rate.clade=rate.dir[[1]], rate = rate.dir[[2]], 
                                         prior=(0.95^(2*length(names1[[1]])*length(names2[[1]]))), from=name1, to=name2, color = species_colors[species_nr[,name1]+1])
  }else{
    new.migration_rate <- data.frame(x1=n.x[,name1]*scaler + median(Ne.cond_true_top[,name1])/2, x2=n.x[,name2]*scaler - median(Ne.cond_true_top[,name2])/2, y1=y, y2=y, 
                                     rate.clade=rate[[1]], rate = rate[[2]], 
                                     prior=(0.95^(length(names1[[1]])*length(names2[[1]]))), from=name1, to=name2, color = species_colors[species_nr[,name1]+1])
  }
  if (i==1){
    migration_rate <- new.migration_rate
  }else{
    migration_rate <- rbind(migration_rate, new.migration_rate)
    if (n.x[,name1]>n.x[,name2]){
      if (dir_exist){
        migration_rate.dir <- rbind(migration_rate.dir, new.migration_rate.dir)
      }else{
        migration_rate.dir <- new.migration_rate.dir
        dir_exist=TRUE
      }
    }
  }
}
#migration_rate$rate <- migration_rate$rate/mean(migration_rate$rate)-1







chimp.rect.data <- rect.data[-c(1,2,10,11), ]
chimp.branch <- branch[-c(1,2,10,11), ]

chimp.rect.data[7,4] <- 1.4*chimp.rect.data[7,3] 
chimp.branch[7,5] <- 1.4*chimp.branch[7,4] 

chimp.edge <- edge[-c(4,5),]

chimp.migration_rate <- migration_rate

offset_y = -0.003
scale_y = 2

chimp.rect.data$ystart = chimp.rect.data$ystart*scale_y-offset_y
chimp.rect.data$yend = chimp.rect.data$yend*scale_y-offset_y
chimp.branch$ystart = chimp.branch$ystart*scale_y-offset_y
chimp.branch$yend = chimp.branch$yend*scale_y-offset_y
chimp.edge$ystart = chimp.edge$ystart*scale_y-offset_y
chimp.edge$yend = chimp.edge$yend*scale_y-offset_y
chimp.migration_rate$y1 = chimp.migration_rate$y1*scale_y-offset_y
chimp.migration_rate$y2 = chimp.migration_rate$y2*scale_y-offset_y


offset_x = -0.0001
scale_x = 2

chimp.rect.data$x1 = chimp.rect.data$x1*scale_x-offset_x
chimp.rect.data$x2 = chimp.rect.data$x2*scale_x-offset_x
chimp.branch$x1 = chimp.branch$x1*scale_x-offset_x
chimp.branch$x2 = chimp.branch$x2*scale_x-offset_x
chimp.edge$x1 = chimp.edge$x1*scale_x-offset_x
chimp.edge$x2 = chimp.edge$x2*scale_x-offset_x
chimp.migration_rate$x1 = chimp.migration_rate$x1*scale_x-offset_x
chimp.migration_rate$x2 = chimp.migration_rate$x2*scale_x-offset_x








color_hex = c("#4137C2", "#4066CF", "#4B8DC2", "#5DA8A3", "#77B67F", "#96BD60", "#B8BC4B", "#D4B13F", "#E59638", "#E4672F", "#DC2F24")

color_hex = color_hex[c(10,8,6,4,1,2,3,5,7,9,11)]

data.axis = data.frame(x1=c(0,0.0017),x2=c(0,0.0022),ystart=c(0,0.001), yend=c(0.01,0.001))
data.y.ticks = data.frame(x1=c(-1,-1,-1,-1,-1)*0.00005,x2=c(0,0,0,0,0),ystart=c(0,0.0025,0.005,0.0075,0.01), yend=c(0,0.0025,0.005,0.0075,0.01))
data.x.ticks = data.frame(x1=c(0,2.5,5)*0.0001+0.0017,x2=c(0,2.5,5)*0.0001+0.0017,ystart=c(1,1,1)*0.00005+0.001, yend=c(1,1,1)*0+0.001)

data.axis.chimp = data.frame(x1=c(0.00185,0.0035),x2=c(0.00185,0.004),ystart=c(0,0.0008)*scale_y-offset_y, yend=c(0.0015,0.0008)*scale_y-offset_y)
data.y.ticks.chimp = data.frame(x1=c(-1,-1,-1,-1)*0.00005+0.00185,x2=c(0,0,0,0)+0.00185,ystart=c(0,0.0005,0.001,0.0015)*scale_y-offset_y, yend=c(0,0.0005,0.001,0.0015)*scale_y-offset_y)
data.x.ticks.chimp = data.frame(x1=c(0,2.5,5)*0.0001+0.0035,x2=c(0,2.5,5)*0.0001+0.0035,ystart=c(1,1,1)*0+.0008*scale_y-offset_y, yend=c(1,1,1)*0.00005+0.0008*scale_y-offset_y)


plot_species <- ggplot() 


plot_species <- plot_species+ 
  geom_rect(data=rect.data, mapping=aes(xmin=x1, xmax=x2, ymin=ystart, ymax=yend, fill=color), alpha=1) +
  scale_fill_manual(values = color_hex) +
  theme(legend.position="none")
plot_species <- plot_species+ 
  geom_rect(data=chimp.rect.data, mapping=aes(xmin=x1, xmax=x2, ymin=ystart, ymax=yend, fill=color), alpha=1)


for (i in seq(1,length(migration_rate.dir$x1))){
  post = migration_rate[i,]$rate.clade
  prio = 1-migration_rate[i,]$prior
  BF = (post/prio)/((1-post)/(1-prio))
  if (BF>2){
    plot_species <- plot_species +
      geom_curve(data = migration_rate.dir[i,],
                 aes(x = x2, y = y1, xend = x1, yend = y2 ), size=BF/4,  arrow = arrow(length = unit(0.01, "npc")), curvature = -0.05)
#       geom_curve(data = migration_rate.dir[i,],
#                aes(x = x1, y = y1, xend = x2, yend = y2), size=BF/4,  arrow = arrow(length = unit(0.01, "npc")), curvature = 0.05)
  }
}

for (i in seq(1,length(chimp.migration_rate$x1))){
  post = chimp.migration_rate[i,]$rate.clade
  prio = 1-chimp.migration_rate[i,]$prior
  BF = (post/prio)/((1-post)/(1-prio))
  
  if (BF>2){
    plot_species <- plot_species +
      geom_curve(data = chimp.migration_rate[i,],
                 aes(x = x2, y = y1, xend = x1, yend = y2), size=BF/4,  arrow = arrow(length = unit(0.01, "npc")), curvature = -0.05)
  }
}




plot_species <- plot_species + 
  geom_segment(data=branch, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend), size=0.25, linetype=b)
# plot_species <- plot_species + 
#   geom_segment(data=chimp.branch, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend), size=0.25, linetype=b)

plot_species <- plot_species + 
  geom_segment(data=data.axis, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend))+
  geom_segment(data=data.y.ticks, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend)) +
  annotate("text", x = -0.00015, y = 0, label = "0.00", size=4) +
  annotate("text", x = -0.00015, y = 0.0025, label = "0.25", size=4) +
  annotate("text", x = -0.00015, y = 0.0050, label = "0.50", size=4) +
  annotate("text", x = -0.00015, y = 0.0075, label = "0.75", size=4) +
  annotate("text", x = -0.00015, y = 0.0100, label = "1.00", size=4) +
  geom_segment(data=data.x.ticks, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend)) +
  annotate("text", x = 0.00195, y = 0.0012, label = "0.5", size=4)



plot_species <- plot_species + 
  geom_segment(data=data.axis.chimp, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend))+
  geom_segment(data=data.y.ticks.chimp, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend)) +
  annotate("text", x = -0.00015+0.00185, y = 0*scale_y-offset_y, label = "0.0", size=4) +
  annotate("text", x = -0.00015+0.00185, y = 0.0005*scale_y-offset_y, label = "0.05", size=4) +
  annotate("text", x = -0.00015+0.00185, y = 0.001*scale_y-offset_y, label = "0.1", size=4) +
  annotate("text", x = -0.00015+0.00185, y = 0.0015*scale_y-offset_y, label = "0.15", size=4) +
  geom_segment(data=data.x.ticks.chimp, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend)) +
  annotate("text", x = 0.00375, y = 0.0008*scale_y-offset_y +0.0002, label = "0.25", size=4) 


plot_species <- plot_species + 
  geom_segment(data=edge, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend), size=0.5)

plot_species <- plot_species + 
  geom_segment(data=chimp.edge, mapping=aes(x=x1, xend=x2, y=ystart, yend=yend), size=0.5) +
  theme(
    axis.line=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    legend.position="none") +
  coord_cartesian(xlim=c(-0.0003, 0.01)) +
  coord_cartesian(ylim=c(-0.0003, 0.01)) +
  theme_void()


plot(plot_species)

ggsave(plot=plot_species,"../../text/figures/greatapes/GreatApes_ratesandtree.pdf",width=10, height=10)














color_hex = c("#4137C2", "#4066CF", "#4B8DC2", "#5DA8A3", "#77B67F", "#96BD60", "#B8BC4B", "#D4B13F", "#E59638", "#E4672F", "#DC2F24")

color_hex = color_hex[c(10,8,6,4,1,2,3,5,7,9,11)]

violin_colors = c("Ggg"=color_hex[1], "Hsa"=color_hex[2], "Ppa"=color_hex[3], "Ptv"=color_hex[4],
                   "Pts"=color_hex[5], "Ptt"=color_hex[6], "Ptts"=color_hex[7],
                   "chimp"=color_hex[8], "cbo"=color_hex[9], "cbh"=color_hex[10],
                   "root"=color_hex[11])


p_violin <- ggplot(data=Ne.cond_true_top) +
  geom_violin(aes(1,Ggg,fill="Ggg")) + 
  geom_violin(aes(2,Hsa,fill="Hsa")) + 
  geom_violin(aes(3,Ppa,fill="Ppa")) + 
  geom_violin(aes(4,Ptv,fill="Ptv")) + 
  geom_violin(aes(5,Ptt,fill="Ptt")) + 
  geom_violin(aes(6,Pts,fill="Pts")) + 
  geom_violin(aes(7,Ptts,fill="Ptts"))+
  geom_violin(aes(8,chimp,fill="chimp"))+
  geom_violin(aes(9,cbo,fill="cbo"))+
  geom_violin(aes(10,cbh,fill="cbh"))+
  geom_violin(aes(11,root,fill="root"))+
   ylab("effective population sizes") +
    scale_fill_manual(values = violin_colors)+
  theme(
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title=element_blank(),
    legend.position="none") 
plot(p_violin)
ggsave(plot=p_violin,"../../text/figures/greatapes/GreatApes_Ne.pdf",width=5, height=2)




p_violin <- ggplot(data=height.cond_true_top) +
  geom_violin(aes(7,Ptts,fill="Ptts"))+
  geom_violin(aes(8,chimp,fill="chimp"))+
  geom_violin(aes(9,cbo,fill="cbo"))+
  geom_violin(aes(10,cbh,fill="cbh"))+
  geom_violin(aes(11,root,fill="root"))+
  ylab("effective population sizes") +
  scale_fill_manual(values = violin_colors)+
  theme(
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title=element_blank(),
    legend.position="none") +
  scale_y_log10()
plot(p_violin)
ggsave(plot=p_violin,"../../text/figures/greatapes/GreatApes_height.pdf",width=2.3, height=2)








# plot the asymmetric clade rates
diff <- c()
for (i in seq(1, length(migration_rate$from))){
  post = migration_rate[i,"rate.clade"]
  prio = 1 - migration_rate[i,"prior"]
  diff[i] = (post/prio)/((1-post)/(1-prio))
}

sorted_indices = sort(diff, decreasing = TRUE, index.return=TRUE)
bar.data = data.frame()
for (j in seq(length(sorted_indices[[2]]), 1)){
  i = sorted_indices[[2]][j]
  new.bar.data.post = data.frame(name = paste(migration_rate[i,"from"], "->", migration_rate[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  migration_rate[i,"rate.clade"], prob="posterior")
  new.bar.data.prior = data.frame(name = paste(migration_rate[i,"from"], "->", migration_rate[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  1-migration_rate[i,"prior"], prob="prior")
  new.bar.data = rbind(new.bar.data.prior, new.bar.data.post)
  if (sorted_indices[[1]][j] > 2){
      bar.data = rbind(bar.data, new.bar.data)
  }
}

plot.Asym.rate.clade <- ggplot(data=bar.data, aes(x=name,y=val,fill=prob)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position="none") + 
  scale_fill_manual(values = c("posterior" = "black", "prior" = "grey"))
plot(plot.Asym.rate.clade)
ggsave(plot=plot.Asym.rate.clade,"../../text/figures/greatapes/GreatApes_rates_asym_clade.pdf",width=4, height=2)

# plot the symmetric clade rates
diff <- c()
for (i in seq(1, length(migration_rate.dir$from))){
  post = migration_rate.dir[i,"rate.clade"]
  prio = 1 - migration_rate.dir[i,"prior"]
  diff[i] = (post/prio)/((1-post)/(1-prio))
  
}

sorted_indices = sort(diff, decreasing = TRUE, index.return=TRUE)
bar.data = data.frame()
for (j in seq(length(sorted_indices[[2]]), 1)){
  i = sorted_indices[[2]][j]
  new.bar.data.post = data.frame(name =  paste(migration_rate.dir[i,"from"], "<->", migration_rate.dir[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  migration_rate.dir[i,"rate.clade"], prob="posterior")
  new.bar.data.prior = data.frame(name =  paste(migration_rate.dir[i,"from"], "<->", migration_rate.dir[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  1-migration_rate.dir[i,"prior"], prob="prior")
  new.bar.data = rbind(new.bar.data.prior, new.bar.data.post)
  if (sorted_indices[[1]][j] > 2){
    bar.data = rbind(bar.data, new.bar.data)
  }
}

plot.Sym.rate.clade <- ggplot(data=bar.data, aes(x=name,y=val,fill=prob)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position="none") + 
  scale_fill_manual(values = c("posterior" = "black", "prior" = "grey"))
plot(plot.Sym.rate.clade)
ggsave(plot=plot.Sym.rate.clade,"../../text/figures/greatapes/GreatApes_rates_sym_clade.pdf",width=4, height=2)





# plot the asymmetric rates
diff <- c()
for (i in seq(1, length(migration_rate$from))){
  post = migration_rate[i,"rate"]
  prio = 0.05
  diff[i] = (post/prio)/((1-post)/(1-prio))
}

sorted_indices = sort(diff, decreasing = TRUE, index.return=TRUE)
bar.data = data.frame()
for (j in seq(length(sorted_indices[[2]]), 1)){
  i = sorted_indices[[2]][j]
  new.bar.data.post = data.frame(name =  paste(migration_rate[i,"from"], "->", migration_rate[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  migration_rate[i,"rate"], prob="posterior")
  new.bar.data.prior = data.frame(name =  paste(migration_rate[i,"from"], "->", migration_rate[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  0.05, prob="prior")
  new.bar.data = rbind(new.bar.data.prior, new.bar.data.post)
  if (sorted_indices[[1]][j] > 2){
    bar.data = rbind(bar.data, new.bar.data)
  }
}

plot.Asym.rate <- ggplot(data=bar.data, aes(x=name,y=val,fill=prob)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position="none") + 
  scale_fill_manual(values = c("posterior" = "black", "prior" = "grey"))
plot(plot.Asym.rate)
ggsave(plot=plot.Asym.rate,"../../text/figures/greatapes/GreatApes_rates_asym.pdf",width=4, height=2)



# plot the symmetric  rates
diff <- c()
for (i in seq(1, length(migration_rate.dir$from))){
  post = migration_rate.dir[i,"rate"]
  prio = 0.0975
  diff[i] = (post/prio)/((1-post)/(1-prio))
  
}

sorted_indices = sort(diff, decreasing = TRUE, index.return=TRUE)
bar.data = data.frame()
for (j in seq(length(sorted_indices[[2]]), 1)){
  i = sorted_indices[[2]][j]
  new.bar.data.post = data.frame(name = paste(migration_rate.dir[i,"from"], "<->", migration_rate.dir[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  migration_rate.dir[i,"rate"], prob="posterior")
  new.bar.data.prior = data.frame(name = paste(migration_rate.dir[i,"from"], "<->", migration_rate.dir[i,"to"], "(BF=",sprintf("%.2f",sorted_indices[[1]][j]),")", sep=""), val =  0.0975, prob="prior")
  new.bar.data = rbind(new.bar.data.prior, new.bar.data.post)
  if (sorted_indices[[1]][j] > 2){
    bar.data = rbind(bar.data, new.bar.data)
  }
}

plot.Sym.rate <- ggplot(data=bar.data, aes(x=name,y=val,fill=prob)) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position="none") + 
  scale_fill_manual(values = c("posterior" = "black", "prior" = "grey"))
plot(plot.Sym.rate)
ggsave(plot=plot.Sym.rate,"../../text/figures/greatapes/GreatApes_rates_sym.pdf",width=4, height=2)











# # print the migration rates and Ne's to ensure ess values
log_vals = data.frame(matrix(ncol = 0, nrow = length(Ne.cond_true_top$Ppa)))
log_vals[,"sample"] = seq(1,length(Ne.cond_true_top$Ppa))
for (i in seq(1, length(possible_routes$from))){
  name1 = as.character(possible_routes[i,1])
  name2 = as.character(possible_routes[i,2])
  new.rates <- migration.cond_true_top[intersect(which(migration.cond_true_top$from==name1),(which(migration.cond_true_top$to==name2))), "rate"]
  
  name = paste(name1, "to", name2, sep=".")
  name2 = paste(name1, "to", name2, "ind", sep=".")
  log_vals[,name] = new.rates
  log_vals[,name2] = as.numeric(new.rates>0)
}

names_Ne = names(Ne.cond_true_top)
for (i in seq(1, length(names_Ne))){
  name = names_Ne[i]
  log_vals[,name] = Ne.cond_true_top[, name]
}

names_heights = names(height.cond_true_top)
for (i in seq(1, length(names_heights))){
  name = paste(names_heights[i], "nodeheight", sep="_")
  log_vals[,name] = height.cond_true_top[, names_heights[i]]
}

write.table(log_vals, file="./combined/greatapes/true_topology_chimps_tmp.log", quote=F, row.names=F)
system(paste("/Applications/BEAST\\ 2.4.7/bin/logcombiner -burnin 0 -log ./combined/greatapes/true_topology_chimps_tmp.log -o ./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.true_topology.log"))
system("rm -r ./combined/greatapes/true_topology_chimps_tmp.log")


# convert the log file to add a topology value
all_log=read.table("./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.log", header=TRUE, sep="\t")
all_log$true.topology = tree_has_true_topology

for (i in seq(1, length(all_log$true.topology))){
  all_log[i,"Sample"] = i
}
write.table(all_log, file="./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.with_topology_tmp.log", quote=F, row.names=F)
system(paste("/Applications/BEAST\\ 2.4.7/bin/logcombiner -burnin 0 -log ./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.with_topology_tmp.log -o ./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.with_topology.log"))
system("rm -r ./combined/greatapes/Chimp_Ggg_Hsa_Ppa_Pts_Ptt_Ptv.with_topology_tmp.log")

