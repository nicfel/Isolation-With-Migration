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
library("OutbreakTools")
library(ggtree)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# define the clock rate used in per million years
clock_rate = 2.8*10^-9*11 *10^6

loci = c("s50", "s200")
mig = c("mig1", "mig2")
chr = c("all", "chrX")

plot_how_many_trees = 1
for (a in seq(2,2)){ #loci
  for (b in seq(1,2)){ # migration priors
    for (c in seq(1,1;)){ # which chromosomes
      tree_files <- list.files(path=paste("./combined/", loci[[a]], sep=""), pattern=paste("*",chr[[c]],"_",loci[[a]], "_", mig[[b]] ,"_species.trees", sep=""), full.names = TRUE)

      rep_nr = 0
      plot_all = list()
      plot_all_sub = list()
      plot_all_collapse = list()
      plot_count=1
      
      for (tf in seq(1,length(tree_files))){
        
        # Read In Data ---------------------------------
        mig_rate <- data.frame()
        counter <- 1
        
        i = 1
        
        # Make the filenames for all possible migration rates
        # trees <- "./combined/s50/anopheles_0_0_50_s50_sub1_mig1_species.trees"
        trees <- tree_files[tf]
        
        # get all the tree topologies in the posterior, if the file isn't ended by END;, this wll cause an error
        est_species_tree_with_nodes <- read.annotated.nexus(trees)
        
        # get all the ossible phylogenies, does however not get the ranked topologies
        unique_topologies <- unique.multiPhylo(est_species_tree_with_nodes,use.edge.length = FALSE, use.tip.label = TRUE)
        
        
        annotations <- vector(, length(est_species_tree_with_nodes))
        
        # get the node order for each of the species trees
        for (i in seq(1,length(est_species_tree_with_nodes))){
          for (j in seq(1,length(unique_topologies))){
            eq <- all.equal(est_species_tree_with_nodes[[i]], unique_topologies[[j]],use.edge.length = FALSE, use.tip.label = TRUE)
            if (eq){
              top_nr <- j
              break;
            }
          }
          # get the node heights and all subtrees
          subs = subtrees(est_species_tree_with_nodes[[i]])
          
          nr_tips = length(est_species_tree_with_nodes[[i]]$tip.label)
          
          # init the annotations data frame
          new.annotations <- data.frame(sample=i)
          
          # get the root heights of all subtrees
          node_heights <- rep(0,2*nr_tips-1)
          node_names <- vector(,2*nr_tips-1)
          for (j in seq(1,nr_tips)){
            node_names[j] <- est_species_tree_with_nodes[[i]]$tip.label[j]
          }
          
          for (j in seq(1,length(subs))){
            all_heights <- nodeHeights(subs[[j]])
            node_heights[nr_tips+j] <- max(all_heights[,2])
            node_names[nr_tips+j] <- paste(sort(subs[[j]]$tip.label), collapse=":")
            
            new_name = paste("height_", node_names[nr_tips+j], sep="")
            new.annotations[,new_name] = max(all_heights[,2])
            
          }
          
          # get the correct mapping of labels to node numbers, make it way too complicated
          tree_edges <- est_species_tree_with_nodes[[i]]$edge
          
          node_species <- rep(-1,2*nr_tips-1)
          node_ne <- rep(-1,2*nr_tips-1)
          
          
          
          # get the nes
          for (j in seq(1, length(tree_edges[,2]))){
            node_label <- tree_edges[j,2]
            node_species[node_label] <- est_species_tree_with_nodes[[i]]$annotations[[j]]$species
            new_name = paste("Ne_", node_names[node_label], sep="")
            new.annotations[,new_name] = est_species_tree_with_nodes[[i]]$annotations[[j]]$Ne
          }
          node_label <- nr_tips+1
          new_name = paste("Ne_", node_names[node_label], sep="")
          new.annotations[,new_name] = est_species_tree_with_nodes[[i]]$root.annotation$Ne
          
          # get the migration rates
          for (j in seq(1, length(tree_edges[,2]))){
            node_label <- tree_edges[j,2]
            for (k in seq(1, length(est_species_tree_with_nodes[[i]]$annotations[[j]]$to))){
              target = which(node_species == est_species_tree_with_nodes[[i]]$annotations[[j]]$to[[k]])
              new_name = paste("mig_", node_names[node_label], "_", node_names[target], sep="")
              new.annotations[,new_name] = est_species_tree_with_nodes[[i]]$annotations[[j]]$rates[[k]]
            }
          }
          
          # get the correct order
          indices <- sort(node_heights,index.return = TRUE)
          
          
          new.dat <- data.frame(topology=top_nr)
          new.heights <- data.frame(topology=top_nr)
          
          annotations[i] <- list(new.annotations)
          
          
          for (j in seq(1,length(indices$ix))){
            new_name = paste("n", j, sep="")
            new.dat[,new_name] = node_names[[indices$ix[[j]]]]
            new.heights[,new_name] = node_heights[[indices$ix[[j]]]]
          }
          
          if (i==1){
            tree_top = new.dat
          }else{
            tree_top = rbind(tree_top,new.dat)
          }
        }
        
        
        # get the unique ranked topologies
        unique_ranked <- unique(tree_top)
        counts = rep(0,length(unique_ranked$topology))
        # get the posterior support for each ranked topology (haven't found a better way for now at least)
        for (i in seq(1,length(unique_ranked$topology))){
          for (j in seq(1, length(tree_top$topology))){
            if (sum(unique_ranked[i,]!=tree_top[j,])==0){
              counts[i] = counts[i]+1
            }
          }
        }
        
        
        # get the sorted count indices
        sorted_counts = sort(counts, index.return=TRUE, decreasing = T)$ix
        unique_ranked = unique_ranked[sorted_counts,]
        
        counts = counts[sorted_counts]
        
        median_node_heights <- vector(, length(unique_ranked$topology))
        lower_node_heights <- vector(, length(unique_ranked$topology))
        upper_node_heights <- vector(, length(unique_ranked$topology))
        
        log <- vector(, length(est_species_tree_with_nodes))
        
        # plot only the first n unique ranked topologies
        plot_nr = min(length(unique_ranked$topology),plot_how_many_trees)
        
        # for each of these trees, make an individual log file
        for (i in seq(1, plot_nr)){
          # get the row and change the row name to compare rows
          ur <- unique_ranked[i,]
          row.names(ur) <- 1
          
          first = TRUE
          
          for (j in seq(1, length(tree_top$n1))){
            tocomp <- tree_top[j,]
            row.names(tocomp) <- 1
            
            if (identical(ur,tocomp)){
              
              if (first){
                log.data <- annotations[[j]]
                all_labels = labels(log.data)
                first <- FALSE
              }else{
                new.dat <- data.frame(sample=j)
                for (k in seq(1,length(all_labels[[2]]))){
                  new.dat[1, all_labels[[2]][[k]]] <- annotations[[j]][,all_labels[[2]][[k]]]
                }
                
                log.data <- rbind(log.data,annotations[[j]])
              }
            }
          }
          
          # get the median node heights of the ranked trees (and the upper and lower limits)
          new.median_node_heights <- vector(, length(node_heights))
          
          new.lower_node_heights <- vector(, length(node_heights))
          new.upper_node_heights <- vector(, length(node_heights))
          
          all_labels = labels(log.data)
          nc = 1
          for (j in seq(1,length(all_labels[[2]]))){
            if (startsWith(all_labels[[2]][[j]],"height_")){
              new.median_node_heights[[nc]] <- median(log.data[,all_labels[[2]][[j]]]/clock_rate)
              hpd <- HPDinterval(as.mcmc(log.data[,all_labels[[2]][[j]]]/clock_rate))
              new.lower_node_heights[[nc]] = hpd[1,"lower"]
              new.upper_node_heights[[nc]] = hpd[1,"upper"]
              nc <- nc+1
            }
          }
          
          median_node_heights[[i]] <- list(new.median_node_heights)
          lower_node_heights[[i]] <- list(new.lower_node_heights)
          upper_node_heights[[i]] <- list(new.upper_node_heights)
          
          log[i] <- list(log.data)
          
          
          fname1 <- gsub("\\.trees", paste("_", i, "tmp.log", sep=""), trees)
          fname2 <- gsub("\\.trees", paste("_", i, ".log", sep=""), trees)
          
          write.table(log.data, file=fname1, quote=F, row.names=F, sep="\t")
          
          system(paste("/Applications/BEAST\\ 2.5.0/bin/logcombiner -renumber -burnin 0 -log", fname1, "-o", fname2, sep=" "))
          system(paste("rm -r ", fname1))
        }
        
        
        # recompute the branch lengths of each unique ranked phylogeny to be the median branch length
        for (i in seq(1, plot_nr)){
          # find the first tree with the same ranked topology
          ur <- unique_ranked[i,]
          row.names(ur) <- 1
          for (j in seq(1, length(tree_top$n1))){
            tocomp <- tree_top[j,]
            row.names(tocomp) <- 1
            if (identical(ur,tocomp)){
              plot_tree <- ladderize(est_species_tree_with_nodes[[j]],right=F)
              break
            }
          }
          
          # remove annotations
          plot_tree$annotations <- NULL
          plot_tree$root.annotation <- NULL
          
          # get the median node heights
          sorted.node.indices = sort(median_node_heights[[i]][[1]], index.return=T) 
          
          node_heights <- median_node_heights[[i]][[1]]
          node_heights.lower <- median_node_heights[[i]][[1]]
          node_heights.upper <- median_node_heights[[i]][[1]]
          for (k in seq(1,length(sorted.node.indices$ix))){
            node_heights[[k]] <- median_node_heights[[i]][[1]][[sorted.node.indices$ix[[k]]]]
            node_heights.lower[[k]] <- lower_node_heights[[i]][[1]][[sorted.node.indices$ix[[k]]]]
            node_heights.upper[[k]] <- upper_node_heights[[i]][[1]][[sorted.node.indices$ix[[k]]]]
          }
          ape_node_edges = plot_tree$edge
          
          # get the ordering of nodes in the ape tree by comparing the heights
          ape_node_heights = abs(nodeHeights(plot_tree)-max(nodeHeights(plot_tree)[,2]))[seq(1,length(ape_node_edges[,1]),2),1]
          # get the corresponding numbers
          ape_node_numbers = ape_node_edges[seq(1,length(ape_node_edges[,1]),2),1]
          
          indices <- order(ape_node_heights)
          
          # get the original indices
          ori_indices <- seq(1,length(indices))
          
          new.node_heights <- node_heights
          new.node_heights.lower <- node_heights
          new.node_heights.upper <- node_heights
          # order the node_heights to match ape noe heights
          for (j in seq(1, length(indices))){
            new.node_heights[indices[j]+length(plot_tree$tip.label)] <- node_heights[j+length(plot_tree$tip.label)]
            new.node_heights.lower[indices[j]+length(plot_tree$tip.label)] <- node_heights.lower[j+length(plot_tree$tip.label)]
            new.node_heights.upper[indices[j]+length(plot_tree$tip.label)] <- node_heights.upper[j+length(plot_tree$tip.label)]
          }
          
          node_heights <- new.node_heights
          new.node_heights <- node_heights
          
          node_heights.lower <- new.node_heights.lower
          new.node_heights.lower <- node_heights.lower
          
          node_heights.upper <- new.node_heights.upper
          new.node_heights.upper <- node_heights.upper
          
          # now also match the node order
          for (j in seq(1,length(ape_node_numbers))){
            new.node_heights[ape_node_numbers[j]] <- node_heights[j+length(plot_tree$tip.label)]
            new.node_heights.lower[ape_node_numbers[j]] <- node_heights.lower[j+length(plot_tree$tip.label)]
            new.node_heights.upper[ape_node_numbers[j]] <- node_heights.upper[j+length(plot_tree$tip.label)]
          }
          
          plot_tree.lower = plot_tree
          plot_tree.upper = plot_tree
          
          # change the edge lengths to match the median heights
          for (j in seq(1, length(plot_tree$edge.length))){
            plot_tree$edge.length[j] <- new.node_heights[ape_node_edges[j,1]] -  new.node_heights[ape_node_edges[j,2]]
            plot_tree.lower$edge.length[j] <- new.node_heights.lower[ape_node_edges[j,1]] -  new.node_heights.lower[ape_node_edges[j,2]]
            plot_tree.upper$edge.length[j] <- new.node_heights.upper[ape_node_edges[j,1]] -  new.node_heights.upper[ape_node_edges[j,2]]
          }
          
          rm(mig_rates)
          
          # get all migration rate with a bayes factor over ...
          prior = 0.05
          log.data = log[[i]]
          all_labels = labels(log.data)
          nc=1;
          for (j in seq(1,length(all_labels[[2]]))){
            # check if the label describes a route of migration
            if (startsWith(all_labels[[2]][[j]],"mig_")){
              # calculate the mean rate
              meanrate = log.data[which(log.data[,all_labels[[2]][[j]]]>0),all_labels[[2]][[j]]]
              # calculate the posterior support for that route to be active
              posterior = length(which(log.data[,all_labels[[2]][[j]]]>0))/length(log.data[,all_labels[[2]][[j]]])
              # caluclate the Bayes Factor for that route being active
              bayes = posterior*(1-prior)/((1-posterior)*prior)
              new.mig_rates <- data.frame(from = strsplit(all_labels[[2]][[j]], "_")[[1]][3], to = strsplit(all_labels[[2]][[j]], "_")[[1]][2], 
                                          BF=bayes, post=posterior, rate=mean(meanrate))
              
              if (nc==1){
                mig_rates = new.mig_rates
                nc <- nc+1
              }else{
                mig_rates = rbind(mig_rates, new.mig_rates)
              }
              
            }
          }
          rm(arrows.data)
          rm(new.arrows)
          
          
          rootheight = max(nodeHeights(plot_tree))
          
          x.width = 0.1
          y.width = 0
          
          # build the arrow dataframe
          if (exists("mig_rates")){
            for (j in seq(1,length(mig_rates$from))){
              from_leaves = strsplit(as.character(mig_rates[j, "from"]), split=":")
              if (length(from_leaves[[1]])==1){
                from_mrca = which(plot_tree$tip.label==from_leaves[[1]][[1]])
              }else{
                from_mrca = findMRCA(plot_tree, tips=from_leaves[[1]], type=c("node","height"))
              }
              
              to_leaves = strsplit(as.character(mig_rates[j, "to"]), split=":")
              if (length(to_leaves[[1]])==1){
                to_mrca = which(plot_tree$tip.label==to_leaves[[1]][[1]])
              }else{
                to_mrca = findMRCA(plot_tree, tips=to_leaves[[1]], type=c("node","height"))
              }
              
              # get the heights of from and to nodes and the heights of their parents
              from_heights = nodeHeights(plot_tree)[which(plot_tree$edge[,2]==from_mrca),]
              to_heights = nodeHeights(plot_tree)[which(plot_tree$edge[,2]==to_mrca),]
              
              # y_coord = (min(from_heights[2], to_heights[2])-max(from_heights[1], to_heights[1]))/2+max(from_heights[1], to_heights[1])
              
              y_coord_max = min(from_heights[2], to_heights[2])
              y_coord_min = max(from_heights[1], to_heights[1])
              
              
              if (node.height(plot_tree)[from_mrca] > node.height(plot_tree)[to_mrca]){
                new.arrows = data.frame(x=node.height(plot_tree)[from_mrca]-x.width, xend=node.height(plot_tree)[to_mrca]+x.width,
                                        y_coord_max=rootheight-y_coord_max, y_coord_min=rootheight-y_coord_min, post=mig_rates[j, "post"], rate=mig_rates[j, "rate"])
              }else{
                new.arrows = data.frame(x=node.height(plot_tree)[from_mrca]+x.width, xend=node.height(plot_tree)[to_mrca]-x.width,
                                        y_coord_max=rootheight-y_coord_max, y_coord_min=rootheight-y_coord_min, post=mig_rates[j, "post"], rate=mig_rates[j, "rate"])
              }
              
              
              if (j==1){
                arrows.data = new.arrows
              }else{
                arrows.data = rbind(arrows.data, new.arrows)
              }
            }
          }
          
          # get the segments to plot the tree
          rm(tree.rectangular.data)
          
          rootheight = max(nodeHeights(plot_tree))
          rootheight.lower = max(nodeHeights(plot_tree.lower))
          rootheight.upper = max(nodeHeights(plot_tree.upper))
          
          edges=plot_tree$edge
          for (j in seq(1,length(edges[,1]))){
            source.height = rootheight-nodeHeights(plot_tree)[j,1]
            source.lower = rootheight.lower-nodeHeights(plot_tree.lower)[j,1]
            source.upper = rootheight.upper-nodeHeights(plot_tree.upper)[j,1]
            dest.height = rootheight-nodeHeights(plot_tree)[j,2]
            
            x.height = node.height(plot_tree)[edges[j,2]]
            
            
            # check if the parent edge is left or right
            parent.x.height = node.height(plot_tree)[edges[j,1]]
            if (parent.x.height<x.height){
              x.end.ver.1 = x.height - x.width
              x.end.ver.2 = x.height + x.width
              
              x.ver.1 = parent.x.height + x.width
              x.ver.2 = parent.x.height
              
              source.height.1 = source.height
              source.height.2 = source.height + y.width
              
              if (edges[j,2]>length(edges[,1])/2){
                dest.height = dest.height+y.width
              }
            }else{
              x.end.ver.1 = x.height -x.width
              x.end.ver.2 = x.height +x.width
              x.ver.1 = parent.x.height
              x.ver.2 = parent.x.height  - x.width
              
              source.height.1 = source.height+y.width
              source.height.2 = source.height
              if (edges[j,2]>length(edges[,1])/2){
                dest.height = dest.height+y.width
              }
            }
            
            new.rectangular = data.frame(xend=x.height-x.width, x=x.height-x.width, y=source.height.1, yend=dest.height)
            new.rectangular = rbind(new.rectangular, data.frame(xend=x.height+x.width, x=x.height+x.width, y=source.height.2, yend=dest.height))
            
            # add vertical lines
            new.rectangular = rbind(new.rectangular, data.frame(xend=x.end.ver.1, x=x.ver.2, y=source.height.1, yend=source.height.1))
            new.rectangular = rbind(new.rectangular, data.frame(xend=x.end.ver.2, x=x.ver.1, y=source.height.2, yend=source.height.2))
            
            # add node height error bars
            new.width = data.frame(xend=parent.x.height, x=parent.x.height, y=source.lower, yend=source.upper, median=source.height.1)
            
            
            if (j==1){
              tree.rectangular.data = new.rectangular
              tree.node.height.width = new.width
            }else{
              tree.rectangular.data = rbind(tree.rectangular.data, new.rectangular)
              tree.node.height.width = rbind(tree.node.height.width, new.width)
            }
          }
          
          rm(plot_arrows)
          # plot_arrows = arrows.data
          plot_arrows = arrows.data[arrows.data$post>0.3448, ]
          # add empty y values
          plot_arrows$y = c()
          # set the y values of the arrows to minimize overlap of arrows
          col_subset = plot_arrows[,c(3,4)]
          unique_from = unique(col_subset)
          for (u in seq(1,length(unique_from$y_coord_max))){
            ind = which(plot_arrows$y_coord_max==unique_from$y_coord_max[[u]] & plot_arrows$y_coord_min==unique_from$y_coord_min[[u]])
            # set the y start and y end values for these arrows
            
            # get the difference between max and min value
            diff.max.min = unique_from$y_coord_max[[u]] - unique_from$y_coord_min[[u]]
            for (j in seq(1,length(ind))){
              plot_arrows[ind[j], "y"] = unique_from$y_coord_min[[u]] + y.width + (diff.max.min-y.width)*j/(length(ind)+1)
            }
          }
          
          plot_arrows$yend = plot_arrows$y
          
          
          
          p <- ggplot()
          p <- p+
            ylim(c(-0.001, 9.5)) +xlim(c(0,9)) +
            geom_segment(data=tree.rectangular.data, aes(x=x,y=y,xend=xend,yend=yend)) +
            geom_segment(data=tree.node.height.width, aes(x=x,y=y,xend=xend,yend=yend), size=2, color="gray48") +
            geom_curve(data=plot_arrows,curvature = 0,
                       aes(x=x, xend=xend,y = y, yend=yend,size=post),
                       arrow = arrow(angle=5,type="closed",ends="last")) +
            scale_size_continuous(limits=c(0,1),range = c(0,1)) + 
            theme_minimal()
          plot(p)
          
          ori_label = c("AmerM2", "AquaS1", "AaraD1", "AmelC2", "AcolM1", "AgamS1", "AchrA1", "AepiE1")
          d = data.frame(label=ori_label, label2 = c("mer", "qua", "ara", "mel", "col", "gam", "chr", "epi"))
          
          # add tip labels
          text_val = c()
          for (j in seq(1,length(ori_label))){
            to_mrca = which(plot_tree$tip.label==ori_label[[j]])
            text_val[[node.height(plot_tree)[to_mrca]]] = as.character(d$label2[[j]])
          }
          p <- p + scale_x_continuous(breaks=seq(1,8),labels=text_val) +
            xlab("") + ylab("Million years") + theme(legend.position="none") +
            ggtitle(sprintf("posterior support = %.2f", counts[[i]]/sum(counts)))
          
          
          # p <- p  
          plot(p)
          
          
          # build the same plot, but only using the 6 closely related species
          p.sub <- p  + scale_x_continuous(breaks=seq(3,8),labels=text_val[seq(3,8)], limits=c(2.5,8.5)) +  ylim(c(-0.001, 0.8))
          plot(p.sub)
          
          
          
          # make the plots with all species, but collapse everything but the outgroup
          tree.rectangular.data.collapse = tree.rectangular.data
          sub.clade = tree.rectangular.data.collapse[which(tree.rectangular.data.collapse$y<2),]
          max.x = max(sub.clade[, 'x'])
          min.x = min(sub.clade[, 'x'])
          # get the row with max y value
          max.y.row = max(sub.clade[, 'y'])
          min.x.end = min(tree.rectangular.data.collapse[which(tree.rectangular.data.collapse$yend==max.y.row &
                                                                 tree.rectangular.data.collapse$y>2),'x'])
          max.x.end = max(tree.rectangular.data.collapse[which(tree.rectangular.data.collapse$yend==max.y.row &
                                                                 tree.rectangular.data.collapse$y>2),'x'])
          
          ## remove subclades
          tree.rectangular.data.collapse = tree.rectangular.data.collapse[-which(tree.rectangular.data.collapse$y<2),]
          ## add collapsed subclades
          col1 = data.frame(x=min.x.end, xend=min.x, yend=0, y=max.y.row)
          col2 = data.frame(x=max.x.end, xend=max.x, yend=0, y=max.y.row)
          
          tree.rectangular.data.collapse = rbind(tree.rectangular.data.collapse, col1, col2)
          # remove the node height error bars
          tree.node.height.width.collapse = tree.node.height.width
          tree.node.height.width.collapse = tree.node.height.width.collapse[-which(tree.node.height.width.collapse$median<2), ]
          # remove all arrows that don't involve the outgroups
          plot_arrows.collapse = plot_arrows
          plot_arrows.collapse = plot_arrows.collapse[-which(plot_arrows.collapse$x>2.5 & plot_arrows.collapse$xend>2.5),]
          
          
          p.collapse <- ggplot()
          p.collapse <- p.collapse+
            ylim(c(-0.001, 9.5)) +xlim(c(0,9)) +
            geom_segment(data=tree.rectangular.data.collapse, aes(x=x,y=y,xend=xend,yend=yend)) +
            geom_segment(data=tree.node.height.width.collapse, aes(x=x,y=y,xend=xend,yend=yend), size=2, color="gray48") +
            geom_curve(data=plot_arrows.collapse,curvature = 0,
                       aes(x=x, xend=xend,y = y, yend=yend,size=post),
                       arrow = arrow(angle=5,type="closed",ends="last")) +
            scale_size_continuous(limits=c(0,1),range = c(0,1)) + 
            theme_minimal() +
            scale_x_continuous(breaks=seq(1,8),labels=text_val) +
            xlab("") + ylab("Million years") + theme(legend.position="none") +
            ggtitle(sprintf("posterior support = %.2f", counts[[i]]/sum(counts)))
          plot(p.collapse)
          
          ggsave(plot=p.sub,paste("/Users/nicmuell/Documents/github/IsolationWithMigration/Figures/Gambia/", chr[[c]], "_", loci[[a]], "_", mig[[b]], "_combined_sub.pdf", sep=""),width=4, height=4)
          ggsave(plot=p.collapse,paste("/Users/nicmuell/Documents/github/IsolationWithMigration/Figures/Gambia/", chr[[c]], "_", loci[[a]], "_", mig[[b]], "_combined_collapse.pdf", sep=""),width=4, height=4)

        }
      }
    }
  }
}


