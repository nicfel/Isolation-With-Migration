######################################################
######################################################
# combine and plot AIM runs
######################################################
######################################################
library(ape) # laod ape
library(ggplot2) # load ggplot

this.dir <- dirname(parent.frame(2)$ofile) # Set get current directory
source(paste(this.dir, "../../Software/plotSpeciesTree.R", sep="/")) # sources the function to plot the species tree with arrows
sample_numbers = c("s50", "s200")# define the sample numbers
prior = 0.695/50 # define the prior probability of an indicator being active (here it's lambda/(nr indicators) )
posterior.threshold = 0.5 # migration has to have a higher posterior support to be plotted
BF.threshold = 0 # migration has to have a higher Bayes Factor to be plotted
forwards.arrows = T # if true, the arrows denote forward in time migration, if false, backwards in time migration

clock_rate = 3.08*10^(-2);
p <- list()
l_count = 1
# for (s in seq(1,1)){
  
for (s in seq(2,length(sample_numbers))){

  system(paste("rm -r ", this.dir, "/combined/", sample_numbers[[s]], sep=""))
  system(paste("mkdir ", this.dir, "/combined/", sample_numbers[[s]], sep=""))
  
  system(paste("rm -r ", this.dir, "/unique/", sample_numbers[[s]], sep=""))
  system(paste("mkdir ", this.dir, "/unique/", sample_numbers[[s]], sep=""))
  
  trees <- list.files(path=paste(this.dir,"/out/", sample_numbers[[s]], "/", sep=""), pattern="*rep0.*\\.trees$", full.names = TRUE)
  for (i in seq(1,length(trees))){
    in_command <- " -b 50 -resample 100000 -log"
    for (j in seq(0,2)){
      in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
    }
    
    out_command = gsub("rep0_", "", trees[i])
    out_command = gsub("out", "combined", out_command)
    
    combined_command = gsub(".trees",".trees", out_command)
    combined_command = paste(" -o ", combined_command, sep="")
    # combine the trees
    system(paste("/Applications/BEAST\\ 2.6.0/bin/logcombiner", in_command, combined_command, "", sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
    system(paste("/Applications/BEAST\\ 2.6.0/bin/logcombiner", gsub(".trees",".log", in_command), gsub(".trees",".log", combined_command), sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
    in_file  = gsub("-o ","", combined_command)
    out_file = gsub(" ", "", gsub("/combined/","/unique/", in_file))
    
    system(paste("java -jar ", this.dir,  "/../../Software/AIMannotator.jar -burnin 0 -userank false ", in_file, " ", out_file, sep=""), intern=TRUE)
    
    species_trees <- read.nexus(file=out_file,force.multi=T) # get all the tree topologies in the posterior
    # compute total number of trees, only works if Minimal tree support for output was 0 in AIM ann
    tot.trees = 0
    for (j in seq(1, length(species_trees))){
      post.occurance = strsplit(names(species_trees)[[j]], split="_")[[1]]
      tot.trees = tot.trees + as.numeric(post.occurance[[length(post.occurance)]])
    }
    print(i)
    # plots the most likely species tree for each run
    for (j in seq(1, min(1,length(species_trees)))){
      plot_tree <- ladderize(species_trees[[j]], right=F) # ladderizes the tree (defines the node orderign on the x-axis for the plot)
      pl <- plotSpeciesTree(out_file, names(species_trees)[[j]], plot_tree, posterior.threshold, BF.threshold, prior, forwards.arrows) # returns a ggplot object that can be modified
      post.occurance = strsplit(names(species_trees)[[j]], split="_")[[1]] # get how often this ranked topology was observed in the posterior
      pl <- pl + 
        ylab("substitutions per site") + # labels the y-axis
        scale_size_continuous(range=c(0.1,1)) + # defines the range of arrow sizes
        ggtitle(paste("posterior support = ", round(as.numeric(post.occurance[[length(post.occurance)]])/tot.trees,2), "", sep="")) + # specify the title
        theme_minimal() + # changes the plotting theme
        # coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) +
        coord_cartesian(xlim=c(0,8.5), ylim=c(0,0.2)) +
        scale_size_continuous(limits=c(0.5,1), range=c(0.1,1), breaks=c(0.5,0.75,1), name="posterior support\nfor gene flow")+
        scale_y_continuous(sec.axis = sec_axis(~ ./(clock_rate), name="Years in Million"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        xlab("")
      plot(pl)
      p[[l_count]] = pl
      l_count = l_count+1
    }
  }
}

library(ggpubr)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

legend <- g_legend(p[[1]]) 

figure <- ggarrange(p[[1]] +theme(legend.position="none"), 
                    p[[2]] +theme(legend.position="none"), 
                    p[[3]] +theme(legend.position="none"), 
                    legend,
                    labels = c("A", "B", "C", ""),
                    ncol = 4, nrow = 1, widths=c(1,1,1,0.27))
ggsave(plot=figure, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles//all_full.pdf", sep=""), height=5, width=15)
figure <- ggarrange(p[[1]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    p[[2]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    p[[3]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    legend,
                    labels = c("A", "B", "C", ""),
                    ncol = 4, nrow = 1, widths=c(1,1,1,0.27))
ggsave(plot=figure, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles/all_zoom.pdf", sep=""), height=5, width=15)

figure <- ggarrange(p[[4]] +theme(legend.position="none"), 
                    p[[5]] +theme(legend.position="none"), 
                    p[[6]] +theme(legend.position="none"), 
                    legend,
                    labels = c("A", "B", "C", ""),
                    ncol = 4, nrow = 1, widths=c(1,1,1,0.27))
ggsave(plot=figure, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles//chrX_full.pdf", sep=""), height=5, width=15)
figure <- ggarrange(p[[4]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    p[[5]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    p[[6]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none"), 
                    legend,
                    labels = c("A", "B", "C", ""),
                    ncol = 4, nrow = 1, widths=c(1,1,1,0.27))
ggsave(plot=figure, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles/chrX_zoom.pdf", sep=""), height=5, width=15)







p_combined <- list()
l_count = 1

for (s in seq(2,length(sample_numbers))){
  trees <- list.files(path=paste(this.dir,"/combined/", sample_numbers[[s]], "/", sep=""), pattern="*sub1.*\\.trees$", full.names = TRUE)
  for (i in seq(1,length(trees))){
    in_command <- " -b 0 -log"
    for (j in seq(1,3)){
      in_command = paste(in_command, " ", gsub("sub1", paste("sub", j,sep=""), trees[i]), sep="")
    }
    
    out_command = gsub("sub1_", "", trees[i])

    combined_command = gsub(".trees",".trees", out_command)
    combined_command = paste(" -o ", combined_command, sep="")
    # combine the trees
    system(paste("/Applications/BEAST\\ 2.6.0/bin/logcombiner", in_command, combined_command, "", sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
    in_file  = gsub("-o ","", combined_command)
    out_file = gsub(" ", "", gsub("/combined/","/unique/", in_file))
    
    system(paste("java -jar ", this.dir,  "/../../Software/AIMannotator.jar -burnin 0 -userank false ", in_file, " ", out_file, sep=""), intern=TRUE)
    
    species_trees <- read.nexus(file=out_file,force.multi=T) # get all the tree topologies in the posterior
    # compute total number of trees, only works if Minimal tree support for output was 0 in AIM ann
    tot.trees = 0
    for (j in seq(1, length(species_trees))){
      post.occurance = strsplit(names(species_trees)[[j]], split="_")[[1]]
      tot.trees = tot.trees + as.numeric(post.occurance[[length(post.occurance)]])
    }
    print(i)
    # plots the most likely species tree for each run
    for (j in seq(1, min(1,length(species_trees)))){
      plot_tree <- ladderize(species_trees[[j]], right=F) # ladderizes the tree (defines the node orderign on the x-axis for the plot)
      pl <- plotSpeciesTree(out_file, names(species_trees)[[j]], plot_tree, posterior.threshold, BF.threshold, prior, forwards.arrows) # returns a ggplot object that can be modified
      post.occurance = strsplit(names(species_trees)[[j]], split="_")[[1]] # get how often this ranked topology was observed in the posterior
      pl <- pl + 
        ylab("substitutions per site") + # labels the y-axis
        scale_size_continuous(range=c(0.1,1)) + # defines the range of arrow sizes
        ggtitle(paste("posterior support = ", round(as.numeric(post.occurance[[length(post.occurance)]])/tot.trees,2), "", sep="")) + # specify the title
        theme_minimal() + # changes the plotting theme
        coord_cartesian(xlim=c(0,8.5), ylim=c(0,0.2)) +
        scale_size_continuous(limits=c(0.5,1), range=c(0.1,1), breaks=c(0.5,0.75,1), name="posterior support\nfor gene flow")+
        scale_y_continuous(sec.axis = sec_axis(~ ./(clock_rate), name="Years in Million"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        xlab("")
      plot(pl)
      p_combined[[l_count]] = pl
      l_count = l_count+1
    }
  }
}

legend <- g_legend(p_combined[[1]]) 

figure <- ggarrange(p_combined[[2]] +theme(legend.position="none") + ggtitle("species history using ChrX only"),
                    p_combined[[2]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none") + ggtitle("ChrX only zoomed in"), 
                    p_combined[[1]] + coord_cartesian(xlim=c(2.5,8.5), ylim=c(0,0.015)) + theme(legend.position="none") + ggtitle("ChrX + Chr3 zoomed in"), 
                    legend,
                    labels = c("A", "B", "C", ""),
                    ncol = 4, nrow = 1, widths=c(1,1,1,0.27))
ggsave(plot=figure, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles/all_combined.pdf", sep=""), height=5, width=15)
