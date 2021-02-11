######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ggtree)
library(coda)
library(grid)
library(gridExtra)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

sample_numbers = c("s50", "s100", "s200")
numbers = c(50,100,200)


# for (s in seq(1,length(sample_numbers))){
for (s in seq(2,2)){
  
  # analyse the loci trees
  locitrees = list.files(path="./yule/out/", pattern="*\\.log$", full.names = TRUE)
  for (i in seq(1,length(locitrees))){
    t = read.table(locitrees[[i]], header=T, sep="\t")
    hpd = HPDinterval(as.mcmc(t$mrca.age.ph_outgroup1.))
    new.height = data.frame(lower=hpd[1,"lower"], upper=hpd[1,"upper"], mean=mean(t$mrca.age.ph_outgroup1.), run=i)
    if (i==1){height = new.height}
    else{height = rbind(height, new.height)}
  }
  attach(height)
  ordered_height = height[order(height[,3]),]
  height1 = 0.0506
  height2 = 0.0815

  print(paste("cov1=", length(which(ordered_height$upper<height1))/length(ordered_height$lower) ))
  print(paste("cov2=", length(which(ordered_height$upper<height2))/length(ordered_height$lower) ))

  p1 <- ggplot(ordered_height)+
    geom_errorbar(aes(x=seq(1,length(ordered_height$lower)), ymin=lower, ymax=upper )) +
    geom_point(aes(x=seq(1,length(ordered_height$lower)), y=mean), color="blue") +
    geom_hline(yintercept = height1, color="red") + 
    geom_hline(yintercept = height2, color="red") + 
    scale_y_log10(limits=c(0.001,0.5)) +
    theme_minimal()+
    ylab("height relative to mrca with An. christyi")+
    xlab("") +
    ggtitle(" 95% HPD intervals of node heights")
  p2 <- ggplot(ordered_height) +
    geom_violin(aes(x=1, y=mean)) + 
    geom_hline(yintercept = height1, color="red") + 
    geom_hline(yintercept = height2, color="red") + 
    scale_y_log10(limits=c(0.001,0.5)) +
    theme_minimal() +
    ylab("height relative to mrca with An. christyi") +
    xlab("") +
    ggtitle("Distribution of mean node heights relative")
  plot.all = grid.arrange(p1, p2,ncol=2)
  ggsave(plot=plot.all,"/Users/nicmuell/Documents/github/IsolationWithMigration/Figures/Gambia/GamAraHeights.pdf",width=8, height=4)
  
  
  for (i in seq(1,length(locitrees))){
    t = read.table(locitrees[[i]], header=T, sep="\t")
    hpd = HPDinterval(as.mcmc(t$mrca.age.ph_outgroup6.))
    new.height = data.frame(lower=hpd[1,"lower"], upper=hpd[1,"upper"], mean=mean(t$mrca.age.ph_outgroup6.), run=i)
    if (i==1){height = new.height}
    else{height = rbind(height, new.height)}
  }
  attach(height)
  ordered_height = height[order(height[,3]),]
  height1 = 0.0506
  height2 = 0.08291139240506329
  
  print(paste("cov1=", length(which(ordered_height$upper<height1))/length(ordered_height$lower) ))
  print(paste("cov2=", length(which(ordered_height$upper<height2))/length(ordered_height$lower) ))
  

  p1 <- ggplot(ordered_height)+
    geom_errorbar(aes(x=seq(1,length(ordered_height$lower)), ymin=lower, ymax=upper )) +
    geom_point(aes(x=seq(1,length(ordered_height$lower)), y=mean), color="blue") +
    geom_hline(yintercept = height1, color="red") + 
    geom_hline(yintercept = height2, color="red") + 
    scale_y_log10(limits=c(0.001,0.5)) +
    theme_minimal()+
    ylab("height relative to mrca with An. christyi")+
    xlab("") +
    ggtitle(" 95% HPD intervals of node heights")
  
  scale_y_log10()
  p2 <- ggplot(ordered_height) +
    geom_violin(aes(x=1, y=mean)) + 
    geom_hline(yintercept = height1, color="red") + 
    geom_hline(yintercept = height2, color="red") + 
    scale_y_log10(limits=c(0.001,0.5)) +
    theme_minimal() +
    ylab("height relative to mrca with An. christyi") +
    xlab("") +
    ggtitle("Distribution of mean node heights relative")
  
  plot.all = grid.arrange(p1, p2,ncol=2)
  ggsave(plot=plot.all,"/Users/nicmuell/Documents/github/IsolationWithMigration/Figures/Gambia/MerQuaHeights.pdf",width=8, height=4)
  
  heights <- list.files(path=paste("./loci/ancestors/", sample_numbers[[s]], "/", sep=""), pattern="*\\.csv$", full.names = TRUE)
  
  for (i in seq(1,length(heights))){
    t = read.table(heights[[i]], header=T, sep=",")
    
    
    t = t[-which(t$species1=="AchrA1_AchrA1" | t$species2=="AchrA1_AchrA1"), ]
    t = t[-which(t$species1=="AepiE1_AepiE1" | t$species2=="AepiE1_AepiE1"), ]
    t = t[-which(t$species1=="AcolM1_AcolM1" | t$species2=="AcolM1_AcolM1"), ]
    t = t[-which(t$species1=="AmelC2_AmelC2" | t$species2=="AmelC2_AmelC2"), ]
    
    t$x = runif(length(t$species1), 0.9,1.1)
    
    gene = t[which(t$tree=="gene"),]
    species = t[which(t$tree=="species"),]
    
    p = ggplot(gene, aes(x=x,y=distance))+
      geom_point(color="black")+
      geom_segment(data=species, aes(x=0.9, xend=1.1, y=distance, yend=distance), color="red")+
      scale_alpha_manual(values = c("gene"=0.1, "species"=1))+
      theme_light()+
      facet_grid(species1~species2)
    plot(p)
    dsa
  }
}



