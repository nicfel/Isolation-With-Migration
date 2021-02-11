######################################################
######################################################
# combine and plot AIM runs
######################################################
######################################################
library(ape) # laod ape
library(ggplot2) # load ggplot

this.dir <- dirname(parent.frame(2)$ofile) # Set get current directory
sample_numbers = c("s50", "s200")# define the sample numbers

p <- list()
l_count = 1
# for (s in seq(1,1)){
  
for (s in seq(2,length(sample_numbers))){

  logs <- list.files(path=paste(this.dir,"/unique/", sample_numbers[[s]], "/", sep=""), pattern="*sub.*occurances.*\\.log$", full.names = TRUE)
  first = T
  for (i in seq(1,length(logs))){
    t = read.table(file=logs[[i]], header=T, sep="\t")
    labs = labels(t)[[2]]
    
    tmp = strsplit(logs[[i]], split="_")[[1]]
    
    if (tmp[[2]]=="all"){
      data = "ChrX + Chr3"
    }else{
      data = "ChrX only"
    }
    
    for (j in seq(1,length(t))){
      firrst.extant = T
      if (startsWith(labs[[j]], "Ne")){
        hpd=HPDinterval(as.mcmc(t[,j]))
        name = gsub("Ne_", "", labs[[j]])
        if (length(strsplit(name, split='\\.')[[1]])>1){
          new.dat = data.frame(Ne = t[,j], data=data, subset=gsub("sub", "", tmp[[4]]), name=gsub('\\.','\n',name))
          if (first){
            plot.dat = new.dat
            first = F
          }else{
            plot.dat = rbind(plot.dat, new.dat)
          }
        }else if(firrst.extant){
          hpd=HPDinterval(as.mcmc(t[,j]))
          name = "extant"
          new.dat = data.frame(Ne = t[,j], data=data, subset=gsub("sub", "", tmp[[4]]), name=gsub('\\.','\n',name))
          if (first){
            plot.dat = new.dat
            first = F
          }else{
            plot.dat = rbind(plot.dat, new.dat)
          }
          firrst.extant=F
        }
      }
    }
  }
  p = ggplot(plot.dat) + 
    geom_violin(aes(x = name, y=Ne, group=interaction(name,subset))) +
    facet_wrap(.~data, ncol=1) +
    scale_y_log10() +
    theme_minimal() +
    xlab("") +
    ylab("Ne in substitutions per site")
  
  plot(p)
  ggsave(plot=p, paste("/Users/nmueller/Documents/github/IsolationWithMigration-Text/Text/figures/anopheles/Ne.pdf", sep=""), height=5, width=15)
  
}
