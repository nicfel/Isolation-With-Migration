######################################################
######################################################
# plot the sampling bias output after checking that 
# every combined has an ESS value of over 200
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
f <- list.files(path="./out", pattern="*.log", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


nenr <- 1
mnr <- 1

# Read In Data ---------------------------------
mig_rate <- data.frame()
counter <- 1
for (i in seq(1,length(f),1)){
  print(i)
  filename1 <- f[i]
  
  bpp_name <- gsub("out","BPP", f[i])
  bpp_name <- gsub("IM_500Genes","BPP", bpp_name)
  bpp_name <- gsub(".log",".txt", bpp_name)
  
  
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
  
  
  # Read in the asco *.logs and combine all 3 runs that had different initial values
  t <- read.table(filename1, header=TRUE, sep="\t")

  # combine the data after a burn in of 20%
  t <- t[-seq(1,ceiling(length(t$N1)/5)), ]


  tmp <- strsplit(filename1,'_')

  ess <- effectiveSize(t)
  if (min(ess[2])<100){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(ess[2:3]),filename1))
  }else{
    for (i in seq(1,length(Ne_true))){
      name <- paste("Ne",i,sep="")
      rates.new <- data.frame(est=median(t[, name]), true=Ne_true[i]/8, col=Ne_color[i], rate="Ne")
      if(nenr==1){
        rates = rates.new
        nenr <- 2
      }else{
        rates <- rbind(rates, rates.new)
      }
    }
    for (i in seq(1,length(m_true))){
      name <- paste("m",i,sep="")
      rates.new <- data.frame(est=median(t[, name]*t[,"eM"])/0.05, true=m_true[i], col=m_color[i], rate="migration")
      rates <- rbind(rates, rates.new)
    }
  }
}


attach(rates)
rates_sorted <- rates[order(rate,col),]
detach(rates)

rates_red_sorted <- rates_sorted[which(rates$true>0.001),]

# make the data frame for the red dotted line
d <- data.frame(x=c(0.1, 10, 0.1, 10, 0.1, 10),
                y=c(0.1, 10, 0.1, 10, 0.1, 10),
                col=c("sampled","sampled","first ancestral","first ancestral","second ancestral","second ancestral"),
                rate=c("migration","migration","migration","migration","migration","migration"))


limits = data.frame(rate=c("Ne","Ne","Ne","Ne","migration","migration","migration"),
                    col=c("sampled","first ancestral","second ancestral","root","sampled","first ancestral","second ancestral"),
                    upper=c(5,5,5,5,10,10,10),lower=c(0.1,0.1,0.1,0.1,0.001,0.001,0.001))

p_mig <- ggplot()+
  geom_point(data=rates_red_sorted, aes(x=true, y=est))+
  scale_x_log10() + scale_y_log10() +
  ylab("estimated") + xlab("true")  + 
  theme(legend.position="none")+
  facet_wrap(rate ~ col, ncol=4, scales="free") 
p_mig <- p_mig +  geom_blank(data=limits,aes(y = upper, x=upper)) +
  geom_blank(data=limits, aes(y = lower, x=lower))


plot(p_mig)


ggsave(plot=p_mig,"../../text/figures/Rates_rates.pdf",width=10, height=5)























# make the data frame for the red dotted line
d <- data.frame(x=c(0.01, 1, 0.01, 1, 0.01, 1, 0.01, 1),
                y=c(0.01, 1, 0.01, 1, 0.01, 1, 0.01, 1),
                facet_posx=c("sampled","sampled","first ancestral","first ancestral","second ancestral","second ancestral", "root", "root"))

p_Ne <- ggplot()+
  geom_point(data=Ne, aes(x=true/8, y=est, color=col))+
  scale_y_log10(limits=c(0.05,0.2)) + scale_x_log10(limits=c(0.4,1.6)) +
  ylab("estimated") + xlab("true") + 
  theme(legend.position="none")+
  facet_grid(col ~ .) 
p_Ne <- p_Ne + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed") +
  facet_grid(facet_posx ~ .)+
  scale_colour_manual("",values = c("sampled" = col0, "first ancestral" = col1, "second ancestral" = col2, "root" = col3)) + guides(colour = guide_legend(override.aes = list(size=2)))


plot(p_Ne)





tree <- read.tree(text = "((A:3,B:3):7,C:10):2;")
#setEPS()
#postscript("../../text/figures/speciesTree.eps",width=5, height=5)
p <- plot(tree, type = "cladogram",root.edge = TRUE,direction = "downwards",show.tip.label=FALSE)
#dev.off()


#p_lab <- ggplot()+
#  geom_violin(data=rates, aes("N_a",N1),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
#  facet_grid(symmetry ~ .) + labs(x=expression(Ne[a]~Ne[b]~Ne[c]~Ne[d]~Ne[e]),
#                                  y=expression(A~B~C)) 
#plot(p_lab)
#ggsave(plot=p_lab,"../../text/figures/Labels.eps",width=5, height=1)

