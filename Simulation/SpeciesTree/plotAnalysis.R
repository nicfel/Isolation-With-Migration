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


# clear workspace
rm(list = ls())

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the mutation rates
mutation_AIM <- read.table(file='AIM_mutation.tsv', sep='\t', header = TRUE, row.names = 1)
mutation_starBeast <- read.table(file='StarBeast_mutation.tsv', header = TRUE, row.names = 1)

d <- data.frame(x=c(0.01,10),
                y=c(0.01,10))

p_AIM_mut <- ggplot(data=mutation_AIM)+
  geom_point(aes(x=true,y=est), color="black", size=0.001) +
  scale_y_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) + scale_x_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) +
  ylab("estimated mutation rate") + xlab("true mutation rate")
p_AIM_mut <- p_AIM_mut + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed")

p_SB_mut <- ggplot(data=mutation_starBeast)+
  geom_point(aes(x=true,y=est), color="black", size=0.001) +
  scale_y_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) + scale_x_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) +
  ylab("estimated mutation rate") + xlab("true mutation rate")
p_SB_mut <- p_SB_mut + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed")

combined  <- data.frame(AIM=mutation_AIM$est, SB=mutation_starBeast$est)
p_CB_mut <- ggplot(data=combined)+
  geom_point(aes(x=AIM,y=SB), color="black", size=0.001) +
  scale_y_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) + scale_x_log10(breaks=c(0.01,0.1,1,10),limits = c(0.01,10)) +
  ylab("estimated mutation rate") + xlab("true mutation rate")
p_CB_mut <- p_CB_mut + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed")

plot(p_AIM_mut)
plot(p_SB_mut)
plot(p_CB_mut)

###################################
###################################
###################################
###################################


#ggsave(plot=p_mut,"../../text/figures/SpeciesTree_mutation.eps",width=3, height=3)
# read in the rest of the data
AIM_param <- read.table(file='AIM_nodeHeights.tsv', sep='\t', header = TRUE, row.names = 1)
StarBeast_param <- read.table(file='StarBeast_param.tsv', header = TRUE, row.names = 1)

top_AIM = order(StarBeast_param$true_topology)
top_SB = order(StarBeast_param$true_topology)

# make dataframes
topology.aim <- data.frame(x = rep(seq(1,length(top_AIM)*2,2), AIM_param$true_topology[top_AIM]), col="AIM")
topology.msc <- data.frame(x = rep(seq(2,length(top_SB)*2,2), StarBeast_param$true_topology[top_SB]), col="MSC")

p_top <- ggplot() + 
  geom_histogram(data=topology.aim,aes(x, fill=col), binwidth=1) +
  geom_histogram(data=topology.msc,aes(x, fill=col), binwidth=1) +
  ylab("support for the true topology") + xlab("replicate")+
  scale_fill_manual("",values = c("AIM" = col0, "MSC" = col1)) + 
  theme(legend.position="none") + 
  scale_x_continuous(breaks=c(0,100,200), labels=c("0", "50", "100")) +
  scale_y_continuous(breaks=c(0,451,902), labels=c("0.0", "0.5", "1.0"))

plot(p_top)

ggsave(plot=p_top,"../../text/figures/SpeciesTree_topology.pdf",width=5, height=3)


p_node <- ggplot()+
  geom_violin(data=AIM_param, aes(1,n1, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=AIM_param, aes(3,n2, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=AIM_param, aes(5,n3, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(2,n3, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(4,n2, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(6,n1, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_segment(data=StarBeast_param, aes(x = 0.6, y = 0.0025, xend = 2.4, yend = 0.0025), color="red", linetype="dashed") +
  geom_segment(data=StarBeast_param, aes(x = 2.6, y = 0.005, xend = 4.4, yend = 0.005), color="red", linetype="dashed") +
  geom_segment(data=StarBeast_param, aes(x = 4.6, y = 0.01, xend = 6.4, yend = 0.01), color="red", linetype="dashed") +
  ylab("estimated node height")+
  scale_x_discrete("",limits=(c(1.5,3.5,5.5)), labels = c(expression("AB"),expression("ABC"),expression("ABCD")))+
  scale_colour_manual("",values = c("AIM" = col0, "StarBeast" = col1),
                      breaks=c("AIM", "StarBeast")) + 
  theme(legend.position="none")

plot(p_node)
ggsave(plot=p_node,"../../text/figures/SpeciesTree_nodeHeight.pdf",width=5, height=3)

tree <- read.tree(text = "(((A:6,B:6):5,C:11):10,D:21):2;")
setEPS()
postscript("../../text/figures/SpeciesTree_speciesTree.eps",width=5, height=5)
p <- plot(tree, type = "cladogram",root.edge = TRUE,direction = "downwards",show.tip.label=FALSE)
dev.off()


# # make the data frame for the red dotted line
# d <- data.frame(x=c(0.1, 10, 0.1, 10, 0.1, 10),
#                 y=c(0.1, 10, 0.1, 10, 0.1, 10),
#                 col=c("sampled","sampled","first ancestral","first ancestral","second ancestral","second ancestral"),
#                 rates=c("migration","migration","migration","migration","migration","migration"))


AIM_rates <- read.table(file='AIM_rates.tsv', sep='\t', header = TRUE, row.names = 1)

AIM_rates$col_f = factor("`first interval`", levels=c('`first interval`','`second interval`','`third interval`','root'))

AIM_rates$rate_f = factor( "N[e]*mu", levels=c("N[e]*mu","migration"))
AIM_rates[which(AIM_rates$rate=='migration'), "rate_f"] = 'migration'


AIM_rates[which(AIM_rates$col=='sampled'), "col_f"] = '`first interval`'
AIM_rates[which(AIM_rates$col=='first ancestral'), "col_f"] = '`second interval`'
AIM_rates[which(AIM_rates$col=='second ancestral'), "col_f"] = '`third interval`'
AIM_rates[which(AIM_rates$col=='root'), "col_f"] = 'root'


AIM_red_rates <- AIM_rates[-intersect(which(AIM_rates$true<0.001),which(AIM_rates$rate=="migration")),]

limits = data.frame(rate_f=c("N[e]*mu","N[e]*mu","N[e]*mu","N[e]*mu","migration","migration","migration"),
                    col_f=c("`first interval`","`second interval`","`third interval`","root","`first interval`","`second interval`","`third interval`"),
                    upper=c(0.001,0.001,0.001,0.001,10,10,10),lower=c(0.0001,0.0001,0.0001,0.0001,0.001,0.001,0.001))

p_mig <- ggplot()+
  geom_point(data=AIM_red_rates, aes(x=true, y=est))+
  scale_x_log10() + scale_y_log10() +
  ylab("estimated") + xlab("true")  + 
  theme(legend.position="none")+
  facet_wrap(rate_f ~ col_f, ncol=4, scales="free", labeller = label_parsed) 
p_mig <- p_mig +  geom_blank(data=limits,aes(y = upper, x=upper)) +
  geom_blank(data=limits, aes(y = lower, x=lower))


plot(p_mig)


ggsave(plot=p_mig,"../../text/figures/SpeciesTree_migration.pdf",width=10, height=5)




p_node <- ggplot()+
  geom_violin(data=AIM_param, aes(1,n1, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=AIM_param, aes(3,n2, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=AIM_param, aes(5,n3, color="AIM"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(2,n3, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(4,n2, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=StarBeast_param, aes(6,n1, color="StarBeast"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_segment(data=StarBeast_param, aes(x = 0.6, y = 0.0025, xend = 2.4, yend = 0.0025), color="red", linetype="dashed") +
  geom_segment(data=StarBeast_param, aes(x = 2.6, y = 0.005, xend = 4.4, yend = 0.005), color="red", linetype="dashed") +
  geom_segment(data=StarBeast_param, aes(x = 4.6, y = 0.01, xend = 6.4, yend = 0.01), color="red", linetype="dashed") +
  ylab("estimated node height")+
  scale_x_discrete("",limits=(c(1.5,3.5,5.5)), labels = c(expression("AB"),expression("ABC"),expression("ABCD")))+
  scale_colour_manual("",values = c("AIM" = col0, "StarBeast" = col1),
                      breaks=c("AIM", "StarBeast")) + 
  theme(legend.position="top")

plot(p_node)
ggsave(plot=p_node,"../../text/figures/SpeciesTree_legend.pdf",width=5, height=3)

