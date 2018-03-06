######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)
library(XML)
library("methods")


# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

trees <- list.files(path="./out/greatapes/", pattern="rep9.*\\.trees$", full.names = TRUE)
system("rm -r combined/greatapes")
system("mkdir combined/greatapes")

system("rm -r mcc/greatapes")
system("mkdir mcc/greatapes")
# run log combiner
for (i in seq(1,length(trees))){
  system("rm tmp.trees")
  in_command <- " -trees -burnin 40000000"
  for (j in seq(0,9)){
    in_command = paste(in_command, " ", gsub("rep9", paste("rep", j,sep=""), trees[i]), sep="")
  }
  combined_out = gsub("_rep9", "", trees[i])
  combined_out = gsub("out", "combined", combined_out)
  
  
  system(paste("/Applications//beast1/bin/logcombiner ", in_command, combined_out, sep=" "))
  out_command = gsub("_rep9", "_mcc", trees[i])
  out_command = gsub("out", "mcc", out_command)
  
  system(paste("/Applications/BEAST\\ 2.4.7/bin/treeannotator -heights mean", combined_out, out_command, sep=" "))
  
}


log <- list.files(path="./out/greatapes/", pattern="rep9.*\\.log$", full.names = TRUE)
# run log combiner
for (i in seq(1,length(log))){
  in_command <- ""
  in_command <- " -burnin 40000000"
  for (j in seq(0,9)){
    in_command = paste(in_command, " ", gsub("rep9", paste("rep", j,sep=""), log[i]), sep="")
  }
  combined_out = gsub("_rep9", "", log[i])
  combined_out = gsub("out", "combined", combined_out)
  
  system(paste("/Applications//beast1/bin/logcombiner ", in_command, combined_out, sep=" "))
}


#system("rm tmp.trees")