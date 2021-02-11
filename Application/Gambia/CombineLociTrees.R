######################################################
######################################################
# combine the gene tree runs and run the mcc trees
######################################################
######################################################
library(ape)
# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

sample_numbers = c("s50", "s100", "s200")

system("rm -r loci/combined")
system("mkdir loci/combined")

system("rm -r loci/mcc")
system("mkdir loci/mcc")


# for (s in seq(1,length(sample_numbers))){
for (s in seq(3,3)){
    
  system(paste("rm -r loci/combined/", sample_numbers[[s]], sep=""))
  system(paste("mkdir loci/combined/", sample_numbers[[s]], sep=""))
  
  system(paste("rm -r loci/mcc/", sample_numbers[[s]], sep=""))
  system(paste("mkdir loci/mcc/", sample_numbers[[s]], sep=""))
  
  
  trees <- list.files(path=paste("./loci/individual/", sample_numbers[[s]], "/", sep=""), pattern="*rep0.*\\.trees$", full.names = TRUE)
  
  for (i in seq(1,length(trees))){
    print(trees[[i]])
    in_command <- " -b 20 -resample 1000000 -log"
    for (j in seq(0,2)){
      in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), trees[i]), sep="")
    }
    
    out_command = gsub("rep0_", "", trees[i])
    out_command = gsub("individual", "combined", out_command)
    
    combined_command = gsub(".trees",".trees", out_command)
    combined_command = paste(" -o ", combined_command, sep="")
    # combine the trees
    system(paste("/Applications/BEAST\\ 2.5.2/bin/logcombiner", in_command, combined_command, "", sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
    
    in_tree = gsub("-o","", combined_command)
    system(paste("/Applications/BEAST\\ 2.5.2/bin/treeannotator -burnin 0 -heights mean ", in_tree, gsub("combined","mcc", in_tree), sep=" "), intern=TRUE, show.output.on.console = F, ignore.stderr=T)
  }
}


