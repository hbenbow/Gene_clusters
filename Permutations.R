library(ggplot2)
require(zoo)
library(plyr)
library(dplyr)
library(tidyr)

# Read in gene coordinate file, containing Chromosome name, gene start posiiton, gene end
# position, gene id, and strand. Do not need strand.

all <- read.csv("../wheat_all.csv", header=T)
colnames(all)<-c( "Chromosome", "start", "end", "GeneID", "Score", "strand")
# Read in expression scores. This is a .csv file in which each row is a gene
# an d order must correspond to the order in the gene coordinate file. Each column after 
# the gene column is 
expression_scores<-read.csv("../expression_scores.csv", header=T)


diseases<-colnames(expression_scores)[-1]
permutations=1000
dir.create("Permutations")
# this section reads in all files and does the window analysis
# test
correlations_consec<-list() 
correlations_density<-list()
stress_together<-list()
conseq<-list()
conseq2<-list()
positives2<-list()
Density<-list()
predictionsc<-list()
predictionsd<-list()

for (d in diseases){
  perm = list() 
  List<-list()
  selection<-as.data.frame(expression_scores[,1])
  selection[[d]]<-expression_scores[[d]]
  colnames(selection)<-c("GeneID", "Disease")
  bed <- join(all, selection, by="GeneID")
  bed[is.na(bed)]<-0
  shuff<-transform(bed, Disease=sample(Disease))
  shuff$density<-rollapply(shuff$Disease, width=10, FUN=mean, by.column=FALSE, fill=0)
  shuff$Consecutive<-sequence(rle(as.character(shuff$Disease))$lengths)
  shuff[is.na(shuff)]<-0
  perm[[length(perm)+1]] = shuff
  
  for(i in seq_along(perm)){
    shuffled<-perm[[i]]
    positives<-shuffled
    positives<-as.data.frame(positives$density)
    colnames(positives)<-"density"
    counts_density<-count(positives, density)
    positives2[[length(positives2)+1]] = counts_density
    consec<-shuffled[!(shuffled$Disease == 0.0 ),]
    counts_consec<-count(consec, Consecutive)
    conseq[[length(conseq)+1]] = counts_consec
  }
  
  positives<-as.data.frame(do.call(rbind.data.frame, positives2))
  counts_density<-aggregate(positives$n, by=list(positives$density), FUN=sum)
  colnames(counts_density)<-c("density", "n")
  consecs<-as.data.frame(do.call(rbind.data.frame, conseq))
  counts_consec<-aggregate(consecs$n, by=list(consecs$Consecutive), FUN=sum)
  colnames(counts_consec)<-c("Consecutive", "n")
  counts_consec$Stress<-paste(d)
  counts_density$Stress<-paste(d)
  conseq2[[length(conseq2)+1]] = counts_consec
  Density[[length(Density)+1]] = counts_density
  
}

conseqs<-do.call(rbind.data.frame, conseq2)
densitys<-do.call(rbind.data.frame, Density)
write.csv(conseqs, file="Permutations/consec.csv", row.names=F)
write.csv(densitys, file="Permutations/density.csv", row.names=F)

