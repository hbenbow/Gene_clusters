library(grid)
library(gridExtra)
setwd("~/Documents/Hotspots/Paper_version_4/QTLs/")
Fg_final <- read.csv("~/Documents/Hotspots/Paper_version_4/Fg_final.csv")
all_markers_positions <- read.csv("~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv", row.names=1)
final_hotspots <- read.csv("~/Documents/Hotspots/Paper_version_4/final_hotspots.csv")
QTL_database <- read.csv("~/Documents/Hotspots/Paper_version_4/QTLs/QTL_database.csv")
QTL_database$Chromosomes<-paste("chr", QTL_database$Chromosomes, sep="")

all_qtls<-list()
qtls<-QTL_database$General.number
for(qtl in qtls){
  data<-QTL_database[(QTL_database$General.number==qtl),]
  marker<-as.character(data$Linked.markers)
  chr<-data$Chromosomes
  position<-subset(all_markers_positions, all_markers_positions$Chromosome %in% chr)
  position<-subset(position, position$Feature %in% marker)
   clusters<-subset(final_hotspots, final_hotspots$Chromosome %in% chr)
  len<-nrow(clusters)
  lenp<-nrow(position)
  if(lenp>0 && len>0){
      clusters$qtl<-paste(qtl)
      clusters$marker<-paste(unique((position$Feature)))
      clusters$marker_position<-as.numeric(paste(unique(position$Start)))
      clusters$difference<-abs(clusters$start-clusters$marker_position)/1000000
      all_qtls[[length(all_qtls)+1]]<-clusters
    }else{}
  } 


qtls_clusters<-do.call(rbind.data.frame, all_qtls)


