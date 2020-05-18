library(grid)
library(gridExtra)
setwd("~/Documents/Hotspots/Paper_version_4/QTLs/")
Fg_final <- read.csv("~/Documents/Hotspots/Paper_version_4/Fg_final.csv")
all_markers_positions <- read.csv("~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv", row.names=1)
final_hotspots <- read.csv("~/Documents/Hotspots/Paper_version_4/final_hotspots.csv")
QTL_database <- read.csv("~/Documents/Hotspots/Paper_version_4/QTLs/QTL_database.csv")

all_qtls<-list()
qtls<-QTL_database$General.number
for(qtl in qtls){
  data<-QTL_database[(QTL_database$General.number==qtl),]
  marker<-as.character(data$Linked.markers)
  position<-all_markers_positions[(all_markers_positions$Feature=marker),]$Start
  chromosome<-as.character(all_markers_positions[(all_markers_positions$Feature=marker),]$Chromosome)
  clusters<-final_hotspots[(final_hotspots$Chromosome==chromosome),]
  clusters$QTL_position<-position
  clusters$difference<-abs(clusters$start-clusters$QTL_position)/1000000
  clusters$qtl<-data$General.number
  all_qtls[[length(all_qtls)+1]]<-clusters
  }
qtls_clusters<-do.call(rbind.data.frame, all_qtls)
qtls_clusters<-na.omit(qtls_clusters)