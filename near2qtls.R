library(grid)
library(gridExtra)
setwd("~/Documents/Hotspots/Paper_version_4/QTLs/")
Fg_final <- read.csv("~/Documents/Hotspots/Paper_version_4/Fg_final.csv")
all_markers_positions <- read.csv("~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv", row.names=1)
final_hotspots <- read.csv("~/Documents/Hotspots/Paper_version_4/final_hotspots.csv")
QTL_database <- read.csv("~/Documents/Hotspots/Paper_version_4/QTLs/QTL_database.csv")
QTL_database$Chromosomes<-paste("chr", QTL_database$Chromosomes, sep="")

number<-unique(subset
               (all_markers_positions, 
                 all_markers_positions$Feature %in% QTL_database$Linked.markers))
table(number$Platform)

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

for(i in qtls_clusters$Hotspot){
  data<-qtls_clusters[(qtls_clusters$Hotspot==i),]
  min<-min(data$marker_position) - 100000
  max<-max(data$marker_position) + 100000
  chr<-as.character(data$Chromosome)


ggplot(Fg_final[(Fg_final$Chromosome==chr),], aes(x=end, y=density)) +
  geom_jitter(size=1.6, aes(colour=Colour), alpha=0.6) +
  xlab("Position (bp)") + 
  ylab("Gene Density") + theme_bw() +
  geom_hline(yintercept=0.7, alpha=0.7) +
  theme(text = element_text(size=16, colour="black")) +
  scale_color_manual( values=c("grey60", "orangered2")) +
  # coord_cartesian(ylim=c(0,1), xlim=c(min, max))+
  geom_vline(data =data, aes(xintercept=marker_position), colour="green")
}
  
  