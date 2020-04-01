library(grid)
library(gridExtra)
setwd("~/Documents/Hotspots/Paper_version_4/Markers/Meta_qtl_maps/")
chroms<-dir(pattern="map")
chromosome<-substr(chroms, 1, 2)

maps<-list()
qtl_list<-list()
all<-list()
for(i in chromosome){
  map<-read.delim(file=paste(i, "_Map.txt.map.txt", sep=""), header=F)
  colnames(map)<-c("order", "Marker", "cM")
  map$Chromosome<-paste(i)
  maps[[length(maps)+1]] = map
  qtl<-read.delim(paste(i, "_QTL.txt.qtl", sep=""), header=F)
  colnames(qtl)<-c("QTL", "Resistance", "Ontology", "Place", "Year", "Chromosome" ,"Linkage_group", "LOD", "r2", "Peak", "cM_Start", "cM_end")
  qtl_list[[length(qtl_list)+1]]= qtl
  qtls<-unique(qtl$QTL)
  List<-list()
  for(q in qtls){
    data<-qtl[(qtl$QTL==q),]
    min<-data$cM_Start
    max<-data$cM_end
    qtl_position<-map[!(map$cM<=min),]
    qtl_position<-qtl_position[!(qtl_position$cM>=max),]
    qtl_position$QTL<-paste(q)
    qtl_position$Chromosome<-paste(i)
    List[[length(List)+1]]=qtl_position
  }
  df<-do.call(rbind.data.frame, List)
  all[[length(all)+1]]<-df
}
maps<-do.call(rbind.data.frame, maps)
maps$Chromosome<-paste("chr", maps$Chromosome, sep="")
all<-do.call(rbind.data.frame, all)
all$Chromosome<-paste("chr", all$Chromosome, sep="")
qtls<-do.call(rbind.data.frame, qtl_list)
qtls$Chromosome<-paste("chr", qtls$Chromosome, sep="")

dir.create("~/Documents/Hotspots/Paper_version_4/QTLs")
setwd("~/Documents/Hotspots/Paper_version_4/QTLs")

physical_maps<-list()
qtls_list<-list()
for(i in 1:nrow(qtls)){
  data<-qtls[i,]
  qtlid<-data$QTL
  min<-min(data$cM_Start)
  max<-max(data$cM_end)
  peak<-data$Peak
  lod<-data$LOD
  chromosome<-data$Chromosome
  data2<-maps[(maps$Chromosome==chromosome),]
  data2<-data2[(data2$cM >=min),]
  data2<-data2[(data2$cM <= max),]
  data2$LOD<-lod
  coords<-merge(data2, all_markers_positions, by.x="Marker", by.y="Feature")
  coords<-as.data.frame(cbind(coords, ifelse(coords$Chromosome.x==coords$Chromosome.y,1,0)))
  coords<-coords[(coords$`ifelse(coords$Chromosome.x == coords$Chromosome.y, 1, 0)` == 1),]
  if(nrow(coords)>0){
    phys_map<-coords[,c(4, 3, 7, 8, 1, 5)]
    phys_map$cM<-paste(data$QTL)
    phys_map<-phys_map[,c(2, 1, 3, 4, 5, 6)]
    physical_maps[[length(physical_maps)+1]]<-phys_map
    newmin<-min(phys_map$Start)
    newmax<-max(phys_map$End)
    lod<-phys_map$LOD
    hotspot_data<-final_hotspots[final_hotspots$Chromosome==chromosome,]
    hotspot_data<-hotspot_data[(hotspot_data$start >= newmin),]
    hotspot_data<-hotspot_data[(hotspot_data$end <= newmax),]
    hotspot_data<-hotspot_data[,c(1, 2, 3, 4, 5, 7)]
    colnames(phys_map)<-colnames(hotspot_data)
    final<-rbind(hotspot_data, phys_map)
    final<-final[order(final$start),]
    final<-final[!duplicated(final$GeneID),]
    final<-merge(final, final_hotspots, by="GeneID", all.x=T, fill=1)
    final$Density.y[is.na(final$Density)] <- 0
    final<-final[order(final$start.x),]
    if(length(paste(unique(is.na(final$Density.y))))>1){
      final<-final[,1:6]
      write.csv(final, file=paste(qtlid, chromosome, ".csv", sep="_"))
      qtls_list[[length(qtls_list)+1]] <- final
    }else{
    }
  }
}

QTL.physical.maps<-do.call(rbind.data.frame, physical_maps)
QTL.positions<-do.call(rbind.data.frame, qtls_list)
for(i in paste(unique(QTL.positions$Chromosome.x))){
  data<-QTL.positions[(QTL.positions$Chromosome.x==i),]
  data<-data[order(data$start.x),]
  data<-data[!duplicated(data$GeneID),]
  write.csv(data, file=paste(i, "_QTLs.csv"))
}

write.csv(QTL.positions, file="All_QTL_and_Hotspots.csv")


physical_map<-merge(maps, all_markers_positions, by="Marker")
physical_map<-as.data.frame(cbind(physical_map, ifelse(physical_map$Chromosome.x==physical_map$Chromosome.y,1,0)))
physical_map<-physical_map[(physical_map$`ifelse(physical_map$Chromosome.x == physical_map$Chromosome.y, `==1),]
physical_map<-physical_map[!(duplicated(physical_map$Marker)),]

ggplot(physical_map, aes(x=cM, y=Start)) + geom_point()+
  facet_wrap(~Chromosome.x, ncol=3) +
  theme(text = element_text(size=20, colour="black")) + 
 theme_bw()

ggplot(Fg_final[!(Fg_final$Chromosome=="chrUn"),], 
       aes(x=end, y=density, group=Colour)) +geom_jitter(size=1.6, aes(colour=Colour),alpha=0.6) + 
  facet_wrap(~Chromosome, ncol=3) + xlab("Position (bp)") + 
  ylab("Gene Density") + theme_bw() +
  geom_hline(yintercept=0.7, alpha=0.7) +
  theme(text = element_text(size=16, colour="black")) + 
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual( values=c("grey60", "orangered2"), labels=c("No cluster", "Cluster"))

ggplot(physical_qtls, aes(x=End.x, y=LOD)) +geom_line() + 
  facet_wrap(~Chromosome.x.x, ncol=3) + xlab("Position (bp)") 
  ylab("Gene Density") + theme_bw() +
  geom_hline(yintercept=0.7, alpha=0.7) +
  theme(text = element_text(size=16, colour="black")) + 
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual( values=c("grey60", "orangered2"), labels=c("No cluster", "Cluster"))

for(i in unique(physical_qtls$Chromosome.x.x)){
  q<-physical_qtls[(physical_qtls$Chromosome.x.x==i),]
  c<-Fg_final[(Fg_final$Chromosome==i),]
  
qtl_plot<-ggplot(q, aes(x=End.x, y=LOD)) +geom_jitter(size=1.6, alpha=0.6) + xlab("Position (bp)") + 
  ylab("LOD Score")  + theme_bw() + 
  geom_hline(aes(yintercept=3))+
  theme(text = element_text(size=16, colour="black"))

cluster_plot<-ggplot(c, aes(x=end, y=density)) +geom_jitter(size=1.6, alpha=0.6, aes(colour=Colour)) + xlab("Position (bp)") + 
    ylab("Density")  + theme_bw() + 
    geom_hline(aes(yintercept=0.7))+
  theme(text = element_text(size=16, colour="black"), legend.position = "none")+
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual( values=c("grey60", "orangered2"), labels=c("No cluster", "Cluster"))

gq <- ggplotGrob(qtl_plot)
gc <- ggplotGrob(cluster_plot)

pdf(file=paste(i, ".pdf", sep=""))
plot(gtable_rbind(gq, gc))
dev.off()
}
