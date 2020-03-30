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


dir.create("~/Documents/Hotspots/Paper_version_4/QTLs")
setwd("~/Documents/Hotspots/Paper_version_4/QTLs")

qtls_list<-list()
for(i in 1:nrow(qtls)){
  data<-qtls[i,]
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
    qtls_list[[length(qtls_list)+1]] <- final
  }else{
  }
}

QTL.positions<-do.call(rbind.data.frame, qtls_list)


ggplot(final, aes(x=end.x, y=Density)) +geom_jitter(size=1.6, aes(colour=Hotspot.x), alpha=0.6) + xlab("Position (bp)") + 
  ylab("Gene Density")  + theme_bw() + 
  geom_hline(yintercept=0.7, alpha=0.7) +
  theme(text = element_text(size=16, colour="black")) + 
  coord_cartesian(ylim=c(0,1))
scale_color_manual( values=c("grey60", "orangered2"))


for(i in unique(final_hotspots$Hotspot)){
  data<-final_hotspots[final_hotspots$Hotspot== i,]
  data<-data[,c(1, 3, 4, 5, 6, 7, 8, 2)]
  min<-min(data$start)
  max<-max(data$end)
  chromosome<-unique(data$Chromosome)
  qtl_data<-QTL_bed[(QTL_bed$Chromosome==chromosome),]
  qtl_data<-qtl_data[(qtl_data$start >= min),]
  qtl_data<-qtl_data[(qtl_data$end <= max),]
  qtl_data<-qtl_data[,c(1, 2, 3, 4, 5, 6, 7, 9)]
  colnames(qtl_data)<-colnames(data)
  hotspot<-rbind(data, qtl_data)
  hotspot<-hotspot[order(hotspot$start),]
  chromosome<-unique(qtl_data$Chromosome)
  qtl<-paste(qtl_data$Fg.response)
  if(length(qtl) > 0){
    all_qtl_info<-qtls[(qtls$Chromosome==chromosome),]
    all_qtl_info<-all_qtl_info[(all_qtl_info$QTL==qtl),]
    write.csv(hotspot, file=paste(i, "hotspot.csv", sep=""), row.names=F)
    write.csv(all_qtl_info, file=paste(chromosome, "hotspot.csv", sep=""), row.names=F)
  }
  # ggplot(data, aes(x=end, y=Density)) +geom_jitter(size=1.6, alpha=0.6) + xlab("Position (bp)") + 
  #   ylab("Gene Density") + ggtitle(paste(chromosome, qtl)) + theme_bw() + 
  #   geom_hline(yintercept=0.7, alpha=0.7) +
  #   geom_vline(xintercept=qtl_data$start) +
  #   theme(text = element_text(size=16, colour="black")) + 
  #   coord_cartesian(ylim=c(0,1))+
  #   scale_color_manual( values=c("grey60", "orangered2"))
  # ggsave(file=paste(chromosome, qtl, '.pdf', sep=""))
}

genes<-list()
for(i in 1:nrow(all_markers)){
  data<-all_markers[i,]
  marker<-data[,1]
  sequence<-as.character(data[,2])
  marker<-paste(">", marker, sep="")
  gene<-rbind(marker, sequence)
  genes[[length(genes)+1]] = gene
}
all_SSRs<-do.call(rbind.data.frame, genes)

