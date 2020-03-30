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
all<-do.call(rbind.data.frame, all)
qtls<-do.call(rbind.data.frame, qtl_list)


setwd("~/Documents/Hotspots/Paper_version_4/Marker_bed/")
chroms<-list()
for(bed in dir(pattern="iSELECT")){
  file<-read.delim(bed, header=F, skip=1)
  chromosome=substr(bed, 12, 13)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file$Platform<-"iSelect"
  chroms[[length(chroms)+1]]<-file
}
iselect<-do.call(rbind.data.frame, chroms)

chroms<-list()
for(bed in dir(pattern="DArT")){
  file<-read.delim(bed, header=F, skip=1)
  chromosome=substr(bed, 14, 15)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file$Platform<-"DArT"
  chroms[[length(chroms)+1]]<-file
}
DarT<-do.call(rbind.data.frame, chroms)


markers<-rbind(iselect, DarT)
write.csv(all_markers_positions, file="~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv")
write.csv(maps, file="~/Documents/Hotspots/Paper_version_4/maps.csv")

