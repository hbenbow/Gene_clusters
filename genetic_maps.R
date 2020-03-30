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

