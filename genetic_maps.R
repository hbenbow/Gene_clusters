setwd("~/Documents/Hotspots/Paper_version_4/Marker_bed/")

chroms<-list()
for(bed in list.files(path="iSelect/")){
  file<-read.delim(paste("iSelect/", bed, sep=""), header=F, skip=1)
  chromosome=substr(bed, 12, 13)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file$Platform<-"iSelect"
  chroms[[length(chroms)+1]]<-file
}
iselect<-do.call(rbind.data.frame, chroms)

chroms<-list()
for(bed in list.files(path="DArT/")){
  file<-read.delim(paste("DArT/", bed, sep=""), header=F, skip=1)
  chromosome=substr(bed, 14, 15)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file$Platform<-"DArT"
  chroms[[length(chroms)+1]]<-file
}
DarT<-do.call(rbind.data.frame, chroms)


chroms<-list()
for(bed in list.files(path="SSR/gff/", pattern="gff3")){
  file<-read.delim(paste("SSR/gff/", bed, sep=""), header=F, comment.char = "#")
  names<-as.character(file$V9)
  split<-do.call('rbind',strsplit(names, split=";"))[,3]
  split<-do.call('rbind',strsplit(split, split="="))[,2]
  file2<-as.data.frame(cbind(as.character(file$V1), file$V4, file$V5, split))
  colnames(file2)<-c("Chromosome", "Start", "End", "Feature")
  file2$Platform<-"SSR"
  file3<-file2
  file3$Feature<-paste("X", tolower(file3$Feature), sep="")
  file4<-rbind(file2, file3)
  chroms[[length(chroms)+1]]<-file4
}
SSRs<-do.call(rbind.data.frame, chroms)

chroms<-list()
for(bed in list.files(path="SSR/", pattern="bed")){
  file<-read.delim(paste("SSR/", bed, sep=""), header=F, skip=1)
  chromosome=substr(bed, 14, 15)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file1<-file
  file1$Feature<-paste("X", tolower(file1$Feature), sep="")
  file2<-rbind(file, file1)
  file2$Platform<-"SSR"
  chroms[[length(chroms)+1]]<-file2
}
barc<-do.call(rbind.data.frame, chroms)



markers<-rbind(iselect, DarT, SSRs, barc)
write.csv(markers, file="~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv")
write.csv(maps, file="~/Documents/Hotspots/Paper_version_4/maps.csv")

