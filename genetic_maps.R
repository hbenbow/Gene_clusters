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
for(bed in list.files(path="Infinium/")){
  file<-read.delim(paste("Infinium/", bed, sep=""), header=F, skip=1)
  chromosome=substr(bed, 14, 15)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file$Platform<-"Infinium"
  chroms[[length(chroms)+1]]<-file
}
infinium<-do.call(rbind.data.frame, chroms)

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
  file5<-file4
  file5$Feature <- gsub('wms', 'gwm', file5$Feature)
  file6<-rbind(file4, file5)
  chroms[[length(chroms)+1]]<-file6
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

chroms<-list()
for(bed in list.files(path="Blast/", pattern=".blast")){
  file<-read.delim(paste("Blast/", bed, sep=""), header=F, skip=1)
  chromosome=substr(bed, 1, nchar(bed)-6)
  file<-file[,c(2,9,10,1)]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file1<-file
  file1$Feature<-paste("X", tolower(file1$Feature), sep="")
  file2<-rbind(file, file1)
  file2$Platform<-paste(chromosome)
  chroms[[length(chroms)+1]]<-file2
  }
blast<-do.call(rbind.data.frame, chroms)
blast$Chromosome<-gsub("Chr", "chr", blast$Chromosome)

chroms<-list()
for(bed in list.files(path="Manual/")){
  file<-read.delim(paste("Manual/", bed, sep=""), header=F, skip=1)
  file<-file[,1:4]
  colnames(file)<-c("Chromosome", "Start", "End", "Feature")
  file1<-file
  file1$Feature<-paste("X", tolower(file1$Feature), sep="")
  file2<-rbind(file, file1)
  file2$Platform<-paste("Manual")
  chroms[[length(chroms)+1]]<-file2
}
manual<-do.call(rbind.data.frame, chroms)

markers<-rbind(iselect, DarT, SSRs, barc, infinium, blast, manual)
write.csv(markers, file="~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv")

