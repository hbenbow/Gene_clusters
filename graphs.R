library(dplyr)
library(ggplot2)
library(extrafont)
library(svglite)

fg_final<-read.csv("Fg_final.csv")
hotspots<-read.csv("final_hotspots.csv")
df<-merge(fg_final, hotspots[,c(5,10, 11)], by="GeneID", all.x=T, fill="Not")
df[is.na(df)]<-"Not"

ggplot(df[!(df$Chrom=="chrUn"),], aes(x=end, y=density)) +
  geom_jitter(size=.4, aes(colour=Density.1), alpha=0.6) + 
  facet_grid(Chromosome~Genome) + xlab("Position (bp)") +
  ylab("Gene Density")  + theme_bw()+
  geom_hline(yintercept=0.7, alpha=0.7, linetype='dashed') +
  theme(text = element_text(size=8, colour="black"), legend.position = "none", axis.text.x = element_text(size=4), 
        panel.background = )+
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(values=c("orangered2", "grey60"))

# write.csv(df, file=paste("Gene_cluster_analysis/", d, "analysed.csv", sep="_"), row.names=F)
ggsave(file="~/Documents/Hotspots/Paper_version_4/graphs/density.svg", width=85, height=100, unit="mm", dpi=300)

df$test<-df$Disease * df$Consecutive
ggplot(df[!(df$Chrom=="chrUn"),], aes(x=end, y=test)) +
  geom_jitter(size=.4, aes(colour=Consecutive.1), alpha=0.6) + 
  facet_grid(Chromosome~Genome) + xlab("Position (bp)") +
  ylab("Gene Consecutiveness")  + theme_bw()+
  geom_hline(yintercept=5, alpha=0.7, linetype='dashed') +
  theme(text = element_text(size=8, colour="black"), legend.position = "none", axis.text.x = element_text(size=4), 
        panel.background = ) +
  scale_color_manual(values=c("orangered2", "grey60"))
ggsave(file="~/Documents/Hotspots/Paper_version_4/graphs/consecutive.svg", width=85, height=100, unit="mm", dpi=300)

