library(dplyr)
library(ggplot2)
Mouse_Fib_pairs<-read.table(file="../Data_needed/Fib_all_reads10",sep="\t",header=TRUE)
names(Mouse_Fib_pairs)<-c("G_1","G_2","G_1_index","G_2_index","Cor","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads")
#correct for the TSS sites
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes$TSS<-ifelse(Mouse_genes$Strand==1,Mouse_genes$Gene.start..bp.,Mouse_genes$Gene.end..bp.)
TSS_match_G1index<-match(Mouse_Fib_pairs$G_1,Mouse_genes$Gene.name)
TSS_match_G2index<-match(Mouse_Fib_pairs$G_2,Mouse_genes$Gene.name)
Mouse_Fib_pairs$G1_TSS<-Mouse_genes$TSS[TSS_match_G1index]
Mouse_Fib_pairs$G2_TSS<-Mouse_genes$TSS[TSS_match_G2index]
Mouse_Fib_pairs$g_1_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G1index]
Mouse_Fib_pairs$g_2_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G2index]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="X",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="3",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="4",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="X",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="3",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="4",]
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr==Mouse_Fib_pairs$g_2_chr,]
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!=Mouse_Fib_pairs$g_2_chr,]
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs_on_samchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs_on_samchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE&is.na(Cor)!=TRUE)
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs_on_difchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs_on_difchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE&is.na(Cor)!=TRUE)
median(Mouse_Fib_pairs_on_difchr$Cor)
#0.184697 all reads, 0.01043835 rpkm
median(Mouse_Fib_pairs_on_samchr$Cor)
Mouse_Fib_pairs_on_samchr$Cor<-Mouse_Fib_pairs_on_samchr$Cor-0.184697
Mouse_Fib_pairs_on_samchr$distance<-abs(Mouse_Fib_pairs_on_samchr$G1_TSS-Mouse_Fib_pairs_on_samchr$G2_TSS)
Bin<-quantile(as.numeric(as.character(Mouse_Fib_pairs_on_samchr$distance)),(0:100)/100)

Mouse_Fib_pairs_on_samchr_bin_by_distance<-split(Mouse_Fib_pairs_on_samchr,cut(Mouse_Fib_pairs_on_samchr$distance,100))
Cor<-numeric(100)
Distance<-numeric(100)
for(i in 1:100){
  Cor[i]<-median(Mouse_Fib_pairs_on_samchr_bin_by_distance[[i]]$Cor)
  Distance[i]<-median(Mouse_Fib_pairs_on_samchr_bin_by_distance[[i]]$distance)
}
plot(Distance,Cor)
Dat<-data.frame(Distance/1000000,Cor)
names(Dat)<-c("TSS_distance","CI")
png(filename="FigS3a",width=10,height=6.18,units="in",res=600)
ggplot(Dat,aes(x=TSS_distance,y=CI))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  xlab(label="TSS distance (Mb)")+
  #ylab(label=expression(paste(italic(Delta)[e]^{"'"})))+
  ylab(label=expression(paste(italic(Delta)[e])))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",angle=0,vjust = 0.6))+
  geom_hline(yintercept =0,color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
#######################################################################################################
#FigS2b
library(dplyr)
library(ggplot2)
Mouse_Fib_pairs<-read.table(file="../Data_needed/Fib_rpkm10",sep="\t",header=TRUE)
names(Mouse_Fib_pairs)<-c("G_1","G_2","G_1_index","G_2_index","Cor","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads")
#correct for the TSS sites
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes$TSS<-ifelse(Mouse_genes$Strand==1,Mouse_genes$Gene.start..bp.,Mouse_genes$Gene.end..bp.)
TSS_match_G1index<-match(Mouse_Fib_pairs$G_1,Mouse_genes$Gene.name)
TSS_match_G2index<-match(Mouse_Fib_pairs$G_2,Mouse_genes$Gene.name)
Mouse_Fib_pairs$G1_TSS<-Mouse_genes$TSS[TSS_match_G1index]
Mouse_Fib_pairs$G2_TSS<-Mouse_genes$TSS[TSS_match_G2index]
Mouse_Fib_pairs$g_1_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G1index]
Mouse_Fib_pairs$g_2_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G2index]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="X",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="3",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!="4",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="X",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="3",]
Mouse_Fib_pairs<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_2_chr!="4",]
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr==Mouse_Fib_pairs$g_2_chr,]
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs[Mouse_Fib_pairs$g_1_chr!=Mouse_Fib_pairs$g_2_chr,]
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs_on_samchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_Fib_pairs_on_samchr<-Mouse_Fib_pairs_on_samchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE&is.na(Cor)!=TRUE)
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs_on_difchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_Fib_pairs_on_difchr<-Mouse_Fib_pairs_on_difchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE&is.na(Cor)!=TRUE)
median(Mouse_Fib_pairs_on_difchr$Cor)
#0.184697 all reads, 0.01043835 rpkm
median(Mouse_Fib_pairs_on_samchr$Cor)
Mouse_Fib_pairs_on_samchr$Cor<-Mouse_Fib_pairs_on_samchr$Cor-0.01043835
Mouse_Fib_pairs_on_samchr$distance<-abs(Mouse_Fib_pairs_on_samchr$G1_TSS-Mouse_Fib_pairs_on_samchr$G2_TSS)
Bin<-quantile(as.numeric(as.character(Mouse_Fib_pairs_on_samchr$distance)),(0:100)/100)

Mouse_Fib_pairs_on_samchr_bin_by_distance<-split(Mouse_Fib_pairs_on_samchr,cut(Mouse_Fib_pairs_on_samchr$distance,100))
Cor<-numeric(100)
Distance<-numeric(100)
for(i in 1:100){
  Cor[i]<-median(Mouse_Fib_pairs_on_samchr_bin_by_distance[[i]]$Cor)
  Distance[i]<-median(Mouse_Fib_pairs_on_samchr_bin_by_distance[[i]]$distance)
}
plot(Distance,Cor)
Dat<-data.frame(Distance/1000000,Cor)
names(Dat)<-c("TSS_distance","CI")
png(filename="FigS3b",width=10,height=6.18,units="in",res=600)
ggplot(Dat,aes(x=TSS_distance,y=CI))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  xlab(label="TSS distance (Mb)")+
  ylab(label=expression(paste(italic(Delta)[e]^{"'"})))+
  #ylab(label=expression(paste(italic(Delta)[e])))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",angle=0,vjust = 0.6))+
  geom_hline(yintercept =0,color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
