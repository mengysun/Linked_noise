library(dplyr)
library(ggplot2)
Mouse_cis_trans<-read.table(file="../Data_needed/cl6_gpairs",sep="\t",header=FALSE)
names(Mouse_cis_trans)<-c("G_1","G_2","G_1_index","G_2_index","cis_c57","cis_cas","c57_cas",
                          "cas_c57","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads")
#correct for the TSS sites
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes$TSS<-ifelse(Mouse_genes$Strand==1,Mouse_genes$Gene.start..bp.,Mouse_genes$Gene.end..bp.)
TSS_match_G1index<-match(Mouse_cis_trans$G_1,Mouse_genes$Gene.name)
TSS_match_G2index<-match(Mouse_cis_trans$G_2,Mouse_genes$Gene.name)
Mouse_cis_trans$G1_TSS<-Mouse_genes$TSS[TSS_match_G1index]
Mouse_cis_trans$G2_TSS<-Mouse_genes$TSS[TSS_match_G2index]
Mouse_cis_trans$g_1_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G1index]
Mouse_cis_trans$g_2_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G2index]
table(Mouse_cis_trans$g_1_chr)

Mouse_cis_trans$cis_sum<-Mouse_cis_trans$cis_c57+Mouse_cis_trans$cis_cas
Mouse_cis_trans$trans_sum<-Mouse_cis_trans$c57_cas+Mouse_cis_trans$cas_c57
Mouse_cis_trans$cis_trans_dif<-(Mouse_cis_trans$cis_sum-Mouse_cis_trans$trans_sum)/2

Mouse_cis_trans_on_samchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr==Mouse_cis_trans$g_2_chr,]
Mouse_cis_trans_on_difchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!=Mouse_cis_trans$g_2_chr,]
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE)
Mouse_cis_trans_on_difchr<-Mouse_cis_trans_on_difchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_cis_trans_on_difchr<-Mouse_cis_trans_on_difchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE)
sum(Mouse_cis_trans_on_samchr$cis_trans_dif>0)
220822/383289
sum(Mouse_cis_trans_on_difchr$cis_trans_dif>0)
3114144/6236052
df <- data.frame(gene_pairs=c("Linked", "Unlinked"),
                 ratio=c(220822/383289,3114144/6236052))
df$bar_order <- factor(df$gene_pairs, as.character(df$gene_pairs))
p_value<-c("p=3e-13","p=0.28")
png(filename="FigS2a",width=10,height=6,units="in",res=600)
ggplot(data=df, aes(x=bar_order, y=ratio)) +
  geom_bar(stat="identity",width=0.3)+
  coord_cartesian(ylim=c(0.4,0.6))+
  xlab(label="")+
  ylab(label=expression(atop("Fraction of gene pairs",paste("with positive ",italic(delta)[e]))))+
  theme(axis.text.x = element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y = element_text(size=30,angle=90,vjust=0.6,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
Mouse_cis_trans_on_samchr$distance<-abs(Mouse_cis_trans_on_samchr$G1_TSS-Mouse_cis_trans_on_samchr$G2_TSS)
Mouse_cis_trans_on_samchr_bin_by_distance<-split(Mouse_cis_trans_on_samchr,cut(Mouse_cis_trans_on_samchr$distance,100))
Cor_increase<-numeric(100)
Distance<-numeric(100)
for(i in 1:100){
  Cor_increase[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$cis_trans_dif)
  Distance[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$distance)
}

Dat<-data.frame(Distance/1000000,Cor_increase)
names(Dat)<-c("TSS_distance","Correlation_increase")
Dat<-read.table(file="DatS2b",sep="\t",header=TRUE)
png(filename="FigS2b",width=10,height=6,units="in",res=600)
ggplot(Dat,aes(x=TSS_distance,y=Correlation_increase))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  xlab(label="TSS distance (Mb)")+
  ylab(label=expression(italic(delta)[e]))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",angle=0,vjust = 0.5))+
  geom_hline(aes(yintercept=0), color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
####################################################################################################
library(dplyr)
library(ggplot2)
Mouse_cis_trans<-read.table(file="../Data_needed/M124170_gpairs",sep="\t",header=FALSE)
names(Mouse_cis_trans)<-c("G_1","G_2","G_1_index","G_2_index","cis_c57","cis_cas","c57_cas",
                          "cas_c57","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads")
#correct for the TSS sites
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes$TSS<-ifelse(Mouse_genes$Strand==1,Mouse_genes$Gene.start..bp.,Mouse_genes$Gene.end..bp.)
TSS_match_G1index<-match(Mouse_cis_trans$G_1,Mouse_genes$Gene.name)
TSS_match_G2index<-match(Mouse_cis_trans$G_2,Mouse_genes$Gene.name)
Mouse_cis_trans$G1_TSS<-Mouse_genes$TSS[TSS_match_G1index]
Mouse_cis_trans$G2_TSS<-Mouse_genes$TSS[TSS_match_G2index]
Mouse_cis_trans$g_1_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G1index]
Mouse_cis_trans$g_2_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G2index]
table(Mouse_cis_trans$g_1_chr)
Mouse_cis_trans<-Mouse_cis_trans%>%
  filter(is.na(g_1_chr)!=TRUE&is.na(g_2_chr)!=TRUE)
Mouse_cis_trans_sam_chr<-Mouse_cis_trans%>%
  filter(g_1_chr==g_2_chr)
Mouse_cis_trans_sam_chr$cis_trans_dif<-(Mouse_cis_trans_sam_chr$cis_c57+Mouse_cis_trans_sam_chr$cis_cas-Mouse_cis_trans_sam_chr$c57_cas-Mouse_cis_trans_sam_chr$cas_c57)/2
table(Mouse_cis_trans_sam_chr$g_1_chr)
Mouse_cis_trans_dif_chr<-Mouse_cis_trans%>%
  filter(g_1_chr!=g_2_chr)
Mouse_cis_trans_dif_chr$cis_trans_dif<-(Mouse_cis_trans_dif_chr$cis_c57+Mouse_cis_trans_dif_chr$cis_cas-Mouse_cis_trans_dif_chr$c57_cas-Mouse_cis_trans_dif_chr$cas_c57)/2
sum(Mouse_cis_trans_sam_chr$cis_trans_dif>0)
276266/433911
sum(Mouse_cis_trans_dif_chr$cis_trans_dif>0)
3543033/7052604

df <- data.frame(gene_pairs=c("Linked", "Unlinked"),
                 ratio=c(276266/433911,3543033/7052604))
df$bar_order <- factor(df$gene_pairs, as.character(df$gene_pairs))
p_value<-c("p=3e-13","p=0.28")
png(filename="FigS2C",width=10,height=6,units="in",res=600)
ggplot(data=df, aes(x=bar_order, y=ratio)) +
  geom_bar(stat="identity",width=0.3)+
  coord_cartesian(ylim=c(0.4,0.65))+
  xlab(label="")+
  ylab(label=expression(atop("Fraction of gene pairs",paste("with positive ",italic(delta)[e]))))+
  theme(axis.text.x = element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y = element_text(size=30,angle=90,vjust=0.6,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
Mouse_cis_trans_sam_chr$distance<-abs(Mouse_cis_trans_sam_chr$G1_TSS-Mouse_cis_trans_sam_chr$G2_TSS)
Mouse_cis_trans_on_samchr_bin_by_distance<-split(Mouse_cis_trans_sam_chr,cut(Mouse_cis_trans_sam_chr$distance,100))
Cor_increase<-numeric(100)
Distance<-numeric(100)
for(i in 1:100){
  Cor_increase[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$cis_trans_dif)
  Distance[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$distance)
}

Dat<-data.frame(Distance/1000000,Cor_increase)
names(Dat)<-c("TSS_distance","Correlation_increase")
Dat<-read.table(file="DatS2d",sep="\t",header=TRUE)
png(filename="FigS2d",width=10,height=6,units="in",res=600)
ggplot(Dat,aes(x=TSS_distance,y=Correlation_increase))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  xlab(label="TSS distance (Mb)")+
  ylab(label=expression(italic(delta)[e]))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",angle=0,vjust = 0.5))+
  geom_hline(aes(yintercept=0), color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
