library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
Mouse_cis_trans<-read.table(file="../Data_needed/Fib_pairs_AllNew",sep="\t",header=TRUE)
names(Mouse_cis_trans)<-c("G_1","G_2","G_1_index","G_2_index","cis_c57","cis_cas","c57_cas",
                          "cas_c57","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads",
                          "g_1_chr","g_2_chr","g_1_start","g_1_end","g_2_start","g_2_end")
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes$TSS<-ifelse(Mouse_genes$Strand==1,Mouse_genes$Gene.start..bp.,Mouse_genes$Gene.end..bp.)
TSS_match_G1index<-match(Mouse_cis_trans$G_1,Mouse_genes$Gene.name)
TSS_match_G2index<-match(Mouse_cis_trans$G_2,Mouse_genes$Gene.name)
Mouse_cis_trans$G1_TSS<-Mouse_genes$TSS[TSS_match_G1index]
Mouse_cis_trans$G2_TSS<-Mouse_genes$TSS[TSS_match_G2index]
Mouse_cis_trans$g_1_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G1index]
Mouse_cis_trans$g_2_chr<-Mouse_genes$Chromosome.scaffold.name[TSS_match_G2index]

Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="X",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="3",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="4",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="X",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="3",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="4",]
table(Mouse_cis_trans$g_2_chr)

Mouse_cis_trans$cis_sum<-Mouse_cis_trans$cis_c57+Mouse_cis_trans$cis_cas
Mouse_cis_trans$trans_sum<-Mouse_cis_trans$c57_cas+Mouse_cis_trans$cas_c57
Mouse_cis_trans$cis_trans_dif<-(Mouse_cis_trans$cis_sum-Mouse_cis_trans$trans_sum)/2

Mouse_cis_trans_on_samchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr==Mouse_cis_trans$g_2_chr,]
Mouse_cis_trans_on_samchr$pair_distance<-abs(Mouse_cis_trans_on_samchr$G1_TSS-Mouse_cis_trans_on_samchr$G2_TSS)
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(is.na(pair_distance)==FALSE)
chr_list<-c(1,2,5:19)
plots<-list()
cor_vector<-numeric(17)
sig_vector<-numeric(17)
for(j in 1:17){
  Mouse_sep_chr<-Mouse_cis_trans_on_samchr%>%
    filter(g_1_chr==chr_list[j])
  Mouse_nearest_bin_by_distance<-split(Mouse_sep_chr,cut(Mouse_sep_chr$pair_distance,100))
  Cor_increase<-numeric(100)
  distance<-numeric(100)
  for(i in 1:100){
    Cor_increase[i]<-median(Mouse_nearest_bin_by_distance[[i]]$cis_trans_dif)
    distance[i]<-median(Mouse_nearest_bin_by_distance[[i]]$pair_distance)
  }
  Dat<-data.frame(distance/1000000,Cor_increase)
  names(Dat)<-c("Genetic_distance","Correlation_increase")
  chr_number<-paste("chr",chr_list[j],sep="")
  plots[[j]]<-ggplot(Dat,aes(x=Genetic_distance,y=Correlation_increase))+
    geom_point(size=0.5)+
    geom_smooth(method="lm",se=FALSE)+
    xlab(label="TSS Distance (Mb)")+
    ylab(label=expression(italic(delta)[e]))+
    theme(axis.title.y=element_text(angle=0,vjust = 0.5))+
    ggtitle(chr_number)+
    theme(plot.title=element_text(angle=0,hjust=0.5,size=12,"TT Arial",color="black"))+
    geom_hline(yintercept =0,color="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_text(size=12,family="TT Arial",color="black"))+
    theme(axis.text.y=element_text(size=12,"TT Arial",color="black"))+
    theme(axis.title.x=element_text(size=12,"TT Arial",color="black"))+
    theme(axis.title.y=element_text(size=12,"TT Arial",color="black",angle=0,vjust = 0.6))
  
  cor_computed<-cor.test(distance[is.na(distance)!=TRUE&is.na(Cor_increase)!=TRUE],Cor_increase[is.na(distance)!=TRUE&is.na(Cor_increase)!=TRUE],method="spearman")
  cor_vector[j]<-cor_computed$estimate
  sig_vector[j]<-cor_computed$p.value
}
cor_vector
sig_vector
png(filename="FigS5.png",width=12,height=16,units="in",res=600)
grid.arrange(grobs=plots,ncol=3)
dev.off()

