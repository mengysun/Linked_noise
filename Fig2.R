library(dplyr)
library(ggplot2)
library(reshape2)
#Fig2b
link_test<-binom.test(230329,377584,p=0.5)
unlink_test<-binom.test(2707864,5417726,p=0.5)
df <- data.frame(gene_pairs=c("Linked", "Unlinked"),
                 ratio=c(230329/377584,2707864/5417726),
                 se=c((link_test$conf.int[2]-link_test$conf.int[1])/2,(unlink_test$conf.int[2]-unlink_test$conf.int[1])/2))
df$bar_order <- factor(df$gene_pairs, as.character(df$gene_pairs))
p_value<-c(expression(paste(italic(P),"<2.2e-16")),expression(paste(italic(P),"=0.25")))
png(file="./Fig2b",units="in",height=6,width=10,res=600)
ggplot(data=df, aes(x=bar_order, y=ratio)) +
  geom_bar(stat="identity",width=0.3)+
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  coord_cartesian(ylim=c(0.4,0.7))+
  xlab(label="")+
  ylab(label=expression(atop("Fraction of gene pairs ",paste("with positive ",italic(delta)[e]))))+
  theme(axis.text.x = element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.y = element_text(size=35,angle=90,vjust = 0.5,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=35,family="TT Arial",color="black"))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
#####################################################################################################
#Fig2c
library(dplyr)
library(ggplot2)
Mouse_cis_trans<-read.table(file="../Data_needed/Fib_pairs_AllNew",sep="\t",header=TRUE)
names(Mouse_cis_trans)<-c("G_1","G_2","G_1_index","G_2_index","cis_c57","cis_cas","c57_cas",
                          "cas_c57","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads",
                          "g_1_chr","g_2_chr","g_1_start","g_1_end","g_2_start","g_2_end")
#correct for the TSS sites
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
Mouse_cis_trans<-Mouse_cis_trans%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE)

Mouse_cis_trans$cis_sum<-Mouse_cis_trans$cis_c57+Mouse_cis_trans$cis_cas
Mouse_cis_trans$trans_sum<-Mouse_cis_trans$c57_cas+Mouse_cis_trans$cas_c57
Mouse_cis_trans$cis_trans_dif<-(Mouse_cis_trans$cis_sum-Mouse_cis_trans$trans_sum)/2
Mouse_cis_trans_on_samchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr==Mouse_cis_trans$g_2_chr,]

chr<-c(1,2,5:19)
numerator<-numeric(17)
denominator<-numeric(17)
for(i in 1:17){
  Mouse_tmp<-Mouse_cis_trans_on_samchr%>%
    filter(g_1_chr==chr[i])
  numerator[i]<-sum(Mouse_tmp$cis_trans_dif>0)
  denominator[i]<-length(as.character(Mouse_tmp$G_1))
}
Ratio<-numerator/denominator
Ratio[18]<-2707864/5417726
chr_vector<-character(18)
for(i in 1:17){
  #chr_vector[i]<-paste("chr",chr[i],sep="")
  chr_vector[i]<-as.character(chr[i])
}
chr_vector[18]<-"Ctl"
#p by independent
Mouse_cis_trans_independently_samchr<-read.table(file="../Data_needed/Mouse_sam_chr_independent_filtered_tab",sep="\t",header=TRUE)
chr_independent<-c(1,2,5:19)
numerator_independent<-numeric(17)
denominator_independent<-numeric(17)
Ratio_independent<-numeric(18)
for(i in 1:17){
  Mouse_tmp<-Mouse_cis_trans_independently_samchr%>%
    filter(g_1_chr==chr[i])
  numerator_independent[i]<-sum(Mouse_tmp$cis_trans_dif>0)
  denominator_independent[i]<-length(as.character(Mouse_tmp$G_1))
}
Ratio_independent[c(1:17)]<-numerator_independent/denominator_independent
Ratio_independent[18]<-2707864/5417726
binomp<-numeric(17)
se<-numeric(18)
se_non_indenpendent<-numeric(18)
for(i in 1:17){
  binom_tmp_origin<-binom.test(numerator[i],denominator[i],p=0.5)
  binom_tmp<-binom.test(numerator_independent[i],denominator_independent[i],p=0.5)
  binomp[i]<-binom_tmp$p.value
  se[i]<-(binom_tmp$conf.int[2]-binom_tmp$conf.int[1])/2
  se_non_indenpendent[i]<-(binom_tmp_origin$conf.int[2]-binom_tmp_origin$conf.int[1])/2
}
se[18]<-(0.5099443-0.4619023)/2
se_non_indenpendent[18]<-(binom.test(2707864,5417726,p=0.5)$conf.int[2]-binom.test(2707864,5417726,p=0.5)$conf.int[1])/2
chr_independent[binomp<0.05]
#chr1,6,9,10,11,8 are significantly higher than 50%.1:***,6:*,9:****,10:****,11:*,18:**
c(1:17)[binomp<0.00005]
df <- data.frame(chr_vector,Ratio,se,Ratio_independent,se_non_indenpendent)
names(df)<-c("chr_vector","Ratio","se","Ratio_independent","se_non_indenpendent")
df$bar_order <- factor(df$chr_vector, as.character(df$chr_vector))
p_value<-character(18)
p_value[1]<-"***"
p_value[c(4,11)]<-"*"
p_value[c(7,8)]<-"****"
p_value[16]<-"**"
png(file="./Fig2c",width=20,height=6,units="in",res=600)
ggplot(data=df, aes(x=chr_vector, y=Ratio)) +
  geom_bar(aes(x=bar_order),stat="identity",width=0.6,data=df)+
  geom_errorbar(aes(ymin=Ratio-se_non_indenpendent, ymax=Ratio+se_non_indenpendent),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  xlab("Chromosomes")+
  theme(axis.text.x = element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  #geom_text(aes(label=p_value),size=5,vjust=0,family="TT Arial",color="black")+
  ylab(label=expression(atop("Fraction of gene pairs",paste("with positive ",delta[e]))))+
  theme(axis.title.y=element_text(size=25,vjust = 0.5,angle=90,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=25,family="TT Arial",color="black"))+
  coord_cartesian(ylim=c(0.4,0.9))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))

dev.off()

# P value
# 
# 
# Wording
# 
# 
# Summary
# 
# < 0.0001
# 
# 
# Extremely significant
# 
# 
# ****
#   
#   0.0001 to 0.001
# 
# 
# Extremely significant
# 
# 
# ***
#   
#   0.001 to 0.01
# 
# 
# Very significant
# 
# 
# **
#   
#   0.01 to 0.05
# 
# 
# Significant
# 
# 
# *
#   
#   ≥ 0.05
# 
# 
# Not significant
# 
# 
# ns

#without p value
ggplot(data=df, aes(x=chr_vector, y=Ratio)) +
  geom_bar(aes(x=bar_order),stat="identity",width=0.3,data=df)+
  theme(axis.text.x = element_text(size=15,face="bold"))+
  theme(axis.text.y=element_text(size=15,face="bold"))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  ylab(label="Fraction of Cis>Trans")+
  theme(axis.title.y=element_text(size=15,face="bold"))+
  theme(axis.title.x=element_blank())+
  coord_cartesian(ylim=c(0.4,0.9))
###################################################################################################
#Fig2d
Mouse_cis_trans<-read.table(file="../Data_needed/Fib_pairs_AllNew",sep="\t",header=TRUE)
names(Mouse_cis_trans)<-c("G_1","G_2","G_1_index","G_2_index","cis_c57","cis_cas","c57_cas",
                          "cas_c57","c57_g1_reads","c57_g2_reads","cas_g1_reads","cas_g2_reads",
                          "g_1_chr","g_2_chr","g_1_start","g_1_end","g_2_start","g_2_end")
#correct for the TSS sites
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

Mouse_cis_trans$cis_sum<-Mouse_cis_trans$cis_c57+Mouse_cis_trans$cis_cas
Mouse_cis_trans$trans_sum<-Mouse_cis_trans$c57_cas+Mouse_cis_trans$cas_c57
Mouse_cis_trans$cis_trans_dif<-(Mouse_cis_trans$cis_sum-Mouse_cis_trans$trans_sum)/2

Mouse_cis_trans_on_samchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr==Mouse_cis_trans$g_2_chr,]
Mouse_cis_trans_on_difchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!=Mouse_cis_trans$g_2_chr,]
#Draw figure for all distance and correlation increase
#filtering higher number of reads.
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE)

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
png(filename="./Fig2d",width=10,height=6,unit="in",res=600)
ggplot(Dat,aes(x=TSS_distance,y=Correlation_increase))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab(label="TSS distance (Mb)")+
  ylab(label=expression(delta[e]))+
  theme(axis.text.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,family="TT Arial",color="black",angle=0,hjust = 0.5,vjust=0.5))+
  geom_hline(aes(yintercept=0), color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
#####################################################################################################
#Fig2e
reads_cutoff<-c(1:10)*5+5
cis_trans_ratio<-numeric(10)
cis_trans_increase<-numeric(10)
for(i in 1:10){
  Mouse_tmp<-Mouse_cis_trans_on_samchr%>%
    filter(c57_g1_reads>=reads_cutoff[i]&c57_g2_reads>=reads_cutoff[i]&cas_g1_reads>=reads_cutoff[i]&cas_g2_reads>=reads_cutoff[i])
  ctd<-as.numeric(as.character(Mouse_tmp$cis_trans_dif))
  cis_trans_increase[i]<-median(ctd)
  cis_trans_ratio[i]<-sum(ctd>0)/length(ctd)
}
dat<-data.frame(reads_cutoff,cis_trans_ratio)
png(filename="Fig2e",width=10,height=6,units="in",res=600)
pd<-position_dodge(1)
ggplot(dat,aes(x=reads_cutoff,y=cis_trans_ratio))+
  #theme(axis.text.x = element_text(size=8,face="bold"))+
  scale_x_continuous(breaks=dat$reads_cutoff,labels=dat$reads_cutoff)+
  xlab(label="Read number cutoff")+ylab(label=expression(atop("Fraction of gene",paste("pairs with positive ",italic(delta)[e]))))+
  geom_point()+
  geom_line()+
  theme(axis.text.y=element_text(size=34,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=34,family="TT Arial",color="black",hjust=0.5))+
  theme(axis.title.y=element_text(size=34,family="TT Arial",color="black",angle=90,vjust=0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks=c(10,20,30,40,50))+
  theme(axis.text.x = element_text(size=34,family="TT Arial",color="black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
######################################################################################################
#Fig2f
dat_f<-data.frame(reads_cutoff,cis_trans_increase)
Mouse_cis_trans_on_samchr$distance<-abs(Mouse_cis_trans_on_samchr$G1_TSS-Mouse_cis_trans_on_samchr$G2_TSS)
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(is.na(distance)!=TRUE)
Bin<-quantile(as.numeric(as.character(Mouse_cis_trans_on_samchr$distance)),(0:100)/100)
Mouse_cis_trans_on_samchr_closely_linked<-Mouse_cis_trans_on_samchr%>%
  filter(distance<Bin[2])
reads_cutoff<-c(1:10)*5+5
cis_trans_ratio<-numeric(10)
cis_trans_increase<-numeric(10)
distance_median<-numeric(10)
for(i in 1:10){
  Mouse_tmp<-Mouse_cis_trans_on_samchr_closely_linked%>%
    filter(c57_g1_reads>=reads_cutoff[i]&c57_g2_reads>=reads_cutoff[i]&cas_g1_reads>=reads_cutoff[i]&cas_g2_reads>=reads_cutoff[i])
  ctd<-as.numeric(as.character(Mouse_tmp$cis_trans_dif))
  cis_trans_increase[i]<-median(ctd)
  cis_trans_ratio[i]<-sum(ctd>0)/length(ctd)
  distance_median[i]<-median(Mouse_tmp$distance)
}
cor.test(c(1:10),cis_trans_increase,method="spearman")$p.value
cor.test(c(1:10),distance_median,method="spearman")
pd<-position_dodge(1)
dat_g<-data.frame(reads_cutoff,cis_trans_increase)
dat_fg<-data.frame(dat_g$reads_cutoff,dat_g$cis_trans_increase,dat_f$cis_trans_increase)
names(dat_fg)<-c("reads_cutoff","closed_linked","all_distance")
cor.test(dat_fg$reads_cutoff,dat_fg$all_distance,method="spearman")
dat_fg_long<-melt(dat_fg,id="reads_cutoff")
dat_fg_long$group<-ifelse(dat_fg_long$variable=="closed_linked","closely linked gene pairs","all linked gene pairs")
png(filename="Fig2f",units="in",width=10,height=6,res=600)
ggplot(dat_fg,aes(x=reads_cutoff))+
  geom_point(aes(y=closed_linked,colour="closely linked gene pairs"))+
  geom_line(aes(y=closed_linked,colour="closely linked gene pairs"))+
  geom_point(aes(y=all_distance,colour="all linked gene pairs"))+
  geom_line(aes(y=all_distance,colour="all linked gene pairs"))+
  scale_colour_manual("", breaks = c("closely linked gene pairs", "all linked gene pairs"),values = c("closely linked gene pairs"="blue", "all linked gene pairs"="red"))+
  xlab(label="Read number cutoff")+ylab(label=expression(delta[e]))+
  theme(axis.text.x = element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",angle=0,vjust=0.5))+
  #theme(legend.position = c(0.6,0.3),legend.text=element_text(size=20,family="TT Arial",color="black"),legend.background = element_rect(color = "black", 
  #fill = "none", size = 0.5, linetype = "solid"),legend.title=element_blank())+
  theme(legend.position = c(0.55,0.25),legend.text=element_text(size=30,family="TT Arial",color="black"),
        legend.background =element_blank(),legend.title=element_blank(),legend.key = element_rect(size = 5),
        legend.key.size = unit(3, 'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))+
  scale_y_continuous(limits = c(0.010, 0.035))
dev.off()
