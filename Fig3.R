library(dplyr)
library(ggplot2)
library(reshape2)
Hic_tab<-read.table(file="../Data_needed/Hic_tab_all",sep="\t",header=TRUE)
Hic_tab_samechr<-Hic_tab%>%
  filter(chr_1==chr_2)
sum(Hic_tab_samechr$Contact_dif<0)
#sam chr: 552261/552736,354/552736,121/552736
Hic_tab_difchr<-Hic_tab%>%
  filter(chr_1!=chr_2)
sum(Hic_tab_difchr$Contact_dif<0)
#dif chr:2112869/9399755,5185752/9399755,2101134/9399755


gene_pairs<-c("Linked","Unlinked")
larger_than_zero<-c(552261/552736,2112869/9399755)
equal_to_zero<-c(354/552736,5185752/9399755)
smaller_than_zero<-c(121/552736,2101134/9399755)
df<-data.frame(gene_pairs,larger_than_zero,equal_to_zero,smaller_than_zero)
names(df)<-c("Gene_pairs","cis-trans>0","cis-trans=0","cis-trans<0")
data.m <- melt(df, id.vars='Gene_pairs')
png(filename="./Fig3b",width=10,height=6,units="in",res=600)
ggplot(data.m, aes(Gene_pairs, value)) + 
  geom_bar(aes(fill = variable), width = 0.5, position = position_dodge(width=0.5), stat="identity") +  
  scale_fill_manual(values=c("#52BE80","#3498DB","#E74C3C"),
                    breaks=c("cis-trans>0", "cis-trans=0", "cis-trans<0"),
                    labels=c(expression(paste(italic(delta)[i]," > 0")), expression(paste(italic(delta)[i]," = 0")), expression(paste(italic(delta)[i]," < 0"))))+
  ylab(label=expression(atop("Fraction of genomic region",paste("pairs with different ",italic(delta)[i]))))+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",vjust=0.6,angle=90))+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=25,family="TT Arial",color="black"))+
  theme(legend.position=c(0.8,0.8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1.1))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
#############################################################################
#Fig3c
library(dplyr)
library(ggplot2)
library(MASS)
library(viridis)
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
Hic_tab_samechr$distance<-abs(Hic_tab_samechr$pos_1-Hic_tab_samechr$pos_2)/1000000
Hic_tab_samechr$Contact_dif_transformed<-Hic_tab_samechr$Contact_dif+5
Hic_tab_samechr$logTransform<-log10(Hic_tab_samechr$Contact_dif_transformed)
png(filename="./Fig3c",units="in",width=10,height=6,res=600)
Hic_tab_samechr$density<-get_density(Hic_tab_samechr$distance,Hic_tab_samechr$Contact_dif_transformed,500)
ggplot(Hic_tab_samechr,aes(x=distance,y=logTransform))+
  geom_point(alpha=0.05)+
  geom_smooth(se=FALSE)+
  xlab(label="Distance (Mb)")+
  ylab(label=expression(paste("log"[10],"(",italic(delta)[i],"+5)")))+
  geom_hline(yintercept = log10(5),color="red")+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=90,vjust=0.5))+
  #scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = c(50, 100, 150))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
#scale_color_viridis()
#stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE) +       
#scale_fill_viridis() +
#coord_cartesian(expand = FALSE) 
dev.off()
##################################################################
#Fig3d
library(dplyr)
library(ggplot2)
ATAC_50_tab<-read.table(file="../Data_needed/ATAC50_gpairs",sep="\t")
names(ATAC_50_tab)<-c("G_1_chr","G_2_chr","G_1_start","G_2_start","G_1_end","G_2_end","G_1_index",
                      "G_2_index","cis_M129","cis_cas","M129_cas","cas_M129","M129_g1_reads",
                      "M129_g2_reads","cas_g1_reads","cas_g2_reads")
ATAC_50_tab_sam_chr<-ATAC_50_tab%>%
  filter(G_1_chr==G_2_chr)
ATAC_50_tab_sam_chr$cis_trans_dif<-((ATAC_50_tab_sam_chr$cis_M129+ATAC_50_tab_sam_chr$cis_cas)-(ATAC_50_tab_sam_chr$M129_cas+ATAC_50_tab_sam_chr$cas_M129))/2
sum(ATAC_50_tab_sam_chr$cis_trans_dif>0)
260015/451684
ATAC_50_tab_dif_chr<-ATAC_50_tab%>%
  filter(G_1_chr!=G_2_chr)
ATAC_50_tab_dif_chr$cis_trans_dif<-((ATAC_50_tab_dif_chr$cis_M129+ATAC_50_tab_dif_chr$cis_cas)-(ATAC_50_tab_dif_chr$M129_cas+ATAC_50_tab_dif_chr$cas_M129))/2
sum(ATAC_50_tab_dif_chr$cis_trans_dif>0)
3047960/5979807
link_test<-binom.test(260015,451684,p=0.5)
unlink_test<-binom.test(3047960,5979807,p=0.5)
df <- data.frame(gene_pairs=c("Linked", "Unlinked"),
                 ratio=c(260015/451684,3047960/5979807),
                 se=c((link_test$conf.int[2]-link_test$conf.int[1])/2,(unlink_test$conf.int[2]-unlink_test$conf.int[1])/2))
df$bar_order <- factor(df$gene_pairs, as.character(df$gene_pairs))
png(filename="Fig3d",units="in",width=10,height=6,res=600)
ggplot(data=df, aes(x=gene_pairs, y=ratio)) +
  geom_bar(aes(x=bar_order),stat="identity",width=0.3)+
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  theme(axis.text.x = element_text(size=25,family="TT Arial",color="black"))+
  coord_cartesian(ylim=c(0.4,0.6))+
  geom_hline(aes(yintercept=0.5,col="red"),show.legend = FALSE)+
  ylab(label=expression(atop("Fraction of ATAC peak",paste("pairs with positive ",italic(delta[a])))))+xlab(label="")+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=90,vjust=0.6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
##########################################################################################
#Fig3e
library(dplyr)
library(ggplot2)
ATAC_50_tab<-read.table(file="../Data_needed/ATAC50_gpairs",sep="\t")
names(ATAC_50_tab)<-c("G_1_chr","G_2_chr","G_1_start","G_2_start","G_1_end","G_2_end","G_1_index",
                      "G_2_index","cis_M129","cis_cas","M129_cas","cas_M129","M129_g1_reads",
                      "M129_g2_reads","cas_g1_reads","cas_g2_reads")
ATAC_50_tab_on_samchr<-ATAC_50_tab%>%
  filter(G_1_chr==G_2_chr)
ATAC_50_tab_on_samchr$g_1_pos<-(ATAC_50_tab_on_samchr$G_1_start+ATAC_50_tab_on_samchr$G_1_end)/2
ATAC_50_tab_on_samchr$g_2_pos<-(ATAC_50_tab_on_samchr$G_2_start+ATAC_50_tab_on_samchr$G_2_end)/2
ATAC_50_tab_on_samchr$distance<-abs(ATAC_50_tab_on_samchr$g_1_pos-ATAC_50_tab_on_samchr$g_2_pos)
ATAC_50_tab_on_samchr$cis_trans_dif<-((ATAC_50_tab_on_samchr$cis_M129+ATAC_50_tab_on_samchr$cis_cas)-(ATAC_50_tab_on_samchr$M129_cas+ATAC_50_tab_on_samchr$cas_M129))/2
ATAC_50_samchr_bin_by_distance<-split(ATAC_50_tab_on_samchr,cut(ATAC_50_tab_on_samchr$distance,100))


Cor_increase<-numeric(100)
Distance<-numeric(100)
dat_points<-numeric(100)
for(i in 1:100){
  Cor_increase[i]<-median(ATAC_50_samchr_bin_by_distance[[i]]$cis_trans_dif)
  Distance[i]<-median(ATAC_50_samchr_bin_by_distance[[i]]$distance)
  dat_points[i]<-length(ATAC_50_samchr_bin_by_distance[[i]]$distance)
}

Dat<-data.frame(Distance[c(1:95,97:100)]/1000000,Cor_increase[c(1:95,97:100)])
Distance[96]
Cor_increase[96]
names(Dat)<-c("peaks_distance","Correlation_increase")
png(filename="Fig3e",width=10,height=6,units="in",res=600)
ggplot(Dat,aes(x=peaks_distance,y=Correlation_increase))+
  geom_point()+
  geom_smooth(method="lm")+xlab(label="Peak distance (Mb)")+ylab(label=expression(italic(delta)[a]))+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=0,vjust=0.6))+
  geom_hline(yintercept = 0,color="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()

##############################################################################################
#Fig3f
library(dplyr)
library(ggplot2)
cis_trans_acceCor_vs_Hic<-read.table(file="../Data_needed/cis_trans_CorTabAllWithP_na_rm_for_sum",sep="\t",header=TRUE)
cis_trans_acceCor_vs_Hic$trans_Hic_sum<-cis_trans_acceCor_vs_Hic$trans_129cas_Hic+cis_trans_acceCor_vs_Hic$trans_cas129_Hic
cis_trans_acceCor_vs_Hic$trans_cor_sum<-cis_trans_acceCor_vs_Hic$V7+cis_trans_acceCor_vs_Hic$V8


cis_trans_acceCor_vs_Hic<-cis_trans_acceCor_vs_Hic%>%
  filter((trans_cor_sum)!="NaN"&(trans_Hic_sum)!="Nan")%>%
  filter(V9>=0&V10>=0&V11>=0&V12>=0)
Hic_trans_sum<-as.numeric(as.character(cis_trans_acceCor_vs_Hic$trans_Hic_sum))
ATAC_trans_corsum<-as.numeric(as.character(cis_trans_acceCor_vs_Hic$trans_cor_sum))
Hic_ATAC_trans<-data.frame(Hic_trans_sum/2,ATAC_trans_corsum/2)

names(Hic_ATAC_trans)<-c("Hic","ATAC")
Hic_ATAC_trans_no_contact<-Hic_ATAC_trans%>%
  filter(Hic==0)
Hic_ATAC_trans_contacted<-Hic_ATAC_trans%>%
  filter(Hic>0)
#just two bins, contact versus no contacted
contacted<-rep("Contacted",358102)
uncontacted<-rep("Uncontacted",916304)
Hic_state<-c(contacted,uncontacted)
co_ATAC<-c(as.numeric(as.character(Hic_ATAC_trans_contacted$ATAC)),as.numeric(as.character(Hic_ATAC_trans_no_contact$ATAC)))
contacted_ATAC<-as.numeric(as.character(Hic_ATAC_trans_contacted$ATAC))
uncontacted_ATAC<-as.numeric(as.character(Hic_ATAC_trans_no_contact$ATAC))
dat<-data.frame(Hic_state,co_ATAC)
png(filename="Fig3f",width=10,height=6,units="in",res=600)
ggplot(dat,aes(factor(Hic_state),co_ATAC))+
  geom_boxplot()+
  xlab(label="")+
  ylab(label="Co-accessibility")+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=90,vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(limits=c(-1,1.3))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
##############################################################################################
#Fig3g
library(dplyr)
Mouse_cis_trans<-read.table(file="../Data_needed/Mouse_cis_trans_for_cor_withATAC_filtered",sep="\t",header=TRUE)
Mouse_ATAC_cor<-read.table(file="../Data_needed/ATAC_cor_gpairs",sep="\t")
names(Mouse_ATAC_cor)<-c("g_1_chr","g_2_chr","g_1_start","g_2_start","g_1_end","g_2_end","g_1_index","g_2_index",
                         "cis_M129","cis_cas","M129_cas","cas_M129","M129_reads_g1","M129_reads_g2","cas_reads_g1",
                         "cas_reads_g2")
Mouse_ATAC_cor$trans_sum<-Mouse_ATAC_cor$M129_cas+Mouse_ATAC_cor$cas_M129
Mouse_ATAC_cor$index_together<-paste(Mouse_ATAC_cor$g_1_index,Mouse_ATAC_cor$g_2_index,sep="_")
Mouse_cis_trans$index_together<-paste(Mouse_cis_trans$g_1_posindex,Mouse_cis_trans$g_2_posindex,sep="_")
Mouse_cis_trans$ATAC_trans_sum<-Mouse_ATAC_cor$trans_sum
Mouse_cis_trans$A_129_g1Reads<-Mouse_ATAC_cor$M129_reads_g1
Mouse_cis_trans$A_129_g2Reads<-Mouse_ATAC_cor$M129_reads_g2
Mouse_cis_trans$A_cas_g1Reads<-Mouse_ATAC_cor$cas_reads_g1
Mouse_cis_trans$A_cas_g2Reads<-Mouse_ATAC_cor$cas_reads_g2
Mouse_cis_trans_group_by_index<-Mouse_cis_trans%>%
  filter(A_129_g1Reads>=10&A_129_g2Reads>=10&A_cas_g1Reads>=10&A_cas_g2Reads>=10)%>%
  group_by(index_together)%>%
  dplyr::summarise(pos_1=mean(g_1_posindex),pos_2=mean(g_2_posindex),cor=mean(trans_sum),cor_ATAC=mean(ATAC_trans_sum))
Bins<-quantile(as.numeric(as.character(Mouse_cis_trans_group_by_index$cor_ATAC)),(0:100)/100)
Mouse_group<-split(Mouse_cis_trans_group_by_index,cut(Mouse_cis_trans_group_by_index$cor_ATAC,Bins))
CG<-numeric(100)
CA<-numeric(100)
for(i in 1:100){
  CG[i]<-median(Mouse_group[[i]]$cor)
  CA[i]<-median(Mouse_group[[i]]$cor_ATAC)
}
Dat<-data.frame(CA/2,CG/2)
names(Dat)<-c("co_accessibility","co_fluctuation")
png(filename="Fig3g",width=10,height=6,units="in",res=600)
ggplot(Dat,aes(x=co_accessibility,y=co_fluctuation))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  xlab(label="Co-accessibility")+
  ylab(label="Co-fluctuation")+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=90,hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))+
  scale_y_continuous(breaks=c(0.10,0.11,0.12))
dev.off()
##############################################################################################
#Fig3H
library(dplyr)
library(ppcor)
library(ggplot2)
Mouse_cor_with_Hic_and_ATAC<-read.table(file="../Data_needed/Mouse_cis_trans_with_Hic_and_ATAC_all",sep="\t",header=TRUE)
Mouse_cor_with_Hic_and_ATAC<-Mouse_cor_with_Hic_and_ATAC%>%
  filter(g_1_chr==g_2_chr)

Mouse_cor_with_Hic_and_ATAC$delta_e<-(Mouse_cor_with_Hic_and_ATAC$Hic_all-2*Mouse_cor_with_Hic_and_ATAC$Hic_trans)/2
Bins<-quantile(as.numeric(as.character(Mouse_cor_with_Hic_and_ATAC$delta_e)),(0:100)/100,na.rm = TRUE)
Mouse_cis_trans_bin_by_Hic<-split(Mouse_cor_with_Hic_and_ATAC,cut(Mouse_cor_with_Hic_and_ATAC$delta_e,unique(Bins)))
length(unique(Bins))

cofluc<-numeric(98)
Hic<-numeric(98)
for(i in 1:98){
  cofluc[i]<-median(Mouse_cis_trans_bin_by_Hic[[i]]$cis_trans_dif)
  Hic[i]<-median(Mouse_cis_trans_bin_by_Hic[[i]]$delta_e)
}
plot(log10(Hic),log10(cofluc))
Dat<-data.frame(log10(Hic),cofluc)
names(Dat)<-c("Hic","co_fluctuation")
summary(lm(cofluc~Hic))
png(filename="Fig3h",width=10,height=6,units="in",res=600)
ggplot(Dat,aes(x=Hic,y=co_fluctuation))+
  geom_point()+
  geom_smooth()+
  xlab(label=expression(paste("log"[10],"(",italic(delta)[i],")")))+
  ylab(label=expression(italic(delta)[e]))+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=25,family="TT Arial",color="black",angle=90,hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
#scale_y_continuous(breaks=c(0.10,0.11,0.12))
dev.off()
