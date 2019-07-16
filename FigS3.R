library(dplyr)
library(ggplot2)
Mouse_cis_trans_samchr_tab<-read.table(file="../Data_needed/Mouse_cis_trans_same_chr_tab",sep="\t",header=TRUE)
Neiborhood_information<-read.table(file="../Data_needed/Neiborhood_pairs_tab",sep="\t",header=TRUE)
Mouse_cis_trans_samchr_tab$label<-paste(Mouse_cis_trans_samchr_tab$G_1,Mouse_cis_trans_samchr_tab$G_2,sep="_")
Neiborhood_information$forward_label<-paste(Neiborhood_information$g_1,Neiborhood_information$g_2,sep="_")
Neiborhood_information$backward_label<-paste(Neiborhood_information$g_2,Neiborhood_information$g_1,sep="_")
forward_match<-match(Mouse_cis_trans_samchr_tab$label,Neiborhood_information$forward_label)
backward_match<-match(Mouse_cis_trans_samchr_tab$label,Neiborhood_information$backward_label)


Mouse_cis_trans_samchr_tab$orientation<-Neiborhood_information$orientation[forward_match]
Mouse_cis_trans_samchr_tab$g_1_start<-Neiborhood_information$g_1_start[forward_match]
Mouse_cis_trans_samchr_tab$g_2_start<-Neiborhood_information$g_2_start[forward_match]
Mouse_cis_trans_samchr_tab$g_1_end<-Neiborhood_information$g_1_end[forward_match]
Mouse_cis_trans_samchr_tab$g_2_end<-Neiborhood_information$g_2_end[forward_match]
Mouse_cis_trans_samchr_tab_filtered<-Mouse_cis_trans_samchr_tab%>%
  filter(is.na(g_1_start)!=TRUE)
Mouse_cis_trans_samchr_tab_filtered$distance<-abs(Mouse_cis_trans_samchr_tab_filtered$g_1_start-Mouse_cis_trans_samchr_tab_filtered$g_2_start)

Convergent_tab<-Mouse_cis_trans_samchr_tab_filtered%>%
  filter(orientation=="convergent")
Divergent_tab<-Mouse_cis_trans_samchr_tab_filtered%>%
  filter(orientation=="divergent")
same_tab<-Mouse_cis_trans_samchr_tab_filtered%>%
  filter(orientation=="same")
median(Convergent_tab$cis_trans_dif)
median(Divergent_tab$cis_trans_dif)
median(same_tab$cis_trans_dif)
wilcox.test(Divergent_tab$distance,Convergent_tab$distance)$p.value
wilcox.test(Divergent_tab$cis_trans_dif,same_tab$cis_trans_dif)
dat<-Mouse_cis_trans_samchr_tab_filtered[,c(23,21,24)]
dat$orientation2<-ifelse(dat$orientation=="convergent","Convergent (n=108)",ifelse(dat$orientation=="divergent","Divergent (n=134)","Same (n=180)"))
png(filename="FigS3",width=10,height=6,units="in",res=600)
ggplot(dat, aes(x = orientation2, y = cis_trans_dif)) +
  geom_boxplot()+
  ylab(label=expression(italic(delta)[e]))+
  xlab(label="Transcription orientations")+
  theme(axis.text.x = element_text(size=18,family="TT Arial",color="black"))+
  theme(axis.text.y = element_text(size=18,family="TT Arial",color="black"))+
  theme(axis.title.x = element_text(size=18,family="TT Arial",color="black",angle=0,vjust=0.5))+
  theme(axis.title.y = element_text(size=18,family="TT Arial",color="black",angle=0,vjust=0.5))+
  
  geom_segment(aes(x="Convergent (n=108)",xend="Divergent (n=134)",y=0.3,yend=0.3),size=0.5)+
  geom_segment(aes(x="Divergent (n=134)",xend="Same (n=180)",y=0.4,yend=0.4),size=0.5)+
  geom_segment(aes(x="Convergent (n=108)",xend="Same (n=180)",y=0.5,yend=0.5),size=0.5)+
  annotate("text",x=1.5,y=0.35,label="NS",size=5)+
  annotate("text",x=2.5,y=0.45,label="NS",size=5)+
  annotate("text",x=2,y=0.55,label="NS",size=5)+
  geom_hline(aes(yintercept=0,col="red"),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
