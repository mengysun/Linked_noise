library(dplyr)
library(ggplot2)
Dat<-read.table(file="../Data_needed/Capture_efficiency_100Single",sep="\t",header=TRUE)
png(filename="FigS6a.png",width=10,height=10,units="in",res=600)
ggplot(Dat, aes(x = cor_assigned, y = cor_generated)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab("True correlation")+
  ylab("Measured correlation across single cells")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black",margin=margin(t=30,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",margin=margin(t=0,r=30,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
##################################################################################################
#FigS5b
Dat<-read.table(file="../Data_needed/Capture_efficiency_100",sep="\t",header=TRUE)
png(filename="FigS6b.png",width=10,height=10,units="in",res=600)
ggplot(Dat, aes(x = cor_assigned, y = cor_estimated)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab("True correlation")+
  ylab("Measured correlation across single cells")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black",margin=margin(t=30,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",margin=margin(t=0,r=30,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()

####################################################################################################
#FigS5c
Dat<-read.table(file="../Data_needed/Capture_efficiency_10Single",sep="\t",header=TRUE)
png(filename="FigS6c.png",width=10,height=10,units="in",res=600)
ggplot(Dat, aes(x = cor_assigned, y = cor_estimated)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab("True correlation")+
  ylab("Measured correlation across single cells")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black",margin=margin(t=30,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",margin=margin(t=0,r=30,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
####################################################################################################
Dat<-read.table(file="../Data_needed/Capture_efficiency_10",sep="\t",header=TRUE)
png(filename="FigS6d.png",width=10,height=10,units="in",res=600)
ggplot(Dat, aes(x = cor_assigned, y = cor_estimated)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  xlab("True correlation")+
  ylab("Measured correlation across single cells")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="black",margin=margin(t=30,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="black",margin=margin(t=0,r=30,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
