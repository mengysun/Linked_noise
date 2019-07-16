library(dplyr)
library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
#Figure 1a component
Sigma <- matrix(c(2.56,2.4,2.4,2.56),2,2)
x<-mvrnorm(n = 1000, c(8,8), Sigma, empirical = TRUE)
Dat<-data.frame(x)
png(filename = "./Fig1a_down",units="in",width=10,height=10,res=600)
ggplot(Dat, aes(x = X1, y = X2)) + 
  geom_point() +
  xlab("Gene A mRNA concentration (a.u.)")+
  ylab("Gene B mRNA concentration (a.u.)")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="blue",margin = margin(t = 30, r = 0, b = 0, l = 0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="red",margin = margin(t = 0, r = 30, b = 0, l = 0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = '#FFF6D5'),axis.line = element_line(colour = "black"),panel.border=element_blank())
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border=element_blank())+
  xlim(1,15)+
  ylim(1,15)+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()

#Figure 1b component
Sigma <- matrix(c(2.56,0,0,2.56),2,2)
x<-mvrnorm(n = 1000, c(8,8), Sigma, empirical = TRUE)
Dat<-data.frame(x)
png(filename = "./Fig1b_down",units="in",width=10,height=10,res=600)
ggplot(Dat, aes(x = X1, y = X2)) + 
  geom_point() +
  xlab("Gene A mRNA concentration (a.u.)")+
  ylab("Gene B mRNA concentration (a.u.)")+
  theme(axis.title.x=element_text(size=30,family="TT Arial",color="blue",margin = margin(t = 30, r = 0, b = 0, l = 0)))+
  theme(axis.title.y=element_text(size=30,family="TT Arial",color="red",margin = margin(t = 0, r = 30, b = 0, l = 0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = '#FFF6D5'),axis.line = element_line(colour = "black"),panel.border=element_blank())
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border=element_blank())+
  xlim(1,15)+
  ylim(1,15)+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
