library(dplyr)
library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
set.seed(7)
#install.packages("gridExtra")
alleles<-c(1,1,1,1)
cv<-0.2
s_c<-0.35
K<-1000000
cor_increase<-0.12
sigma<-matrix(0.1*0.04,nrow=4,ncol=4)
cell_number=10000
for(i in 1:4){
  sigma[i,i]<-0.2^2
}
sigma[1,2]<-0.04*s_c
sigma[2,1]<-0.04*s_c
sigma[3,4]<-0.04*s_c
sigma[4,3]<-0.04*s_c
sigma_1<-sigma
sigma_1[1,3]<-0.04*cor_increase
sigma_1[3,1]<-0.04*cor_increase
sigma_1[2,4]<-0.04*cor_increase
sigma_1[4,2]<-0.04*cor_increase

effective_concentration<-function(E_correlated){
  EA1<-E_correlated[,1]
  EA2<-E_correlated[,2]
  EB1<-E_correlated[,3]
  EB2<-E_correlated[,4]
  EA1[EA1<0]<-0
  EA2[EA2<0]<-0
  EB1[EB1<0]<-0
  EB2[EB2<0]<-0
  E_1<-EA1+EA2
  E_2<-EB1+EB2
  concentration_vector<-numeric(cell_number)
  for(j in 1:length(concentration_vector)){
    f1c<-c(K,K*(E_2[j]-E_1[j])+1,-E_1[j])
    f2c<-c(K,K*(E_1[j]-E_2[j])+1,-E_2[j])
    E1f<-(-f1c[2]+sqrt(f1c[2]^2-4*f1c[1]*f1c[3]))/(2*f1c[1])
    E2f<-(-f2c[2]+sqrt(f2c[2]^2-4*f2c[1]*f2c[3]))/(2*f2c[1])
    concentration_vector[j]<-K*E1f*E2f
  }
  return(concentration_vector)
}
unlinked<-numeric(1000)
linked<-numeric(1000)
for(i in 1:1000){
  E_pair_correlated<-mvrnorm(n = cell_number, alleles, sigma, empirical = TRUE)
  E_pair_correlated1<-mvrnorm(n = cell_number, alleles, sigma_1, empirical = TRUE)
  unlinked[i]<-mean(effective_concentration(E_pair_correlated))
  linked[i]<-mean(effective_concentration(E_pair_correlated1))
}
dat_simulated<-data.frame((linked-unlinked)/unlinked*100)
names(dat_simulated)<-c("relative_concentration_increase")
png(file="FigS8c.png",width=10,height=10,units="in",res=600)
ggplot(dat_simulated, aes(x=relative_concentration_increase)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)))+
  ylab(label="Frequency")+
  xlab(label=expression(atop("Percentage increase in concentration", "(linked vs unlinked)")))+
  theme(axis.text.y=element_text(size=30,color="black"))+
  theme(axis.title.y=element_text(size=30,angle=90,vjust=0.5,color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,color="black"))+
  theme(axis.title.x=element_text(size=30,color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
dev.off()
