library("MCMCpack")
library(ggplot2)
#install.packages('lmodel2')
library('lmodel2')
cell_num<-50000
population_num<-16
Capturing_efficiency<-0.15*0.55
#biological parameters.
probv<-rdirichlet(10000, c(1,1,1,1) )
r<-probv[,1]
p<-probv[,2]
q<-probv[,3]
#random process parameters

E_x<-p+q
E_y<-p+r
EXY<-p
sdx<-(E_x*(1-E_x))^(1/2)
sdy<-(E_y*(1-E_y))^(1/2)
cor_assigned1<-(EXY-(E_x*E_y))/(sdx*sdy)
cor_estimated1<-numeric(10000)
for(i in 1:10000){
  total_reads_1<-numeric(population_num)
  total_reads_2<-numeric(population_num)
  for(j in 1:population_num){
    cells<-sample(sample(x=c("01","11","10","00"),cell_num,prob=c(r[i],p[i],q[i],1-r[i]-p[i]-q[i]),replace=TRUE))
    total_reads_1[j]<-rbinom(n=1,size=sum(cells=="11")+sum(cells=="10"),prob=Capturing_efficiency)
    total_reads_2[j]<-rbinom(n=1,size=sum(cells=="11")+sum(cells=="01"),prob=Capturing_efficiency)
  }
  cor_estimated1[i]<-cor(total_reads_1,total_reads_2)
}

probv<-rdirichlet(10000, c(1,1,1,1) )
r<-probv[,1]
p<-probv[,2]
q<-probv[,3]
#random process parameters

E_x<-p+q
E_y<-p+r
EXY<-p
sdx<-(E_x*(1-E_x))^(1/2)
sdy<-(E_y*(1-E_y))^(1/2)
cor_assigned2<-(EXY-(E_x*E_y))/(sdx*sdy)
cor_estimated2<-numeric(10000)
for(i in 1:10000){
  total_reads_1<-numeric(population_num)
  total_reads_2<-numeric(population_num)
  for(j in 1:population_num){
    cells<-sample(sample(x=c("01","11","10","00"),cell_num,prob=c(r[i],p[i],q[i],1-r[i]-p[i]-q[i]),replace=TRUE))
    total_reads_1[j]<-rbinom(n=1,size=sum(cells=="11")+sum(cells=="10"),prob=Capturing_efficiency)
    total_reads_2[j]<-rbinom(n=1,size=sum(cells=="11")+sum(cells=="01"),prob=Capturing_efficiency)
  }
  cor_estimated2[i]<-cor(total_reads_1,total_reads_2)
}
True_delta<-(cor_assigned1-cor_assigned2)
estimated_delta<-(cor_estimated1-cor_estimated2)
cor.test(True_delta,estimated_delta)
dat<-data.frame(True_delta,estimated_delta)
png(filename="./FigS6a",width=10,height=10,unit="in",res=600)
ggplot(dat,aes(x=True_delta,y=estimated_delta))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab(label=expression(paste("True ",delta[a])))+
  ylab(label=expression(paste("Estimated ",delta[a])))+
  theme(axis.text.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,family="TT Arial",color="black",angle=90,hjust = 0.5,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
summary(lm(estimated_delta~True_delta))
summary(lm(True_delta~estimated_delta))
#######################################################################################################
#FigS6b
library(MASS)
library(ggplot2)
expression<-readRDS("../Data_needed/mouse.rpkm.3digits.rds") 
C57_expression<-readRDS("../Data_needed/mouse.c57.counts.rds")
Cas_expression<-readRDS("../Data_needed/mouse.cast.counts.rds")
total_reads<-readRDS("../Data_needed/mouse.counts.rds")
mouse_genes<-readRDS("../Data_needed/mouse.gene.annotation.rds")
C57_expression<-C57_expression[,61:120]
Cas_expression<-Cas_expression[,61:120]
total_reads<-total_reads[,61:120]
cl7<-expression[,61:120]
C57_high<-rowMeans(C57_expression)
Cas_high<-rowMeans(Cas_expression)
C57_expression<-C57_expression[C57_high>=10&Cas_high>=10,]
Cas_expression<-Cas_expression[C57_high>=10&Cas_high>=10,]
mouse_genes<-mouse_genes[C57_high>=10&Cas_high>=10,]
total_reads<-total_reads[C57_high>=10&Cas_high>=10,]
cl7<-cl7[C57_high>=10&Cas_high>=10,]
non_aneuploid_index<-(mouse_genes$chrom!="chr3"&mouse_genes$chrom!="chr4"&mouse_genes$chrom!="X")
cl7<-cl7[non_aneuploid_index,]
expression_for_sam<-as.numeric(round(rowMeans(cl7)))

cell_num=60
c_1<-0.5
c_2<-0.5
true_delta_e<-numeric(10000)
estimate_delta_e<-numeric(10000)
for(j in 1:10000){
  e_1<-round(sample(expression_for_sam,size=1)/2)
  e_2<-round(sample(expression_for_sam,size=1)/2)
  cor1<-runif(n=1,min=-1,max=1)
  cor2<-runif(n=1,min=-1,max=1)
  true_delta_e[j]<-cor1-cor2
  #correlation 1
  cor_assigned=cor1
  Sigma <- matrix(c((e_1*c_1)^2,e_2*e_1*c_1*c_2*cor_assigned,e_2*e_1*c_1*c_2*cor_assigned,(e_2*c_2)^2),2,2)
  Gene_pair_correlated<-mvrnorm(n = cell_num,c(e_1,e_2),Sigma,empirical = TRUE)
  Gene_pair_correlated<-round(Gene_pair_correlated)
  g_1_reads<-numeric(cell_num)
  g_2_reads<-numeric(cell_num)
  i<-1
  for(i in 1:cell_num){
    if(Gene_pair_correlated[i,1]<0){
      Gene_pair_correlated[i,1]<-0
    }
    if(Gene_pair_correlated[i,2]<0){
      Gene_pair_correlated[i,2]<-0
    }
    g_1_reads[i]<-rbinom(n=1,size=Gene_pair_correlated[i,1],prob=0.15*0.17)
    g_2_reads[i]<-rbinom(n=1,size=Gene_pair_correlated[i,2],prob=0.15*0.17)
  }
  estimated_cor1<-cor(g_1_reads,g_2_reads,method="pearson")
  #correlation 2
  cor_assigned=cor2
  Sigma <- matrix(c((e_1*c_1)^2,e_2*e_1*c_1*c_2*cor_assigned,e_2*e_1*c_1*c_2*cor_assigned,(e_2*c_2)^2),2,2)
  Gene_pair_correlated<-mvrnorm(n = cell_num,c(e_1,e_2),Sigma,empirical = TRUE)
  Gene_pair_correlated<-round(Gene_pair_correlated)
  g_1_reads<-numeric(cell_num)
  g_2_reads<-numeric(cell_num)
  i<-1
  for(i in 1:cell_num){
    if(Gene_pair_correlated[i,1]<0){
      Gene_pair_correlated[i,1]<-0
    }
    if(Gene_pair_correlated[i,2]<0){
      Gene_pair_correlated[i,2]<-0
    }
    g_1_reads[i]<-rbinom(n=1,size=Gene_pair_correlated[i,1],prob=0.0255)
    g_2_reads[i]<-rbinom(n=1,size=Gene_pair_correlated[i,2],prob=0.0255)
  }
  estimated_cor2<-cor(g_1_reads,g_2_reads,method="pearson")
  #estimate_delta_e
  estimate_delta_e[j]<-estimated_cor1-estimated_cor2
}

dat<-data.frame(true_delta_e,estimate_delta_e)
png(filename="./FigS6b",width=10,height=10,unit="in",res=600)
ggplot(dat,aes(x=true_delta_e,y=estimate_delta_e))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab(label=expression(paste("True ",delta[e])))+
  ylab(label=expression(paste("Estimated ",delta[e])))+
  theme(axis.text.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,family="TT Arial",color="black",angle=90,hjust = 0.5,vjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=2))+
  theme(axis.ticks.length=unit(5,"mm"))
dev.off()
