library(dplyr)
library(ggplot2)
library(reshape2)
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
Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
Mouse_cis_trans_on_difchr<-Mouse_cis_trans_on_difchr%>%
  filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
median(Mouse_cis_trans_on_difchr$c57_cas)
#allele correlation
C57_genes<-data.frame(readRDS("../Data_needed/mouse.c57.counts.rds"))
Cas_genes<-data.frame(readRDS("../Data_needed/mouse.cast.counts.rds"))
Mouse_genes_all<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_expression<-data.frame(readRDS("../Data_needed/mouse.rpkm.3digits.rds"))
C57_genes<-C57_genes[,61:120]
Cas_genes<-Cas_genes[,61:120]
C57_high<-rowMeans(C57_genes)
Cas_high<-rowMeans(Cas_genes)
focus_gene_index<-((C57_high>=10)&(Cas_high>=10))
C57_focus_genes<-C57_genes[focus_gene_index,]
Cas_focus_genes<-Cas_genes[focus_gene_index,]
Focus_expression<-rowMeans(Mouse_expression[focus_gene_index,])
Focus_genes<-names(Focus_expression)
allele_correlation<-numeric(3985)
for(i in 1:3985){
  allele_correlation[i]<-cor(as.numeric(as.character(C57_focus_genes[i,])),as.numeric(as.character(Cas_focus_genes[i,])))
}
Chr_index<-match(Focus_genes,Mouse_genes_all$Gene.name)
chr<-as.character(Mouse_genes_all$Chromosome.scaffold.name[Chr_index])
median(allele_correlation[chr!="3"&chr!="4"&chr!="X"],na.rm = TRUE)
median(allele_correlation)
median(c(as.numeric(as.character(Mouse_cis_trans_on_difchr$c57_cas)),as.numeric(as.character(Mouse_cis_trans_on_difchr$cas_c57))))
median(Mouse_cis_trans_on_samchr$cis_trans_dif)/(0.3753041-0.1105693)
median(Mouse_cis_trans_on_samchr$cis_trans_dif)
#sample allele
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
true_cor<-numeric(10000)
estimate_cor<-numeric(10000)
for(j in 1:10000){
  e_1<-round(sample(expression_for_sam,size=1)/2)
  #e_2<-round(sample(expression_for_sam,size=1)/2)
  e_2<-e_1
  cor1<-runif(n=1,min=-1,max=1)
  true_cor[j]<-cor1-cor2
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
  estimate_cor[j]<-cor(g_1_reads,g_2_reads,method="pearson")
  
}
summary(lm(estimate_cor~true_cor))
median(estimate_cor/true_cor,na.rm = TRUE)

#density plot
#data
correlation_label<-c(rep("c57-c57",377584),rep("cas-cas",377584),rep("c57-cas",377584),rep("cas-c57",377584))
correlation_type_label<-c(rep("cis",377584*2),rep("trans",377584*2))
correlation<-c(as.numeric(as.character(Mouse_cis_trans_on_samchr$cis_c57)),
               as.numeric(as.character(Mouse_cis_trans_on_samchr$cis_cas)),
               as.numeric(as.character(Mouse_cis_trans_on_samchr$c57_cas)),
               as.numeric(as.character(Mouse_cis_trans_on_samchr$cas_c57)))
dat_for_density<-data.frame(correlation,correlation_label,correlation_type_label)
ggplot(dat_for_density[c(1:377584),],aes(x=correlation,fill=correlation_label)) + 
  geom_density(alpha=.3,fill="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank())+
  ggtitle(label = "c57-c57")
ggplot(dat_for_density[c(377585:(377584*2)),],aes(x=correlation,fill=correlation_label)) + 
  geom_density(alpha=.3,fill="green")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank())+
  ggtitle(label = "cas-cas")
ggplot(dat_for_density[c((377584*2+1):(377584*3)),],aes(x=correlation,fill=correlation_label)) + 
  geom_density(alpha=.3,fill="blue")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank())+
  ggtitle(label = "c57-cas")
ggplot(dat_for_density[c((377584*3+1):(377584*4)),],aes(x=correlation,fill=correlation_label)) + 
  geom_density(alpha=.3,fill="yellow")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_blank())+
  ggtitle(label = "cas-c57")

