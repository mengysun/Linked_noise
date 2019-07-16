#Fig2 tests
#Test for Fig2b and Fig2c
library(dplyr)
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
Mouse_cis_trans_on_difchr<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!=Mouse_cis_trans$g_2_chr,]
sum(Mouse_cis_trans_on_samchr$cis_trans_dif>0)
sum(Mouse_cis_trans_on_difchr$cis_trans_dif>0)
#select the independent pairs:
length(unique(c(as.character(Mouse_cis_trans$G_1),as.character(Mouse_cis_trans$G_2))))
3405*3404/2

3405/2
set.seed(7)
#sam chr independent
Mouse_sam_chr_independent_filtered<-Mouse_cis_trans_on_samchr
Mouse_sam_chr_independent_filtered<-Mouse_cis_trans_on_samchr[sample(c(1:377584)),]
sam_chr_pairsnum<-377584
k<-1
while((k<sam_chr_pairsnum)){
  G1<-as.character(Mouse_sam_chr_independent_filtered$G_1[k])
  G2<-as.character(Mouse_sam_chr_independent_filtered$G_2[k])
  subset_vector<-c(rep(TRUE,k),(Mouse_sam_chr_independent_filtered$G_1[-c(1:k)]!=G1)&(Mouse_sam_chr_independent_filtered$G_1[-c(1:k)]!=G2)&
                     (Mouse_sam_chr_independent_filtered$G_2[-c(1:k)]!=G1)&(Mouse_sam_chr_independent_filtered$G_2[-c(1:k)]!=G2))
  Mouse_sam_chr_independent_filtered<-Mouse_sam_chr_independent_filtered[subset_vector,]
  sam_chr_pairsnum<-length(as.character(Mouse_sam_chr_independent_filtered$G_1))
  k<-k+1
}
sum(Mouse_sam_chr_independent_filtered$cis_trans_dif>0)
1018/1698
length(unique(c(as.character(Mouse_sam_chr_independent_filtered$G_1),as.character(Mouse_sam_chr_independent_filtered$G_2))))
#done for sam chr pairs
3396/2
x<-binom.test(1018,1698,p=0.5)
x$p.value

#dif chr independent
set.seed(7)
Mouse_dif_chr_independent_filtered<-Mouse_cis_trans_on_difchr
Mouse_dif_chr_independent_filtered<-Mouse_dif_chr_independent_filtered[sample(c(1:5417726)),]
dif_chr_pairsnum<-5417726
k<-1
while((k<dif_chr_pairsnum)){
  G1<-as.character(Mouse_dif_chr_independent_filtered$G_1[k])
  G2<-as.character(Mouse_dif_chr_independent_filtered$G_2[k])
  subset_vector<-c(rep(TRUE,k),(Mouse_dif_chr_independent_filtered$G_1[-c(1:k)]!=G1)&(Mouse_dif_chr_independent_filtered$G_1[-c(1:k)]!=G2)&
                     (Mouse_dif_chr_independent_filtered$G_2[-c(1:k)]!=G1)&(Mouse_dif_chr_independent_filtered$G_2[-c(1:k)]!=G2))
  Mouse_dif_chr_independent_filtered<-Mouse_dif_chr_independent_filtered[subset_vector,]
  dif_chr_pairsnum<-length(as.character(Mouse_dif_chr_independent_filtered$G_1))
  k<-k+1
}
sum(Mouse_dif_chr_independent_filtered$cis_trans_dif>0)
length(unique(c(as.character(Mouse_dif_chr_independent_filtered$G_1),as.character(Mouse_dif_chr_independent_filtered$G_2))))
y<-binom.test(827,1702,p=0.5)
y$p.value
write.table(Mouse_sam_chr_independent_filtered,file="Mouse_sam_chr_independent_filtered_tab",row.names=FALSE,col.names=TRUE,sep="\t")
write.table(Mouse_dif_chr_independent_filtered,file="Mouse_dif_chr_independent_filtered_tab",row.names=FALSE,col.names=TRUE,sep="\t")

##########################################################################################################
#Fig2d
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

Mouse_cis_trans<-Mouse_cis_trans%>%
  filter(is.na(G1_TSS)!=TRUE&is.na(G2_TSS)!=TRUE)
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="X",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="3",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr!="4",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="X",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="3",]
Mouse_cis_trans<-Mouse_cis_trans[Mouse_cis_trans$g_2_chr!="4",]

#all_gene_Shuffling
all_genes_S<-unique(c(as.character(Mouse_cis_trans$G_1),as.character(Mouse_cis_trans$G_2)))
gene_index<-match(all_genes_S,Mouse_genes$Gene.name)
chr_all<-Mouse_genes$Chromosome.scaffold.name[gene_index]
TSS_all<-Mouse_genes$TSS[gene_index]
index_G1_back<-match(Mouse_cis_trans$G_1,all_genes_S)
index_G2_back<-match(Mouse_cis_trans$G_2,all_genes_S)
cor_shuffled<-numeric(999)
for(j in 1:999){
  index_shuffled<-sample(c(1:3405))
  chr_all_shuffled<-chr_all[index_shuffled]
  TSS_all_shuffled<-TSS_all[index_shuffled]
  
  Mouse_cis_trans$shuffled_G1_chr<-chr_all_shuffled[index_G1_back]
  Mouse_cis_trans$shuffled_G2_chr<-chr_all_shuffled[index_G2_back]
  Mouse_cis_trans$shuffled_G1_TSS<-TSS_all_shuffled[index_G1_back]
  Mouse_cis_trans$shuffled_G2_TSS<-TSS_all_shuffled[index_G2_back]
  
  #correlation after shuffling
  Mouse_cis_trans_on_samchr<-Mouse_cis_trans[Mouse_cis_trans$shuffled_G1_chr==Mouse_cis_trans$shuffled_G2_chr,]
  Mouse_cis_trans_on_samchr<-Mouse_cis_trans_on_samchr%>%
    filter(c57_g1_reads>=10&c57_g2_reads>=10&cas_g1_reads>=10&cas_g2_reads>=10)
  Mouse_cis_trans_on_samchr$cis_trans_dif<-((Mouse_cis_trans_on_samchr$cis_c57+Mouse_cis_trans_on_samchr$cis_cas)-(Mouse_cis_trans_on_samchr$c57_cas+Mouse_cis_trans_on_samchr$cas_c57))/2
  Mouse_cis_trans_on_samchr$distance<-abs(Mouse_cis_trans_on_samchr$shuffled_G1_TSS-Mouse_cis_trans_on_samchr$shuffled_G2_TSS)
  #Bin<-quantile(as.numeric(as.character(Mouse_cis_trans_on_samchr$distance)),(0:100)/100)
  Mouse_cis_trans_on_samchr_bin_by_distance<-split(Mouse_cis_trans_on_samchr,cut(Mouse_cis_trans_on_samchr$distance,100))
  Cor_increase<-numeric(100)
  Distance<-numeric(100)
  for(i in 1:100){
    Cor_increase[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$cis_trans_dif)
    Distance[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance[[i]]$distance)
  }
  cor_shuffled[j]<-cor(Distance,Cor_increase,method="spearman")
}
hist(cor_shuffled)
sum(cor_shuffled<=(-0.60))
mean(cor_shuffled)
Mouse_cis_trans_on_samchr_observed<-Mouse_cis_trans[Mouse_cis_trans$g_1_chr==Mouse_cis_trans$g_2_chr,]
Mouse_cis_trans_on_samchr_observed$cis_trans_dif<-((Mouse_cis_trans_on_samchr_observed$cis_c57+Mouse_cis_trans_on_samchr_observed$cis_cas)-(Mouse_cis_trans_on_samchr_observed$c57_cas+Mouse_cis_trans_on_samchr_observed$cas_c57))/2
Mouse_cis_trans_on_samchr_observed$distance<-abs(Mouse_cis_trans_on_samchr_observed$G1_TSS-Mouse_cis_trans_on_samchr_observed$G2_TSS)
Mouse_cis_trans_on_samchr_bin_by_distance_observed<-split(Mouse_cis_trans_on_samchr_observed,cut(Mouse_cis_trans_on_samchr_observed$distance,100))
O_Cor_increase<-numeric(100)
ODistance<-numeric(100)
for(i in 1:100){
  O_Cor_increase[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance_observed[[i]]$cis_trans_dif)
  ODistance[i]<-median(Mouse_cis_trans_on_samchr_bin_by_distance_observed[[i]]$distance)
}
cor(O_Cor_increase,ODistance,method="spearman")
#observed=-0.6004242,p<0.001
#spearman:cor=-0.7628683,p<0.001
cor_simulated<-data.frame(cor_shuffled)
names(cor_simulated)<-c("spcor")
write.table(cor_simulated,file="simulated_spcor_cf_vs_distance",sep="\t",col.names=TRUE,row.names=FALSE)

