#Fig3b test
library(dplyr)
library(ggplot2)
library(reshape2)
Hic_tab<-read.table(file="../Data_needed/Hic_tab_all",sep="\t",header=TRUE)
Hic_tab_samechr<-Hic_tab%>%
  filter(chr_1==chr_2)
Hic_tab_difchr<-Hic_tab%>%
  filter(chr_1!=chr_2)

#same chr
Hic_tab_samechr_independent_filtered<-Hic_tab_samechr[sample(c(1:552736)),]
samechr_pairs_remaining<-552736
same_chr_g1<-as.numeric(as.character(Hic_tab_samechr_independent_filtered$g_1))
same_chr_g2<-as.numeric(as.character(Hic_tab_samechr_independent_filtered$g_2))
cis_trans_difv<-as.numeric(as.character(Hic_tab_samechr_independent_filtered$Contact_dif))
samechr_genes_1<-list()
samechr_genes_2<-list()
cis_trans_dif<-list()
c<-1
while(samechr_pairs_remaining>0){
  samechr_genes_1[c]<-same_chr_g1[1]
  samechr_genes_2[c]<-same_chr_g2[1]
  cis_trans_dif[c]<-cis_trans_difv[1]
  filter_vector<-(same_chr_g1!=samechr_genes_1[c]&
                    same_chr_g2!=samechr_genes_1[c]&
                    same_chr_g1!=samechr_genes_2[c]&
                    same_chr_g2!=samechr_genes_2[c])
  samechr_pairs_remaining<-sum(filter_vector)
  same_chr_g1<-same_chr_g1[filter_vector]
  same_chr_g2<-same_chr_g2[filter_vector]
  cis_trans_difv<-cis_trans_difv[filter_vector]
  c<-c+1
}
cis_trans_difr<-as.numeric(cis_trans_dif)
sum(cis_trans_difr==0)
#same chr: 2222,0,2(>0,==0,<0)
binom.test(2222,2224,p=0.5)$p.value
binom.test(2222,2224,p=0.5)

# Exact binomial test
# 
# data:  2222 and 2224
# number of successes = 2222, number of trials = 2224, p-value < 2.2e-16
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#   0.9967553 0.9998911
# sample estimates:
#   probability of success 
# 0.9991007 
#dif chr
Hic_tab_difchr_independent_filtered<-Hic_tab_difchr[sample(c(1:9399755)),]
difchr_pairs_remaining<-9399755
dif_chr_g1<-as.numeric(as.character(Hic_tab_difchr_independent_filtered$g_1))
dif_chr_g2<-as.numeric(as.character(Hic_tab_difchr_independent_filtered$g_2))
dcis_trans_difv<-as.numeric(as.character(Hic_tab_difchr_independent_filtered$Contact_dif))
difchr_genes_1<-list()
difchr_genes_2<-list()
dcis_trans_dif<-list()
c<-1
while(difchr_pairs_remaining>0){
  difchr_genes_1[c]<-dif_chr_g1[1]
  difchr_genes_2[c]<-dif_chr_g2[1]
  dcis_trans_dif[c]<-dcis_trans_difv[1]
  filter_vector<-(dif_chr_g1!=difchr_genes_1[c]&
                    dif_chr_g2!=difchr_genes_1[c]&
                    dif_chr_g1!=difchr_genes_2[c]&
                    dif_chr_g2!=difchr_genes_2[c])
  difchr_pairs_remaining<-sum(filter_vector)
  dif_chr_g1<-dif_chr_g1[filter_vector]
  dif_chr_g2<-dif_chr_g2[filter_vector]
  dcis_trans_difv<-dcis_trans_difv[filter_vector]
  c<-c+1
}
dcis_trans_difr<-as.numeric(dcis_trans_dif)
sum(dcis_trans_difr==0)
#dif chr: 499,1228,504(>0,==0,<0)
binom.test(499,499+504,p=0.5)
#Exact binomial test

# data:  499 and 499 + 504
# number of successes = 499, number of trials = 1003, p-value = 0.8995
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#   0.4661121 0.5289176
# sample estimates:
#   probability of success 
# 0.4975075 

###############################################################################################
#Fig3c test
library(dplyr)
cis_129_Hic<-read.table(file="../Data_needed/chrAll_cis_129.matrix",sep="\t")
cis_129_HicMat<-as.matrix(cis_129_Hic[2:4956,2:4956])
cis_Cas_Hic<-read.table(file="../Data_needed/chrAll_cis_cas.matrix",sep="\t")
cis_Cas_HicMat<-as.matrix(cis_Cas_Hic[2:4956,2:4956])
trans_129Cas_Hic<-read.table(file="../Data_needed/chrAll_trans_129cas.matrix",sep="\t")
trans_129Cas_HicMat<-as.matrix(trans_129Cas_Hic[2:4956,2:4956])

bin_infor<-as.character(cis_129_Hic$V1)[2:4956]
bin_chr_infor<-as.character(lapply(bin_infor,function(x){unlist(strsplit(x,split="\\|"))[3]}))
chr_vector<-as.character(lapply(bin_chr_infor,function(x){unlist(strsplit(x,split="-"))[1]}))
pos_vector<-as.character(lapply(bin_chr_infor,function(x){unlist(strsplit(x,split=":"))[2]}))
pos_start<-as.numeric(lapply(pos_vector,function(x){unlist(strsplit(x,split="-"))[[1]]}))
pos_end<-as.numeric(lapply(pos_vector,function(x){unlist(strsplit(x,split="-"))[[2]]}))
pos_mid<-(pos_start+pos_end)/2
filter_idx<-(rowSums(cis_129_HicMat=="nan")!=4955)&(rowSums(cis_Cas_HicMat=="nan")!=4955)
chr_vector<-chr_vector[filter_idx]
pos_mid<-pos_mid[filter_idx]
cis_129_HicMat<-cis_129_HicMat[filter_idx,filter_idx]
cis_Cas_HicMat<-cis_Cas_HicMat[filter_idx,filter_idx]
trans_129Cas_HicMat<-trans_129Cas_HicMat[filter_idx,filter_idx]

cis_129_HicMat<-apply(cis_129_HicMat,1,as.numeric)
cis_Cas_HicMat<-apply(cis_Cas_HicMat,1,as.numeric)
trans_129Cas_HicMat<-apply(trans_129Cas_HicMat,1,as.numeric)

HiC_dif_all<-(cis_129_HicMat+cis_Cas_HicMat)-(trans_129Cas_HicMat+t(trans_129Cas_HicMat))
#shuffling
shuffled_cor<-numeric(999)
shuffled_index<-c(1:4462)
for(j in 1:999){
  shuffled_index<-sample(shuffled_index)
  chr_vector_sampled<-chr_vector[shuffled_index]
  pos_mid_sampled<-pos_mid[shuffled_index]
  pair_num<-10000
  k<-0
  gene_1_index<-numeric(pair_num)
  gene_2_index<-numeric(pair_num)
  while(k<pair_num){
    gene_pair<-sample(4462,2,replace=FALSE)
    if(chr_vector_sampled[gene_pair[1]]==chr_vector_sampled[gene_pair[2]]){
      k<-k+1
      gene_1_index[k]<-gene_pair[1]
      gene_2_index[k]<-gene_pair[2]
    }
  }
  
  pos_gene_1<-pos_mid_sampled[gene_1_index]
  pos_gene_2<-pos_mid_sampled[gene_2_index]
  cis_trans_dif<-numeric(pair_num)
  for(i in 1:pair_num){
    cis_trans_dif[i]<-as.numeric(HiC_dif_all[gene_1_index[i],gene_2_index[i]])
  }
  bin_distance<-abs(pos_gene_1-pos_gene_2)
  Dat<-data.frame(bin_distance,cis_trans_dif)
  names(Dat)<-c("distance","Contact_frequency_difference")
  #Bin<-quantile(as.numeric(as.character(Dat$distance)),(0:100)/100)
  Contact_difference_bin_by_distance<-split(Dat,cut(Dat$distance,100))
  Contact_difference<-numeric(100)
  Distance<-numeric(100)
  for(l in 1:100){
    Contact_difference[l]<-median(Contact_difference_bin_by_distance[[l]]$Contact_frequency_difference)
    Distance[l]<-median(Contact_difference_bin_by_distance[[l]]$distance)
  }
  
  shuffled_cor[j]<-cor(Distance,Contact_difference,method="spearman")
  
}
cor(c(0,1),c(0,0),method="spearman")
min(shuffled_cor[is.na(shuffled_cor)!=TRUE])
#p<0.001

#Fig3d test needed
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
#same chr independent:

ATAC_sam_chr_independent_filtered<-ATAC_50_tab_sam_chr[sample(c(1:451684)),]
sam_chr_pairsnum<-451684
k<-1
while((k<sam_chr_pairsnum)){
  G1<-as.numeric(as.character(ATAC_sam_chr_independent_filtered$G_1_index[k]))
  G2<-as.numeric(as.character(ATAC_sam_chr_independent_filtered$G_2_index[k]))
  subset_vector<-c(rep(TRUE,k),(ATAC_sam_chr_independent_filtered$G_1_index[-c(1:k)]!=G1)&(ATAC_sam_chr_independent_filtered$G_1_index[-c(1:k)]!=G2)&
                     (ATAC_sam_chr_independent_filtered$G_2_index[-c(1:k)]!=G1)&(ATAC_sam_chr_independent_filtered$G_2_index[-c(1:k)]!=G2))
  ATAC_sam_chr_independent_filtered<-ATAC_sam_chr_independent_filtered[subset_vector,]
  sam_chr_pairsnum<-length(as.character(ATAC_sam_chr_independent_filtered$G_1_index))
  k<-k+1
}
sum(ATAC_sam_chr_independent_filtered$cis_trans_dif>0)
1020/1788
binom.test(1020,1788,p=0.5)
#Exact binomial test

# data:  1020 and 1788
# number of successes = 1020, number of trials = 1788, p-value = 2.756e-09
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#   0.5471457 0.5935623
# sample estimates:
#   probability of success 
# 0.5704698 
ATAC_dif_chr_independent_filtered<-ATAC_50_tab_dif_chr[sample(c(1:5979807)),]
dif_chr_pairsnum<-5979807
k<-1
while((k<dif_chr_pairsnum)){
  G1<-as.numeric(as.character(ATAC_dif_chr_independent_filtered$G_1_index[k]))
  G2<-as.numeric(as.character(ATAC_dif_chr_independent_filtered$G_2_index[k]))
  subset_vector<-c(rep(TRUE,k),(ATAC_dif_chr_independent_filtered$G_1_index[-c(1:k)]!=G1)&(ATAC_dif_chr_independent_filtered$G_1_index[-c(1:k)]!=G2)&
                     (ATAC_dif_chr_independent_filtered$G_2_index[-c(1:k)]!=G1)&(ATAC_dif_chr_independent_filtered$G_2_index[-c(1:k)]!=G2))
  ATAC_dif_chr_independent_filtered<-ATAC_dif_chr_independent_filtered[subset_vector,]
  dif_chr_pairsnum<-length(as.character(ATAC_dif_chr_independent_filtered$G_1_index))
  k<-k+1
}
sum(ATAC_dif_chr_independent_filtered$cis_trans_dif>0)
911/1793
binom.test(911,1793,p=0.5)
#binom.test(911,1793,p=0.5)

# Exact binomial test
# 
# data:  911 and 1793
# number of successes = 911, number of trials = 1793, p-value = 0.5085
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#   0.4846705 0.5314770
# sample estimates:
#   probability of success 
# 0.508087 
write.table(ATAC_dif_chr_independent_filtered,sep="\t",col.names=TRUE,row.names=FALSE,file="ATAC_50_tab_dif_chr_independent_control")
write.table(ATAC_sam_chr_independent_filtered,sep="\t",col.names=TRUE,row.names=FALSE,file="ATAC_50_tab_sam_chr_independent_control")
##########################################################################################################################################
#Fig3e
library(dplyr)
library(ggplot2)
ATAC_50_tab<-read.table(file="../Data_needed/ATAC50_gpairs",sep="\t")
names(ATAC_50_tab)<-c("G_1_chr","G_2_chr","G_1_start","G_2_start","G_1_end","G_2_end","G_1_index",
                      "G_2_index","cis_M129","cis_cas","M129_cas","cas_M129","M129_g1_reads",
                      "M129_g2_reads","cas_g1_reads","cas_g2_reads")
ATAC_50_tab$g_1_pos<-(ATAC_50_tab$G_1_start+ATAC_50_tab$G_1_end)/2
ATAC_50_tab$g_2_pos<-(ATAC_50_tab$G_2_start+ATAC_50_tab$G_2_end)/2
all_bins<-c(as.numeric(as.character(ATAC_50_tab$G_1_index)),
            as.numeric(as.character(ATAC_50_tab$G_2_index)))
all_pos<-c(as.numeric(as.character(ATAC_50_tab$g_1_pos)),
           as.numeric(as.character(ATAC_50_tab$g_2_pos)))
all_chr<-c(as.character(ATAC_50_tab$G_1_chr),
           as.character(ATAC_50_tab$G_2_chr))
dat_for_shuffled<-data.frame(all_bins,all_chr,all_pos)
dat_simplified<-dat_for_shuffled%>%
  group_by(all_bins)%>%
  dplyr::summarise(chr=all_chr[1],pos=mean(all_pos))
all_index<-as.numeric(as.character(dat_simplified$all_bins))
all_pos_forsam<-as.numeric(as.character(dat_simplified$pos))
all_chr_forsam<-as.character(dat_simplified$chr)
g1_back_index<-match(ATAC_50_tab$G_1_index,all_index)
g2_back_index<-match(ATAC_50_tab$G_2_index,all_index)
cor_shuffled<-numeric(999)
for(i in 1:999){
  index_shuffled<-sample(c(1:3587))
  chr_shuffled<-all_chr_forsam[index_shuffled]
  pos_shuffled<-all_pos_forsam[index_shuffled]
  
  ATAC_50_tab$g_1_chrsam<-chr_shuffled[g1_back_index]
  ATAC_50_tab$g_2_chrsam<-chr_shuffled[g2_back_index]
  ATAC_50_tab$g_1_possam<-pos_shuffled[g1_back_index]
  ATAC_50_tab$g_2_possam<-pos_shuffled[g2_back_index]
  
  ATAC_50_tab_on_sam_chr<-ATAC_50_tab%>%
    filter(g_1_chrsam==g_2_chrsam)
  ATAC_50_tab_on_sam_chr$cis_trans_dif<-((ATAC_50_tab_on_sam_chr$cis_M129+ATAC_50_tab_on_sam_chr$cis_cas)-(ATAC_50_tab_on_sam_chr$M129_cas+ATAC_50_tab_on_sam_chr$cas_M129))/2
  ATAC_50_tab_on_sam_chr$distance<-abs(ATAC_50_tab_on_sam_chr$g_1_possam-ATAC_50_tab_on_sam_chr$g_2_possam)
  ATAC_50_tab_on_sam_chr_bin_by_distance<-split(ATAC_50_tab_on_sam_chr,cut(ATAC_50_tab_on_sam_chr$distance,100))
  cor_increase<-numeric(100)
  Distance<-numeric(100)
  for(j in 1:100){
    cor_increase[j]<-median(ATAC_50_tab_on_sam_chr_bin_by_distance[[j]]$cis_trans_dif)
    Distance[j]<-median(ATAC_50_tab_on_sam_chr_bin_by_distance[[j]]$distance)
  }
  cor_shuffled[i]<-cor(Distance,cor_increase,method="spearman")
}
hist(cor_shuffled)
#p<0.001
ATAC_50_tab_observed_on_samchr<-ATAC_50_tab%>%
  filter(G_1_chr==G_2_chr)
ATAC_50_tab_observed_on_samchr$cis_trans_dif<-(ATAC_50_tab_observed_on_samchr$cis_M129+ATAC_50_tab_observed_on_samchr$cis_cas-ATAC_50_tab_observed_on_samchr$M129_cas-ATAC_50_tab_observed_on_samchr$cas_M129)/2
ATAC_50_tab_observed_on_samchr$distance<-abs(ATAC_50_tab_observed_on_samchr$g_1_pos-ATAC_50_tab_observed_on_samchr$g_2_pos)
ATAC_50_tab_observed_on_samchr_bin<-split(ATAC_50_tab_observed_on_samchr,cut(ATAC_50_tab_observed_on_samchr$distance,100))
OCI<-numeric(100)
OD<-numeric(100)
for(i in 1:100){
  OCI[i]<-median(ATAC_50_tab_observed_on_samchr_bin[[i]]$cis_trans_dif)
  OD[i]<-median(ATAC_50_tab_observed_on_samchr_bin[[i]]$distance)
}
cor(OD,OCI,method="spearman")
#spearman:-0.8963456
cor_simulated<-data.frame(cor_shuffled)
names(cor_simulated)<-c("spearmanATAC")
write.table(cor_simulated,sep="\t",file="Spcor_simulated_ATAC",col.names=TRUE,row.names=FALSE)

##########################################################################################################################################
#Fig3f
library(dplyr)
library(ggplot2)
library(ade4)
library(vegan)
cis_trans_acceCor_vs_Hic<-read.table(file="../Data_needed/cis_trans_CorTabAllWithP_na_rm_for_sum",sep="\t",header=TRUE)
cis_trans_acceCor_vs_Hic$trans_Hic_sum<-cis_trans_acceCor_vs_Hic$trans_129cas_Hic+cis_trans_acceCor_vs_Hic$trans_cas129_Hic
cis_trans_acceCor_vs_Hic$trans_cor_sum<-cis_trans_acceCor_vs_Hic$V7+cis_trans_acceCor_vs_Hic$V8
cis_trans_acceCor_vs_Hic<-cis_trans_acceCor_vs_Hic%>%
  filter(is.na(trans_cor_sum)!=TRUE)
old_index<-unique(c(as.numeric(as.character(cis_trans_acceCor_vs_Hic$V3)),
                    as.numeric(as.character(cis_trans_acceCor_vs_Hic$V4))))
old_index<-sort(old_index)
new_index<-c(1:1597)
bin_1_matchi<-match(cis_trans_acceCor_vs_Hic$V3,old_index)
bin_2_matchi<-match(cis_trans_acceCor_vs_Hic$V4,old_index)
bin_1_index<-new_index[bin_1_matchi]
bin_2_index<-new_index[bin_2_matchi]
Hic_trans_sum<-as.numeric(as.character(cis_trans_acceCor_vs_Hic$trans_Hic_sum))
ATAC_trans_corsum<-as.numeric(as.character(cis_trans_acceCor_vs_Hic$trans_cor_sum))

Hic_trans_matrix<-matrix(NA,1597,1597)
ATAC_trans_cormatrix<-matrix(NA,1597,1597)
for(i in 1:1274406){
  Hic_trans_matrix[bin_1_index[i],bin_2_index[i]]<-Hic_trans_sum[i]
  Hic_trans_matrix[bin_2_index[i],bin_1_index[i]]<-Hic_trans_sum[i]
  ATAC_trans_cormatrix[bin_1_index[i],bin_2_index[i]]<-ATAC_trans_corsum[i]
  ATAC_trans_cormatrix[bin_2_index[i],bin_1_index[i]]<-ATAC_trans_corsum[i]
  
}


Hic_trans_dist<-as.dist(Hic_trans_matrix)
ATAC_trans_cordist<-as.dist(ATAC_trans_cormatrix)
spcor_simulated<-mantel(Hic_trans_dist,ATAC_trans_cordist,method="spearman",permutations = 9999)
hist(spcor_simulated$perm,breaks=50)
hist(spcor_simulated$perm,col="grey",xlab="Simulated correlation",breaks=50,main=NULL,xlim=c(-0.02,0.08))
arrows(x0=spcor_simulated$statistic,y0=200,x1=spcor_simulated$statistic,y1=0,col="red",angle=10)
text(x=spcor_simulated$statistic,y=260,labels="observed")
text(x=spcor_simulated$statistic,y=310,labels="p<0.0001")
spcor_simulated
#Mantel statistic based on Spearman's rank correlation rho 

# Call:
#   mantel(xdis = Hic_trans_dist, ydis = ATAC_trans_cordist, method = "spearman",      permutations = 9999) 
# 
# Mantel statistic r: 0.0691 
# Significance: 1e-04 
# 
# Upper quantiles of permutations (null model):
#   90%     95%   97.5%     99% 
#   0.00434 0.00565 0.00672 0.00816 
# Permutation: free
# Number of permutations: 9999
#p<0.0001
Hic_coATAC_simulated<-as.numeric(spcor_simulated$perm)
Hic_coATAC_Mantel_tab<-data.frame(Hic_coATAC_simulated)
names(Hic_coATAC_Mantel_tab)<-c("simulated correlation")
write.table(Hic_coATAC_Mantel_tab,sep="\t",row.names=FALSE,col.names=TRUE,file="Hic_coATAC_Mantel")
###################################################################################################
#Fig3g test
library(dplyr)
library(ade4)
library(vegan)
Mouse_cis_trans<-read.table(file="../Data_needed/Mouse_cis_trans_for_cor_withATAC_filtered",sep="\t",header=TRUE)
Mouse_ATAC_cor<-read.table(file="../Data_needed/ATAC_cor_gpairs",sep="\t")
names(Mouse_ATAC_cor)<-c("g_1_chr","g_2_chr","g_1_start","g_2_start","g_1_end","g_2_end","g_1_index","g_2_index",
                         "cis_M129","cis_cas","M129_cas","cas_M129","M129_reads_g1","M129_reads_g2","cas_reads_g1",
                         "cas_reads_g2")

1562*1561/2
sum(Mouse_cis_trans$g_1_posindex==Mouse_ATAC_cor$g_1_index)
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


cor.test(Mouse_cis_trans_group_by_index$cor,Mouse_cis_trans_group_by_index$cor_ATAC,method="spearman")

all_index<-unique(c(as.numeric(as.character(Mouse_cis_trans_group_by_index$pos_1)),as.numeric(as.character(Mouse_cis_trans_group_by_index$pos_2))))
ordered_index<-all_index[order(all_index)]
index_vector<-c(1:978)
p_1_indexback<-match(Mouse_cis_trans_group_by_index$pos_1,ordered_index)
p_2_indexback<-match(Mouse_cis_trans_group_by_index$pos_2,ordered_index)
matrix_index_1<-index_vector[p_1_indexback]
matrix_index_2<-index_vector[p_2_indexback]
cor_trans_matrix<-matrix(NA,978,978)
ATAC_trans_matrix<-matrix(NA,978,978)
cor_vector<-as.numeric(as.character(Mouse_cis_trans_group_by_index$cor))
ATAC_vector<-as.numeric(as.character(Mouse_cis_trans_group_by_index$cor_ATAC))
for(i in 1:477998){
  cor_trans_matrix[matrix_index_1[i],matrix_index_2[i]]<-cor_vector[i]
  cor_trans_matrix[matrix_index_2[i],matrix_index_1[i]]<-cor_vector[i]
  ATAC_trans_matrix[matrix_index_1[i],matrix_index_2[i]]<-ATAC_vector[i]
  ATAC_trans_matrix[matrix_index_2[i],matrix_index_1[i]]<-ATAC_vector[i]
}
cor_dist<-as.dist(cor_trans_matrix)
ATAC_dist<-as.dist(ATAC_trans_matrix)
cor_simulated<-mantel.randtest(cor_dist,ATAC_dist,nrepet = 999)
cor_simulated
spcor_simulated<-mantel(cor_dist,ATAC_dist,method="spearman",permutations = 999)
spcor_simulated
spcor_simulated_large<-mantel(cor_dist,ATAC_dist,method="spearman",permutations = 9999)
spcor_simulated_large$statistic
spcor_simulated_large$signif
cor_simulated_large<-mantel(cor_dist,ATAC_dist,method="pearson",permutations = 9999)
Cofluctuation_coATAC_tab<-data.frame(cor_simulated_large$perm,spcor_simulated_large$perm)
names(Cofluctuation_coATAC_tab)<-c("pearson","spearman")
write.table(Cofluctuation_coATAC_tab,file="Cofluctuation_coATAC_tab",sep="\t",row.names=FALSE,col.names=TRUE)
cor_simulated_large
#Mantel statistic based on Pearson's product-moment correlation 

# Call:
#   mantel(xdis = cor_dist, ydis = ATAC_dist, method = "pearson",      permutations = 9999) 
# 
# Mantel statistic r: 0.02466 
# Significance: 0.012 
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0138 0.0181 0.0216 0.0254 
# Permutation: free
# Number of permutations: 9999
spcor_simulated_large
#Mantel statistic based on Spearman's rank correlation rho 

# Call:
#   mantel(xdis = cor_dist, ydis = ATAC_dist, method = "spearman",      permutations = 9999) 
# 
# Mantel statistic r: 0.02184 
# Significance: 0.027 
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0148 0.0186 0.0222 0.0268 
# Permutation: free
# Number of permutations: 9999
hist(spcor_simulated_large$perm,col="grey",xlab="Simulated correlation",breaks=50,main=NULL)
arrows(x0=spcor_simulated_large$statistic,y0=300,x1=spcor_simulated$statistic,y1=150,col="red",angle=10)
text(x=spcor_simulated$statistic,y=330,labels="observed")
text(x=spcor_simulated$statistic,y=380,labels="p=0.027")
