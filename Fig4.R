library(dplyr)
library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
#Fig4a
A<-1
B<-1
K_v<-c(10^5,10^2,10^3)
cell_number<-10000
covar_vector<-c(0:100)/2500
plots<-list()
Mean_concentration<-numeric(101)
K<-K_v[1]
for(i in 1:101){
  Sigma <- matrix(c(0.04,covar_vector[i],covar_vector[i],0.04),2,2)
  E_pair_correlated<-mvrnorm(n = cell_number, rep(A,2), Sigma, empirical = TRUE)
  E_1<-E_pair_correlated[,1]
  E_2<-E_pair_correlated[,2]
  E_1[E_1<0]<-0
  E_2[E_2<0]<-0
  concentration_vector<-numeric(cell_number)
  for(j in 1:length(concentration_vector)){
    f1c<-c(K,K*(E_2[j]-E_1[j])+1,-E_1[j])
    f2c<-c(K,K*(E_1[j]-E_2[j])+1,-E_2[j])
    E1f<-(-f1c[2]+sqrt(f1c[2]^2-4*f1c[1]*f1c[3]))/(2*f1c[1])
    E2f<-(-f2c[2]+sqrt(f2c[2]^2-4*f2c[1]*f2c[3]))/(2*f2c[1])
    concentration_vector[j]<-K*E1f*E2f
  }
  Mean_concentration[i]<-mean(concentration_vector)
}
Dat<-data.frame(covar_vector*25,Mean_concentration)
names(Dat)<-c("correlation","Mean_concentration")
plot(Dat$correlation,Dat$Mean_concentration)
Dat$Mean_concentration[100]/Dat$Mean_concentration[1]

kt<-paste("K=",K,sep="")
png(filename="Fig4a",width=10,height=10,units="in",res=600)
ggplot(Dat, aes(x = correlation, y = Mean_concentration)) + 
  geom_point() +
  xlab(expression(paste("Correlation between ","[A]"["t"]," and ","[B]"["t"])))+
  ylab("Mean [AB]")+
  theme(axis.title.x=element_text(size=33,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=35,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=33,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=35,family="TT Arial",color="black"))+
  theme(plot.title=element_text(hjust=0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits=c(0,1.1))+
  scale_y_continuous(breaks=c(0.90,0.95,1.00))+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
dev.off()
####################################################################################################
#Fig4b
dat_simulated<-read.table(file="../Data_needed/simulated_complex_protein_coding_only",sep="\t",header=TRUE)
median(dat_simulated$samchr_pair_number)
png(file="Fig4b",width=10,height=10,units="in",res=600)
ggplot(dat_simulated, aes(x=samchr_pair_number)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)))+
  geom_segment(aes(x = 200, y = 240/10000, xend = 200, yend = 50/10000), size=2,arrow = arrow(length = unit(0.5, "cm")))+
  ylab(expression(atop("Frequency",paste("(same complex)"))))+
  xlab(label=expression(atop("Number of pairs", "of linked genes")))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,angle=90,vjust=0.5,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
dev.off()
####################################################################################################
#Fig4c
dat_simulated<-read.table(file="../Data_needed/simulated_complex_from_all_protein_coding_controlled",sep="\t",header=TRUE)
png(file="Fig4c",width=10,height=10,units="in",res=600)
ggplot(dat_simulated, aes(x=samchr_pair_number)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)))+
  geom_segment(aes(x = 14101, y = 910/10000, xend = 14101, yend = 760/10000), size=2,arrow = arrow(length = unit(0.5, "cm")))+
  ylab(expression(atop("Frequency",paste("(different complexes)"))))+
  xlab(label=expression(atop("Number of pairs", "of linked genes")))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,angle=90,vjust=0,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#panel.background = element_blank(), axis.line = element_blank())
dev.off()
####################################################################################################
#Fig4d
dat_simulated<-read.table(file="../Data_needed/simulated_complex_protein_coding_only",sep="\t",header=TRUE)
png(file="Fig4d",width=10,height=10,units="in",res=600)
ggplot(dat_simulated, aes(x=simuated_distance/1000000)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)))+
  geom_segment(aes(x = 33267366/1000000, y = 1000/10000, xend = 33267366/1000000, yend = 850/10000), size=2,arrow = arrow(length = unit(0.5, "cm")))+
  ylab(expression(atop("Frequency",paste("(same complex)"))))+xlab(label=expression(atop("Distance between", "linked genes (Mb)")))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,angle=90,vjust=0.5,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour ="black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
dev.off()
#####################################################################################################
#Fig4e
dat_simulated<-read.table(file="../Data_needed/simulated_complex_from_all_protein_coding_controlled",sep="\t",header=TRUE)
png(file="Fig4e",width=10,height=10,units="in",res=600)
ggplot(dat_simulated, aes(x=simuated_distance/1000000)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)))+
  geom_segment(aes(x = 34114963/1000000, y = 500/10000, xend = 34114963/1000000, yend = 300/10000), size=2,arrow = arrow(length = unit(0.5, "cm")))+
  ylab(expression(atop("Frequency",paste("(different complexes)"))))+xlab(label=expression(atop("Distance between", "linked genes (Mb)")))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,angle=90,vjust=0.5,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))
dev.off()
######################################################################################################
#Fig4f
library(dplyr)
library(ggplot2)
dat_simulated<-read.table(file="..Data_needed/sampled_linked_tab500",sep="\t",header=TRUE)
png(file="Fig4f",width=10,height=6,units="in",res=600)
ggplot(dat_simulated, aes(x=sampled_linked_Mouse)) +
  theme_bw() +
  geom_histogram(colour="black", fill="grey70",aes(y=..count../sum(..count..)),bins=30)+
  geom_segment(aes(x = 25, y = 30/1000, xend = 25, yend = 10/1000), size=2,arrow = arrow(length = unit(0.5, "cm")))+
  ylab(label="Frequency")+
  xlab(label=expression(atop("Number of evolved", "linked gene pairs")))+
  theme(axis.text.y=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=35,angle=90,vjust=0,family="TT Arial",color="black",margin=margin(t=0,r=35,b=0,l=0)))+
  theme(axis.text.x=element_text(size=30,family="TT Arial",color="black"))+
  theme(axis.title.x=element_text(size=35,family="TT Arial",color="black",margin=margin(t=35,r=0,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(axis.ticks.length = unit(5,"mm"))+
  theme(axis.ticks=element_line(size=2))+
  xlim(c(1,32))
#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#panel.background = element_blank(), axis.line = element_blank())
dev.off()
###################################################################################################
#Fig4g
library(dplyr)
library(ggplot2)
Mouse_tissue_expressionRaw<-read.table(file="../Data_needed/mouse_rpkm.txt",sep="\t",header=TRUE)
Mouse_genes<-read.table(file="../Data_needed/Mouse_genes_all.txt",sep="\t",header=TRUE)
H1_family<-c("Hist1h1a","Hist1h1b","Hist1h1c",
             "Hist1h1d","Hist1h1e","H1f0","H1foo","Hist1h1t","H1fnt","H1fx", "Hils1")
gene_names_index<-match(row.names(Mouse_tissue_expressionRaw),Mouse_genes$Gene.stable.ID)
Mouse_tissue_expressionRaw$geneName<-Mouse_genes$Gene.name[gene_names_index]
Mouse_tissue_expression<-data.frame(as.character(Mouse_tissue_expressionRaw$geneName))
names(Mouse_tissue_expression)<-c("geneName")
#Pancreas: 1, 14, 27
#liver: 2, 15,28
#Stomach: 3, 16,29
#DuoDenum: 4, 17, 30
#Jejunum: 5, 18, 31
#Ileum: 6,19,32,
#Colon:7, 20, 33
#Kidney: 8,21,34
#Quadriceps: 9,22,35
#Thymus:10,23,36
#Heart: 11, 24,37
#Esophagus: 12, 25,38
#Brain: 13, 26, 39
Mouse_tissue_expression$Pancreas<-rowMeans(Mouse_tissue_expressionRaw[,c(1,14,27)])
Mouse_tissue_expression$liver<-rowMeans(Mouse_tissue_expressionRaw[,c(2,15,28)])
Mouse_tissue_expression$Stomach<-rowMeans(Mouse_tissue_expressionRaw[,c(3,16,29)])
Mouse_tissue_expression$DuoDenum<-rowMeans(Mouse_tissue_expressionRaw[,c(4,17,30)])
Mouse_tissue_expression$Jejunum<-rowMeans(Mouse_tissue_expressionRaw[,c(5,18,31)])
Mouse_tissue_expression$Ileum<-rowMeans(Mouse_tissue_expressionRaw[,c(6,19,32)])
Mouse_tissue_expression$Colon<-rowMeans(Mouse_tissue_expressionRaw[,c(7,20,33)])
Mouse_tissue_expression$Kidney<-rowMeans(Mouse_tissue_expressionRaw[,c(8,21,34)])
Mouse_tissue_expression$Qudriceps<-rowMeans(Mouse_tissue_expressionRaw[,c(9,22,35)])
Mouse_tissue_expression$Thymus<-rowMeans(Mouse_tissue_expressionRaw[,c(10,23,36)])
Mouse_tissue_expression$Heart<-rowMeans(Mouse_tissue_expressionRaw[,c(11,24,37)])
Mouse_tissue_expression$Esophagus<-rowMeans(Mouse_tissue_expressionRaw[,c(12,25,38)])
Mouse_tissue_expression$Brain<-rowMeans(Mouse_tissue_expressionRaw[,c(13,26,39)])
#histone data
H1_gene_index<-match(H1_family,Mouse_tissue_expression$geneName)
Mouse_H1_expression<-Mouse_tissue_expression[H1_gene_index,]
Histone_expression_pattern<-colSums(Mouse_H1_expression[,c(2:14)])
#selected complex data
Mouse_complex_pairs<-read.table(file="../Data_needed/Mouse_complex_pairs_without_removing_duplicated",sep="\t",header=TRUE)
Rat_complex_pairs<-read.table(file="../Data_needed/p_pairs_rat_in_sam_com",sep="\t",header=TRUE)
Human_complex_pairs<-read.table(file="../Data_needed/p_pairs_Human_in_sam_com",sep="\t",header=TRUE)
Human_pair_vector<-c(as.character(paste(Human_complex_pairs$g_1,Human_complex_pairs$g_2,sep="_")),as.character(paste(Human_complex_pairs$g_2,Human_complex_pairs$g_1,sep="_")))
Mouse_pair_vector<-c(as.character(paste(Mouse_complex_pairs$g_1,Mouse_complex_pairs$g_2,sep="_")),as.character(paste(Mouse_complex_pairs$g_2,Mouse_complex_pairs$g_1,sep="_")))
Rat_pair_vector<-c(as.character(paste(Rat_complex_pairs$g_1,Rat_complex_pairs$g_2,sep="_")),as.character(paste(Rat_complex_pairs$g_2,Rat_complex_pairs$g_1,sep="_")))
Human_Mouse_rat_homolog<-read.table(file="../Data_needed/Mouse_Rat_Human_homologs.txt",sep="\t",header=TRUE)
Human_Mouse_rat_homolog_one2one<-Human_Mouse_rat_homolog%>%
  filter(Rat.homology.type=="ortholog_one2one")%>%
  filter(Human.homology.type=="ortholog_one2one")%>%
  filter(Chromosome.scaffold.name!="MT"&Rat.chromosome.scaffold.name!="MT"&Human.chromosome.scaffold.name!="MT")
Mouse_g1_one2one_index<-match(Mouse_complex_pairs$g_1,Human_Mouse_rat_homolog_one2one$Gene.name)
Mouse_g2_one2one_index<-match(Mouse_complex_pairs$g_2,Human_Mouse_rat_homolog_one2one$Gene.name)
Mouse_complex_pairs$g_1_one2one_rat<-Human_Mouse_rat_homolog_one2one$Rat.gene.name[Mouse_g1_one2one_index]
Mouse_complex_pairs$g_2_one2one_rat<-Human_Mouse_rat_homolog_one2one$Rat.gene.name[Mouse_g2_one2one_index]
Mouse_complex_pairs$g_1_one2one_Human<-Human_Mouse_rat_homolog_one2one$Human.gene.name[Mouse_g1_one2one_index]
Mouse_complex_pairs$g_2_one2one_Human<-Human_Mouse_rat_homolog_one2one$Human.gene.name[Mouse_g2_one2one_index]
Mouse_complex_pairs$g_1_chr<-Human_Mouse_rat_homolog_one2one$Chromosome.scaffold.name[Mouse_g1_one2one_index]
Mouse_complex_pairs$g_2_chr<-Human_Mouse_rat_homolog_one2one$Chromosome.scaffold.name[Mouse_g2_one2one_index]
Mouse_complex_pairs$Rat_g1_chr<-Human_Mouse_rat_homolog_one2one$Rat.chromosome.scaffold.name[Mouse_g1_one2one_index]
Mouse_complex_pairs$Rat_g2_chr<-Human_Mouse_rat_homolog_one2one$Rat.chromosome.scaffold.name[Mouse_g2_one2one_index]
Mouse_complex_pairs$Human_g1_chr<-Human_Mouse_rat_homolog_one2one$Human.chromosome.scaffold.name[Mouse_g1_one2one_index]
Mouse_complex_pairs$Human_g2_chr<-Human_Mouse_rat_homolog_one2one$Human.chromosome.scaffold.name[Mouse_g2_one2one_index]
Mouse_complex_pairs<-Mouse_complex_pairs%>%
  filter(is.na(g_1_one2one_rat)!=TRUE&is.na(g_2_one2one_rat)!=TRUE)

Mouse_complex_pairs$linked<-(Mouse_complex_pairs$g_1_chr==Mouse_complex_pairs$g_2_chr)
Mouse_complex_pairs$Human_linked<-(Mouse_complex_pairs$Human_g1_chr==Mouse_complex_pairs$Human_g2_chr)
Mouse_complex_pairs$rat_linked<-(Mouse_complex_pairs$Rat_g1_chr==Mouse_complex_pairs$Rat_g2_chr)
Mouse_complex_pairs$same_complex_rat<-as.character(paste(Mouse_complex_pairs$g_1_one2one_rat,Mouse_complex_pairs$g_2_one2one_rat,sep="_"))%in%Rat_pair_vector
Mouse_complex_pairs$same_complex_human<-as.character(paste(Mouse_complex_pairs$g_1_one2one_Human,Mouse_complex_pairs$g_2_one2one_Human,sep="_"))%in%Human_pair_vector


Conserved_complex<-Mouse_complex_pairs%>%
  filter(same_complex_human)
Conserved_complex_unlinked_Human_rat<-Conserved_complex%>%
  filter((!Human_linked)&(!rat_linked))%>%
  filter((g_1_chr!="MT"&g_2_chr!="MT"&Human_g1_chr!="MT"&Human_g2_chr!="MT"&Rat_g1_chr!="MT"&Rat_g2_chr!="MT"))
Selected_complex_candidate<-Conserved_complex_unlinked_Human_rat%>%
  filter(linked==TRUE)
complex_candidate<-unique(c(as.character(Selected_complex_candidate$g_1),as.character(Selected_complex_candidate$g_2)))
complex_gene_index<-match(complex_candidate,Mouse_tissue_expression$geneName)
Mouse_complex_expression<-Mouse_tissue_expression[complex_gene_index,]
#all complex pairs
Mouse_complex_pairsAll<-read.table(file="../Data_needed/Mouse_in_same_complex_pairs",sep="\t",header=TRUE)
Mouse_all_linkedCom<-Mouse_complex_pairsAll%>%
  filter(g_1_chr==g_2_chr)
Mouse_unlinkedCom<-Mouse_complex_pairsAll%>%
  filter(g_1_chr!=g_2_chr)
all_complex_can<-unique(c(as.character(Mouse_all_linkedCom$g_1),as.character(Mouse_all_linkedCom$g_2)))
all_unlinkedCom<-unique(c(as.character(Mouse_unlinkedCom$g_1),as.character(Mouse_unlinkedCom$g_2)))
abs_unlinked<-all_unlinkedCom[!all_unlinkedCom%in%all_complex_can]
linked_com_index<-match(all_complex_can,Mouse_tissue_expression$geneName)
Mouse_all_linkedcomplex<-Mouse_tissue_expression[linked_com_index,]
Mouse_all_linkedcomplex<-Mouse_all_linkedcomplex%>%
  filter(is.na(geneName)!=TRUE)
unlinked_com_index<-match(abs_unlinked,Mouse_tissue_expression$geneName)
Mouse_unlinked_ComExp<-Mouse_tissue_expression[unlinked_com_index,]
Mouse_unlinked_ComExp<-Mouse_unlinked_ComExp%>%
  filter(is.na(geneName)!=TRUE)
linked_com_expPattern<-colSums(Mouse_all_linkedcomplex[,c(2:14)])
unlinked_com_pattern<-colSums(Mouse_unlinked_ComExp[,c(2:14)])

#stratified & comparison
is_complex<-c(rep(1,length(all_complex_can)),rep(0,length(abs_unlinked)))
complex_all_dat<-data.frame(c(all_complex_can,abs_unlinked),is_complex)
names(complex_all_dat)<-c("genes","is_complex")
exp_pattern_index<-match(complex_all_dat$genes,Mouse_tissue_expression$geneName)
complex_all_dat[,3:15]<-Mouse_tissue_expression[exp_pattern_index,c(2:14)]
complex_all_dat$mean_exp<-rowMeans(complex_all_dat[,3:15])
complex_all_dat<-complex_all_dat%>%
  filter(is.na(mean_exp)!=TRUE)
bins<-quantile(complex_all_dat$mean_exp,prob=c(1:21)/21)
complex_all_stratified<-split(complex_all_dat,cut(complex_all_dat$mean_exp,breaks=bins))
stratified_complex_tab<-complex_all_dat[FALSE,]
set.seed(7)
for(i in 1:20){
  complex_stratified_i<-complex_all_stratified[[i]]
  complex_stratified_linked<-complex_stratified_i%>%
    filter(is_complex==1)
  complex_stratified_unlinked<-complex_stratified_i%>%
    filter(is_complex==0)
  linked_per_bin<-sum(complex_stratified_i$is_complex==1)
  unlinked_per_bin<-sum(complex_stratified_i$is_complex==0)
  if(min(linked_per_bin,unlinked_per_bin)>0){
    print(linked_per_bin)
    print(unlinked_per_bin)
    if(linked_per_bin<=unlinked_per_bin){
      sub_tab_index<-sample(c(1:unlinked_per_bin),linked_per_bin)
      sub_tab<-complex_stratified_unlinked[sub_tab_index,]
      stratified_complex_tab<-rbind(stratified_complex_tab,sub_tab)
      stratified_complex_tab<-rbind(stratified_complex_tab,complex_stratified_linked)
    }else{
      sub_tab_index<-sample(c(1:linked_per_bin),unlinked_per_bin)
      sub_tab<-complex_stratified_linked[sub_tab_index,]
      stratified_complex_tab<-rbind(stratified_complex_tab,sub_tab)
      stratified_complex_tab<-rbind(stratified_complex_tab,complex_stratified_unlinked)
    }
  }
}
wilcox.test(stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex==0],
            stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex==1])
# Wilcoxon rank sum test with continuity correction
# 
# data:  stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex ==  and stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex ==     0] and     1]
# W = 18289, p-value = 0.8958
# alternative hypothesis: true location shift is not equal to 0
# Wilcoxon rank sum test with continuity correction
# 
# data:  stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex ==  and stratified_complex_tab$mean_exp[stratified_complex_tab$is_complex ==     0] and     1]
# W = 18307, p-value = 0.9089
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(stratified_complex_tab$mean_exp,rowMeans(Mouse_complex_expression[,2:14]))
#Wilcoxon rank sum test with continuity correction

# data:  stratified_complex_tab$mean_exp and rowMeans(Mouse_complex_expression[, 2:14])
# W = 8960.5, p-value = 0.6843
# alternative hypothesis: true location shift is not equal to 0

cor_vector<-numeric(384)
for(i in 1:384){
  cor_vector[i]<-cor.test(as.numeric(as.character(stratified_complex_tab[i,c(3:15)])),Histone_expression_pattern,method="pearson")$estimate
}
median(Histone_expression_pattern)


complex_exp_pattern<-stratified_complex_tab[,c(3:15)]


mean(Histone_expression_pattern)
hist(Histone_expression_pattern)

selected_complex_expPattern<-Mouse_complex_expression[,2:14]

selected_complex_expPattern$mean_exp<-rowMeans(selected_complex_expPattern[,1:13])
wilcox.test(stratified_complex_tab$mean_exp,selected_complex_expPattern$mean_exp)

stratified_complex_tab$cor<-cor_vector
wilcox.test(stratified_complex_tab$cor[stratified_complex_tab$is_complex==0],
            stratified_complex_tab$cor[stratified_complex_tab$is_complex==1])
wilcox.test(histone_cor_selected,
            stratified_complex_tab$cor[stratified_complex_tab$is_complex==0])
histone_cor_selected<-numeric(45)
for(i in 1:45){
  histone_cor_selected[i]<-cor.test(as.numeric(as.character(Mouse_complex_expression[i,c(2:14)])),Histone_expression_pattern,method="pearson")$estimate
}
median(histone_cor_selected)
#Figure
selected_complexes<-as.character(Mouse_complex_expression$geneName)
tissue_specific_dat<-data.frame(c(as.character(stratified_complex_tab$genes),selected_complexes),
                                c(as.character(stratified_complex_tab$is_complex),rep("2",45)),
                                c(as.numeric(as.character(stratified_complex_tab$cor)),histone_cor_selected))
names(tissue_specific_dat)<-c("genes","linked","cor")
sum(tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==2]<1)
png(filename="Fig4g",width=10,height=6,units="in",res=600)
colorder<-c("0","1","2")
ggplot(tissue_specific_dat,aes(factor(linked),cor))+
  geom_boxplot(fill="gray")+
  xlab(label="")+
  ylab(label="Correlation in expression level\n with linker histone genes")+
  theme(axis.text.x=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.text.y=element_text(size=25,family="TT Arial",color="black"))+
  theme(axis.title.y=element_text(size=23,family="TT Arial",color="black",angle=90,vjust = 0.5,margin=margin(t=0,r=20,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks=element_line(size=1.5))+
  theme(axis.ticks.length=unit(5,"mm"))+
  scale_x_discrete(limits=colorder,labels=c("Unlinked","Linked", "Evolved linked\nin mouse lineage"))
dev.off()
wilcox.test(tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==2],
            tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==0],alternative = "greater")
wilcox.test(tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==1],
            tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==0],alternative = "greater")
wilcox.test(tissue_specific_dat$cor[tissue_specific_dat$linked==2],
            tissue_specific_dat$cor[tissue_specific_dat$linked==0],alternative = "less")
#one-side wilcox test:
# Wilcoxon rank sum test with continuity correction
# 
# data:  tissue_specific_dat$cor[tissue_specific_dat$linked == 2] and tissue_specific_dat$cor[tissue_specific_dat$linked == 0]
# W = 2993.5, p-value = 0.0006795
# alternative hypothesis: true location shift is less than 0
wilcox.test(tissue_specific_dat$cor[tissue_specific_dat$linked==2],
            tissue_specific_dat$cor[tissue_specific_dat$linked==0],alternative = "less")
#one-side Wilcoxon rank sum test with continuity correction

# data:  tissue_specific_dat$cor[tissue_specific_dat$linked == 1] and tissue_specific_dat$cor[tissue_specific_dat$linked == 0]
# W = 15979, p-value = 0.01206
# alternative hypothesis: true location shift is less than 0

#linked vs unlinked:0.02339
#evolve vs unlinked:0.009178
m1<-matrix(c(192,192-86,192,192-68),nrow=2,byrow=TRUE)
fisher.test(m1)
m2<-matrix(c(45,45-9,192,192-86),nrow=2,byrow=TRUE)
median(tissue_specific_dat$low_histone_enrichment[tissue_specific_dat$linked==1])
fisher.test(m2)
wilcox.test(tissue_specific_dat$cor[tissue_specific_dat$linked==2],
            tissue_specific_dat$cor[tissue_specific_dat$linked==0])
# #Wilcoxon rank sum test with continuity correction
# 
# data:  tissue_specific_dat$cor[tissue_specific_dat$linked == 2] and tissue_specific_dat$cor[tissue_specific_dat$linked == 0]
# W = 2958.5, p-value = 0.00101
# alternative hypothesis: true location shift is not equal to 0
wilcox.test(tissue_specific_dat$cor[tissue_specific_dat$linked==1],
            tissue_specific_dat$cor[tissue_specific_dat$linked==0])
# Wilcoxon rank sum test with continuity correction
# 
# data:  tissue_specific_dat$cor[tissue_specific_dat$linked == 1] and tissue_specific_dat$cor[tissue_specific_dat$linked == 0]
# W = 15896, p-value = 0.01973
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(tissue_specific_dat$cor[tissue_specific_dat$linked==1],
            tissue_specific_dat$cor[tissue_specific_dat$linked==2])
#Wilcoxon rank sum test with continuity correction

# data:  tissue_specific_dat$cor[tissue_specific_dat$linked == 1] and tissue_specific_dat$cor[tissue_specific_dat$linked == 2]
# W = 5181.5, p-value = 0.03753
# alternative hypothesis: true location shift is not equal to 0