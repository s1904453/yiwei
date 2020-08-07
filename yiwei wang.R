library(tibble)
library(dplyr)
library(tidyr)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(amap)
library(rio)
library(ggplot2)
library(ggbiplot)
library(devtools)
library("factoextra")
library(glmnet)
install_github("vqv/ggbiplot")
gene_expression <- import("~/cluster/CW-kallisto-abundance-foldchange-long-bygene.csv")
k_mer <-import("~/cluster/H99_allorfs_promoter500nt_5mercounts 9.txt")
k_mer = as_tibble(k_mer)
gene_expression = as_tibble(gene_expression)
# we only consider coding genes
gene_coding <- subset(gene_expression, grepl("^CNAG_0", gene_expression$Gene))
# remove the records containing low expression genes

# first find the gene labels whose ESTCounts and TPM is 0 in some conditions, 
# and these genes are considered as low expression genes
Gene_labels <- gene_coding %>%
  filter(gene_coding$EstCounts ==0) %>%
  select(Gene) %>%
  distinct()
# get the vector which contains lables of low expression genes
Gene_labels <- as.vector(unlist(Gene_labels[,]))
# create a new column to identify whether a gene is low expressed
gene_coding$low_nonlow <- ifelse(gene_coding$Gene %in% Gene_labels,
                                 "low", "nonlow")
# filter out low expression genes
gene_coding <- gene_coding %>%
  filter(gene_coding$low_nonlow == "nonlow") %>%
  select(-low_nonlow)

# similarity between different conditions
code_tpm <- gene_coding %>%
  select(-EstCounts,-FoldChange)
code_tpm$TPM <- log(code_tpm$TPM+1, base = 2)

condition_differ <- spread(code_tpm, Gene, TPM)
condition_differ <- as.data.frame(condition_differ)
condition_differ <- condition_differ %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Code')
# scale value of log(tpm+1) in each sample
condition_differ_scaled <- scale(condition_differ)
# here we compute the correlation based distance
condition.dist.cor <- get_dist(condition_differ_scaled, method = "pearson")
as.matrix(condition.dist.cor)
# visualize distance matrices
fviz_dist(condition.dist.cor)

# The color level is proportional to value of similarity between sample conditions. Pure red if distance between samples is 0.
# and pure blue if distance is 2, i.e., a large negative correlation between two samples. And samples belonging to the
# same cluster are displayed in consecutive order.thus we can find that our sample conditions can be divided into there groups roughly.

# The expression level of all genes in YPD rich media in time point 2,3,4 in both low and high temperature is similar which can be
# denoted as the first smaple cluster. Also, sample for extra time point 9 in YPD media can also be included in cluster 1 and the 
# gene expression level in this sample is very similar to sample of YPD media in time point 2,3 and 4. thus, we will omit extra sample
# YM9A and YM9B.

# Also, sample conditions in time point 1 for both YPD and RPMI media and both low and high tempature which includes 8 different
# samples are clustered together. It indicates that there is not a big difference between cells grown in different media
# and tempature for only 10 minutes.

# The third contains RPMI media with both low and high temperature in time points 2,3 and 4. Also, the similarity between conditions
# under the same level of tempature is higher. And we can find that in "RC" condition,time point 1 and 0 are similar and there is a 
# between time point 1 and 2,3,4. Thus, we will omit "YC0A" and "YC0B" because the first ten minutes would not cause much 
# difference on cells.

# Also, it is clear that the distance between replicate A and B is always smaller than distance between different conditions.
# It can be concluded from the distance matrix plot that most of replicates samples are adjacent and the red color between
# replicates A and B are deeper than beween different conditions even through these conditions are in the same cluster.


# remove sample "RC0A", "RC0B","YM9A" and "YM9B"
condition_differ <- condition_differ[row.names(condition_differ) !=
                                       c("YC0A","YC0B","YM9A","YM9B"),]
condition_differ_scaled <- scale(condition_differ)
sample_dist <- get_dist(condition_differ_scaled, method = "pearson")
sample_hc <- hclust(d = sample_dist, method = "ward.D")
fviz_dend(sample_hc)

# create a new data frame gene_cluster,
# whose columns are samples(Code) and rows are genes, 
# and values are from fold change.
code_foldchange <- gene_coding %>%
  select(-EstCounts,-TPM)
code_foldchange$FoldChange <- log(code_foldchange$FoldChange,base = 2)
gene_cluster <- spread(code_foldchange, Code, FoldChange)
gene_cluster <- gene_cluster %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')

# scatter plots for RH1A and RH1B, RH2A and RH2B, RH3A and RH3B and RH4A and RH4B
p1 <- ggplot(gene_cluster, aes(x=gene_cluster$RH1A, y=gene_cluster$RH1B,xlab("Dose (mg)"))+ geom_point()+
               geom_abline(intercept=0, slope=1,color = "red"))
p1 <-p1+labs(x ="RH1A", y = "RH1B")
p2 <- ggplot(gene_cluster, aes(x=gene_cluster$RH2A, y=gene_cluster$RH2B)) + geom_point()+
  geom_abline(intercept=0, slope=1,color = "red")
p2 <- p2+labs(x ="RH2A", y = "RH2B")
p3 <- ggplot(gene_cluster, aes(x=gene_cluster$RH3A, y=gene_cluster$RH3B)) + geom_point()+
  geom_abline(intercept=0, slope=1,color = "red")
p3 <- p3+labs(x ="RH3A", y = "RH3B")
p4 <- ggplot(gene_cluster, aes(x=gene_cluster$RH4A, y=gene_cluster$RH4B)) + geom_point()+
  geom_abline(intercept=0, slope=1,color = "red")
p4 <- p4+labs(x ="RH4A", y = "RH4B")
ggarrange(p1, p2, p3,p4,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
# remove the bias between RH1A and RH1B, RH2A and RH2B, RH3A and RH3B and RH4A and RH4B
ss1 = smooth.spline(gene_cluster$RH1A,
                    gene_cluster$RH1B-gene_cluster$RH1A,
                    cv=TRUE)
bias1 <- rep(1,nrow(gene_cluster))
for (i in 1:nrow(gene_cluster)) {
  bias1[i] <-as.numeric(predict(ss1,x=gene_cluster$RH1A[i])$y)
}
gene_cluster$RH1B <- gene_cluster$RH1B-bias1
ss2 = smooth.spline(gene_cluster$RH2A,
                    gene_cluster$RH2B-gene_cluster$RH2A,
                    cv=TRUE)
bias2 <- rep(1,nrow(gene_cluster))
for (i in 1:nrow(gene_cluster)) {
  bias2[i] <-as.numeric(predict(ss1,x=gene_cluster$RH2A[i])$y)
}
gene_cluster$RH2B <- gene_cluster$RH2B-bias2
ss3 = smooth.spline(gene_cluster$RH3A,
                    gene_cluster$RH3B-gene_cluster$RH3A,
                    cv=TRUE)
bias3 <- rep(1,nrow(gene_cluster))
for (i in 1:nrow(gene_cluster)) {
  bias3[i] <-as.numeric(predict(ss3,x=gene_cluster$RH3A[i])$y)
}
gene_cluster$RH3B <- gene_cluster$RH3B-bias3
ss4 = smooth.spline(gene_cluster$RH4A,
                    gene_cluster$RH4B-gene_cluster$RH4A,
                    cv=TRUE)
bias4 <- rep(1,nrow(gene_cluster))
for (i in 1:nrow(gene_cluster)) {
  bias4[i] <-as.numeric(predict(ss4,x=gene_cluster$RH4A[i])$y)
}
gene_cluster$RH4B <- gene_cluster$RH4B-bias4
gene_cluster$RC1 <- (gene_cluster$RC1A+gene_cluster$RC1B)/2
gene_cluster$RC2 <- (gene_cluster$RC2A+gene_cluster$RC2B)/2
gene_cluster$RC3 <- (gene_cluster$RC3A+gene_cluster$RC2B)/2
gene_cluster$RC4 <- (gene_cluster$RC4A+gene_cluster$RC4B)/2
gene_cluster$RH1 <- (gene_cluster$RH1A+gene_cluster$RH1B)/2
gene_cluster$RH2 <- (gene_cluster$RH2A+gene_cluster$RH2B)/2
gene_cluster$RH3 <- (gene_cluster$RH3A+gene_cluster$RH3B)/2
gene_cluster$RH4 <- (gene_cluster$RH4A+gene_cluster$RH4B)/2
gene_cluster$YC1 <- (gene_cluster$YC1A+gene_cluster$YC1B)/2
gene_cluster$YC2 <- (gene_cluster$YC2A+gene_cluster$YC2B)/2
gene_cluster$YC3 <- (gene_cluster$YC3A+gene_cluster$YC3B)/2
gene_cluster$YC4 <- (gene_cluster$YC4A+gene_cluster$YC4B)/2
gene_cluster$YH1 <- (gene_cluster$YH1A+gene_cluster$YH1B)/2
gene_cluster$YH2 <- (gene_cluster$YH2A+gene_cluster$YH2B)/2
gene_cluster$YH3 <- (gene_cluster$YH3A+gene_cluster$YH3B)/2
gene_cluster$YH4 <- (gene_cluster$YH4A+gene_cluster$YH4B)/2
gene_cluster <- gene_cluster %>%
  select(-RC1A,-RC1B,-RC2A,-RC2B,-RC3A,-RC3B,-RC4A,-RC4B,
         -RH1A,-RH1B,-RH2A,-RH2B,-RH3A,-RH3B,-RH4A,-RH4B,-YC0A,-YC0B,
         -YC1A,-YC1B,-YC2A,-YC2B,-YC3A,-YC3B,-YC4A,-YC4B,
         -YH1A,-YH1B,-YH2A,-YH2B,-YH3A,-YH3B,,-YH4A,-YH4B,-YM9A,-YM9B)
gene_cluster_scaled <- scale(gene_cluster)
# get the optimal number of clusters
fviz_nbclust(gene_cluster_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# k means clustering
set.seed(3)
final <- kmeans(gene_cluster_scaled, 4, nstart = 25)
fviz_cluster(final, data = gene_cluster_scaled)
# assign cluster number for each gene
gene_cluster <-cbind(gene_cluster, cluster = final$cluster)
# compute the mean of foldchange each condition in each cluster 
cluster_result <- aggregate(2^(gene_cluster),by=list(cluster=final$cluster),mean)
cluster1 <- c()
cluster2 <- c()
cluster3 <- c()
cluster4 <- c()
for (i in 2:17) {
  cluster1 <- append(cluster1,cluster_result[[i]][1])
  cluster2 <- append(cluster2,cluster_result[[i]][2])
  cluster3 <- append(cluster3,cluster_result[[i]][3])
  cluster4 <- append(cluster4,cluster_result[[i]][4])
}
# rearrange dataframe cluster_result
clustered_mean <- data.frame(cluster1,cluster2,cluster3,cluster4)
clustered_mean$Sample <- c("RC1","RC2","RC3","RC4","RH1","RH2","RH3","RH4",
                           "YC1","YC2","YC3","YC4","YH1","YH2","YH3","YH4")
clustered_mean <- gather(clustered_mean,'cluster1','cluster2','cluster3','cluster4',
                         key = "cluster",value = "foldchange")
# divide clustered_mean into four dataframes in terms of conditions
clustered_RC <- clustered_mean %>%
  filter(Sample==c("RC1","RC2","RC3","RC4"))
clustered_RH <- clustered_mean %>%
  filter(Sample==c("RH1","RH2","RH3","RH4"))
clustered_YC <- clustered_mean %>%
  filter(Sample==c("YC1","YC2","YC3","YC4"))
clustered_YH <- clustered_mean %>%
  filter(Sample==c("YH1","YH2","YH3","YH4"))
# plot mean foldchange value in each cluster in each condition
s1 <- ggplot(data=clustered_RC, aes(x=Sample, y=foldchange, group=cluster,colour=cluster)) +
  geom_line()+geom_point()
s1 <-s1+labs(x ="Sample RC", y = "mean foldchange")
s2 <- ggplot(data=clustered_RH, aes(x=Sample, y=foldchange, group=cluster,colour=cluster)) +
  geom_line()+geom_point()
s2 <-s2+labs(x ="Sample RH", y = "mean foldchange")
s3 <- ggplot(data=clustered_YC, aes(x=Sample, y=foldchange, group=cluster,colour=cluster)) +
  geom_line()+geom_point()
s3 <-s3+labs(x ="Sample YC", y = "mean foldchange")
s4 <- ggplot(data=clustered_YH, aes(x=Sample, y=foldchange, group=cluster,colour=cluster)) +
  geom_line()+geom_point()
s4 <-s4+labs(x ="Sample YH", y = "mean foldchange")
ggarrange(s1, s2, s3, s4,ncol = 2, nrow = 2)

# filter out low expression genes
# create a new column to identify whether a gene is low expressed
k_mer$low_nonlow <- ifelse(k_mer$Gene %in% Gene_labels,
                           "low", "nonlow")
k_mer <- k_mer%>%
  filter(k_mer$low_nonlow == "nonlow") %>%
  select(-low_nonlow)
k_mer <- k_mer %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')
# filter out genes in gene_cluster that do not have kmers 
gene_cluster <- gene_cluster[rownames(gene_cluster) %in% rownames(k_mer),]
k_mer <- k_mer[rownames(k_mer) %in% rownames(gene_cluster),]
motif_label <- colnames(k_mer)

letter <- rep(0,5)
reverse_letter <- rep(0,5)
j <-1
# create a function get_reverse_motif
# for each motif, create its biologically equivalent label
get_reverse_motif <- function(label) {
  for (j in 1:5) {
    if (label[j] == "A") {
      reverse_letter[j] <- "T"
    }
    else if (label[j] == "T") {
      reverse_letter[j] <- "A"
    }
    else if (label[j] == "C") {
      reverse_letter[j] <- "G"
    }
    else if (label[j] == "G") {
      reverse_letter[j] <- "C"
    }}
  k <- 2
  inter <- 0
  for (i in 4:5) {
    inter <- reverse_letter[i]
    reverse_letter[i] <- reverse_letter[k]
    reverse_letter[k] <- inter
    k <- k-1
  }
  return(reverse_letter)
}
# check get_reverse_motif function
get_reverse_motif(unlist(strsplit(motif_label[1], split = "")))

# find the biologically equivalent motif for all of motif in k-mer dataset
equal_index <- c()
for (i in 1:length(motif_label)) {
  j <- i+1
  while (j<=length(motif_label)){
    if (get_reverse_motif(unlist(strsplit(motif_label[i], split = "")))[1]==
        unlist(strsplit(motif_label[j], split = ""))[1]&
        get_reverse_motif(unlist(strsplit(motif_label[i], split = "")))[2]==
        unlist(strsplit(motif_label[j], split = ""))[2] &
        get_reverse_motif(unlist(strsplit(motif_label[i], split = "")))[3]==
        unlist(strsplit(motif_label[j], split = ""))[3] &
        get_reverse_motif(unlist(strsplit(motif_label[i], split = "")))[4]==
        unlist(strsplit(motif_label[j], split = ""))[4] &
        get_reverse_motif(unlist(strsplit(motif_label[i], split = "")))[5]==
        unlist(strsplit(motif_label[j], split = ""))[5]
    ) {
      # replace number of a 5-mer with the mean number of it and number of its biologically equivalent motif
      k_mer[,i]=(k_mer[,i]+k_mer[,j])/2
      # record the higher index of each pair of biologically equivalent motifs
      equal_index <- append(equal_index,j)
    }
    j <- j+1
  }
}
# remove the higher index of each pair of biologically equivalent motifs
k_mer<- k_mer[,-equal_index]
# there are 512 motifs after combining
ncol(k_mer)
# there are 34 different motifs containing "AAA"
sub_AAA <- c()
for (i in 1:512) {
  if ( grepl("AAA", colnames(k_mer)[i], fixed = TRUE)){
    sub_AAA <- append(sub_AAA ,colnames(k_mer)[i])
  }
}
#there are 26 different motifs containing "CCC"
sub_CCC <- c()
for (i in 1:512) {
  if ( grepl("CCC", colnames(k_mer)[i], fixed = TRUE)){
    sub_CCC <- append(sub_CCC,colnames(k_mer)[i])
  }
}


X1 <- as.matrix(k_mer)
y1 <- as.vector(gene_cluster$cluster)
lasso1<- cv.glmnet(X1,y1,family = "multinomial")
fit_lasso1 <- glmnet(X1,y1,lambda = lasso1$lambda.min,family = "multinomial")
# optimal lambda value is 0.007090064
coefficients <- coef(fit_lasso1)
# get motifs with larger than absolute value of 0.02 
# remove the intercept term
sig_kmer1pos <- colnames(X1)[which(unlist(coefficients[[1]][-1]>=0.02))]
sig_kmer1neg <- colnames(X)[which(unlist(coefficients[[1]][-1]<=-0.02))]
sig_kmer2pos <- colnames(X)[which(unlist(coefficients[[2]][-1]>=0.02))]
sig_kmer2neg <- colnames(X)[which(unlist(coefficients[[2]][-1]<=-0.02))]
sig_kmer3pos <- colnames(X)[which(unlist(coefficients[[3]][-1]>=0.02))]
sig_kmer3neg <- colnames(X)[which(unlist(coefficients[[3]][-1]<=-0.02))]
sig_kmer4pos <- colnames(X)[which(unlist(coefficients[[4]][-1]>=0.02))]
sig_kmer4neg <- colnames(X)[which(unlist(coefficients[[4]][-1]<=-0.02))]
# output mitif table
five_mer11 <- t(data.frame(sig_kmer1pos))
five_mer12 <- t(data.frame(sig_kmer1neg))
xtable(five_mer11)
xtable(five_mer12)
five_mer21 <- t(data.frame(sig_kmer2pos))
xtable(five_mer21)
five_mer22 <- t(data.frame(sig_kmer2neg))
xtable(five_mer22)
five_mer31 <- t(data.frame(sig_kmer3pos))
five_mer32 <- t(data.frame(sig_kmer3neg))
xtable(five_mer31)
xtable(five_mer32)
five_mer41 <- t(data.frame(sig_kmer4pos))
five_mer42 <- t(data.frame(sig_kmer4neg))
xtable(five_mer41)
xtable(five_mer42)

# boxplot of coefficients of important motifs
cluster1_5 <- unlist(coefficients[[1]])[unlist(abs(coefficients[[1]])>=0.02)][-1]
cluster2_5 <- unlist(coefficients[[2]])[unlist(abs(coefficients[[2]])>=0.02)][-1]
cluster3_5 <- unlist(coefficients[[3]])[unlist(abs(coefficients[[3]])>=0.02)][-1]
cluster4_5 <- unlist(coefficients[[4]])[unlist(abs(coefficients[[4]])>=0.02)][-1]
box1_5<- data.frame(rep("cluster1",length(cluster1_5)),c(cluster1_5))
colnames(box1_5) <- c("cluster1","coeffi")
box2_5<- data.frame(rep("cluster2",length(cluster2_5)),c(cluster2_5))
colnames(box2_5) <- c("cluster2","coeffi")
box3_5 <- data.frame(rep("cluster3",length(cluster3_5)),c(cluster3_5))
colnames(box3_5) <- c("cluster3","coeffi")
box4_5<- data.frame(rep("cluster4",length(cluster4_5)),c(cluster4_5))
colnames(box4_5) <- c("cluster4","coeffi")
library(ggplot2)
# Basic box plot
p1 <- ggplot(box1_5, aes(x=cluster1,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p2 <- ggplot(box2_5, aes(x=cluster2,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p3 <- ggplot(box3_5, aes(x=cluster3,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p4 <- ggplot(box4_5, aes(x=cluster4,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
ggarrange(p1, p2, p3, p4,ncol = 2, nrow = 2)

# get the mean coeffi of positive and negative motifs in each cluster
box1_5pos <- box1_5 %>%
  filter(box1_5$coeffi>0)
2*mean(box1_5pos[,"coeffi"]) #0.275167

box2_5pos <- box2_5 %>%
  filter(box2_5$coeffi>0)
2*mean(box2_5pos[,"coeffi"]) #0.1258698

box3_5pos <- box3_5 %>%
  filter(box3_5$coeffi>0)
2*mean(box3_5pos[,"coeffi"]) #0.1321469

box4_5pos <- box4_5 %>%
  filter(box4_5$coeffi>0)
2*mean(box4_5pos[,"coeffi"]) # 0.1661267

box1_5neg <- box1_5 %>%
  filter(box1_5$coeffi<0)
2*mean(box1_5neg[,"coeffi"]) #-0.1956674

box2_5neg <- box2_5 %>%
  filter(box2_5$coeffi<0)
2*mean(box2_5neg[,"coeffi"]) #-0.1269633

box3_5neg <- box3_5 %>%
  filter(box3_5$coeffi<0)
2*mean(box3_5neg[,"coeffi"]) #-0.09466986

box4_5neg <- box4_5 %>%
  filter(box4_5$coeffi<0)
2*mean(box4_5neg[,"coeffi"]) #-0.1242249

# distribution of non zero coefficients
c1 <- unlist(coefficients[[1]])[unlist(coefficients[[1]])!=0]
c1 <- data.frame(c1[-1])
p1 <- ggplot(c1,aes(x=c1..1.)) +
  geom_density()+labs(x ="non-zero coefficients for cluster 1")
c2 <- unlist(coefficients[[1]])[unlist(coefficients[[2]])!=0]
c2 <- data.frame(c2[-1])
p2 <- ggplot(c2,aes(x=c2..1.)) +
  geom_density()+labs(x ="non-zero coefficients for cluster 2")
c3 <- unlist(coefficients[[1]])[unlist(coefficients[[3]])!=0]
c3 <- data.frame(c3[-1])
p3 <- ggplot(c3,aes(x=c3..1.)) +
  geom_density()+labs(x ="non-zero coefficients for cluster 3")
c4 <- unlist(coefficients[[1]])[unlist(coefficients[[4]])!=0]
c4 <- data.frame(c4[-1])
p4 <- ggplot(c4,aes(x=c4..1.)) +
  geom_density()+labs(x ="non-zero coefficients for cluster 4")
library(ggpubr)
ggarrange(p1, p2, p3, p4,ncol = 2, nrow = 2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
# r setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE, cache=TRUE,
                      cache.path = "cache/CountKmers-")
library(Biostrings)
# r kmer_functions
make_kmer_list <- function(k=2L, alphabet = c("A","C","G","T")) {
  list(alphabet) %>%
    base::rep(times = k) %>%
    purrr::set_names(LETTERS[1:k]) %>%
    base::expand.grid() %>%
    tidyr::unite("kmer",LETTERS[1:k],sep="") %>%
    dplyr::pull(kmer)
}
count_kmer_one <- function(string,kmer_list,k=2L) {
  stopifnot( nchar(kmer_list[1]) == k )
  kCount <- rep(0L,length(kmer_list)) %>% purrr::set_names(kmer_list)
  kminus1 <- k - 1L
  for(i in 1:( nchar(string) - kminus1 ) ) {
    kmer <- stringr::str_sub(string, start=i, end = i + kminus1)
    if(kmer %in% kmer_list) {
      kCount[kmer] <- kCount[kmer] + 1L
    }
  }
  return(kCount)
}
count_kmer_tibble <- function(stringset,kmer_list,k=2L) {
  lapply(stringset,count_kmer_one,kmer_list=kmer_list,k=k) %>%
    lapply(as.list) %>%
    dplyr::bind_rows()
}
library(here)
# r count_6mers_H99promoters
promoterseqs <- Biostrings::readDNAStringSet("/Users/wangyiwei/Desktop/H99_allorfs_p500 (4).fasta")
promoterseqids <- names(promoterseqs) %>%
  stringr::str_extract(pattern="\\w+")
all_6mers <- make_kmer_list(k=6L)
promoter_6mer_countsonly <- promoterseqs %>%
  as.character() %>%
  count_kmer_tibble(kmer_list = all_6mers, k=6L)
promoter_6mer_counts <- bind_cols(tibble(Gene=promoterseqids),
                                  promoter_6mer_countsonly)
promoter_6mer_counts
promoter_6mer_counts <- promoter_6mer_counts%>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')
motif_label2 <- colnames(promoter_6mer_counts)

letter <- rep(0,6)
reverse_letter <- rep(0,6)
j <-1
# create a function get_reverse_motif2
# for each motif, create its biologically equivalent label
get_reverse_motif2 <- function(label) {
  for (j in 1:6) {
    if (label[j] == "A") {
      reverse_letter[j] <- "T"
    }
    else if (label[j] == "T") {
      reverse_letter[j] <- "A"
    }
    else if (label[j] == "C") {
      reverse_letter[j] <- "G"
    }
    else if (label[j] == "G") {
      reverse_letter[j] <- "C"
    }}
  k <- 3
  inter <- 0
  for (i in 4:6) {
    inter <- reverse_letter[i]
    reverse_letter[i] <- reverse_letter[k]
    reverse_letter[k] <- inter
    k <- k-1
  }
  return(reverse_letter)
}
# check get_reverse_motif function
get_reverse_motif2(unlist(strsplit(motif_label2[1], split = "")))

# find the biologically equivalent motif for all of motif in k-mer dataset
equal_index <- c()
for (i in 1:length(motif_label2)) {
  j <- i+1
  while (j<=length(motif_label2)){
    if (get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[1]==
        unlist(strsplit(motif_label2[j], split = ""))[1]&
        get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[2]==
        unlist(strsplit(motif_label2[j], split = ""))[2] &
        get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[3]==
        unlist(strsplit(motif_label2[j], split = ""))[3] &
        get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[4]==
        unlist(strsplit(motif_label2[j], split = ""))[4] &
        get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[5]==
        unlist(strsplit(motif_label2[j], split = ""))[5] &
        get_reverse_motif2(unlist(strsplit(motif_label2[i], split = "")))[6]==
        unlist(strsplit(motif_label2[j], split = ""))[6]
    ) {
      # replace number of a 5-mer with the mean number of it and number of its biologically equivalent motif
      promoter_6mer_counts[,1]=(promoter_6mer_counts[,i]+promoter_6mer_counts[,j])/2
      # record the higher index of each pair of biologically equivalent motifs
      equal_index <- append(equal_index,j)
    }
    j <- j+1
  }
}
# remove the higher index of each pair of biologically equivalent motifs
promoter_6mer_counts<- promoter_6mer_counts[,-equal_index]
# there are 2080 motifs after combining
ncol(promoter_6mer_counts)


# filter out genes in gene_cluster that do not have kmers 
gene_cluster <- gene_cluster[rownames(gene_cluster) %in% rownames(promoter_6mer_counts),]
promoter_6mer_counts <- promoter_6mer_counts[rownames(promoter_6mer_counts) %in% rownames(gene_cluster),]
nrow(promoter_6mer_counts)

# lasso regression with 6 mers
X2 <- as.matrix(promoter_6mer_counts)
y2 <- as.vector(gene_cluster$cluster)
lasso2 <- cv.glmnet(X2,y2,family = "multinomial")
fit_lasso2 <- glmnet(X2,y2,lambda = lasso2$lambda.min,family = "multinomial")
# optimal lambda value is 0.01134407
coefficients2 <- coef(fit_lasso2)
# get motifs with larger than absolute value of 0.02 
# remove the intercept term
sig_6mer1pos <- colnames(X2)[which(unlist(coefficients2[[1]][-1]>=0.02))]
sig_6mer1neg <- colnames(X2)[which(unlist(coefficients2[[1]][-1]<=-0.02))]
sig_6mer2pos <- colnames(X2)[which(unlist(coefficients2[[2]][-1]>=0.02))]
sig_6mer2neg <- colnames(X2)[which(unlist(coefficients2[[2]][-1]<=-0.02))]
sig_6mer3pos <- colnames(X2)[which(unlist(coefficients2[[3]][-1]>=0.02))]
sig_6mer3neg <- colnames(X2)[which(unlist(coefficients2[[3]][-1]<=-0.02))]
sig_6mer4pos <- colnames(X2)[which(unlist(coefficients2[[4]][-1]>=0.02))]
sig_6mer4neg <- colnames(X2)[which(unlist(coefficients2[[4]][-1]<=-0.02))]
cluster1_6 <- unlist(coefficients2[[1]])[unlist(abs(coefficients2[[1]])>=0.02)][-1]
cluster2_6 <- unlist(coefficients2[[2]])[unlist(abs(coefficients2[[2]])>=0.02)][-1]
cluster3_6 <- unlist(coefficients2[[3]])[unlist(abs(coefficients2[[3]])>=0.02)][-1]
cluster4_6 <- unlist(coefficients2[[4]])[unlist(abs(coefficients2[[4]])>=0.02)][-1]
box1_6 <- data.frame(rep("cluster1",length(cluster1_6)),c(cluster1_6))
colnames(box1_6) <- c("cluster1","coeffi")
box2_6 <- data.frame(rep("cluster2",length(cluster2_6)),c(cluster2_6))
colnames(box2_6) <- c("cluster2","coeffi")
box3_6 <- data.frame(rep("cluster3",length(cluster3_6)),c(cluster3_6))
colnames(box3_6) <- c("cluster3","coeffi")
box4_6 <- data.frame(rep("cluster4",length(cluster4_6)),c(cluster4_6))
colnames(box4_6) <- c("cluster4","coeffi")
library(ggplot2)
# Basic box plot
p1 <- ggplot(box1_6, aes(x=cluster1,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p2 <- ggplot(box2_6, aes(x=cluster2,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p3 <- ggplot(box3_6, aes(x=cluster3,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
p4 <- ggplot(box4_6, aes(x=cluster4,y = coeffi)) + 
  geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))
library(ggpubr)
ggarrange(p1, p2, p3, p4,ncol = 2, nrow = 2)

# output mitif table
six_mer11 <- t(data.frame(sig_6mer1pos))
six_mer12 <- t(data.frame(sig_6mer1neg))
xtable(six_mer11)
xtable(six_mer12)
six_mer21 <- t(data.frame(sig_6mer2pos))
six_mer22 <- t(data.frame(sig_6mer2neg))
xtable(six_mer21)
xtable(six_mer22)
six_mer31 <- t(data.frame(sig_6mer3pos))
six_mer32 <- t(data.frame(sig_6mer3neg))
xtable(six_mer31)
xtable(six_mer32)
six_mer41 <- t(data.frame(sig_6mer4pos))
six_mer42 <- t(data.frame(sig_6mer4neg))
xtable(six_mer41)
xtable(six_mer42)































