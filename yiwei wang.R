library(stringr)
Gene <- rep(1,nrow(gene_cluster))
for (i in as.numeric(rownames(gene_cluster))) {
  Gene[i] <- str_c(cbind(gene_cluster$Gene[i],"B"), collapse = "_")
}

RC1 <- gene_cluster$RC1B
RC2 <- gene_cluster$RC2B
RC3 <- gene_cluster$RC3B
RC4 <- gene_cluster$RC4B
RH1 <- gene_cluster$RH1B
RH2 <- gene_cluster$RH2B
RH3 <- gene_cluster$RH3B
RH4 <- gene_cluster$RH4B
YC1 <- gene_cluster$YC1B
YC2 <- gene_cluster$YC2B
YC3 <- gene_cluster$YC3B
YC4 <- gene_cluster$YC4B
YH1 <- gene_cluster$YH1B
YH2 <- gene_cluster$YH2B
YH3 <- gene_cluster$YH3B
YH4 <- gene_cluster$YH4B

gene_cluster_B <- data.frame(Gene, RC1,RC2,RC3,RC4,RH1,RH2,RH3,RH4,
                             YC1,YC2,YC3,YC4,YH1,YH2,YH3,YH4)
gene_cluster_A <- gene_cluster %>%
  select(-RC1B,-RC2B,-RC3B,-RC4B,-RH1B,-RH2B,-RH3B,-RH4B,
         -YC0A,-YC0B,-YC1B,-YC2B,-YC3B,-YC4B,-YH1B,
         -YH2B,-YH3B,-YH4B,-YM9A,-YM9B)
library(data.table)
setnames(gene_cluster_A, 
         old = c("RC1A","RC2A","RC3A","RC4A","RH1A","RH2A","RH3A",
                 "RH4A","YC1A","YC2A","YC3A","YC4A","YH1A","YH2A",
                 "YH3A","YH4A"), 
         new = c("RC1","RC2","RC3","RC4","RH1","RH2","RH3",
                 "RH4","YC1","YC2","YC3","YC4","YH1","YH2",
                 "YH3","YH4"))
gene_coreplicate <- rbind(gene_cluster_A,gene_cluster_B)
gene_coreplicate <- gene_coreplicate %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')

gene_coreplicate_scaled <- scale(gene_coreplicate)
# hierarcical clustering
gene_dist <- dist(gene_coreplicate_scaled, method = "euclidean")
gene_hc <- hclust(d = gene_dist, method = "ward.D")
fviz_dend(gene_hc)
heatmap(gene_coreplicate_scaled, scale = "none")

fviz_nbclust(gene_coreplicate_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
fviz_nbclust(gene_coreplicate_scaled, kmeans, method = "wss")+
  labs(subtitle = "Silhouette method")

coreplicate_final <- kmeans(gene_coreplicate_scaled, 4, nstart = 25)
fviz_cluster(coreplicate_final, data = gene_coreplicate_scaled)
gene_coreplicate <-cbind(gene_coreplicate, cluster = coreplicate_final$cluster)
for (i in 1:nrow(gene_cluster_A)) {
  gene_coreplicate$noise[i] <- ifelse(gene_coreplicate$cluster[i]==gene_coreplicate$cluster[i+4781],
                                      "not_noise","noise")
  gene_coreplicate$noise[i+4781] <- ifelse(gene_coreplicate$cluster[i]==gene_coreplicate$cluster[i+4781],
                                           "not_noise","noise")
}
noise_num <-0
for (i in 1:nrow(gene_cluster_A)) {
  if (gene_coreplicate$noise[i]=="noise") {
    noise_num <- noise_num+1
  }
}
noise_index <- 1
noise_gene <- rep(0,noise_num)
for (i in 1:nrow(gene_cluster_A)) {
  if (gene_coreplicate$noise[i]=="noise") {
    noise_gene[noise_index] <- rownames(gene_coreplicate)[i]
    noise_index <- noise_index+1
  }
}
gene_cluster <- gene_cluster %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')
RC1A <- gene_cluster %>%
  select(RC1A)
RC1A_final <- kmeans(scale(RC1A), 4, nstart = 25)
RC1B <- gene_cluster %>%
  select(RC1B)
RC1B_final <- kmeans(scale(RC1B), 4, nstart = 25)
RC1 <-cbind(RC1A, clusterA = RC1A_final$cluster,clusterB=RC1B_final$cluster)

heatmap(gene_coreplicate_scaled, scale = "none")
# hierarcical clustering
gene_cluster <- gene_cluster %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Gene')
gene_cluster_scaled <- scale(gene_cluster)
gene_dist <- dist(gene_cluster_scaled, method = "euclidean")
gene_hc <- hclust(d = gene_dist, method = "ward.D")
fviz_dend(gene_hc)

# Cut in 4 groups and color by groups
fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

fviz_nbclust(gene_cluster_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(gene_cluster_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

final <- kmeans(gene_cluster_scaled, 4, nstart = 25)
fviz_cluster(final, data = gene_cluster_scaled)
fviz_nbclust(sample_cluster_sclaed, hcut, method = "wss")
sample_final <- kmeans(sample_cluster_sclaed, 5, nstart = 25)
sample_cluster <-cbind(sample_cluster, cluster = sample_final$cluster)



sample_cluster$cluster <- as.factor(sample_cluster$cluster)
ggparcoord(sample_cluster, columns = 1:4781, groupColumn = 4782)+
  scale_color_viridis(discrete=TRUE) 
fviz_cluster(sample_final, data = sample_cluster[1:4781])
sample_final$cluster

condition_differ_scaled <- scale(condition_differ)
fviz_nbclust(condition_differ_scaled, kmeans, method = "wss",iter.max = 15)

sample_dist <- dist(condition_differ_scaled, method = "euclidean")
sample_hc <- hclust(d = sample_dist, method = "ward.D")
fviz_dend(sample_hc)

# Cut in 4 groups and color by groups
fviz_dend(sample_hc, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
library(Rmisc)
2^(CI(cluster1$RH3, ci=0.95))
2^(final$centers)

# Display parallel coordinates plots, one for each cluster
library(plotly)

display_parallel_coordinates(gene_cluster, 3)
#gene_dist <- dist(gene_cluster_scaled, method = "euclidean")
gene_hr <- hclust(as.dist(1-cor(t(gene_cluster_scaled), method="spearman")), method="ward.D")
gene_hc <- hclust(as.dist(1-cor(gene_cluster_scaled, method="pearson")), method="ward.D")
mycl <- cutree(gene_hr, h=max(gene_hr$height)/3)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)] 
pdf("heatmap3.pdf",width = 10,height = 30)
## Plot heatmap 
heatmap.2(gene_cluster_scaled, Rowv=as.dendrogram(gene_hr), Colv=as.dendrogram(gene_hc), col=redgreen(75), scale="row", 
          density.info="none", trace="none", RowSideColors=mycolhc) 
#heatmap()
dev.off()
