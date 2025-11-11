# Creating the dendrogram from the hclust created in hclustering_script 
# plot dendrogram for whole data
hclust <- readRDS("results/hclust.RDS")

png(file = "cluster.png", width = 600, height = 600)
cluster <- plot(hclust, labels = FALSE)
dev.off()

# Reading more information about dendrogram clusters and implementing it 
# Every datapoint eventually joins into one big cluster, so you can cut the dendrogram at certain points to create the desired number of clusters
clustercut <- rect.hclust(hclust , h = 50, border = 2:6)
abline(h = 50, col = 'red')
# colouring the cluster into the clusters
install.packages('dendextend')
install.packages('colorspace')
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust)
avg_col_dend <- color_branches(avg_dend_obj, h = 50)
plot(avg_col_dend)
# cutting the cluster at h=40 to see how many clusters it produces. 
cluster40 <- cutree(hclust, h = 40)
# clusters is a named vector, where the names are the proteins and the values are the cluster IDs
cluster40[1:10]
range(cluster40)
# 21 clusters
plot(hclust, labels = FALSE)

# NCOA should be clustered together as they are nuclear receptor coactivator proteins for steroid receptors - GR
# NCOR should also be clustered together- nuclear corepressor protein. 
# Protein IDs: 
# NCOA1   Q15788
# NCOA2   Q15596
# NCOA3   Q9Y6Q9
# Creating a vector with just those in 
NCOAs <- c("Q15788", "Q15596", "Q9Y6Q9")
cluster40[NCOAs]
# This tells us that Q15788 and Q15596, NCOA1 and NCOA2, are clustered together in cluster 21, and Q9Y6Q9, NCOA3, is clustered alone in cluster 7. 
length(cluster40[cluster40 == 21])
# There are 207 proteins in this cluster

# Are they still clustered together when you cut it at 25?
cluster25 <- cutree(hclust, h = 25)
range(clusters_h25)
# 63 clusters
cluster25[NCOAs]
# NCOA1 and 2 are still in the same cluster (cluster 39)
length(cluster25[cluster25 == 39])
# 139 proteins in this cluster
# Investigating this cluster further: 
# What are the other proteins in this cluster? 
cluster39h25 <- names(cluster25[cluster25 == 39])
write.csv(cluster39h25, file = "results/NCOA1-2_cluster.csv")

# Where is the glucocorticoid receptor clustered?
GR