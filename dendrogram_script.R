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
