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
range(cluster25)
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
GR <- c("P04150")
cluster25[GR]
# It is in cluster 32. What other proteins are in this cluster?
clusterGR <- names(cluster25[cluster25 == 32])
write.csv(clusterGR, file = "results/GR_cluster_32.csv")
# how many proteins in the cluster 
length(cluster25[cluster25 == 32])
# 253 proteins in the cluster. 

# Create a diagram with coloured clusters - NCOA1
dend <- as.dendrogram(hclust)
# plot uncoloured dend
plot(dend)
# Matching the order of the cluster labels in cluster25 and dend object
clusters_ordered <- cluster25[order.dendrogram(dend)]
# assigning target cluster
target_cluster <- clusters_ordered["Q15788"]
target_cluster 
# this is cluster 39. correct. will include both Q15788 and Q15596
# each cluster ID must be assigned a colour
cluster_ids <- unique(clusters_ordered) # get cluster IDs
cluster_ids
# these are all the cluster ids 
# set it so that if it is not the target cluster it will be black, if it is red. 
cluster_colours <- ifelse(cluster_ids == target_cluster, "red", "black")
cluster_colours
# everything is black except 39 which is red. 
# cluster_colours is currently an unnamed vector which contains just colours
names(cluster_colours) <- (cluster_ids)
# check again 
cluster_colours 
# create coloured dendrogram
dend_coloured <- color_branches(dend, h = 25, groups = clusters_ordered,
                                col = cluster_colours)
plot(dend_coloured, leaflab = "none")
# add custom labels to this plot
poi <- NCOAs
custom_labels <- ifelse(labels(dend) %in% poi, labels(dend), "")
labels(dend_coloured) <- custom_labels
plot(dend_coloured)

# save as png
png(filename = "cluster_NCOAs.png", width = 600, height = 600)
plot(dend_coloured, cex = 0.1, main = "NCOA1/2/3")
dev.off()