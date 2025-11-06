##############
# Mbiol Code #
##############
# Caleb summary ------------
# RIME data shows proteins which interact with the glucocorticoid receptor 
# vs interactions with IgG (control). log fold change for each protein can be found
# in DE_results.csv (this data was first processed in fragpipe)
# roughly shows whether or not proteins interacts with GR

# changes
# Dex activates the glucocorticoid receptor. RNA seq data of cells +/- Dex is found
# in Dex_RNA-seq.rds.

# aim is to compare LFC of RNA-seq data to LFC of RIME data (LFC of RNA vs LFC of 
# proteins when GR is active) and see which proteins and transcripts correlate
# then create a matrix of the different R values and cluster
###################################
# Installing and Loading Packages #
###################################
options(repos = c(CRAN = "http://cran.rstudio.com"))
library(tidyverse) 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

#################################
# Loading and manipulating data #
#################################
#### RIME data ####
# Import
rime_data <- read_csv("data-raw/DE_results.csv")
# Make into smaller dataset which only contains cellines also in RNAseq dataset (Jurkat, MCF7, KMBC, M231)
rime_lfcs <- rime_data[,(names(rime_data) %in% 
                           c("Protein ID",
                             "Gene Name",
                             "Jurkat_IgG_vs_Jurkat_GR_log2 fold change",
                             "KMBC2_IgG_vs_KMBC2_GR_log2 fold change",
                             "MCF7_IgG_vs_MCF7_GR_log2 fold change",
                             "X231_IgG_vs_X231_GR_log2 fold change"))]
# Rename columns
names(rime_lfcs)[names(rime_lfcs) == "Jurkat_IgG_vs_Jurkat_GR_log2 fold change"] <- "Jurkat"
names(rime_lfcs)[names(rime_lfcs) == "KMBC2_IgG_vs_KMBC2_GR_log2 fold change"] <- "KMBC2"
names(rime_lfcs)[names(rime_lfcs) == "MCF7_IgG_vs_MCF7_GR_log2 fold change"] <- "MCF7"
names(rime_lfcs)[names(rime_lfcs) == "X231_IgG_vs_X231_GR_log2 fold change"] <- "231"

# get rid of the gene ID
# just keeping protein ID
rime_lfcs_p <- rime_lfcs[,(names(rime_lfcs) %in% 
                             c("Protein ID",
                               "Jurkat",
                               "KMBC2",
                               "MCF7",
                               "231"))]
# set protein IDs as row names
rime_lfcs_p <- as.data.frame(rime_lfcs_p) # no longer a tibble
rownames(rime_lfcs_p) <- rime_lfcs_p$`Protein ID`
rime_lfcs_p$`Protein ID` <- NULL

# give this data frame a simple name and convert to a matrix for later
rime <- as.matrix(rime_lfcs_p)

write.csv(rime, file = "data-tidy/rime.csv")

#### RNAseq data ####
# make sure DESeq2 is installed above
dds <- readRDS(file="data-raw/Dex_RNA-seq.rds")

#############
# RIME data #
#############


################
# RNA seq data #
################
# vst normalises gene variance for the PCA plot. 
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "cellline"))

# add new column  to dds which combines condition and cellline (necessary to 
# compare lfcs of individual cell types)
colData(dds)$condition_cellline <- factor(paste0(dds$condition, "_", dds$cellline))
# update DESeq design to include new combined factor
design(dds) <- ~ condition_cellline

# pre-filtering to remove low expression counts
dds <- dds[rowSums(counts(dds)) > 1,]
# run DESeq
dds <- DESeq(dds)
# use contrast to compare DMSO vs Dex treatment in Jurkat cells
res_Jurkat <- results(dds, contrast = c("condition_cellline", "DMSO_Jurkat", "Dex_Jurkat"))
# output Jurkat results as data frame
results_Jurkat <- as.data.frame(res_Jurkat)
# add a column with the cell line
results_Jurkat$cellline <-c("Jurkat")
# delete all columns except  cellline, pvalue (padj) and log2FoldChange
results_Jurkat <- results_Jurkat |> 
  select(cellline, padj, log2FoldChange)
# rename log2foldchange to specify this is the RNA LFC
names(results_Jurkat)[names(results_Jurkat) == "log2FoldChange"] <- "rna_lfc"

# do the same for other cell lines
# Jurkat, MCF7, KMBC2 and 231 are the four cell types which appear in both the RNASeq
# and RIME datasets

# MCF7 
res_MCF7 <- results(dds, contrast = c("condition_cellline", "DMSO_MCF7", "Dex_MCF7"))
results_MCF7 <- as.data.frame(res_MCF7)
results_MCF7$cellline <- c("MCF7")
results_MCF7 <- results_MCF7 |> 
  select(cellline, padj, log2FoldChange)
names(results_MCF7)[names(results_MCF7) == "log2FoldChange"] <- "rna_lfc"

# KMBC2
res_KMBC2 <- results(dds, contrast = c("condition_cellline", "DMSO_KMBC2", "Dex_KMBC2"))
results_KMBC2 <- as.data.frame(res_KMBC2)
results_KMBC2$cellline <- c("KMBC2")
results_KMBC2 <- results_KMBC2 |> 
  select(cellline, padj, log2FoldChange)
names(results_KMBC2)[names(results_KMBC2) == "log2FoldChange"] <- "rna_lfc"

# 231
res_231 <- results(dds, contrast = c("condition_cellline", "DMSO_M231", "Dex_M231"))
results_231 <- as.data.frame(res_231)
results_231$cellline <- c("231")
results_231 <- results_231 |> 
  select(cellline, padj, log2FoldChange)
names(results_231)[names(results_231) == "log2FoldChange"] <- "rna_lfc"



# merge all rna seq datasets into one tidy dataset
rnaseq_lfcs <- bind_rows(results_Jurkat, results_KMBC2, results_231, results_MCF7)

# add column with gene names
rnaseq_lfcs$gene <- rownames(rnaseq_lfcs)
rownames(rnaseq_lfcs) <- NULL
# reorder columns
rnaseq_lfcs <- rnaseq_lfcs[, c("gene", "cellline", "rna_lfc")]

# make it wide format, instead of tidy. each column is a different cell type

# ensembl ID is used for gene names, which contains a gene version at the end of the name
# the genes are named different versions
# remove version from ensembl IDs
rnaseq_lfcs <- rnaseq_lfcs |> 
  mutate(gene = sub("\\..*", "", gene))

# now put into wide format
rnaseq_wide <- rnaseq_lfcs |> 
  pivot_wider(names_from = cellline, values_from = rna_lfc)

# put gene names back into the column name
rnaseq_wide <- as.data.frame(rnaseq_wide) # no longer a tibble
rownames(rnaseq_wide) <- rnaseq_wide$gene
rnaseq_wide$gene <- NULL

# convert this data frame to a matrix so it's numeric and give it a simple name
rnaseq <- as.matrix(rnaseq_wide)

#############################
# Running Correlation Tests #
#############################
# compare all proteins and all transcripts just to see what happens!
# create matrix large enough to compare all proteins and all transcripts
mat1 <- matrix(ncol=nrow(rime), nrow=nrow(rnaseq))

colnames(mat1) <- rownames(rime)
rownames(mat1) <- rownames(rnaseq)

mat2 <- mat1

# create a for loop which runs Pearson test for each combination of protein and transcript
# p value is put into mat1
# R squared is put into mat2 (the correlation statistic)
for (xsamples in rownames(rime)) {
  
  # convert rime row to numeric vector
  x <- rime[xsamples, ]
  
  cor_mat <- apply(rnaseq, 1, function(y) { 
    test <- cor.test(x, y)
    c(cor = test$estimate, p.value = test$p.value)
  } )
  
  cor_mat <- t(cor_mat)
  
  # putting names here will slow up a bit, but protects against reordering.
  mat1[names(cor_mat[,'cor.cor']), xsamples] <- cor_mat[,'cor.cor']
  mat2[names(cor_mat[,'p.value']), xsamples] <- cor_mat[,'p.value']
  
  
}

# save results so they can be accessed after viking is finished with them
saveRDS(mat1, "results/mat1.rds")
saveRDS(mat2, "results/mat2.rds")