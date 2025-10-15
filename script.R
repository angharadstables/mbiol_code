##############
# Mbiol Code #
##############
# Caleb summary ------------
# RIME data shows proteins which interact with the glucocorticoid receptor 
# vs interactions with IgG (control). log fold change for each protein can be found
# in DE_results.csv (this data was first processed in fragpipe)
# roughly shows whether or not proteins interacts with GR

# Dex activates the glucocorticoid receptor. RNA seq data of cells +/- Dex is found
# in Dex_RNA-seq.rds.

# aim is to compare LFC of RNA-seq data to LFC of RIME data (LFC of RNA vs LFC of 
# proteins when GR is active) and see which proteins and transcripts correlate
# then create a matrix of the different R values and cluster
###################################
# Installing and Loading Packages #
###################################
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