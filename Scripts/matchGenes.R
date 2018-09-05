library(edgeR)
library(limma)
library(magrittr)
library(dplyr)
library(readr)
library(tibble)
library(RColorBrewer)
library(DESeq)
library(splitstackshape)
library(pheatmap)
library(ggdendro)
library(plotly)
library(ggthemes)

# My plan here is to match up the gene names between the results I got vs the results our friend Azar got
# I'm going to left join the ensembl gene IDs together and then see which ones I got turn up in Azars. I'm hoping
# to find all of my genes

countsDGE <- readRDS("~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Results/DGElist_72hrsvs0hrs_forskolin.rds")

# First up, we need to load in the data
topGenesGLM <- readRDS("Results/topTags_72hrs_vs_0hrs.rds")
topGeneslim <- readRDS("Results/topTable_72hrs_vs_0hrs.rds")

# Now download the biomart object
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

# And pull out your GO terms
GO_Limma <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id",
                  values = topGeneslim$GeneID, mart = mart)

colnames(GO_Limma) <- c("GeneID", "Gene")

GO_GLM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id",
                values = topGenesGLM$GeneID, mart = mart)

colnames(GO_GLM) <- c("GeneID", "Gene")

# Now for the first left join
TopGLM_Names <- left_join(topGenesGLM, GO_GLM, by = "GeneID")
TopLimma_Names <- left_join(topGeneslim, GO_Limma, by = "GeneID")

# Now to join them to Azars results
AzarResults <- read_csv("Files/edgeR_results.csv")

GLM_in_Azar <- left_join(TopGLM_Names, AzarResults, by = "Gene")
Limma_in_Azar <- left_join(TopLimma_Names, AzarResults, by = "Gene")
