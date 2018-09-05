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

countsDGE <- readRDS("Results/DGElist_72hrsvs0hrs_forskolin.rds")

countsDGEnorm <- calcNormFactors(countsDGE)

countsDGEdisp <- estimateDisp(countsDGEnorm, design = NULL)

# Generalised linear model fit
fit <- glmFit(countsDGEdisp)

# First we can check for significant DE genes using genewise negative binomial generalized linear models
lrt <- glmLRT(fit, coef = 2)

topGenes <- topTags(lrt, n = 50)

# Since the gene IDs had decimal places we should take those off here
rownames(topGenes$table) <- str_remove(rownames(topGenes$table), "\\..*")

# Now download the biomart object
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

# And pull out your GO terms
GO_LRT <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                values = rownames(topGenes$table), mart = mart)
