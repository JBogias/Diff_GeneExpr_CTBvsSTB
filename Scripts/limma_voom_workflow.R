# First, we need to load all of our packages!

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

# Load in the data
countsDGE <- readRDS("~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Results/DGElist_72hrsvs0hrs_forskolin.rds")

# Define the model matrix, this will be needed for the limma workflow
variable <- "time"
numericVar <- is.numeric(BeWo_meta[[variable]])
int <- dplyr::if_else(numericVar, 1, 0)
designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
  model.matrix(data = BeWo_meta)
designMatrix

# Calculate normal factors
countsDGEnorm <- calcNormFactors(countsDGE)


cpm <- cpm(countsDGEnorm)
# Filter
table(rowSums(countsDGEnorm$counts == 0) == 3)

# We want to adjust our threshold to incorporate data that had 20 or more reads, we set a cpm cut-off of 1 which
# is a log-cpm of 0. This ensures that there are at least 20 reads recorded for each gene. We select 3 as our
# sample cut-off as we would need at least 3 samples to have above 20 read counts each to be used in downstream
# analysis on limma.
# Let 'keep.expr' represent our cut-off
keep.exprs <- rowSums(cpm > 1) >= 3

# Check dimensions to see how many genes we have now (to compare later, save it as an object if you want)
before_filter <- dim(countsDGEnorm)

## Now we will subset the data here
countsDGEfilter <- countsDGEnorm[
  keep.exprs,
  keep.lib.sizes = FALSE
  ]

# Now check to see how the filtering worked
after_filter <- dim(countsDGEfilter)

# Estimate dispersions for convoluted reasons
countsDGEdisp <- estimateDisp(countsDGEfilter, design = designMatrix)

# Work that limma-voom magic
fitGeneExpr <- voom(countsDGEdisp, designMatrix, plot = FALSE) %>%
  lmFit(design = designMatrix) %>%
  eBayes()

# Now I think we're trying to make the contrast between gene expression
nGenesTotal <- nrow(fitGeneExpr)

slopeContrast  <- colnames(fitGeneExpr$design)[2]

# This toptable will contain all of the differentially expressed genes and the p value associated with them
# from the linear model fit
topGeneslim <- topTable(fitGeneExpr,
                     coef = slopeContrast, number = nGenesTotal) %>%
  rownames_to_column("GeneID")

# Since the gene IDs had decimal places we should take those off here
topGeneslim$GeneID <- str_remove(topGenes$GeneID, "\\..*")

# Now download the biomart object
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

# And pull out your GO terms
GO_Limma <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                values = topGeneslim$GeneID, mart = mart)
