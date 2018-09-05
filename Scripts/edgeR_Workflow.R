library(edgeR)
library(limma)
library(stringr)
library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(DESeq)
library(splitstackshape)
library(biomaRt)
library(GO.db)
library(pheatmap)

# Define the sub-directory
subdir <- "~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Data/newCounts"

#Run through each of the count files
countFiles <- list.files(
  subdir,
  pattern=".counts",
  full.names = TRUE
)

# A little nifty way to make a data frame from all of your files
counts <- countFiles %>%
  lapply(function(x){
    read.delim(x, 
               header = FALSE,
               row.names = 1) %>%
      set_names(x)}) %>%
  as.data.frame()

#Dont worry about this bad code, it was me being a fool. I wanted to easily make the new column names for my data frame
#new_colnames <- counts[2,]

# The first two lines don't contain any proper data, but they do contain info on the colnames
counts_df <- counts[-c(1:2),]

manual_colnames_CTBvsSTB <- c("72hrs_no_forskolin_Chr", "72hrs_no_forskolin_Start", "72hrs_no_forskolin_End",
                              "72hrs_no_forskolin_Strand", "72hrs_no_forskolin_Length", "72hrs_no_forskolin_Counts",
                              "24hrs_with_forskolin_Chr", "24hrs_with_forskolin_Start", "24hrs_with_forskolin_End",
                              "24hrs_with_forskolin_Strand", "24hrs_with_forskolin_Length", "24hrs_with_forskolin_Counts",
                              "48hrs_with_forskolin_Chr", "48hrs_with_forskolin_Start", "48hrs_with_forskolin_End",
                              "48hrs_with_forskolin_Strand", "48hrs_with_forskolin_Length", "48hrs_with_forskolin_Counts",
                              "72hrs_with_forskolin_Chr", "72hrs_with_forskolin_Start", "72hrs_with_forskolin_End",
                              "72hrs_with_forskolin_Strand", "72hrs_with_forskolin_Length", "72hrs_with_forskolin_Counts")

# Make meta data
samples <- colnames(countsDGE)
exposure <- c(0, 1, 1, 1)
time <- c(0, 24, 48, 72)

BeWo_meta <- data.frame(x = exposure, y = time)  
rownames(BeWo_meta) <- samples
colnames(BeWo_meta) <- c("exposure", "time")

#saveRDS(BeWo_meta, file = "~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Data/BeWo_meta.rds")

variable <- "time"
numericVar <- is.numeric(BeWo_meta[[variable]])
int <- dplyr::if_else(numericVar, 1, 0)
designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
  model.matrix(data = BeWo_meta)
designMatrix


# Add in the new colnames, first for larger counts
colnames(counts_df) <- manual_colnames_CTBvsSTB


# This code selects for the columns which contain the gene counts, the DGElist only wants gene count data, doesn't like anything non-numeric
counts_matrix <- dplyr::select(counts_df, ends_with("Counts"))
# Ive already got the counts for BeWo


# Changes the inputs from factors to numeric
indx <- sapply(counts_matrix, is.factor)
counts_matrix[indx] <- lapply(counts_matrix[indx], function(x) as.numeric(as.character(x)))

# Now create the countsDGE!
countsDGE <- DGEList(counts_matrix)


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

## We won't use those filtering objects again, so feel free to get rid of them

## Record the sample names here
samplenames <- substring(colnames(DGE_RNAseq),
                         0,
                         nchar(colnames(DGE_RNAseq)
                         )
)

# Estimate dispersions
countsDGEdisp <- estimateDisp(countsDGEfilter, design = designMatrix)


# Generalised linear model fit
fit <- glmFit(countsDGEdisp, designMatrix)

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

GO_ID_and_ontology <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                            filters = "go",
                            values = ens_and_GO_ID$go_id,
                            mart = mart)


# Next we can try filtering by logFC and adjusting the FDR appropriately using glmTreat
fcFilter <- glmTreat(fit, coef = ncol(fit$design), lfc = log2(1.2))

topFCgenes <- topTags(fcFilter, n = 50)

rownames(topFCgenes$table) <- str_remove(rownames(topFCgenes$table), "\\..*")

# our mart should still be loaded so just run the GO analysis for this now and compare results
GO_Treat<- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                       values = rownames(topFCgenes$table), mart = mart)


### Create heatmap ###
rownames(cpm) <- str_remove(rownames(cpm), "\\..*")
cpmdf <- as.data.frame(cpm) %>% rownames_to_column(var = "IDs")

ensID_df <- dplyr::select(ens_and_GO_ID, .data$ensembl_gene_id, .data$external_gene_name) %>% as.data.frame()
colnames(ensID_df) <- c("IDs", "genes")

cpm2plot <- left_join(ensID_df, cpmdf)

genenames <- cpm2plot$genes

plotcpm <- dplyr::select(cpm2plot, ends_with("Counts")) %>% as.matrix
rownames(plotcpm) <- genenames

crude <- plotcpm[-c(1:34), ]

pheatmap(head(unique(plotcpm), n = 10))

pheatmap(head(crude, n = 10))

#### Workflow for only 72hrs no forskolin vs 72hrs with forskolin ####

BeWo <- data.frame(counts_matrix$`72hrs_no_forskolin_Counts`, counts_matrix$`72hrs_with_forskolin_Counts`)
rownames(BeWo) <- rownames(countsDGE$counts)

samples <- colnames(BeWo)
exposure <- c(0, 1)
time <- c(0, 72)

BeWo_meta <- data.frame(exposure, time)
rownames(BeWo_meta) <- samples
colnames(BeWo_meta) <- c("exposure", "time")

variable <- "exposure"
numericVar <- is.numeric(BeWo_meta[[variable]])
int <- dplyr::if_else(numericVar, 1, 0)
designMatrixBeWo <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
  model.matrix(data = BeWo_meta)
designMatrixBeWo

# Then for the smaller counts
colnames(BeWo) <- c("72hrs_no_forskolin", "72hrs_with_forskolin")
rownames(designMatrixBeWo) <- colnames(BeWo)


indx <- sapply(BeWo, is.factor)
BeWo[indx] <- lapply(BeWo[indx], function(x) as.numeric(as.character(x)))

BeWoDGE <- DGEList(BeWo)

countsBeWonorm <- calcNormFactors(BeWoDGE)


BeWoDGEdisp <- estimateDisp(countsBeWonorm, design = NULL)

# Now for BeWo
fit <- glmFit(BeWoDGEdisp, design = designMatrixBeWo)

lrt <- glmLRT(fit, coef = 2)

topGenes <- topTags(lrt)

rownames(topGenes$table) <- str_remove(rownames(topGenes$table), "\\..*")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

ens_and_GO_ID_BeWo <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                       values = rownames(topGenes$table), mart = mart)

GO_ID_and_ontology_BeWo <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                            filters = "go",
                            values = ens_and_GO_ID$go_id,
                            mart = mart)

ensIDs <- ens_and_GO_ID$ensembl_gene_id

