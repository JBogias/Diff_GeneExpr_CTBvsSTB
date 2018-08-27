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

# Define the sub-directory
subdir <- "~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Data/countsBeWo"

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

saveRDS(BeWo_meta, file = "~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Data/BeWo_meta.rds")

variable <- "exposure"
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



# Estimate dispersions
countsDGEdisp <- estimateDisp(countsDGEnorm, design = designMatrix)


# Generalised linear model fit
fit <- glmFit(countsDGEdisp, designMatrix)

lrt <- glmLRT(fit, coef = 2)

topGenes <- topTags(lrt)

rownames(topGenes$table) <- str_remove(rownames(topGenes$table), "\\..*")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

ens_and_GO_ID <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                       values = rownames(topGenes$table), mart = mart)

GO_ID_and_ontology <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                            filters = "go",
                            values = ens_and_GO_ID$go_id,
                            mart = mart)


#### Workflow for only 72hrs no forskolin vs 72hrs with forskolin ####

BeWo <- data.frame(counts_matrix$`72hrs_no_forskolin_Counts`, counts_matrix$`72hrs_with_forskolin_Counts`)

samples <- colnames(BeWo)
exposure <- c(0, 1)
time <- c(0, 72)

BeWo_meta <- data.frame(exposure, time)
rownames(BeWo_meta) <- samples
colnames(BeWo_meta) <- c("exposure", "time")

variable <- "exposure"
numericVar <- is.numeric(BeWo_meta[[variable]])
int <- dplyr::if_else(numericVar, 1, 0)
designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
  model.matrix(data = BeWo_meta)
designMatrix

# Then for the smaller counts
colnames(BeWo) <- c("72hrs_no_forskolin", "72hrs_with_forskolin")
rownames(designMatrix) <- colnames(BeWo)


indx <- sapply(BeWo, is.factor)
BeWo[indx] <- lapply(BeWo[indx], function(x) as.numeric(as.character(x)))

BeWoDGE <- DGEList(BeWo)

countsBeWonorm <- calcNormFactors(BeWoDGE)


BeWoDGEdisp <- estimateDisp(countsBeWonorm, design = designMatrix)

# Now for BeWo
fit <- glmFit(BeWoDGEdisp, design = NULL)
