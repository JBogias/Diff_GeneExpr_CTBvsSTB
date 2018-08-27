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

# Add in the new colnames
colnames(counts_df) <- manual_colnames_CTBvsSTB

# This code selects for the columns which contain the gene counts, the DGElist only wants gene count data, doesn't like anything non-numeric
counts_matrix <- dplyr::select(counts_df, ends_with("Counts"))

# Changes the inputs from factors to numeric
indx <- sapply(counts_matrix, is.factor)
counts_matrix[indx] <- lapply(counts_matrix[indx], function(x) as.numeric(as.character(x)))

# Now create the countsDGE!
countsDGE <- DGEList(counts_matrix)

# Get the cpm that will be required for filtering out lowly expressed genes
cpm <- cpm(countsDGE, log = FALSE)

# So lets check our read quantities
# Here I'm checking how many genes have no reads across all samples
table(rowSums(countsDGE$counts == 0) == 3)
table(rowSums(lcpm == 0) >= 3)

# We want to adjust our threshold to incorporate data that had 20 or more reads, we set a cpm cut-off of 1 which
# is a log-cpm of 0. This ensures that there are at least 20 reads recorded for each gene. We select 3 as our
# sample cut-off as we would need at least 3 samples to have above 20 read counts each to be used in downstream
# analysis on limma.
# Let 'keep.expr' represent our cut-off
keep.exprs <- rowSums(cpm > 1) >= 2

# Check dimensions to see how many genes we have now (to compare later, save it as an object if you want)
before_filter <- dim(countsDGE)

## Now we will subset the data here
countsDGE <- countsDGE[
  keep.exprs,
  keep.lib.sizes = FALSE
  ]

# Now check to see how the filtering worked
after_filter <- dim(countsDGE)

## We won't use those filtering objects again, so feel free to get rid of them
## Record the sample names here
samplenames <- substring(colnames(countsDGE),
                         0,
                         nchar(colnames(countsDGE)
                         )
)


### Multi-dimensional scaling ###
# First set the variable of interest, in this case it is "pH"
currentVar <- "pH"

# Now is the variable numeric or categorical
numericVar <- is.numeric(appmetaG[[currentVar]])
int <- if_else(numericVar, 1, 0)

# Create the model matrix here using a formula from our predefined parameters
# designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
#   model.matrix(data = appmetaG)

## Now we can fit the gene expression to a linear model
# This is essentially what I was going on about about, the voom
# method forces it to be linear data and from there we can run a lmFit
# on it, then run some baysien magic on it with eBayes()

# There is no metadata that I could find from the study, so I've set the paramter to NULL
# this means that each array is treated as a replicates

fitGeneExpr <- voom(countsDGE, design = NULL, plot = FALSE) %>%
  lmFit(design = NULL) %>%
  eBayes()

# Now I think we're trying to make the contrast between gene expression
nGenesTotal <- nrow(fitGeneExpr)

slopeContrast  <- colnames(fitGeneExpr$design)[1]

# This toptable will contain all of the differentially expressed genes and the p value associated with them
# from the linear model fit
topGenes <- topTable(fitGeneExpr,
                     coef = slopeContrast, number = nGenesTotal) %>%
  rownames_to_column("GeneID")

GeneOrder <- order(topGenes$GeneID)

topGenes <- topGenes[GeneOrder, ] %>%
  as_data_frame()

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

topGenes$GeneID <- str_remove(topGenes$GeneID, "\\..*")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

ens_and_GO_ID <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                       values = topGenes$GeneID, mart = mart)

GO_ID_and_ontology <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                            filters = "go",
                            values = ens_and_GO_ID$go_id,
                            mart = mart)


