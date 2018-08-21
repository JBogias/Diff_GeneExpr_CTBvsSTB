library(edgeR)
library(limma)
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

table(rowSums(countsDGE$counts == 0) == 3)

currentVar <- "Cell_type"






