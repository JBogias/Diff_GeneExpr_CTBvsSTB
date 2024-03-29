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

# Manually make the new colnames, if theres and easier way, I'm all ears
manual_colnames_mVSf <- c("FemCTB_Chr", "FemCTB_Start", "FemCTB_End", "FemCTB_Strand", "FemCTB_Length", "FemCTB_Counts",
                     "FemSTB_Chr", "FemSTB_Start", "FemSTB_End", "FemSTB_Strand", "FemSTB_Length", "FemSTB_Counts",
                     "MCTB_Chr", "MCTB_Start", "MCTB_End", "MCTB_Strand", "MCTB_Length", "MCTB_Counts",
                     "MSTB_Chr", "MSTB_Start", "MSTB_End", "MSTB_Strand", "MSTB_Length", "MSTB_Counts")

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

# fit the gene expression to a linear model
names <- c("Sex", "Cell_type")
samples <- c("FemCTB", "FemSTB", "MCTB", "MSTB")
FemCTB_meta <- c("Female", "Cytotrophoblast")
FemSTB_meta <- c("Female", "Syncytiotrophoblast")
MCTB <- c("Male", "Cytotrophoblast")
MSTB <- c("Male", "Syncytiotrophoblast")

# Create the metadata
metadata <- rbind(FemCTB_meta, FemSTB_meta, MCTB, MSTB)
rownames(metadata) <- samples
colnames(metadata) <- names
metadata <- as.data.frame(metadata)


# This ias a little convenient function to get the design matrix
getDesignMatrix <- function(variable, metadata) {
  numericVar <- is.numeric(metadata[[currentVar]])
  int <- dplyr::if_else(numericVar, 1, 0)
  as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
             model.matrix(data = metadata)
}

designMatrix <- getDesignMatrix(currentVar, metadata = metadata)

fitGeneExpr <- limma::voom(countsDGE, designMatrix, plot = FALSE) %>%
  limma::lmFit(design = designMatrix) %>%
  limma::eBayes()
nGenesTotal <- nrow(fitGeneExpr)
slopeContrast <- colnames(fitGeneExpr$design)[2]
topGenes <- limma::topTable(fitGeneExpr,
                            coef = slopeContrast, number = nGenesTotal) %>%
  tibble::rownames_to_column("GeneID")

## Define gene and sample order to subset for in the topGenes object
inputVals <- metadata[[currentVar]]
sampleOrder <- order(inputVals)

GeneOrder <- order(topGenes$GeneID)
DE_Genes <- topGenes[GeneOrder, ] %>%
  tibble::as_data_frame()

### Select only for FDR and gene names

geneExprPValFull <- dplyr::select(DE_Genes, .data$GeneID, .data$P.Value)

## Retireve Gene IDs for all samples
#First we want to get the gene names and FDR for all of the samples, this way we can make an FDR bar that will give a colour indication of the FDRcutoff

topGenesTotal <- rownames(limma::topTable(fitGeneExpr, coef = slopeContrast, number = nGenesTotal))
inputValsTotal <- metadata[[currentVar]]
sampleOrderTotal <- order(inputValsTotal)

### Get the cpm for all values
#This is the step required for the FDR

geneExprCPM2plot <- cpm(countsDGE, log = TRUE)[topGenesTotal, sampleOrder]
colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)

## Scale by rows
#In this function we are calculating the mean and standard deviation for every row of the matrix. Then we subtract the mean from the raw value and divide by the calculated standard deviation to generate our scaled values. Dividing by the standard deviation means that

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)

saveRDS(scaledCPMgeneExpr, file = "~/Documents/PhD/Data/scaledCPMgeneExpr.rds")
saveRDS(topGenes, file = "~/Documents/PhD/Data/topGenes.rds")

geneExprCPM2ggplot <- reshape2::melt(scaledCPMgeneExpr, value.name = "scaledCPM")
colnames(geneExprCPM2ggplot) <- c("GeneID", "Sample", "scaledCPM") #%>%
  #as.data.frame()


geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
GeneOrder <- levels(geneExprCPM2ggplot$GeneID)
GeneOrder <- data.frame(GeneOrder, 1)
rownames(GeneOrder) <- GeneOrder$GeneOrder
colnames(GeneOrder) <- c("GeneID", "boo")
geneExprPVal <- dplyr::left_join(GeneOrder, geneExprPVal, by = "GeneID")
geneExprPVal$boo <- NULL
geneExprPVal$booo <- NULL

geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = 'GeneID')
geneExprCPM2ggplotPVal$GeneID <- as.factor(geneExprCPM2ggplotPVal$GeneID)

geneExprCPM2ggplotPVal %>% head

saveRDS(geneExprCPM2ggplotPVal, file = "~/Documents/PhD/Data/geneExprCPM2ggplotPVal.rds")


# Start on the heatmap

## Set colour information

nColours <- 101
varRange <- range(metadata[[currentVar]])
varGradient <- floor(seq(varRange[1], varRange[2], length.out = nColours))
varPositions <- findInterval(appmetaG[[currentVar]], varGradient)
varPalette <- grDevices::colorRampPalette(c("brown", "yellow"))(nColours)
varColours <- varPalette[varPositions]

saveRDS(varPalette, "~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/varPalette.rds")
nColoursPVal <- nGenesTotal
PValrange <- range(geneExprPValFull$P.Value)
PValgradient <- seq(from = PValrange[1], to = PValrange[2], length.out = nColoursPVal)
PValpositions <- findInterval(geneExprPVal$P.Value, PValgradient)
PValpalette <- grDevices::colorRampPalette(c("green", "light green", "yellow", "pink", "red", "dark red"))(nColoursPVal)
PValcolours <- PValpalette[PValpositions]

heatmapColours <- grDevices::colorRampPalette(c("dark blue", "blue", "cyan", "white", "pink", "red", "dark red"))(540)

geneExprPVal <- tibble::as.tibble(geneExprPVal)
geneExprPVal$GeneID <- factor(geneExprPVal$GeneID)
geneExprCPM2ggplot <- tibble::as.tibble(geneExprCPM2ggplot)
levels(geneExprPVal$GeneID) <- levels(geneExprCPM2ggplot$GeneID)

GeneExprheatmap <- ggplot2::ggplot(data = head(geneExprCPM2ggplotPVal, 100),
                                   aes(x = Sample, y = GeneID, fill = scaledCPM)) +
  ggplot2::geom_raster(stat = "identity", position = "identity") +
  ggplot2::xlab("Sample") +
  ggplot2::scale_fill_gradientn(colours = heatmapColours) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

GeneExprheatmap

# Now for the GO term stuff

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

ens_and_GO_ID <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id",
                       values = topGenes$GeneID, mart = mart)

GO_ID_and_ontology <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                            filters = "go",
                            values = ens_and_GO_ID$go_id,
                            mart = mart)






