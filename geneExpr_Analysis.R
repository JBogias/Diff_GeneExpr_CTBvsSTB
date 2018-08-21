# A revised workflow for the gene expression analysis using EdgeR and limm from my honours in 2017

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

### Part 1: Load in all data ###

## Preparation and loading of data required for analysis

# Metadata containing information about our samples will need to be loaded in
# It is worth being mindful that the metadata should have the same dimensions as our
# experimental data

# Let 'appmetaG' represent our metadata
appmetaG <- readRDS("~/Documents/PhD/R_Projects/InteGRAPE/InteGRAPE_App_earlyVersion/vitisData/appmetaG_ALL.rds")

# If there isn't a DGElist object already loaded, we'll need to generate one from a counts matris
# If you don't have a counts matrix, there are plenty of ways of making one from dile types, primarily BAM files
counts <- readRDS(file ="~/Documents/PhD/R_Projects/InteGRAPE/InteGRAPE_App_earlyVersion/vitisData/counts.rds")

# This file here contains information about the samples that are in the counts data, in my case I have a csv file which
# has the vineyards in which my grapevine samples have been sourced from
phenoData <- read_csv(file = "~/Documents/PhD/R_Projects/InteGRAPE/InteGRAPE_App_earlyVersion/vitisData/phenoData.csv") %>%
	lapply(as.factor) %>%
	as_data_frame()

# And from that, we can now use the standard edgeR function to create our DGElist
DGE_RNAseq <- DGEList(counts = counts, group = factor(phenoData$Region))

# Alternatively, you can just load your own data
DGE_RNAseq <- readRDS("~/Documents/PhD/R_Projects/InteGRAPE/InteGRAPE_App_earlyVersion/vitisData/DGE_RNAseqALL.rds")

### Part 2: Data exploration and processing ###

# Now we will want to convert to a standard counts format, in this case I will convert to
# gene counts per million

cpm <- cpm(DGE_RNAseq)

# we can also  convert to log-cpm
lcpm <- cpm(DGE_RNAseq, log = TRUE)


## Remove any lowly expressed genes
# An important step is to remove any noise that is in our data - which is lowly expressed molecules that
# happen to be in our transcriptomic snapshot. We are only looking for gene expression here, so removing
# lowly expressed transcripts won't do us any harm

# So lets check our read quantities
# Here I'm checking how many genes have no reads across all samples
table(rowSums(DGE_RNAseq$counts == 0) == 3)

# We want to adjust our threshold to incorporate data that had 20 or more reads, we set a cpm cut-off of 1 which
# is a log-cpm of 0. This ensures that there are at least 20 reads recorded for each gene. We select 3 as our
# sample cut-off as we would need at least 3 samples to have above 20 read counts each to be used in downstream
# analysis on limma.
# Let 'keep.expr' represent our cut-off
keep.exprs <- rowSums(cpm > 1) >= 3

# Check dimensions to see how many genes we have now (to compare later, save it as an object if you want)
before_filter <- dim(DGE_RNAseq)

## Now we will subset the data here
DGE_RNAseq <- DGE_RNAseq[
	keep.exprs,
	keep.lib.sizes = FALSE
]

# Now check to see how the filtering worked
after_filter <- dim(DGE_RNAseq)

## We won't use those filtering objects again, so feel free to get rid of them

## Record the sample names here
samplenames <- substring(colnames(DGE_RNAseq),
	0,
	nchar(colnames(DGE_RNAseq)
		)
	)



### Multi-dimensional scaling ###

# First set the variable of interest, in this case it is "pH"
currentVar <- "pH"

# Now is the variable numeric or categorical
numericVar <- is.numeric(appmetaG[[currentVar]])
int <- if_else(numericVar, 1, 0)

# Create the model matrix here using a formula from our predefined parameters
designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
	model.matrix(data = appmetaG)

# Now we want to generate a forced linear version of our cpm counts using the voom method
# because a lot of the statistical methods we will be using will be from limma, which is
# pretty much linear statistics
voomGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE)

# Now we can begin to plot the MDS
mdsGeneExpr <- plotMDS(voomGeneExpr$E,
			labels = colnames(DGE_RNAseq$counts))

# We want to add a "tool tip" which will allow us to retrieve values from our plot
# by just hovering over a certain data point, we require this function here to make it work
all_values <- function(x) {
	if(is.null(x)) return(NULL)
	paste0(names(x), ": ", format(x), collapse = "<br />")
	}

# Now lets set our variable
Variable <- appmetaG[[currentVar]]

# Set the regions info
Region <- DGE_RNAseq$samples$group

# Set the colour palette we will use
MDScolour <- colorRampPalette(c("brown", "yellow"))(100)


## Here we will put it all together to visualise our plot!
# I won't go over what every little thing does, the package documentation
# exists for that hun <3 xo
mdsGeneExpr[[3]] %>%
	set_colnames(c("Dim1", "Dim2")) %>%
	as.data.frame() %>%
	ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
	layer_points() %>%
	add_legend("fill", title = paste("Variable:", currentVar),
		orient = "left") %>%
	add_legend("shape", title = "Regions",
		orient = "right") %>%
	add_tooltip(all_values, "hover")

# If you will be running this bad boy in shiny, then you merely need to have it inside of a
# reactive function then you need to bind it to the package, an example:

MDStoPlot <- reactive({
		code for MDS plot (mine above essentially
		})

MDStoPlot %>% bind_shiny("MDS", "p_ui")


### Differential Gene expression analysis ###

# Now we're going to get into the actual analysis

# Make another design matrix (model matrix)
currentVar <- "Planting_Date"
numericVar <- is.numeric(appmetaG[[currentVar]])
int <- if_else(numericVar, 1, 0)

designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
	model.matrix(data = appmetaG)


## Now we can fit the gene expression to a linear model
# This is essentially what I was going on about about, the voom
# method forces it to be linear data and from there we can run a lmFit
# on it, then run some baysien magic on it with eBayes()

fitGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) %>%
	lmFit(design = designMatrix) %>%
	eBayes()

# Now I think we're trying to make the contrast between gene expression
nGenesTotal <- nrow(fitGeneExpr)

slopeContrast  <- colnames(fitGeneExpr$design)[2]

# This toptable will contain all of the differentially expressed genes and the p value associated with them
# from the linear model fit
topGenes <- topTable(fitGeneExpr,
	coef = slopeContrast, number = nGenesTotal) %>%
	rownames_to_column("GeneID")

GeneOrder <- order(topGenes$GeneID)

topGenes <- topGenes[GeneOrder, ] %>%
	as_data_frame()

# Now we want to retrieve all of the genes, this is important for plotting
geneExprPValFull <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)

# This code here will be important for making our FDR bar
# This will help us in setting our cutoff
topGenesTotal <- rownames(topTable(fitGeneExpr, coef = slopeContrast, number = nGenesTotal))

inputVals <- appmetaG[[currentVar]]

sampleOrder <- order(inputVals)


# Now we will select for only FDR and the gene names
geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)


### Bonus Section: Visualise the linear representation of Genes ###
# This can be used to make your heatmap even more spiffy than ever, we can add it in later to it
# First create our voom object
voomObj <- voom(DGE_RNAseq, designMatrix, plot = FALSE)

# Select a gene, any random one will do (selecting a row here which will have a gene in it)
selectGene <- 1234

# Pull the row out
Gene <- topGenes[[selectGene, 1]]

# Get its pvalue from the data
pvalue <- topGenes[[selectGene, 2]]

# Subset for the counts of the selected gene
GeneCounts <- voomObj$E[Gene, ]

# Now we can plot the linear trend of the gene expression
LinearRep <- tibble(factor(names(GeneCounts[sampleOrder]),
	levels = unique(names(GeneCounts[sampleOrder]))),
	GeneCounts[sampleOrder])
LinearRep$Variable <- inputVals[sampleOrder]

colnames(LinearRep) <- c("Region", "Counts", "Variable")

# From here we can now plot our gene expression
plotLinearRep <- ggplot(data = LinearRep, aes(x = Variable, y = Counts)) +
	geom_point(stat = "identity")

plotLinearRep + geom_smooth(method = "lm") +
	ggtitle(paste(Gene, ": p =", round(pvalue, 11))) +
	theme_bw()


### Part 5: The heatmap ###
geneExprCPM2plot <- cpm(DGE_RNAseq, log = TRUE)[topGenesTotal, sampleOrder]

colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)

## Define scaling function to make colour values not be silly
scale_rows = function(x) {
	m = apply(x, 1, mean, na.rm = TRUE)
	s = apply(x, 1, sd, na.rm = TRUE)
	return((x - m) / s)
	}

## Generate CPM dataframe for the heatmap
scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)

geneExprCPM2ggplot <- reshape2::melt(scaledCPMgeneExpr, value.name = "scaledCPM")

colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")

geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPValFull, by = "GeneID")

head(geneExprCPM2ggplotPVal)


# Think of this as our little heatmap control pad, we can control the genes we view here by changing nGenes
nGenes <- 30

topGenes <- topTable(fitGeneExpr, coef = slopeContrast, number = nGenes)

topGenes <- rownames_to_column(topGenes, "GeneID")

inputVals <- appmetaG[[currentVar]]

sampleOrder <- order(inputVals)


# Now we run the same code as before to get the FDR and gene names, this one for the
# new data frame we've made that contains the top x genes
geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)



# Now we will make and object which contains all of the cpm data
geneExprCPM2plot <- cpm(DGE_RNAseq, log = TRUE)[topGenes$GeneID, sampleOrder]

# Also set the colnames here
colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)


### Generate clustered and sorted CPM data frame ###
# A necessity for the data visualisation in the heatmap is to cluster the reads, to
# keep all of the colours nicely grouped, all the data will be left-joined so it isn't
# all jumbled up either
# It's also worth mentioning that this step is a high memory process, it takes a LONG time
# and is one of the main reasons that the app is kinda bad at the moment
scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)

dendogramObject <- as.dendrogram(hclust(dist(scaledCPMgeneExpr), method = "ward.D2"))

clusterOrder <- order.dendrogram(dendogramObject)

clusteredMatrix <- scaledCPMgeneExpr[clusterOrder, ]

geneExprCPM2ggplot <- reshape2::melt(clusteredMatrix, value.name = "scaledCPM")

colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")

# now run the left_join
geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = "GeneID")

geneExprCPM2ggplotPVal$GeneID <- as.factor(geneExprCPM2ggplotPVal$GeneID)

head(geneExprCPM2ggplotPVal)


## From here we can set aesthetics for the heatmap, no point replicating it I don't think...





