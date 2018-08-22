
plotHeatmap <- function(DGElist, variable, metadata, nGenes = 30) {
  ### Differential Gene expression analysis ###

  # Now we're going to get into the actual analysis

  # Make another design matrix (model matrix)
  currentVar <- variable
  appmetaG <- metadata
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
  nGenes <- nGenes

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

  ### This should now all work downstream to the heatmap

  ### Will chuck in heatmap code here


  # Purely data visualisation, this will require p-value information and metadata information,
  # most of it should be prepared by this stage but there area few extra steps as well.
  ## Set colour information
  nColours <- 101
  varRange <- range(appmetaG[[currentVar]])
  varGradient <- floor(seq(varRange[1], varRange[2], length.out = nColours))
  varPositions <- findInterval(appmetaG[[currentVar]], varGradient)
  varPalette <- colorRampPalette(c("brown", "yellow"))(nColours)
  varColours <- varPalette[varPositions]

  nColoursPVal <- nGenesTotal
  PValrange <- range(geneExprPValFull$P.Value)
  PValgradient <- seq(from = PValrange[1], to = PValrange[2], length.out = nColoursPVal)
  PValpositions <- findInterval(geneExprPVal$P.Value, PValgradient)
  PValpalette <- colorRampPalette(c("green", "light green", "yellow", "pink", "red", "dark red"))(nColoursPVal)
  PValcolours <- PValpalette[PValpositions]

  #For the heatmap
  heatmapColours <- colorRampPalette(c("dark blue", "blue", "cyan", "white", "pink", "red", "dark red"))(540)

  ## Create colour bars
  ### Create data for colour bar
  i <- nrow(unique(appmetaG[currentVar]))

  if(int == 1) {
    x <- seq(from = min(appmetaG[currentVar]), to = max(appmetaG[currentVar]), length.out = nColours)
    y <- 1
  } else {
    x <- 1:i
    y <- 1
  }

  xy <- data_frame(x, y)


  ### Create data for FDR bar
  a <- 1:101
  b <- 1
  c <- data_frame(a, b)

  ### Create data for FDR bar title
  texta <- 1
  textb <- 1
  text <- data_frame(texta, textb)

  ### Make data for variables ordered
  orderedVariables <- inputVals[sampleOrder]
  yAxis <- 1
  variableAnnotation <- data_frame(orderedVariables, yAxis)


  ### Organise title
  title <- if(int == 1) {
    paste("Top", nGenes, "Genes - Gene Expression Analysis of", currentVar, ":", round(min(appmetaG[currentVar]), 2), "to",
          round(max(appmetaG[currentVar]), 2))
  } else {
    paste("Analysis by", currentVar)
  }

  ### Plot colour bar
  # Sorry this is a big chunk of code, typically this code will design the aesthetics of the heatmap by calling back upon your set objects,
  # it should run on its own. If not you may need to back track to check your names and data dimensions.

  colourBarXaxis <- levels(geneExprCPM2ggplotPVal$Vineyard)
  colourBarXaxis <- as.factor(colourBarXaxis)
  colourBarDF <- data_frame(colourBarXaxis, orderedVariables)
  colourBarDF$colourBarXaxis <-  factor(colourBarDF$colourBarXaxis, levels = colourBarDF$colourBarXaxis)

  varBar <- ggplot(data = colourBarDF, aes(x = colourBarXaxis, y = orderedVariables, fill = orderedVariables)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_gradientn(colours = varPalette) +
    ggtitle(paste(currentVar)) +
    coord_cartesian(ylim = c(min(colourBarDF$orderedVariables),
                             (max(colourBarDF$orderedVariables)))) +
    guides(fill = FALSE) +
    theme_hc()

  colourBar <- ggplot(data = xy, aes(
    x = seq(from = 0, to = 18.99, length.out = 101),
    y = y, fill = x)) +
    geom_raster(stat = "identity", position = "identity") +
    scale_fill_gradientn(colours = varPalette) +
    ggtitle(paste("Gene Expression Analysis of", currentVar, "from",
                  round(min(appmetaG[[currentVar]]), 2), "to",
                  round(max(appmetaG[[currentVar]]), 2))) +
    guides(fill = FALSE) +
    expand_limits(FALSE) +
    theme_void()

  forpal <- seq(from = 0, to = 12, by = 0.1)

  # Define Pvals
  roundedPvals <- round(log10(geneExprPVal$P.Value), 1)

  # find interval
  palRange <- findInterval(range(-roundedPvals), forpal)

  pPalette <- colorRampPalette(c(rgb(0, 0, 0), rgb(0, 1, 0)))(length(forpal))[seq(palRange[2], palRange[1])]

  geneExprPVal <- as.tibble(geneExprPVal)
  geneExprPVal$GeneID <- factor(geneExprPVal$GeneID)
  geneExprCPM2ggplot <- as.tibble(geneExprCPM2ggplot)
  levels(geneExprPVal$GeneID) <- levels(geneExprCPM2ggplot$GeneID)

  PValbarWhole <- geneExprPVal %>%
    mutate(log10P = round(log10(P.Value), 1),
           P.Value = sprintf("%.2e", P.Value)) %>%
    ggplot(aes(x = currentVar, y = GeneID, fill = log10P, label = P.Value)) +
    geom_raster(stat = "identity", position = "identity") +
    scale_fill_gradientn(colours = pPalette) +
    guides(fill = FALSE) +
    theme_void()

  PValbarWhole <- suppressMessages(plotly::ggplotly(PValbarWhole,
                                                    hoverinfo = c(
                                                      "x", "y", "label"
                                                    ))
  )


  PValtitle <- ggplot(data = text, aes(x = texta, y = textb)) +
    geom_text(label = "P Values") +
    theme_void()

  ## Create interactive heatmap
  ### Prepare hide object used to remove labels and lines
  # This object will allow us to remove any excess information which will clutter our plot. In doing so, information regarding
  # the samples and genes will be accessible through the hover function. Turning on the sparklines in our heatmap will allow
  # us to keep track of what sample and vineyard it is from
  # (In the shiny app, I plan to add a statement that will allow the user to turn axes on and off)
  # Don't really need to do this for methylation, but I just made the code and it looked pretty
 hide <- list(
   title = "",
   zeroline = FALSE,
   showline = FALSE,
   showticklabels = FALSE,
   showgrid = FALSE
   )

 ### Create the interactive heatmap
 # The values are in log-cpm, which means that the 0 baseline equals 1cpm and any counts with less than 1 cpm are represented as negative.
 GeneExprheatmap <- ggplot(data = geneExprCPM2ggplot,
                           aes(x = Vineyard, y = GeneID, fill = scaledCPM)) +
   geom_raster(stat = "identity", position = "identity") +
   xlab("Vineyard") +
   scale_fill_gradientn(colours = heatmapColours) +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

 GeneExprheatmap <- GeneExprheatmap %>% ggplotly()


 GeneExprheatmapSubplot <- suppressMessages(subplot(varBar, PValtitle, GeneExprheatmap,
                                                    PValbarWhole,
                                                    nrows = 2, heights = c(0.1, 0.75),
                                                    widths = c(0.8, 0.1),
                                                    shareX = TRUE, shareY = TRUE, titleX = FALSE,
                                                    titleY = FALSE))

 layout(GeneExprheatmapSubplot, title = title, showlegend = TRUE,
        margin = list(l = 140,
                      r = 50,
                      b = 60,
                      t = 40),
        plot_bgcolor = "white")

}
