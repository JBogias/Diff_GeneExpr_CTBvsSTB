subdir <- "~/Documents/PhD/Data/counts"

countFiles <- list.files(
  subdir,
  pattern=".counts",
  full.names = TRUE
)

counts <- countFiles %>%
  lapply(function(x){
    read.delim(x, 
               header = FALSE,
               row.names = 1) %>%
      set_names(x)}) %>%
  as.data.frame()

#new_colnames <- counts[2,]

counts_df <- counts[-c(1:2),]

manual_colnames <- c("FemCTB_Chr", "FemCTB_Start", "FemCTB_End", "FemCTB_Strand", "FemCTB_Length", "FemCTB_Counts",
                     "FemSTB_Chr", "FemSTB_Start", "FemSTB_End", "FemSTB_Strand", "FemSTB_Length", "FemSTB_Counts",
                     "MCTB_Chr", "MCTB_Start", "MCTB_End", "MCTB_Strand", "MCTB_Length", "MCTB_Counts",
                     "MSTB_Chr", "MSTB_Start", "MSTB_End", "MSTB_Strand", "MSTB_Length", "MSTB_Counts")

colnames(counts_df) <- manual_colnames

counts_matrix <- dplyr::select(counts_df, ends_with("Counts"))

countsDGE <- DGEList(counts_matrix)

indx <- sapply(counts_matrix, is.factor)
counts_matrix[indx] <- lapply(counts_matrix[indx], function(x) as.numeric(as.character(x)))

