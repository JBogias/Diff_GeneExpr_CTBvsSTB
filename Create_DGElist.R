library(dplyr)


# Define the sub-directory in which you can find sorted bam counts files (the output of featureCounts)
subdir <- "~/Path/to/files"

#Run through each of the count files. You shouldn't need to change the pattern as ".counts" will most likely be in
# the suffix of the file.
countFiles <- list.files(
  subdir,
  pattern=".counts",
  full.names = TRUE
)

# A little nifty way to make a data frame from all of your files
# This is just a function that will create the data frame of your data
counts <- countFiles %>%
  lapply(function(x){
    read.delim(x, 
               header = FALSE,
               row.names = 1) %>%
      set_names(x)}) %>%
  as.data.frame()

#Dont worry about this bad code, it was me being a fool. I wanted to easily make the new column names for my data frame
#new_colnames <- counts[2,]

# So when you have your counts dataframe, there may be a few rows that do not contain any numerical data. These
# rows can be deleted by using `-c()` and stating the row numbers you want to remove (eg. 1:5 means "1 to 5")
# This is just an example
counts_df <- counts[-c(1:2),]


# Manually make the new colnames, if theres and easier way, I'm all ears
# You'll typically make these colnames to give meaningful titles to your data so you can identify which sample it is
# Below is an example from my own work
manual_colnames <- c("FemCTB_Chr", "FemCTB_Start", "FemCTB_End", "FemCTB_Strand", "FemCTB_Length", "FemCTB_Counts",
                     "FemSTB_Chr", "FemSTB_Start", "FemSTB_End", "FemSTB_Strand", "FemSTB_Length", "FemSTB_Counts",
                     "MCTB_Chr", "MCTB_Start", "MCTB_End", "MCTB_Strand", "MCTB_Length", "MCTB_Counts",
                     "MSTB_Chr", "MSTB_Start", "MSTB_End", "MSTB_Strand", "MSTB_Length", "MSTB_Counts")

# Then when you've made it, add in the new colnames
colnames(counts_df) <- manual_colnames

# So you may have some extra information such as chromosome number, start and end position, length and strand information.
# You will notice in my `manual_colnames` I have specific column names for count data. For a DGElist to be created
# you need a counts matrix with just the numerical counts in then (trust me, I got errors when I tried the other way)
# You can just use `dplyr::select()` and target "Counts" as shown below
# Note it's case senseitive to the column name I set previously
# This code selects for the columns which contain the gene counts, the DGElist only wants gene count data, doesn't like anything non-numeric
counts_matrix <- dplyr::select(counts_df, ends_with("Counts"))

# Before you turn your data into a DGElist, you'll need to change the inputs,
# from factors to numeric as shown below
indx <- sapply(counts_matrix, is.factor)
counts_matrix[indx] <- lapply(counts_matrix[indx], function(x) as.numeric(as.character(x)))

# Now create the countsDGE simply by running `DGEList` on your counts object!
countsDGE <- DGEList(counts_matrix)