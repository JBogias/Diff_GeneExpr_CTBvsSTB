library(sleuth)
library(cowplot)
library(biomaRt)
library(tximport)

files <- file.path(dir, "kallisto_boot", samples$run, "abundance.h5")
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)


BeWometa <- readRDS("~/Documents/PhD/R_Projects/Differential_GeneExpr_CTBvsSTB/Data/BeWo_meta.rds")


mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

ens_and_GO_ID_BeWo <- getBM(attributes = c("ensembl_gene_id", "go_id", "external_gene_name", "external_transcript_name", "description", "name_1006", "definition_1006"), filters = "ensembl_gene_id",
                            values = rownames(topGenes$table), mart = mart)
