# Salmon to counts table

# Set WD
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-5-f12-f9/f9/salmon/")
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-5-f12-f9/f12/salmon/")
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-1/salmon/")
setwd("//research.files.med.harvard.edu/ImStor/sorger/data/rnaseq/ajit_johnson/covid/alignment/final")
setwd("D:/gw/final")
setwd("C:/Users/ajn16/Desktop/cell_line/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/for_others/Ran/pickseq_scrnaseq/15july2021_1st_trial/ran and lydie samples/kallisto")





# Library
library("tximport")
library("tximportData")
library("EnsDb.Hsapiens.v86")
library("stringr")
library("rhdf5")


# Load data
#files <- file.path(getwd(), list.files(getwd()))
files <- file.path(getwd(), list.files(pattern = ".sf", recursive = TRUE))
files <- file.path(getwd(), list.files(pattern = ".h5", recursive = TRUE))

# If needed
tx2gene <- read.csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/iris/100_celllines/tx2gene.csv", header = F)
files <- gsub("\\\\","//",files)
files <- substring(files, 3)


# Prepare to convert to gene
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id", "gene_id")])
tx2gene <- tx2gene[complete.cases(tx2gene),]

# Run trimport
txi.tx <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE) # transcript level abundance
txi.tx <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE) # transcript level abundance

head(txi.tx$counts)

# add column name
file_names <- list.files(pattern = ".sf", recursive = TRUE)
file_names <- list.files(pattern = ".h5", recursive = TRUE)
file_names <- str_remove(file_names, "/quant.sf")
file_names <- str_remove(file_names, "/salmon/quant.sf")
file_names <- str_remove(file_names, "/quant/abundance.h5")

# add columname
data <- txi.tx$counts
colnames(data) <- file_names
data <- round(data, digits = 0)

# Write the data out
write.csv(data, file = "counts_tpm.csv")
write.csv(data, file = "kallisto_tpm.csv")

# Convert transcripts to genes
d <- arseq.ensembl2genename (data)
write.csv(d, file = "kallisto_tpm_genenames.csv")

setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/for_others/Ran/pickseq_scrnaseq/15july2021_1st_trial/ran and lydie samples/")

# Merge a lot of CSV files
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP/pickseq/fastq processing/viral/htseq-count/")
library(dplyr)
library(readr)
library(tidyverse)
d <- list.files(path=getwd(), full.names = TRUE) %>% lapply(read.table) %>% bind_cols
row.names(d) <- d[,1]
c <- 1:ncol(d)
d <- d[, c%%2==0]
# rename the columns
file_names <- list.files(path=getwd(), full.names = TRUE)
file_names <- substring(file_names, 101)
file_names <- str_remove(file_names, ".counts")
colnames(d) <- file_names

write.csv(d, file = "counts.csv")

# Functions to convert to gene names
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/for_others/Ran/pickseq_scrnaseq/15july2021_1st_trial/ran and lydie samples/")
data <- read.csv("kallisto_tpm.csv", header = T, row.names = 1)

library("biomaRt")
arseq.ensembl2genename <- function(data,ensemblmirror=NULL){
  print("Converting ENSEMBL ID's to gene names")
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror=ensemblmirror)
  genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)
  data_m <- merge(data, genes, by.x="row.names", by.y= "ensembl_gene_id")[,-1]
  data_m <- data_m[!(data_m$hgnc_symbol==""), ]
  # Reduce multiple transcripts into a single gene by taking the sum of related transcripts
  data_m$hgnc_symbol <- as.factor(data_m$hgnc_symbol)
  data <- aggregate(.~hgnc_symbol, data=data_m, sum)
  rownames(data) <- data[,1]
  data <- data[,-1]
  return(data)
}
