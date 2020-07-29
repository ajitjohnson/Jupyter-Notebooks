# Salmon to counts table

# Set WD
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-5-f12-f9/f9/salmon/")
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-5-f12-f9/f12/salmon/")

# Library
library("tximport")
library("tximportData")
library("EnsDb.Hsapiens.v86")
library("stringr")

# Load data
#files <- file.path(getwd(), list.files(getwd()))
files <- file.path(getwd(), list.files(pattern = ".sf", recursive = TRUE))
  
  

# Prepare to convert to gene
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id", "gene_id")])
tx2gene <- tx2gene[complete.cases(tx2gene),]

# Run trimport
txi.tx <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE) # transcript level abundance
head(txi.tx$counts)

# add column name
file_names <- list.files(pattern = ".sf", recursive = TRUE)
file_names <- str_remove(file_names, "/quant.sf")

# add columname
data <- txi.tx$counts
colnames(data) <- file_names
data <- round(data, digits = 0)

# Write the data out
write.csv(data, file = "counts.csv")

# Convert transcripts to genes
data <- arseq.ensembl2genename (data)

