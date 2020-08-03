# COMBAT

# Set WD
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/Figures for manuscript/Fig 6/")

# Library
library(sva)

# Load data
d1 <-  read.csv(file= "f8_counts.csv", row.names = 1, header = T, check.names = F)
d2 <-  read.csv(file= "f12_counts.csv", row.names = 1, header = T, check.names = F)


# Counts COMBAT
common_genes <- intersect(row.names(d1), row.names(d2))
d1_c <- d1 [row.names(d1) %in% common_genes,]
d2_c <- d2 [row.names(d2) %in% common_genes,]
combined_1 <- cbind(d1_c,d2_c)
batch_1 <- c(rep(1, ncol(d1_c)), rep(2, ncol(d2_c)))
adjusted <- ComBat_seq(as.matrix(combined_1), batch=batch_1, group=NULL)
adjusted = data.frame(adjusted)

plotba (log1p(combined_1),log1p(adjusted))

write.csv(adjusted, file='combat_data.csv')


# Normal COMBAT
common_genes <- intersect(row.names(d1_n), row.names(d3_n))
d1_c <- d1_n [row.names(d1_n) %in% common_genes,]
d3_c <- d3_n [row.names(d3_n) %in% common_genes,]
combined_1 <- cbind(d1_c,d3_c)
batch_1 <- c(rep(1, ncol(d1_c)), rep(2, ncol(d3_c)))
# combat
combat_edata1 = ComBat(dat=combined_1, batch=batch_1, mod=NULL, par.prior=TRUE, prior.plots=T)