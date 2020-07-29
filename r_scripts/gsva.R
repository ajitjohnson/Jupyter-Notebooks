# GSVA setup

# Library
library(GSEABase)
library(GSVA)


# Convert dataframe colums to named list
named_list <- function(df){
  final_list <- c()
  for (i in 1: ncol(df)){
    tmp <- as.character(df[,i])
    tmp <- list(tmp[nchar(tmp) > 1])
    names(tmp) <- colnames(df[,i,drop=F])
    final_list <- c(final_list, tmp)
  }
  return(final_list)
}
x = named_list (custom_sig)

# Run the enrichment
gbm_es <- gsva(as.matrix(ndata), x,method='ssgsea',ssgsea.norm=T)

# Viz using heatmap
heatmap(gbm_es)