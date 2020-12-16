# GSEA
library(clusterProfiler)

# load data
deg <- read.csv(file= "ARSeq_pickseq/IB vs IT_ET/Differential expression/IB vs IT_ET.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq_pickseq/ET vs IT/Differential expression/ET vs IT.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq_pickseq/IT vs ET/Differential expression/IT vs ET.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq_pickseq/ET vs MIS/Differential expression/ET vs MIS.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq_pickseq/IT vs MIS/Differential expression/IT vs MIS.csv", header = T, check.names=F, row.names = 1)


# convert genes to ID
arseq.gsea.preprocess <- function(deg){
  # Create a ranked list of genes from deg
  ranked_list <- data.frame(deg[!is.na(deg$padj),])
  ranked_list$score <- -log(ranked_list$padj)*sign(ranked_list$log2FoldChange)
  ranked_list <- ranked_list[order(-ranked_list$score),]
  r_list <- ranked_list$score
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  ranked.list = r_list[is.finite(r_list)]
  return(ranked.list)
} # adj p
arseq.gsea.preprocess <- function(deg){
  # Create a ranked list of genes from deg
  ranked_list <- data.frame(deg[!is.na(deg$pvalue),])
  ranked_list$score <- -log(ranked_list$pvalue)*sign(ranked_list$log2FoldChange)
  ranked_list <- ranked_list[order(-ranked_list$score),]
  r_list <- ranked_list$score
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  ranked.list = r_list[is.finite(r_list)]
  return(ranked.list)
} # just P
ranked_deg <- arseq.gsea.preprocess(deg)

# Process signature
custom_sig <- read.csv(file= "/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/signatures/all_sig.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/nfkb.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/hallmark.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/apoptosis_supression.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/EMT_all.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/caf.csv", row.names=1, header = F)
custom_sig <- read.csv(file= "signatures/SENESCENCE.csv", row.names=1, header = F)



custom_sig <- custom_sig[,-1] # new method add this line
h <- data.frame(t(custom_sig))
hh <- NA
for (i in 1: ncol(h)){
  print(i)
  tmp <- h[,i,drop=F]
  tmp <- tmp[!apply(is.na(tmp) | tmp == "", 1, all),,drop=F]
  colnames(tmp) <- c("gene")
  tmp$ont <- names(h[i])
  hh <- rbind(hh, tmp)
}
hh <- hh[-1,]
h <- hh[,c(2,1)]

# Run GSEA
egmt2 <- GSEA(ranked_deg, TERM2GENE=h, verbose=FALSE,pvalueCutoff = 0.5)
egmt2$ID
geneset = 1
gseaplot(egmt2, geneSetID = geneset, by = "runningScore", color ="#000000", color.line = "red", title=egmt2$ID[geneset])
s# Volcano