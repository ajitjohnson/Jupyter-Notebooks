# GSEA analysis of COVID samples

# Lib
library(clusterProfiler)
library("msigdbr")
library(enrichplot)
#library(DOSE)
library(ReactomePA)
library(tidyr)

# Set WD
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP/archive/ARSeq_7/")
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/proteomics/")

# Load the differentially expressed genes
deg <- read.csv(file= "covid_d2 vs uninfected/Differential expression/covid_d2 vs uninfected.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "covid_high_viral vs covid_d2/Differential expression/covid_high_viral vs covid_d2.csv", header = T, check.names=F, row.names = 1)
# Human
deg <- read.csv(file= "ARSeq/Hu_p1 vs Hu_c/Differential expression/Hu_p1 vs Hu_c.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq/Hu_p2 vs Hu_c/Differential expression/Hu_p2 vs Hu_c.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq/Hu_p3 vs Hu_c/Differential expression/Hu_p3 vs Hu_c.csv", header = T, check.names=F, row.names = 1)
deg <- read.csv(file= "ARSeq/Patient_2 vs control_Patient_1_Patient_3/Differential expression/Patient_2 vs control_Patient_1_Patient_3.csv", header = T, check.names=F, row.names = 1)

# Load custom GSEA signature
h  <- read.csv(file= "/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/signatures/lm22.csv", header = T)
h  <- read.csv(file= "/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/signatures/all_sig.csv", header=F, row.names = 1)
h <- data.frame(t(h))

custom_sig <- read.csv(file= "/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/signatures/all_sig.csv", row.names=1, header = F)
custom_sig <- data.frame(t(custom_sig))


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

# Order the genes for passing it into GSEA
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
}
ranked_deg <- arseq.gsea.preprocess(deg)

# corvert gene name to enterzid
gene2entrez <- function(genes){
  require(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  genes <- names(genes)
  ENT <- select (hs, keys = genes,columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  ENT <- ENT[complete.cases(ENT), ]
  entid <- ENT$ENTREZID
  return(entid)
}
names(ranked_deg)  <- gene2entrez(ranked_deg)

# mSigDB list
h = msigdbr(species = "Homo sapiens", category = "H")
h = h[,c('gs_name','gene_symbol')]
colnames(h) <- c("ont","gene")

# Enrichment analysis
#egmt <- enricher(names(ranked_deg), TERM2GENE=h)
#head(egmt)
egmt2 <- GSEA(ranked_deg, TERM2GENE=h, verbose=FALSE,pvalueCutoff = 0.05)
egmt2$ID

# plot
#gseaplot2(egmt2, geneSetID = 3:4, pvalue_table = TRUE, color = c("#E495A5", "#86B875"))
geneset =7 # # 5, 9 # 21, 10 # 13, 7
gseaplot(egmt2, geneSetID = geneset, by = "runningScore", color ="#000000", color.line = "red", title=egmt2$ID[geneset])

#4X4

# Gene-Concept Network
go <- enrichGO(names(ranked_deg),'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
goplot(go)

# Reactome pathway analysis
# Subset genes that are downregulated
deg <- deg[deg$padj <= 0.05, ]
deg <- deg[complete.cases(deg),]
arseq.gsea.preprocess.mod <- function(deg){
  require(org.Hs.eg.db)
  # Create a ranked list of genes from deg
  ranked_list <- data.frame(deg[!is.na(deg$padj),])
  r_list <- ranked_list$log2FoldChange
  names(r_list) <- row.names(ranked_list)
  # Remove any infinity values
  ranked.list = r_list[is.finite(r_list)]
  return(ranked.list)
}
ranked_deg <- arseq.gsea.preprocess.mod(deg)
names(ranked_deg)  <- gene2entrez(ranked_deg)

x <- enrichPathway(gene=names(ranked_deg),pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory=25)
emapplot(x)
cnetplot(x, categorySize="pvalue", foldChange=ranked_deg)


