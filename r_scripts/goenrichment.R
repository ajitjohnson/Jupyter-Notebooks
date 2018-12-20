# Functio for GO enrichment analysis right after differential expression analysis by DESeq
# Function requires output fromDESeq with a column named "padj"

goenrichment <- function(x) {
	require(topGO)
	require(org.Hs.eg.db)
	allgenes <- toTable(org.Hs.egSYMBOL)
	# Subset genes with a significant P value (<0.05)
	sig_genes <- x[which(x$padj <= 0.05), ]
	listEG= merge(allgenes,sig_genes, by.x="symbol", by.y= "row.names")
	universeEG <- Lkeys(org.Hs.egGENENAME)
	geneList <- factor(as.integer(universeEG %in% listEG$gene_id))
	names(geneList) <- universeEG
	# Create GO object
	sampleGOdata <- new("topGOdata", ontology ="BP", allGenes= geneList, nodeSize= 5, annot = annFUN.org, 
						mapping = "org.Hs.eg.db", ID = "entrez")
	# Run enrichment analysis
	resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
	# Cummlulate into a table
	allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = 5)
	return(allRes)
	}