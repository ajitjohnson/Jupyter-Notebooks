# Function to collapse multiple transcripts into one gene
# To use name the column containing gene names as "gene"

genesummary <- function(x) {
	x$gene <- as.factor(x$gene)
	y=aggregate(.~gene, data=x, max)
	rownames(y)=y[,1]
	y=y[,-1]
	return(y)
}