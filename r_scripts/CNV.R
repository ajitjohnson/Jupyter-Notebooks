


variant_priority <- function(vcf_file, SNPEff= T, COSMIC=T, DEG= NULL, PASS=T, minAF=T, clean=T){
  require(vcfR) # Load the necessary packages
  require(stringr)
  require(tidyr)
  require(dplyr)
  require(biomaRt)
  
  # Load data
  vcf <- read.vcfR( vcf_file, verbose = FALSE)
  
  # Create a dataframe
  vcf_info = data.frame(vcf@fix)
  w = subset(vcf_info_gene, select = -c(INFO))
  
  # Add a new column with SNPEff
  ann_extact <- function(data, k=5){
    split <- str_split(as.character(data), ";")
    ANN <- grep('ANN', split[[1]], value=TRUE)
    ANN <- str_remove(ANN, "ANN=")
    ANN_SUB <- str_split(ANN, ",")
    ANN_SUB <- do.call(paste, c(read.table(text = ANN_SUB[[1]], sep = "|")[1:k], sep = "|"))
    return(ANN_SUB)}
  # SNPEff Impact scoring
  impact_scoring <- function(data){
    score = c()
    if(data=="HIGH"){score = c(score, 1)}
    else if(data=="MODERATE"){score = c(score, 0.5)}
    else if(data=="LOW"){score = c(score, 0.25)}
    else if(data=="MODIFIER"){score = c(score, 0.125)}
    else {score = c(score, NA)}
    return(score)
  }
  
  # Extracting foldchange
  fc_extact <- function(data, pattern, numeric = TRUE){
    split <- str_split(as.character(data), ";")
    ANN <- grep(pattern, split[[1]], value=TRUE)
    ANN <- str_remove(ANN, pattern)
    if (numeric == TRUE){ANN <- as.numeric(ANN)}
    return(ANN)}
  
  # Add FC
  vcf_info$FC <- mapply(fc_extact, vcf_info$INFO, pattern='FOLD_CHANGE_LOG=')
  vcf_info$PROBES <- mapply(fc_extact, vcf_info$INFO, pattern='PROBES=')
  vcf_info$SVTYPE <- mapply(fc_extact, vcf_info$INFO, pattern='SVTYPE=', numeric=FALSE)
  vcf_info$SVLEN <- mapply(fc_extact, vcf_info$INFO, pattern='SVLEN=')
  vcf_info$END <- mapply(fc_extact, vcf_info$INFO, pattern='END=')
  vcf_info$CN <- mapply(fc_extact, vcf_info$INFO, pattern='CN=')
  
  # Apply the function
  if (SNPEff==T){
    vcf_info$SNPEFF <- lapply(vcf_info$INFO, function(x) ann_extact(x))
    vcf_info = tidyr::unnest(vcf_info, SNPEFF)
    vcf_info$EFFECT <- vcf_info$SNPEFF
    vcf_info = vcf_info %>% separate(SNPEFF, c("mutation", "mutation_type", "impact", "Gene","ENSEMBL_ID"), "\\|")
    x <- gsub("&.*","",vcf_info$ENSEMBL_ID) # remove anything after &
    x <- gsub("-.*","",x) # remove anything after -
    vcf_info$ENSEMBL_ID <- x
    vcf_info = subset(vcf_info, select = -c(mutation))
  }
  
  # Convert ensembl ID to genename
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)
  vcf_info_gene <- merge(vcf_info, genes, by.x="ENSEMBL_ID", by.y= "ensembl_gene_id")[,-1]

  length(unique(vcf_info_gene$hgnc_symbol))
  
  geneA <- unique(vcf_info_gene$hgnc_symbol)
  geneB <- unique(vcf_info_gene$hgnc_symbol)
  
  length(intersect(geneA,geneB))
  
  #
  v = vcf_info[vcf_info$impact == 'HIGH',]
  
  # Collapse variants by the maximun score for gene
  vcf_info = vcf_info %>% group_by(Gene) %>% slice(which.max(totalscore))
  
  
  # Clean the dataframe
  col = colnames(vcf_info[,(grep("EFFECT", colnames(vcf_info)) + 1): ncol(vcf_info)])
  columns = c("Gene", "impact", "mutation_type", "ID", "CHROM", "POS", "REF", "ALT", col)
  if (clean == T){vcf_info = vcf_info[,columns]}
  
  # Return results
  return(vcf_info)
  
}





# Load data
vcf_file <- "~/Desktop/cnv/normal/MSK70/Tumor70-gatk-cnv.vcf"
vcf_file <- "~/Desktop/cnv/normal/MSK76/Tumor76-gatk-cnv.vcf"
vcf_file <- "~/Desktop/cnv/tumor/MSK76/Tumor76-gatk-cnv.vcf"
