# VCF Variant prioritizer

variant_priority <- function(vcf_file, SNPEff= T, COSMIC=T, DEG= NULL, PASS=T, minAF=T, clean=T){
  require(vcfR) # Load the necessary packages
  require(stringr)
  require(tidyr)
  require(dplyr)

  # Load data
  vcf <- read.vcfR( vcf_file, verbose = FALSE)

  # Create a dataframe
  vcf_info = data.frame(vcf@fix)

  # Add a new column with SNPEff
  ann_extact <- function(data, k=11){
    split <- str_split(as.character(data), ";")
    ANN <- grep('ANN', split[[1]], value=TRUE)
    ANN <- str_remove(ANN, "ANN=")
    #ANN_SUB <- str_split(ANN, ";")
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

  # Apply the function
  if (SNPEff==T){
    vcf_info$SNPEFF <- lapply(vcf_info$INFO, function(x) ann_extact(x))
    vcf_info = tidyr::unnest(vcf_info, SNPEFF)
    vcf_info$EFFECT <- vcf_info$SNPEFF
    vcf_info = vcf_info %>% separate(SNPEFF, c("mutation", "mutation_type", "impact", "Gene",
                                               "ENSEMBL_ID", "TRANSCRIPT", "ENSEMBL_ID_1",
                                               "CODING","SOMETHING","CLINVAR","CODON_CHANGE"), "\\|")
    vcf_info = subset(vcf_info, select = -c(mutation,TRANSCRIPT,ENSEMBL_ID_1,CODING,SOMETHING))
    # Score the impact
    impactscore <- lapply(vcf_info$impact, function(x) impact_scoring(x))
    vcf_info$impactscore <- as.numeric(impactscore)
    }

  # Scoring COSMIC FILTER
  if (COSMIC==T){vcf_info$cosmicscore <- ifelse(grepl("COS", vcf_info$ID), 1, 0)}

  # Scoring PASS FILTER
  if (PASS==T){vcf_info$passscore <- ifelse(grepl("PASS", vcf_info$FILTER), 1, 0)}

  # Scoring mimimum allele frequency
  if (minAF==T){vcf_info$minaf <- ifelse(grepl("MinAF", vcf_info$FILTER), -1, 0)}

  # Scoring Differentially expressed genes
  if (!is.null(DEG)){vcf_info$degscore <- ifelse((vcf_info$Gene %in% DEG), 1, 0)}

  # Total score
  score_columns <- grep("score", colnames(vcf_info))
  vcf_info$totalscore = rowSums(vcf_info[,score_columns])

  # Order the final column
  vcf_info <- vcf_info[order(-vcf_info$totalscore),]

  # Collapse variants by the maximun score for gene
  vcf_info = vcf_info %>% group_by(Gene) %>% slice(which.max(totalscore))

  # Clean the dataframe
  col = colnames(vcf_info[,(grep("EFFECT", colnames(vcf_info)) + 1): ncol(vcf_info)])
  columns = c("Gene", "impact", "mutation_type", "ID", "CHROM", "POS", "REF", "ALT", "CLINVAR","CODON_CHANGE",col)
  if (clean == T){vcf_info = vcf_info[,columns]}

  # Return results
  return(vcf_info)

}
cnv_vcf_splitter <- function (vcf_file){

  require(vcfR) # Load the necessary packages
  require(stringr)
  require(tidyr)
  require(dplyr)

  # Load data
  vcf <- read.vcfR( vcf_file, verbose = FALSE)

  # Create a dataframe
  vcf_info = data.frame(vcf@fix)

  # Add a new column with SNPEff
  ann_extact <- function(data){
    split <- str_split(as.character(data), ";")
    s = sub("\\ANN.*", "", split[[1]])
    return(s)}

  # Apply the function
  vcf_info$Other <- lapply(vcf_info$INFO, function(x) ann_extact(x))
  keeps <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","Other")
  vcf_info <- vcf_info[keeps]
  vcf_info = tidyr::unnest_wider(vcf_info, Other)

  #return
  return(vcf_info)
}


# Load VCF file
vcf_file <- "C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/merged VCF files/two or more variants/0000.vcf"
vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/merged VCF files/two or more variants/0000.vcf"
vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/exome_analysis_trial_2_12172020/mutations/Da-mutect2-annotated.vcf"
vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/tmp/cnv_kristen/D-gatk-cnv.vcf"
vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/tmp/cnv_kristen/B-gatk-cnv.vcf"

# run vcf priority scoring
Da <- variant_priority(vcf_file, SNPEff= T, COSMIC=T, DEG= NULL, PASS=T, minAF=T, clean=T)
Da <- cnv_vcf_splitter (vcf_file)

# Combining multiple vcf files for oncoprint
oncoprint <- function(cleaned_vcf){

  combined_df <- data.frame()

  for (i in 1: length(cleaned_vcf)){
    columns <- c("sample","Gene", "impact", "mutation_type")
    temp_df <- cleaned_vcf[[i]]
    temp_df$sample <- i
    temp_df <- temp_df %>% select(columns)
    combined_df <- bind_rows(combined_df, temp_df)
  }

  # change 5 categories

  return(combined_df)

}
cleaned_vcf <- list(c1,c2,c3)
for_oncoprint <- oncoprint(cleaned_vcf)
write.csv(for_oncoprint, file = "Oncoprint.csv")

oncoprint <- function(cleaned_vcf){
  combined_df <- data.frame()
  for (i in 1: length(cleaned_vcf)){
    temp_df <- cleaned_vcf[[i]]
    temp_df$sample <- i
    combined_df <- bind_rows(combined_df, temp_df)
  }

  # change 5 categories

  return(combined_df)

}
cleaned_vcf <- list(c1,c2,c3)
for_oncoprint <- oncoprint(cleaned_vcf)
write.csv(for_oncoprint, file = "Oncoprint_all.csv")


# write
setwd('/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/')
write.csv(c1, 'c1_vsf.csv')
write.csv(c2, 'c2_vsf.csv')
write.csv(Da, '/Users/aj/Dropbox (Partners HealthCare)/Data/tmp/cnv_kristen/D_vcf_info.csv')

dim(c3)



# Load the differentially expressed genes
DEG = read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/RNAseq/Survival_timepoint/Survival_timepoint_clean/Duvelisib vs Vehicle/Differential expression/Duvelisib vs Vehicle.csv", header = T, row.names = 1)
DEG = read.csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/RNAseq/Survival_timepoint/Survival_timepoint_clean/Duvelisib vs Vehicle/Differential expression/Duvelisib vs Vehicle.csv", header = T, row.names = 1)

require(stringr)
DEG = row.names(DEG[DEG$padj <= 0.05,])
DEG = DEG[!str_detect(DEG,pattern="NA")]

# Load each VCF file (D1,D2,D3) independently
D1 <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/D1.vcf"
D2 <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/D2.vcf"
D3 <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/D3.vcf"

results_D1 = variant_priority(D3, DEG=DEG)

# Find the variants with >1.5
r = results_D1[(results_D1$impact == "HIGH" | results_D1$impact == "MODERATE") & results_D1$totalscore >= 2,]

# Save the result
write.csv(r, file="/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/D3_mutated_genes.csv")
