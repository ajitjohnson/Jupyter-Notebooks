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
  ann_extact <- function(data, k=4){
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

  # Apply the function
  if (SNPEff==T){
    vcf_info$SNPEFF <- lapply(vcf_info$INFO, function(x) ann_extact(x))
    vcf_info = tidyr::unnest(vcf_info, SNPEFF)
    vcf_info$EFFECT <- vcf_info$SNPEFF
    vcf_info = vcf_info %>% separate(SNPEFF, c("mutation", "mutation_type", "impact", "Gene"), "\\|")
    vcf_info = subset(vcf_info, select = -c(mutation))
    # Score the impact
    impactscore <- lapply(vcf_info$impact, function(x) impact_scoring(x))
    vcf_info$impactscore <- as.numeric(impactscore)
    }

  # Scoring COSMIC FILTER
  if (COSMIC==T){vcf_info$cosmicscore <- ifelse(grepl("COS", vcf_info$ID), 1, 0)}

  # Scoring PASS FILTER
  if (PASS==T){vcf_info$passscore <- ifelse(grepl("PASS", vcf_info$FILTER), 1, 0)}

  # Scoring mimimum allele frequency
  if (PASS==T){vcf_info$minaf <- ifelse(grepl("MinAF", vcf_info$FILTER), 1, 0)}

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
  columns = c("Gene", "impact", "mutation_type", "ID", "CHROM", "POS", "REF", "ALT", col)
  if (clean == T){vcf_info = vcf_info[,columns]}

  # Return results
  return(vcf_info)

}


# Load VCF file
vcf_file <- "C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/merged VCF files/two or more variants/0000.vcf"
vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/merged VCF files/two or more variants/0000.vcf"
# Variants in all three samples
vcf_file <- "C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/duvelisib_study/pdx/exome/vechile_merged/merged VCF files/variants in all three/0000.vcf"

vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/1-mutect2-annotated.vcf"
results_C1 = variant_priority(vcf_file)
c1 = results_C1[(results_C1$passscore == 1),]

vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/2-mutect2-annotated.vcf"
results_C2 = variant_priority(vcf_file)
c2 = results_C2[(results_C2$passscore == 1),]

vcf_file <- "/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/3-mutect2-annotated.vcf"
results_C3 = variant_priority(vcf_file)
c3 = results_C3[(results_C3$passscore == 1),]

# Drop the low impact variants
c1 = c1[(c1$impact == "HIGH" | c1$impact == "MODERATE"),]
c2 = c2[(c2$impact == "HIGH" | c2$impact == "MODERATE"),]
c3 = c3[(c3$impact == "HIGH" | c3$impact == "MODERATE"),]

# Figire for distribution of variants
type <- data.frame(table(c3$mutation_type))

ggplot(data=type, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity")+ coord_flip() + theme_minimal()+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(face="bold",size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  annotate(geom="text", x=5, y=20, label=paste(sum(type$Freq),'variants'), fontface="bold")


MISSENSE: a missense mutation
INFRAME: a inframe mutation
TRUNC: a truncation mutation
PROMOTER: a promoter mutation
OTHER: any other kind of mutation

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



setwd('/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/')

write.csv(c1, 'c1_vsf.csv')
write.csv(c2, 'c2_vsf.csv')
write.csv(c3, 'c3_vsf.csv')

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
