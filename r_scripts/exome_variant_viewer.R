# MAF viewer
# LIB
library(maftools)

#set wd
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/duvelisib_study/clinical samples/p1_exome/")

# Merge maf
merged <- maftools:::merge_mafs(maf=c('c1.maf','c2.maf','c3.maf'), verbose = TRUE)
# rename
merged@data$Tumor_Sample_Barcode=merged@data$Source_MAF
merged@data$Tumor_Sample_Barcode = as.factor(merged@data$Tumor_Sample_Barcode)
#merged@maf.silent$Tumor_Sample_Barcode=merged@data$Source_MAF
write.mafSummary(maf = merged, basename = 'merged')
# subset merged to included only passed filter
merged_pass = subsetMaf(merged, query = "FILTER == 'PASS'")
write.mafSummary(maf = merged_pass, basename = 'merged_pass')



# Loaddata
merged = read.maf(maf = 'merged_maftools.maf', vc_nonSyn=NULL)
merged_pass = read.maf(maf = 'merged_pass_maftools.maf', vc_nonSyn=NULL)

# Summary Plot
plotmafSummary(maf = merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes=T)
plotmafSummary(maf = merged_pass, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes=T)

# Onco plot
oncoplot(maf = merged, top = 30)
oncoplot(maf = merged_pass, top = 30, altered=T, showTitle=F)

oncostrip(maf = merged, top = 30, altered=T, showTitle=F)
oncostrip(maf = merged_pass, top = 50, altered=T, showTitle=F)






laml.titv = titv(maf = merged_subset, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)

lollipopPlot(maf = merged_subset, gene = 'TET2', AACol = 'Protein_position', showMutationRate = TRUE)

plotVaf(maf = merged_subset, vafCol = 'Allele')

somaticInteractions(maf = merged_subset, top = 25, pvalue = c(0.05, 0.1))

laml.sig = oncodrive(maf = merged, AACol = 'Protein_position', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.2, useFraction = TRUE)


