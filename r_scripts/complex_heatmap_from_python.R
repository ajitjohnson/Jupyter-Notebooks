# Created on Tue Feb 25 11:43:59 2020
# @author: Ajit Johnson Nirmal
# PTCL TMA restaining analysis

# Library
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/AITL/tumor heterogenity/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/ALCL/tumor heterogenity/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/ALCL/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/ALCL/tumor_heterogenity/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/AITL/cell_phenotyping/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/AITL/tumor_heterogenity/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/NOS/phenotype/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/NOS/tumor_het/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/NOS/fht_het/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/ALI/tumor_heterogenity/")
setwd("/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/for R/all_whole_slide/stroma/")

# tumor heterogenity plot
# Import data
median_expression <- read.csv("median_expression.csv", row.names = 1, header = T)
median_cellfeatures <- read.csv("median_cellfeatures.csv", row.names = 1, header = T)
pheno_prop <- read.csv("pheno_prop.csv", row.names = 1, header = T, check.names = F)
pheno_group <- read.csv("pheno_group.csv", row.names = 1, header = T, check.names = F)

# rename phenogroup colum
colnames(pheno_group) <- 'group'

# remove stroma
remove_pheno = c('Stroma')
pheno_group <- pheno_group[!pheno_group$group %in% remove_pheno,,drop=F]
median_expression <- median_expression[row.names(pheno_group),]
median_cellfeatures <-  median_cellfeatures[row.names(pheno_group),]
pheno_prop <-  pheno_prop[row.names(pheno_group),]

# scale data if needed
#median_expression <- scale(median_expression)
#median_expression <- t(scale(t(median_expression)))
#median_cellfeatures <- t(scale(t(median_cellfeatures)))

h1_col = colorRamp2(c(0.20, 0.4, 0.5, 0.6, 0.8), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
h2_col = colorRamp2(c(0.20, 0.4, 0.5, 0.6, 0.8), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B"))
#h1_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
#h2_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B"))

cols <- colorRampPalette(brewer.pal(length(colnames(pheno_prop)), "Spectral"))(length(colnames(pheno_prop))) # extend colors
pheno_p = rowAnnotation(axis_reverse = anno_barplot(pheno_prop,
                                                   border = F,
                                                   gp = gpar(fill = cols),
                                                   width = unit(4, "cm"),
                                                   bar_width = 0.8,
                                                   axis_param = list(direction = "reverse"),
                                                   labels_gp = gpar(col = "white", fontsize = 10)
                                                   ))
pheno_lgd = Legend(title = "barplot", labels = colnames(pheno_prop), legend_gp = gpar(fill = cols))

h1 = Heatmap(as.matrix(median_expression), col = h1_col,
             left_annotation = pheno_p,
             show_row_dend = FALSE, name='scale',
             row_title_rot = 0,
             #row_split = factor(pheno_group$group, levels = unique(pheno_group$group)),
             row_title_gp = gpar(fontsize = 8, fontface = "bold")
             )
h2 = Heatmap(as.matrix(median_cellfeatures), col = h1_col,name='scale')

ht_list = h1 + h2

draw(ht_list, heatmap_legend_list = list(pheno_lgd)
     )


# sample clustering -------------------------------------------------------

# WD
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/AITL/patient_clustering/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/ALCL/patient_clustering/")

# import data
phenotype_prop <- read.csv("phenotype_proportion.csv", row.names = 1, header = T, check.names = F)
phenotype_prop = scale(phenotype_prop)
# heatmap
Heatmap(as.matrix(phenotype_prop))

# Aftersubgrouping heatmap
phenotype_group <- read.csv("phenotype_grouping.csv", row.names = 1, header = T, check.names = F)
# heatmap
Heatmap(as.matrix(phenotype_prop),row_title_rot = 0,
        row_split = factor(phenotype_group$group, levels = unique(phenotype_group$group)))


# sample clustering heatmap with markers ----------------------------------

setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/AITL/patient_clustering_heatmap/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/ALCL/patient_clustering_heatmap/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/NOS/immune_prop/")

# Import data
median_expression <- read.csv("median_expression.csv", row.names = 1, header = T)
median_cellfeatures <- read.csv("median_cellfeatures.csv", row.names = 1, header = T)
pheno_prop <- read.csv("pheno_prop.csv", row.names = 1, header = T, check.names = F)
pheno_group <- read.csv("pheno_group.csv", row.names = 1, header = T, check.names = F)

# rename phenogroup colum
colnames(pheno_group) <- 'group'

# scale data if needed
#median_expression <- scale(median_expression)
median_expression <- t(scale(t(median_expression)))
median_cellfeatures <- t(scale(t(median_cellfeatures)))

# remove stroma
remove_pheno = c('Stroma')
pheno_group <- pheno_group[!pheno_group$group %in% remove_pheno,,drop=F]
median_expression <- median_expression[row.names(pheno_group),]
median_cellfeatures <-  median_cellfeatures[row.names(pheno_group),]
pheno_prop <-  pheno_prop[row.names(pheno_group),]

#h1_col = colorRamp2(c(0.20, 0.4, 0.5, 0.6, 0.8), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
#h2_col = colorRamp2(c(0.25, 0.25, 0.5, 0.75, 1), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B"))
h1_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
h2_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B"))


cols <- colorRampPalette(brewer.pal(length(colnames(pheno_prop)), "Spectral"))(length(colnames(pheno_prop))) # extend colors
pheno_p = rowAnnotation(axis_reverse = anno_barplot(pheno_prop,
                                                    border = F,
                                                    gp = gpar(fill = cols),
                                                    width = unit(4, "cm"),
                                                    bar_width = 0.8,
                                                    axis_param = list(direction = "reverse"),
                                                    labels_gp = gpar(col = "white", fontsize = 10)
))
pheno_lgd = Legend(title = "barplot", labels = colnames(pheno_prop), legend_gp = gpar(fill = cols))

h1 = Heatmap(as.matrix(median_expression), col = h1_col,
             left_annotation = pheno_p,
             show_row_dend = FALSE, name='scale',
             row_title_rot = 0,
             row_title_gp = gpar(fontsize = 8, fontface = "bold"),
             row_split = factor(pheno_group$group, levels = unique(pheno_group$group)))
h2 = Heatmap(as.matrix(median_cellfeatures), name='scale', col = h2_col) # name='scale', col = h2_col

ht_list = h1 + h2

draw(ht_list, heatmap_legend_list = list(pheno_lgd))



# Interaction map ---------------------------------------------------------

setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/AITL/tumor heterogenity/")
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA/for R/ALCL/tumor_interaction/")


# Import data
median_expression <- read.csv("interaction_data.csv", row.names = 1, header = T, check.names = F)
median_cellfeatures <- read.csv("median_cellfeatures.csv", row.names = 1, header = T)
pheno_prop <- read.csv("pheno_prop.csv", row.names = 1, header = T, check.names = F)
pheno_group <- read.csv("pheno_group.csv", row.names = 1, header = T, check.names = F)

# rename phenogroup colum
colnames(pheno_group) <- 'group'

# remove rows
remove_pheno = c('Unknown','Immune cells','Non-immune cells','T cells')
pheno_group <- pheno_group[!pheno_group$group %in% remove_pheno,,drop=F]
median_expression <- median_expression[row.names(pheno_group),]
median_cellfeatures <-  median_cellfeatures[row.names(pheno_group),]
pheno_prop <-  pheno_prop[row.names(pheno_group),]

# remov columsn in interaction
remove_pheno = c('Unknown','Immune cells','Non-immune cells','T cells')
median_expression <- median_expression[ , -which(names(median_expression) %in% remove_pheno)]

# scale data if needed
median_cellfeatures <- t(scale(t(median_cellfeatures)))

#h1_col = colorRamp2(c(0.20, 0.4, 0.5, 0.6, 0.8), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
#h2_col = colorRamp2(c(0.25, 0.25, 0.5, 0.75, 1), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B"))
h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8
h2_col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#4B6AAF",  '#55B0AE', "#FFFFFF","#F5A15B","#B11E4B")) #F8F6B8

cols <- colorRampPalette(brewer.pal(length(colnames(pheno_prop)), "Spectral"))(length(colnames(pheno_prop))) # extend colors
pheno_p = rowAnnotation(axis_reverse = anno_barplot(pheno_prop,
                                                    border = F,
                                                    gp = gpar(fill = cols),
                                                    width = unit(4, "cm"),
                                                    bar_width = 0.8,
                                                    axis_param = list(direction = "reverse"),
                                                    labels_gp = gpar(col = "white", fontsize = 10)
))
pheno_lgd = Legend(title = "barplot", labels = colnames(pheno_prop), legend_gp = gpar(fill = cols))

h1 = Heatmap(as.matrix(median_expression), col = h1_col,
             left_annotation = pheno_p,
             show_row_dend = FALSE, name='scale',
             row_title_rot = 0,border='black',
             rect_gp = gpar(col = "#22223b", lwd = 1),
             row_title_gp = gpar(fontsize = 8, fontface = "bold"),
             row_split = factor(pheno_group$group, levels = unique(pheno_group$group)))
h2 = Heatmap(as.matrix(median_cellfeatures), col = h2_col,name='cell_features')

ht_list = h1 + h2

draw(ht_list, heatmap_legend_list = list(pheno_lgd))


# Stromal patient grouping ------------------------------------------------\
setwd("/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/for R/AITL/stromal_grouping/")
pheno_prop <- read.csv("immune_prop.csv", row.names = 1, header = T, check.names = F)
pheno_group <- read.csv("pheno_group.csv", row.names = 1, header = T, check.names = F)

# rename phenogroup colum
colnames(pheno_group) <- 'group'

# scale
#pheno_prop <- t(scale(t(pheno_prop)))
pheno_prop <- scale(pheno_prop)

# Heaatmap
Heatmap(as.matrix(pheno_prop),
        show_row_dend = T,
        row_title_rot = 0,border='black',
        rect_gp = gpar(col = "#22223b", lwd = 1),
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        #row_split = factor(pheno_group$group, levels = unique(pheno_group$group))
        )





