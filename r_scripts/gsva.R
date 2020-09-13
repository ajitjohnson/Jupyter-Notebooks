# GSVA setup


# load signature
custom_sig <- read.csv(file= "/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/signatures/stat3_sig.csv", row.names=1, header = F)
custom_sig <- data.frame(t(custom_sig))

# Load data
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP")
data <- read.csv(file= "exp_cleaned.csv", header = T, check.names=F,row.names = 1) #row.names = 1
meta <- read.csv(file= "meta_cleaned.csv", row.names = 1, header = T,check.names=F)
ndata <- read.csv(file= "ARSeq/normalized_data.csv", header = T, check.names=F, row.names = 1)

# ssGSEA analysis 
ssgsea_analysis <- function(data,signature,meta=NULL,ssgsea.norm=T, intgroup=NULL, custom_color=NULL, show_column_names=F, padding = unit(c(2, 2, 2, 60), "mm")){
  
  # load required libraries
  require(GSEABase)
  require(GSVA)
  require(circlize)
  require(ComplexHeatmap)
  require(RColorBrewer)
  
  # Make a copy of the group type
  group <- signature[,1,drop=FALSE]
  colnames(group) <- c('cluster')
  # Drop the groups from signature
  signature <- signature[,-1]
  # transpose the signature
  signature <- data.frame(t(signature))
  
  # Convert signature into a named list
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
  x = named_list (signature)
  
  # Run the enrichment
  gbm_es <- gsva(as.matrix(data), x,method='ssgsea',ssgsea.norm=ssgsea.norm)
  
  # Viz using heatmap with metadata
  ssgsea_scaled = t(scale(t(gbm_es)))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#4B6AAF",  '#55B0AE', "#F8F6B8","#F5A15B","#B11E4B"))
  h1_col = colorRamp2(c(-1, 0, 1), c("#F8F6B8","#F5A15B","#B11E4B"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#f2e9e4",  '#c9ada7', "#9a8c98","#4a4e69","#22223b"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#e0e1dd",  '#778da9', "#415a77","#1b263b","#0d1b2a"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#f8f9fa",  '#dee2e6', "#ced4da","#6c757d","#343a40"))
  
  # converge grouping with ssgsea
  group <- group[row.names(group) %in% row.names(gbm_es), , drop = FALSE]
  
  if(!is.null(meta)){
    # Add column annotation (Sample group)
    col_ann <- data.frame(as.character(meta[,intgroup]))
    colnames(col_ann) <- c('Groups')
    # Create color pallete
    colourCount <- length(unique(col_ann[,1]))
    getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
    col_color <- getPalette(colourCount)
    if (!is.null(custom_color)){col_color <- custom_color}
    names(col_color) <- unique(col_ann[,1])
    col_color <- list(Groups = col_color)
    col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
  }
  
  # Heatmap
  if(is.null(meta)){
    h1 <- Heatmap(as.matrix(ssgsea_scaled), col = h1_col, border=T, cluster_columns = T,  
                  rect_gp = gpar(col = "#22223b", lwd = 1), show_column_names = show_column_names,
                  row_names_gp = gpar(fontsize = 12), cluster_rows = T,
                  row_split = factor(group$cluster, levels = unique(group$cluster)),
                  show_row_names = TRUE, name='ssGSEA score')
    draw(h1, heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)
    
  } else {
    h1 <- Heatmap(as.matrix(ssgsea_scaled), col = h1_col, border=T, cluster_columns = F,  
                  rect_gp = gpar(col = "#22223b", lwd = 1), show_column_names = show_column_names,
                  row_names_gp = gpar(fontsize = 12), cluster_rows = T, top_annotation=col_Ann,
                  row_split = factor(group$cluster, levels = unique(group$cluster)),
                  show_row_names = TRUE, name='ssGSEA score')
    draw(h1, heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)
    
  }
  
}

# with meta
ssgsea_analysis (data,meta,signature=custom_sig,intgroup="group")
# without meta
ssgsea_analysis (data,signature=custom_sig)




