# Script to figure out the genes contibuting towards Principal components
# Created on 5/28/2019
# Author: Ajit Johnson Nirmal


pc_genes <- function(data){
  require(Seurat)
  #data <- (data/rowSums(data))*median(unlist(as.list(as.data.frame(t(data)))))
  data <- log1p(data)
  exp_seurat <- CreateSeuratObject(counts = data)
  # Scale data
  exp_seurat <- ScaleData(exp_seurat)
  # PCA
  exp_seurat <- RunPCA(exp_seurat)
  
  VizDimLoadings(exp_seurat, dims = 1:5, reduction = "pca")
  
  DimHeatmap(exp_seurat, dims = 1:10, cells = 500, balanced = TRUE)

  exp_seurat <- FindVariableFeatures(exp_seurat, selection.method = "vst", nfeatures = 100)
  
}

