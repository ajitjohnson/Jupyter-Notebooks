# -*- coding: utf-8 -*-
"""
Created on Sun May 29 22:02:02 2022
@author: ajn16
"""

# lib
import scimap as sm
import scanpy as sc


# load data
adata = sc.read('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/PCA/F8 image analysis/Z147_1_750.h5ad')

# add rois

image_path = '//research.files.med.harvard.edu/ImStor/sorger/data/RareCyte/Connor/Z147_PCAc/1_750/registration/1_750_reformat.ome.tif'
adata = sm.pl.addROI_image (image_path, adata, label="ROI")

#umap
#subsample
bdata = sc.pp.subsample(adata, fraction=0.01, copy=True)

bdata = sm.tl.umap(bdata, use_raw=True, log=True, n_neighbors=30, min_dist=0.1)

#plot umap
sm.pl.umap (bdata, color='phenotype', tight_layout=True, figsize=(9,5))

#animate
color_mapping = {
    'APCs': '#34a0a4', 
    'Blood Vessels': '#eb5e28', 
    'CD11C+ PDL1+ cells': '#34a0a4', 
    'Cytotoxic T cells': '#212529',
    'Keratinocytes': '#2b2d42', 
    'Langerhan cells': '#99d98c', 
    'Macrophages': '#3c6e71', 
    'Mast cells': '#fff3b0',
    'Melanocytes': '#9b2226', 
    'Myeloid Lineage': '#99d98c', 
    'Myofibroblast': '#d88c9a',
    'Patially Exhausted T cells': '#ffba08', 
    'Regulatory T cells': '#8e7dbe', 
    'T cells': '#00b4d8',
    'Terminally Exhausted T cells': '#ffba08', 
    'Tumor': '#9b2226', 
    'Unknown': '#f7ede2'
    }

sm.hl.animate(bdata,color='phenotype', s=1, n_frames=80, alpha=0.8,
              interval=50, final_frame=20, palette=color_mapping,
              save_animation="C:/Users/ajn16/Downloads/animate")

np.unique(adata.obs['phenotype'])





['APCs', 'Blood Vessels', 'CD11C+ PDL1+ cells', 'Cytotoxic T cells',
       'Keratinocytes', 'Langerhan cells', 'Macrophages', 'Mast cells',
       'Melanocytes', 'Myeloid Lineage', 'Myofibroblast',
       'Patially Exhausted T cells', 'Regulatory T cells', 'T cells',
       'Terminally Exhausted T cells', 'Tumor', 'Unknown']