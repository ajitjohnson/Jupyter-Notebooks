# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 11:45:52 2022
@author: Ajit Johnson Nirmal
The very first 3D data from LSP 
"""

# import lib
import pandas as pd
import glob
import os
import anndata as ad
import scanpy as sc
import numpy as np
import seaborn as sns
import scimap as sm

# wd
os.chdir('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/images/')

# import data
adata = sc.read('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/data.h5ad')

# initial processing of data
all_files = glob.glob("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/mean/*.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename,skiprows=3)
    # make ID
    df.index = 'cell_' + df['ID'].astype(str)
    # remove unwanted columns
    df.drop(['Unit', 'Category', 'Channel', 'Image', 'Time', 'Unnamed: 7', 'ID'], axis=1, inplace=True)
    li.append(df)

data = pd.concat(li, axis=1, ignore_index=True)

# add marker names and remove DNA
# load markers.csv
markers = pd.read_csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/markers.csv")

data.columns = markers['marker_name'].tolist()
data = data.loc[:,~data.columns.str.contains('Hoechst', case=False)]
data.drop(['backgroundPM', 'contoursPM', 'nucleiPM','cytochromeC','PDL1', 'aSMA_1'], axis=1, inplace=True)

# initial processing of meta
all_meta = glob.glob("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/meta/*.csv")

li = []

for filename in all_meta:
    df = pd.read_csv(filename,skiprows=3)
    # make ID
    df.index = 'cell_' + df['ID'].astype(str)
    # remove unwanted columns
    df.drop(['Unit', 'Category', 'Collection', 'Time', 'Unnamed: 6', 'ID'], axis=1, inplace=True)
    li.append(df)

meta = pd.concat(li, axis=1, ignore_index=True)
meta.columns = ['X', 'Y', 'Z']

# Add volume
df = pd.read_csv('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/volume\\cycle1_Volume.csv',skiprows=3)
df.index = 'cell_' + df['ID'].astype(str)
df.drop(['Unit', 'Category', 'Time', 'Unnamed: 5', 'ID'], axis=1, inplace=True)
meta = pd.concat([meta,df], axis=1, ignore_index=True)
meta.columns = ['X', 'Y', 'Z', 'volume']

# Add cellid column
meta['CellID'] = 'image_1'
meta['imageid'] = 'imageid'

# create anndata object
adata = ad.AnnData (data)
adata.obs = meta

# add image meta
markers_image = pd.read_csv("C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/markers_in_image.csv")
adata.uns['all_markers'] = markers_image['marker_name'].tolist()

# copy raw data and log the data
adata.raw = adata
adata.X = np.log1p(adata.X)

# scale the x and y coordinates
adata.obs['X'] = adata.obs['X'] * 3.33333333
adata.obs['Y'] = adata.obs['Y'] * 3.33333333
adata.obs['Z'] = adata.obs['Z'] * 2

# look at distribution of expression
logdata = np.log1p(data)
sns.displot(logdata, x="SOX10", kind="kde")
sns.displot(logdata, x="laminB1", kind="kde")

# convert volume to log scale
adata.obs['volume'] = np.log1p(adata.obs['volume'])

# gate 
#['CD11c', 'CD45', 'laminB1', 'S100B', 'aSMA_1', 'nucleoporin', 'CD163', 'collagen', 'btubulin', 'vimentin', 'aSMA', 'SOX10', 'actin', 'MART1']
image_path = "Y:/lsp-data/cycif-techdev/confocal/40xMET15microns/outputb/cycle1_3D_MERGED_PRJ.tif"
sm.pl.gate_finder(image_path, adata, marker_of_interest='MART1', from_gate=2, to_gate=5, increment=0.05, x_coordinate='X', y_coordinate='Y', point_size=10, imageid='imageid', subset=None, seg_mask=None)

# rescale based on manual gates
manual_gate = pd.read_csv('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/manual_gates.csv')
adata = sm.pp.rescale (adata, gate=manual_gate)

# phenotype cells
# Load the gating workflow
phenotype = pd.read_csv('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/phenotype_workflow.csv')
adata = sm.tl.phenotype_cells (adata, phenotype=phenotype, label="phenotype") 
adata.obs['phenotype'].value_counts()

# view image
image_path = "Y:/lsp-data/cycif-techdev/confocal/40xMET15microns/outputb/cycle1_3D_MERGED_PRJ.tif"
sm.pl.image_viewer(image_path, adata, overlay='phenotype',  x_coordinate='X', y_coordinate='Y')


# UMAP
data.columns
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata) # Build a UMAP to visualize the neighbourhood graph
# cluster data
sc.tl.leiden(adata, resolution = 0.3)
# viz umap
sc.pl.umap(adata, color=['SOX10', 'S100B','MART1',
                         'CD45', 'CD11c', 'CD163', 'aSMA',
                         'laminB1','actin', 'volume','leiden','phenotype'], cmap= 'vlag', ncols=3,
           use_raw=False, s=30) # Plot the UMAP

# scatter plot
plotly3D (adata,phenotype='phenotype',image_id=None,x='X',y='Y',z='Z', size=6)

# find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# heataps
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype', dendrogram=True, use_raw=False, cmap="vlag") #standard_scale='var'


# interaction analysis
# spatial_interaction heatmap for a single image
adata = spatial_interaction(adata, method='radius', radius=70,pval_method='zscore', x_coordinate='X',y_coordinate='Y',z_coordinate='Z')
adata = spatial_interaction(adata, method='radius', radius=70,pval_method='zscore', x_coordinate='X',y_coordinate='Y',z_coordinate=None)

sm.pl.spatial_interaction(adata, summarize_plot=True, row_cluster=False, linewidths=0.75, linecolor='black')
sm.pl.spatial_interaction(adata, summarize_plot=True, binary_view=True,row_cluster=False, linewidths=0.75, linecolor='black')

# spatial_interaction heatmap for multiple images
sns.set(font_scale=0.6)
sm.pl.spatial_interaction(adata, summarize_plot=False, row_cluster=True, col_cluster=True, yticklabels=True)

# voronoi plor
sm.pl.voronoi(adata, color_by='phenotype', colors=None, 
         x_coordinate='X', y_coordinate='Y',
         imageid='ImageId',subset=None,
         voronoi_edge_color = 'black',voronoi_line_width = 0.2, 
         voronoi_alpha = 0.5, size_max=np.inf,
         overlay_points='phenotype', overlay_points_categories=None, 
         overlay_drop_categories=None,
         overlay_point_size = 5, overlay_point_alpha= 1, 
         overlay_point_shape=".", 
         plot_legend=True, 
         legend_size=6)


# spatial distance
adata = spatial_distance (adata,x_coordinate='X', y_coordinate='Y', z_coordinate='Z', phenotype='phenotype')
adata = spatial_distance (adata,x_coordinate='X', y_coordinate='Y', z_coordinate=None, phenotype='phenotype')
#plot
sm.pl.spatial_distance (adata, method='numeric', distance_from='Tumor', distance_to = ['Macrophages','Myeloid','Other Immune','Blood Vessel'])

# density
sc.tl.embedding_density(adata, basis='umap', groupby='phenotype')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_phenotype',ncols=3)

# density of all cells
adata.obsm['X_umap']
adata.obsm['X_cells'] = meta[['X', 'Y']].to_numpy()
sc.tl.embedding_density(adata, basis='umap', groupby='phenotype')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_phenotype',ncols=2)


# save result
adata.write('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/data.h5ad')
output = sm.hl.scimap_to_csv(adata, data_type='raw', output_dir=None)
output.to_csv('C:/Users/ajn16/Dropbox (Partners HealthCare)/Data/3D_imaging/data/first_15micron/processed_data.csv')
