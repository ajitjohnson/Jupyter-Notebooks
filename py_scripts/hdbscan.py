#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:57:48 2021

@author: aj
"""

import sys, os
import anndata as ad
import hdbscan
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


os.chdir ("/Users/aj/Dropbox (Partners HealthCare)/packages/scimap")
image_path = os. getcwd() + '/scimap/tests/_data/example_data.h5ad'
adata = ad.read(image_path)

# data
data = pd.DataFrame(np.log1p(adata.raw.X))

#hdbscan
clusterer = hdbscan.HDBSCAN(min_cluster_size=5,metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False)
cluster_labels = clusterer.fit_predict(data)
cluster_labels = clusterer.fit_predict(d)


x = adata.obs['X_position'].values
y = adata.obs['Y_position'].values
c = np.log1p(adata.raw.X[:,4])

#plot
fig, axes = plt.subplots()
#axes.scatter(x, y, c=c, edgecolors='none', s=25, cmap='vlag')
axes.scatter(x, y, c=cluster_labels, edgecolors='none', s=25)
axes.invert_yaxis()
plt.xticks([]) ; plt.yticks([]);


#plot
fig, axes = plt.subplots()
#axes.scatter(x, y, c=c, edgecolors='none', s=25, cmap='vlag')
axes.scatter(data[:,0], data[:,1], c=cluster_labels, edgecolors='none', s=25)
axes.invert_yaxis()
plt.xticks([]) ; plt.yticks([]);


# 
sc.pp.neighbors(adata, n_neighbors=5)
sc.tl.umap(adata)
sc.tl.leiden(adata)

x = adata.obsm['X_umap'][:,0]
y = adata.obsm['X_umap'][:,1]


data = pd.DataFrame(adata.obsm['X_umap'])
data = data.todense()


np.unique(cluster_labels)
np.unique(adata.obs['leiden'])
d = c.reshape(-1, 1)
