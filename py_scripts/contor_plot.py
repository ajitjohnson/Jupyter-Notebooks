#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 19:15:50 2021
@author: aj
Contor Map for a ROI
"""

# Import library
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import scipy

x_coordinate='X_centroid';y_coordinate='Y_centroid'

bdata = adata.copy()
bdata = bdata[bdata.obs['phenotype_2'] == 'Tumor']

# artificial grid
x = np.arange(0, 36810, 100)
y = np.arange(0, 25760, 100)
xx, yy = np.meshgrid(x, y)
XYpairs = np.vstack([ xx.reshape(-1), yy.reshape(-1) ])
pseudo_x = XYpairs[0]
pseudo_y = XYpairs[1]

pseudo_data = pd.DataFrame({'x': pseudo_x, 'y': pseudo_y})
real_data = pd.DataFrame({'x': bdata.obs[x_coordinate], 'y': bdata.obs[y_coordinate]})

# build the tree
tree = BallTree(real_data, leaf_size= 2)
dist, ind = tree.query(pseudo_data, k=500, return_distance= True)

# replace the dist matrix with weights between 0-1 (shotest:1; farthest:0)
temp = np.reshape(dist, dist.shape[0] * dist.shape[1])
temp_scaled = (temp - np.min(temp)) / (np.max(temp) - np.min(temp))
temp_scaled = 1-temp_scaled # make lowest distance have the highest weight
dist_scaled = np.reshape(temp_scaled, dist.shape)

# replace the indeces with value of a marker of interest
goi = 'CD63'
goi_loc = bdata.var.index.get_loc(goi)
#maker_data = bdata.X[:,goi_loc]
maker_data = np.log1p(bdata.raw.X[:,goi_loc])

#reshap ind into one long string
A = np.reshape(ind, ind.shape[0] * ind.shape[1])
A = A.tolist()
B = maker_data.tolist()
C = [B[i] for i in A]
ind_mapped = np.reshape(C, ind.shape)


# multiply the arrays to get a weighted matrix
wighted_distance = ind_mapped * dist_scaled
# spatial lag
spatial_lag = wighted_distance.mean(axis = 1)

# spatial lag
#spatial_lag = ind_mapped.mean(axis = 1)


# data
x = pseudo_data['x'].values
y = pseudo_data['y'].values
z = spatial_lag

# real data
a = adata.obs['X_centroid'].values
b = adata.obs['Y_centroid'].values
c = adata.X[:,goi_loc]

# contor plot
fig, ax = plt.subplots(figsize=(8,6))
#ax.scatter(pseudo_x, pseudo_y, s=0.2, c = 'b')
ax.set_facecolor('#adb5bd')
ax.scatter(a, b, s=0.2, c = c, cmap='RdYlBu_r') # '#83c5be'
ax.tricontour(x, y, z,  linewidths=1.2, colors='white', levels=150) #levels=15,
#ax.tricontourf(x, y, z,  cmap="cividis") # levels=15,
ax.invert_yaxis()
plt.xticks([]) ; plt.yticks([]);


# plotyl
plotly (IM,phenotype='phenotype_2',size=2, opacity = 1)
plotly (adata,phenotype='general_roi',size=2, opacity = 1)
