#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 20:15:45 2021
@author: Ajit Johnson Nirmal
KNN classifier
"""

# Lib
import numpy as np
import pandas as pd

import scipy
import matplotlib.pyplot as plt
import sklearn
from sklearn.neighbors import KNeighborsClassifier
from sklearn import neighbors
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn import metrics
import mpl_scatter_density
from sklearn.ensemble import RandomForestClassifier


# train data
train_adata = adata[adata.obs['general_roi'] == 'aj_IT']
train_adata = sc.pp.subsample(train_adata, fraction=0.3, copy=True)
train_data = pd.DataFrame(train_adata.X, index = train_adata.obs.index, columns = train_adata.var.index)
train_data = train_data[['SOX10', 'S100B', 'MITF','KI67','PCNA','HLADPB1', 'S100A','CCNA2', 'CCND1', 'CD63']]
train_data = train_data.values

# Import data
data = pd.DataFrame(adata.X, index = adata.obs.index, columns = adata.var.index)

# Subset data with required markers
# ['SOX10', 'S100B','MITF','KI67','PCNA','S100A','CCNA2', 'CCND1', 'CD63']
#data = data[['S100B','MITF','KI67','PCNA','CCNA2']]
data = data[['SOX10', 'S100B', 'MITF','KI67','PCNA','HLADPB1', 'S100A','CCNA2', 'CCND1', 'CD63']]
data = data[['SOX10','S100B','MITF','KI67','PCNA','CCNA2', 'CCND1']]
#data = data[['SOX10', 'S100B','MITF','KI67','PCNA','S100A','CCNA2', 'CCND1']]
d = data.values

# Classification
classification = adata.obs
classification = classification['general_roi'].tolist()

# convert to binary classification
c = [1 if x == 'aj_IT' else 0 for x in classification]


# pre-processing
X = preprocessing.scale (d)
train = preprocessing.scale (train_data)

# split into train and test
x_train, x_test, y_train, y_test = train_test_split (X, c, test_size = .70)

# build the model
#clf = neighbors.KNeighborsClassifier()
#clf.fit(x_train, y_train)
#print(clf)

# RF
clf=RandomForestClassifier(class_weight='balanced')
clf=RandomForestClassifier()
clf.fit(x_train, y_train)

# Evaluate model
#y_expect = y_test
#y_pred = clf.predict(x_test)
#metrics
#x = metrics.classification_report(y_expect, y_pred)


# probability
prob = clf.predict_proba(X)
# extract the first element of the probabilities
#lst2 = [item[0] for item in prob]
lst2 = [item[1] for item in prob]
# sns plot
x = adata.obs['X_centroid']
y = adata.obs['Y_centroid']

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
density = ax.scatter_density(x, y, dpi=70, c=lst2, cmap='nipy_spectral_r') #nipy_spectral gist_stern_r gist_ncar
ax.invert_yaxis()
fig.colorbar(density, label='Probability Prediction')
plt.xticks([]) ; plt.yticks([]);


# normal scatter plot
fig = plt.figure()
plt.scatter(x, y,  c = prediction, s=0.1, cmap='nipy_spectral_r') #nipy_spectral gist_stern_r gist_ncar
plt.invert_yaxis()
plt.colorbar(density, label='Probability Prediction')
plt.xticks([]) ; plt.yticks([]);


  
%matplotlib inline
%matplotlib qt

for i in maps:
    print(i)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    ax.scatter_density(x, y, dpi=300, c=lst2, cmap=i)
    ax.invert_yaxis()
    plt.xticks([]) ; plt.yticks([]);
    plt.savefig("/Users/aj/Downloads/fig/" + str(i))

/Users/aj/Downloads/fig

maps = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r']
