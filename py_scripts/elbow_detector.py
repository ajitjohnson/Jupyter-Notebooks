#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 12:04:45 2020
@author: Ajit Johnson Nirmal
Elbow finder and K means elbow finder
"""

# Import library
import os
from kneed import DataGenerator, KneeLocator
import pandas as pd
import numpy as np
# For Kmenas clustering
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
from sklearn import metrics
from scipy.spatial.distance import cdist

# Set working direstory
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP/ARSeq/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP/ARSeq_Phospho/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/proteomics/ARSeq/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP/geomx/ARSeq/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/Figures for manuscript/Fig 3 F8 analysis (Round 4)/ARSeq")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/proteomics/ARSeq/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/atlas/ARSeq/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/ARSeq_geomx")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/ARSeq_pickseq")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/ARSeq_geomx/ET/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/atlas/")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/endotoxin_study")
os.chdir("/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/MVG")



# Import data
data = pd.read_csv("normalized_data.csv", index_col=0)
data = pd.read_csv("normalized_data_stable.csv", index_col=0)
data = pd.read_csv("et.csv", index_col=0)
data = pd.read_csv("combat_normalized_data.csv", index_col=0)
data = pd.read_csv("data.csv", index_col=0)
data = pd.read_csv("tumor_only.csv", index_col=0)
data = pd.read_csv("exophytic_geo.csv", index_col=0) #


# calculate variance
variance = data.var(axis=1).sort_values( ascending=False)

# x and y values
x = np.array(list(range(len(variance))))
y = variance.values

# Elbow finder
kneedle = KneeLocator(x, y, S=5, curve='convex', direction='decreasing')

kneedle = KneeLocator(x, y, S=2, curve='convex', direction='decreasing')


# Plot
kneedle.plot_knee_normalized()
sns.set_style("white")
kneedle.plot_knee()

# Print results
kneedle.elbow
kneedle.knee_y

# Subset the genes
mvg = variance[0:kneedle.elbow].index
mvg = data.loc[mvg,:]

# write out the data
mvg.to_csv("mvg_knee.csv")


##############################################################################
# Kmeans clustering of the most variable genes
sns.heatmap(mvg, cmap='vlag')

# k means determine k
distortions = []
K = range(1,20)
for k in K:
    kmeanModel = KMeans(n_clusters=k).fit(mvg)
    kmeanModel.fit(mvg)
    distortions.append(sum(np.min(cdist(mvg, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / mvg.shape[0])

# Knee locator
kn = KneeLocator(list(K), distortions, S=10, curve='convex', direction='decreasing')

# Plot the elbow
plt.xlabel('KMeans clusters')
plt.ylabel('Distortion')
plt.title('K means clustering of the Most Variable Genes')
plt.plot(K, distortions, 'bx-')
plt.vlines(10, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xticks(np.arange(1, 20, step=1))

# Knee
kn.knee


# Alternate method
def calculate_WSS(points, kmax):
  sse = []
  for k in range(1, kmax+1):
    kmeans = KMeans(n_clusters = k).fit(points)
    centroids = kmeans.cluster_centers_
    pred_clusters = kmeans.predict(points)
    curr_sse = 0
    
    # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
    for i in range(len(points)):
      curr_center = centroids[pred_clusters[i]]
      curr_sse += (points[i, 0] - curr_center[0]) ** 2 + (points[i, 1] - curr_center[1]) ** 2
      
    sse.append(curr_sse)
  return sse

K = calculate_WSS(mvg.values, 20)

# Plot the elbow
plt.xlabel('KMeans clusters')
plt.ylabel('Within-Cluster-Sum of Squared Errors (WSS)')
plt.title('K means clustering of the Most Variable Genes')
plt.plot(list(range(20)), K, 'bx-')
plt.vlines(10, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
plt.xticks(np.arange(1, 20, step=1))




















