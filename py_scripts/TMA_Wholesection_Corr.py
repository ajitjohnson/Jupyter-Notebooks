#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 20:14:31 2021
@author: Ajit Johnson Nirmal
Function to calculate correlation between whole section and corresponding TMA's for each marker
"""

# Lib
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import BallTree
from sklearn.preprocessing import MinMaxScaler

# Function

def marker_corr (adata, corr_from, corr_to, threshold=0.5, imageid='imageid'):
    
    # EXAMPLE
    #corr_from = ['PTCL6_n']
    #corr_to = ['unmicst-86','unmicst-88']
    #marker_corr (adata, corr_from, corr_to, threshold=0.5, imageid='imageid')
    
    # create a copy of the anndata object
    adata = adata.copy()
    
    # sanitycheck
    if isinstance(corr_from, str):
        corr_from = [corr_from]
    if isinstance(corr_to, str):
        corr_to = [corr_to]
    
    # subset the anndata object
    # identify all the samples needed
    adata = adata[adata.obs[imageid].isin(corr_from+corr_to)]
    
    # all data
    data = pd.DataFrame(adata.X, index = adata.obs.index, columns=adata.var.index)
    
    # theshold data
    data[data >= threshold] = 1
    data[data < threshold] = 0
    
    # divide the data into two
    data_corr_from = data.loc[adata[adata.obs[imageid].isin(corr_from)].obs.index]
    # compute the proportion
    data_corr_from = data_corr_from.mean()
    
    data_tmp = pd.DataFrame(index= data_corr_from.index)
    final_corr = []
    for i in corr_to:
        print('processing: ' + str(i))
        d = data.loc[adata[adata.obs[imageid].isin([i])].obs.index].mean()
        data_tmp = pd.concat([d,data_tmp],axis=1)
        corr = pearsonr(data_corr_from, d)
        final_corr = final_corr + [corr[0]]

    # add columns to the dataframe
    data_tmp.columns = corr_to
    
    # Merge with the from data column
    data_tmp = data_tmp.merge(pd.DataFrame(data_corr_from, index = data_corr_from.index, columns=['corr_from']), left_index=True, right_index=True)
    
    # coorelation df
    final_corr = pd.DataFrame(final_corr, index = corr_to, columns=['Pearson Corr'] )
    
    # return
    return [data_tmp, final_corr]


def marker_corr_random (adata, subset=None, number_vtma=10, 
                        x_coordinate='X_centroid',y_coordinate='Y_centroid',
                        threshold=0.5, min_cells=1000, vtma_radius=500, custom_radius=None,
                        imageid='imageid', label='vtma'):
    
    # EXAMPLE
    #corr_from = ['PTCL6_n']
    #corr_to = ['unmicst-86','unmicst-88']
    #marker_corr (adata, corr_from, corr_to, threshold=0.5, imageid='imageid')
    
    # create a copy of the anndata object
    bdata = adata.copy()
    
    # subset data if needed
    if subset is not None:
        if isinstance(subset, str):
            subset = [subset]
        bdata = bdata[bdata.obs[imageid].isin(subset)]
        
    # print image id
    print (np.unique(bdata.obs[imageid]))
        
    
    # create data
    data = pd.DataFrame(bdata.obs)
    
    # expression data
    exp_data = pd.DataFrame(bdata.X, index = bdata.obs.index, columns=bdata.var.index)
        
    # theshold data
    exp_data[exp_data >= threshold] = 1
    exp_data[exp_data < threshold] = 0
    
    # compute the proportion of pos cells all data
    data_corr_from = exp_data.mean()
    #dcf = exp_data.mean()
    #x = MinMaxScaler().fit_transform(dcf.values.reshape(-1,1))
    #data_corr_from = pd.DataFrame(x, index = dcf.index )

    # build the neighhood map
    print("Identifying neighbours within " + str(vtma_radius) + " pixels of every cell")
    kdt = BallTree(data[[x_coordinate,y_coordinate]], metric='euclidean')
    
    
    # create vTMA based on the size provided by user
    
    if isinstance(vtma_radius, int):
        vtma_radius = [vtma_radius]
        
    # Process:
    # for each TMA size - calculate X permuations of correlation values
    # pick a point - increase the size gradually - compute the corr
    # repeat for the next permutation
    
    # Check if there is a custom radius given to adjust for number of cells within the vtma incase a range of vtma is given
    if custom_radius is not None:
        custom_radius = custom_radius
    else:
        custom_radius = vtma_radius[-1]   
    
    # create columns for all vTMA sizes to later modify
    all_labels = []
    for i in vtma_radius:
        data[label + '_' + str(i)] = 'rest'
        all_labels = all_labels + [label + '_' + str(i)]
        
    # permutation & vTMA range
    final_df = pd.DataFrame(index = vtma_radius)
    permutation = 0
    while permutation <= number_vtma:
        random_point = data.sample(n=1) # create a random sample       
        ind = kdt.query_radius(random_point[[x_coordinate,y_coordinate]], r=custom_radius, return_distance=False)
        if len(ind[0]) > min_cells:
            permutation = permutation+1
            print("Permutation: " + str(permutation))
            final_corr = []
            for j in vtma_radius:
                print ('VTMA radius: ' + str(j))
                inde = kdt.query_radius(random_point[[x_coordinate,y_coordinate]], r=j, return_distance=False)
                # for each TMA size create a seperate column
                column = label + '_' + str(j)
                data.loc[data.iloc[inde[0].tolist()].index, column] = 'vTMA-'+ str(permutation)
                
                # Perform the correlation analysis
                d = exp_data.loc[data.loc[data.iloc[inde[0].tolist()].index, column].index].mean().fillna(0)
                corr = pearsonr(data_corr_from, d)
                final_corr = final_corr + [corr[0]]   
            # merge the output correlation into a larger DF
            fc = pd.DataFrame(final_corr, index = vtma_radius, columns= ['perm_'+ str(permutation)] )
            final_df = final_df.merge(fc, left_index=True, right_index=True)        
                
        if i == number_vtma:
            continue
    
    # add to adata file if needed
    #bdata.obs[label] = data['vtma_40000']    
    # return
    #return [data_tmp, final_corr]
    return final_df


# Figure function
def scatter_text(x, y, text_column, data):

    # Create the scatter plot
    p1 = sns.scatterplot(x=x, y=y, data=data, legend=False)
    # Add text besides each point
    for line in range(0,data.shape[0]):
         p1.text(data[x][line]+0.01, data[y][line], 
                 data[text_column][line], horizontalalignment='left', 
                 size='medium', color='black', weight='semibold')
    # Set title and axis labels
    return p1




# %% Example

adata = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/ptcl_nos_aj.h5ad')
adata = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/TMA_and_wholeslide.h5ad')


np.unique(adata.obs['imageid'])

# Remove nan's from data
adata.X[np.isnan(adata.X)] = 0
np.isnan(adata.X).any()



corr_from = ['PTCL6_n']; corr_to = ['unmicst-86','unmicst-88']
corr_from = ['unmicst-86','unmicst-88']; corr_to = ['PTCL6_n']
corr_from = ['PTCL7_n']; corr_to = ['unmicst-18','unmicst-21', 'unmicst-25']
corr_from = ['unmicst-18','unmicst-21', 'unmicst-25']; corr_to = ['PTCL7_n']
result =  marker_corr (adata, corr_from, corr_to, threshold=0.5, imageid='imageid')

vtma_radius = [2, 10, 20, 50,100,500,1000,2000,8000,15000,20000,40000]
result = marker_corr_random (adata, number_vtma=1000, subset='PTCL7_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)


result = result.T
result = result.add_prefix('size_')
x = pd.melt(result)
y = x[x['value'].notnull()]
#x = x.fillna(0)
sns.lineplot(data=y, x="variable", y="value")
plt.xticks(rotation=90)
#y['variable'].value_counts()



# Figure
scatter_text(x = 'corr_from', y = 'unmicst-25', text_column = 'names', data = r)
scatter_text(x = 'corr_from', y = 'PTCL6_n', text_column = 'names', data = r)
scatter_text(x = 'corr_from', y = 'vTMA-45', text_column = 'names', data = r)

bdata = adata.copy()
bdata.obs[label] = data [label]
plotly (bdata,phenotype='vtma',size=2, opacity = 1)
plotly (bdata,phenotype='vtma',size=2, opacity = 1, x='X_position', y='Y_position')


# %% Bulk Run

# Remove nan's from data
ptcl_nos.X[np.isnan(ptcl_nos.X)] = 0
np.isnan(tonsil.X).any()
np.unique(alcl.obs['imageid'])

vtma_radius = [2, 10, 20, 50,100,500,1000,2000,8000,15000,20000,40000]

# PTCL1

# ALK- ALCL
#alcl = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/ALCL_aj.h5ad')
#alcl = alcl[alcl.obs['imageid'].isin(['PTCL2_n', 'PTCL3_n', 'PTCL4_n'])]
# PTCL2
PTCL2_n = marker_corr_random (alcl, number_vtma=1000, subset='PTCL2_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)

# PTCL3
PTCL3_n = marker_corr_random (alcl, number_vtma=1000, subset='PTCL3_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)

# PTCL4
PTCL4_n = marker_corr_random (alcl, number_vtma=1000, subset='PTCL4_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)

# PTCL5


# PTCL NOS
#ptcl = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/ptcl_nos_aj.h5ad')
#ptcl = ptcl[ptcl.obs['imageid'].isin(['PTCL6_n', 'PTCL7_n'])]
# PTCL6
PTCL6_n = marker_corr_random (ptcl, number_vtma=1000, subset='PTCL6_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)


# PTCL7
PTCL7_n = marker_corr_random (ptcl, number_vtma=1000, subset='PTCL7_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)


# ALCL
#aitl = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/AITL_aj.h5ad')
#aitl = aitl[aitl.obs['imageid'].isin(['PTCL8_n', 'PTCL9_n'])]
# PTCL8
PTCL8_n = marker_corr_random (aitl, number_vtma=1000, subset='PTCL8_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)

# PTCL9
PTCL9_n = marker_corr_random (aitl, number_vtma=1000, subset='PTCL9_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)

# Tonsil
#tonsil = sc.read('/Users/aj/Dropbox (Partners HealthCare)/PTCL Analysis/H5AD Files/TMA_and_wholeslide.h5ad')
#tonsil = tonsil[tonsil.obs['imageid'].isin(['Ton_192_n'])]
# PTCL10
PTCL10_n = marker_corr_random (tonsil, number_vtma=1000, subset='Ton_192_n',
                             x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)


# Melanoma
mel = sc.read('/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/F8 image analysis/Z147_1_750.h5ad')
#tonsil = tonsil[tonsil.obs['imageid'].isin(['Ton_192_n'])]
# PTCL10
mel_results = marker_corr_random (mel, number_vtma=1000, #subset='Ton_192_n',
                             #x_coordinate='X_position',y_coordinate='Y_position',
                             threshold=0.5, min_cells=1000, vtma_radius=vtma_radius,
                             imageid='imageid', label='vtma', custom_radius=500)


# write out the correlation results
mel_results.to_csv('/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/10-correlation_analysis/mel_results.csv')


# plot
result = mel_results.T
result = result.add_prefix('size_')
x = pd.melt(result)
y = x[x['value'].notnull()]
#x = x.fillna(0)
sns.lineplot(data=y, x="variable", y="value")
plt.xticks(rotation=90)
plt.yticks(np.arange(0, 1.2, 0.1))
#y['variable'].value_counts()
pt = result['size_500'].mean() # connection pint
x1, y1 = [0,5], [pt,pt]
x2, y2 = [5,5], [0,pt]
plt.plot(x1, y1,'--',color='grey')
plt.plot(x2, y2,'--',color='grey')
plt.text(0.1, pt+0.01, str(round(pt,2)), style='italic', bbox={'facecolor': 'red', 'alpha': 0.6, 'pad': 10})





# %% Left over

# =============================================================================
#     # calculate the correlation for all vTMA range
#     #data_tmp = pd.DataFrame(index= data_corr_from.index)
#     # identify the vTMA id's
#     corr_to = np.unique(data[all_labels[0]]).tolist()
#     if 'rest' in corr_to:
#         corr_to.remove('rest')
#     # intialte a df to hold all correlation values
#     final_df = pd.DataFrame(index = corr_to)
#     for i in all_labels:
#         print('processing VTMA -------- : ' + str(i))
#         final_corr = []
#         for j in corr_to:
#             print('processing: ' + str(j))
#             d = exp_data.loc[data[data[i].isin([j])].index].mean().fillna(0)
#             #dd = exp_data.loc[data[data[i].isin([j])].index].mean()
#             #y = MinMaxScaler().fit_transform(dd.values.reshape(-1,1))
#             #d = pd.DataFrame(y, index = dd.index )
# 
#             #data_tmp = pd.concat([d,data_tmp],axis=1)
#             corr = pearsonr(data_corr_from, d)
#             #corr = spearmanr(data_corr_from, d)
#             final_corr = final_corr + [corr[0]]
#             
#             #error = abs(data_corr_from - d)
#             #error = (error/data_corr_from) * 100
#             #final_corr = final_corr + [np.mean(error)]
#             #corr = pearsonr(data_corr_from[0], d[0])
#             
#         # add column names
#         fc = pd.DataFrame(final_corr, index = corr_to, columns= [str(i)] )
#         final_df = final_df.merge(fc, left_index=True, right_index=True)
# =============================================================================
