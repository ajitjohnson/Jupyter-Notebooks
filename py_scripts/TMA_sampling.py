# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 11:58:57 2022
@author: Ajit Johnson Nirmal
Sampling TMA from WSI data
"""
# Lib
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
import seaborn as sns
import matplotlib.pyplot as plt


# Load data
#adata = sc.read('Y:\lsp-analysis\spatialpower\data\MEL1.h5ad')

# write data
#adata.write("Y:\lsp-analysis\spatialpower\data\MEL1.h5ad")



# Function to generate random virtual TMA's
def vtma (adata, vtma_radius=500, number_of_tma=10,
          layer = 'markers',phenotype='phenotype',
          min_cells=100, min_cells_radius=None,
          randomisation_tries = 100, method = 'non-overlapping',
          overlap_sensitivity = 10,
          x_coordinate='X_centroid',y_coordinate='Y_centroid',
          label='vtma'):
    """
    

    Parameters
    ----------
    adata : object
        AnnData Object.
    vtma_radius : list, optional
        A list of vtma sizes that needs to be constructed (e.g. [100,200,500]). The default is 500.
    number_of_tma : int, optional
        How many vTMAs are needed for each size requested. The default is 10.
    min_cells : int, optional
        The minimum number of cells within the largest vTMA or custom chosen size by min_cells_radius. 
        The default is 100.
    min_cells_radius : int, optional
        This parameter affects two things
        a) The size of the vTMA used for satisfying the min_cells parameter
        b) The size of the vTMA used to check for overlap between two TMAs when
        method = 'non-overlapping'. The default is the largest vTMA size provided.
    randomisation_tries : int, optional
        When method = 'non-overlapping', the number of times the algorithm tries to identify
        non-overlapping vTMA regions before the program gives up. The default is 100.
    method : string, optional
        Two options are available
        a) method = 'non-overlapping': The vTMAs are chosen such as they do not overlap.
        b) method = 'overlapping': Random vTMAs are sampled without consideration for overlap.
        The default is 'non-overlapping'.
    overlap_sensitivity : int, optional
        When method = 'non-overlapping', this number is the max number of cells that are allowed to
        to overlap between any two vTMA regions. The default is 10.
    x_coordinate : string, optional
        Column name containing x_coordinate. The default is 'X_centroid'.
    y_coordinate : string, optional
        Column name containing y_coordinate. The default is 'Y_centroid'.
    label : string, optional
        Name of the resultant column. The default is 'vtma'.
        The results will be saved in adata.uns['label']

    Returns
    -------
    Modified adata object.
    
    Example
    -------
    adata = vtma (adata, vtma_radius=[200,500,1000], number_of_tma=100,
              min_cells=1000, min_cells_radius=500,
              randomisation_tries = 100, 
              method = 'non-overlapping',
              overlap_sensitivity = 10,
              x_coordinate='X_centroid',y_coordinate='Y_centroid',
              label='vtma')

    """
    
    # create data with obs data
    data = pd.DataFrame(adata.obs)
    
    if layer == 'phenotype':
        # Do this
        exp_data = pd.get_dummies(adata.obs[phenotype])
        exp_data.reindex(adata.obs.index)
        true_proportion = exp_data.mean()
        #true_proportion = exp_data.sum()
        
    elif layer == 'markers':
        # extract the adata object
        exp_data = pd.DataFrame(adata.X, index = adata.obs.index, columns=adata.var.index)
        # theshold data
        exp_data[exp_data >= 0.5] = 1
        exp_data[exp_data < 0.5] = 0
        # compute the true proportion of + cells for each marker
        true_proportion = exp_data.mean()
    
    # Build vTMA
    if isinstance(vtma_radius, int):
        vtma_radius = [vtma_radius]
    
    # create columns for each tma size in the obs for vizualization purposes
    all_labels = []
    for i in vtma_radius:
        data[label + '_' + str(i)] = 'rest'
        all_labels = all_labels + [label + '_' + str(i)]
    

    # build the neighhood map
    print("Identifying neighbours within " + str(vtma_radius) + " pixels of every cell")
    kdt = BallTree(data[[x_coordinate,y_coordinate]], metric='euclidean')
    
    # Check if there is a custom radius given to adjust for number of cells within the vtma incase a range of vtma is given
    if min_cells_radius is not None:
        min_cells_radius = min_cells_radius
    else:
        min_cells_radius = vtma_radius[-1] 
        
    # permutation & vTMA range
    final_df = pd.DataFrame(index = exp_data.columns)
    permutation = 0
    stored_indeces = np.array([])
    
    while permutation <= number_of_tma:
        
        if method == 'non-overlapping':
            # check if any of these indeces are in the previously profiled region
            for i in range(randomisation_tries):
                print ("Finding a random seed point" + " " + str(i))
                random_point = data.sample(n=1) # create a random sample
                ind = kdt.query_radius(random_point[[x_coordinate,y_coordinate]], r=min_cells_radius, return_distance=False)
                if len(np.intersect1d(ind[0], stored_indeces)) < overlap_sensitivity: 
                    stored_indeces = np.append(stored_indeces,ind[0])
                    break
                elif i == randomisation_tries:
                    print ("Reached maximum number of randomisation_tries")
                    break
        
        if method == 'overlapping':
            random_point = data.sample(n=1) # create a random sample       
            ind = kdt.query_radius(random_point[[x_coordinate,y_coordinate]], r=min_cells_radius, return_distance=False)
        

        if len(ind[0]) > min_cells:
            permutation = permutation+1
            print("Permutation: " + str(permutation))
            final_corr = []
            for j in vtma_radius:
                print ('VTMA radius: ' + str(j))
                inde = kdt.query_radius(random_point[[x_coordinate,y_coordinate]], r=j, return_distance=False)
                
                # for each TMA size create a seperate column
                column = label + '_' + str(j)
                data.loc[data.iloc[inde[0].tolist()].index, column] = 'vTMA-'#+ str(permutation)
                
                # Calculate the proportion of cells within vTMA
                d = exp_data.loc[data.loc[data.iloc[inde[0].tolist()].index, column].index].mean().fillna(0)
                #d = exp_data.loc[data.loc[data.iloc[inde[0].tolist()].index, column].index].sum().fillna(0)
                # cumulate all the results into a list
                final_corr = final_corr + [d]   
            # merge the output vTMA's for a single permutation
            mystring = '_'+ 'perm_'+ str(permutation)
            prepare_index = [str(s) + mystring for s in vtma_radius]
            # DF
            fc = pd.DataFrame(final_corr, index = prepare_index ).T
            final_df = final_df.merge(fc, left_index=True, right_index=True)        
                
        if i == number_of_tma:
            continue    
    
    # Sort the such that they are by vTMA radius
    final_df = final_df.sort_index(axis = 1)
    # Add WSI proportion to the resultant matrix
    final_df.insert(0,'WSI_proportion',true_proportion)
    
    # save to adata
    adata.uns[label] = final_df
    adata.obs = data
    
    # return
    return adata
    

    
    
    
