# Preprocess the imaging data: Convert normalized data into AnnData
# Author: Ajit Johnson Nirmal
# Last updated: 1/16/19

import anndata as ad
import numpy as np
def mi_pp_anndata (data):
    print ("Converting dataframe into anndata - Annotated Data...")
    # Extract array of expression values
    X = data.values
    # Convert the cell names/ index into observations (obs)
    obs = pd.DataFrame()
    obs['cells'] = data.index
    # Convert the protein Markers/ genenames into variables (var)
    var_names = list(data)
    n_vars = len(var_names)
    var = pd.DataFrame(index=var_names)
    # Create AnnData
    adata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
    return (adata)
