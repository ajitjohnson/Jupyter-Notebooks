# Plotting the user defined top expressed genes/ protein markers
# Author: Ajit Johnson Nirmal
# Last updated: 1/17/19

import seaborn as sns; sns.set(style="white", color_codes=True)
import matplotlib.pyplot as plt
import anndata as ad
def mi_pl_topgenes (adata, n_top: int = 10):
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
    adata = ad.AnnData(X, obs=obs, var=var)
    return (adata)


d = {'avg': adata.var['avg_exp']}
df = pd.DataFrame(data=d)
df = df.sort_values(by=['avg'], ascending=False)
