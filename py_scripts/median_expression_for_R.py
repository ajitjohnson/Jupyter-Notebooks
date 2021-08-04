#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 15:51:35 2020
@author: Ajit Johnson Nirmal
"""

# %% median expression of cell phenotyping for R

adata = alcl_tumor.copy()

# set wd
import os
os.chdir('/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA')
os.chdir('/Volumes/SSD/Dropbox (Partners HealthCare)/PTCL Analysis/')
os.chdir('/Users/aj/Dropbox (Partners HealthCare)/PCA_Atlas1-2020/Figures/figures_aj_pdf')
os.chdir('/Users/aj/Dropbox (Partners HealthCare)/PTCL_Jackson/tumor_subgroups/ALCL/heatmaps') # ALCL AITL PTCL

# median_expression
comb_data = pd.DataFrame(np.log1p(adata.raw.X), columns=adata.var.index)
comb_data['phenotype'] = adata.obs['tumor_het_phenograph_scaled'].values
median_expression = comb_data.groupby('phenotype').median()
# Gene normalize the data
median_expression = median_expression.apply(lambda x: (x - x.min()), axis=0).apply(lambda x: (x / x.max()), axis=0)
median_expression.to_csv('for R/median_expression.csv', index=True)

# median_cellfeatures
median_cellfeatures = pd.DataFrame(adata.obs)[['Area', 'Eccentricity', 'Extent', 'Solidity']]
median_cellfeatures = median_cellfeatures.apply(pd.to_numeric)
median_cellfeatures['phenotype'] = adata.obs['tumor_het_phenograph_scaled'].values
median_cellfeatures = median_cellfeatures.groupby('phenotype').median()
# normalize the data
median_cellfeatures = median_cellfeatures.apply(lambda x: (x - x.min()), axis=0).apply(lambda x: (x / x.max()), axis=0)
median_cellfeatures.to_csv('for R/median_cellfeatures.csv', index=True)

# phenotype distribution for ROI
#subset_yaxis = [ 'aj_IR', 'aj_BTIL', 'aj_IB', 'aj_MIS', 'aj_IT', 'aj_ET']
pheno_prop = sm.pl.stacked_barplot (adata,x_axis='tumor_het_phenograph_scaled',y_axis='patientid',method='percent',return_data=True)
pheno_prop.to_csv('for R/pheno_prop.csv', index=True)

# Pheno group
pheno_group = pd.DataFrame(adata.obs[['tumor_het_phenograph_scaled','tumor_het_phenograph_scaled_consolidated']])
pheno_group = pheno_group.drop_duplicates(subset=['tumor_het_phenograph_scaled'])
pheno_group.to_csv('for R/pheno_group.csv', index=False)


sc.pl.matrixplot(adata, var_names=adata.var.index, groupby='Kmeans', use_raw=False,
                         cmap='RdBu_r', dendrogram=True
                         )


# %% Using scaled data

adata = aitl.copy()
adata = adata[adata.obs['tumor_stroma'] == 'Tumor']


# input
each_row = 'tumor_expression_TC' # tumor_het_phenograph_scaled patientid
within_barplot = 'general_roi' #patientid phenotype
row_split ='tumor_het_phenograph_scaled_consolidated' # tumor_het_phenograph_scaled_consolidated patient_group

# median_expression
comb_data = pd.DataFrame(adata.X, columns=adata.var.index)
comb_data['phenotype'] = adata.obs[each_row].values
median_expression = comb_data.groupby('phenotype').median()
median_expression.to_csv('for R/median_expression.csv', index=True)

# median_cellfeatures
median_cellfeatures = pd.DataFrame(adata.obs)[['Area', 'Eccentricity', 'Extent', 'Solidity']]
median_cellfeatures = median_cellfeatures.apply(pd.to_numeric)
median_cellfeatures['phenotype'] = adata.obs[each_row].values
median_cellfeatures = median_cellfeatures.groupby('phenotype').median()
# normalize the data
median_cellfeatures = median_cellfeatures.apply(lambda x: (x - x.min()), axis=0).apply(lambda x: (x / x.max()), axis=0)
median_cellfeatures.to_csv('for R/median_cellfeatures.csv', index=True)

# phenotype distribution for ROI
#subset_yaxis = [ 'aj_IB', 'aj_IT', 'aj_ET']
pheno_prop = sm.pl.stacked_barplot (adata,x_axis=each_row,y_axis=within_barplot,method='percent',return_data=True) #subset_yaxis=subset_yaxis
pheno_prop.to_csv('for R/pheno_prop.csv', index=True)

# Pheno group
pheno_group = pd.DataFrame(adata.obs[[each_row,row_split]])
pheno_group = pheno_group.drop_duplicates(subset=[each_row])
pheno_group.to_csv('for R/pheno_group.csv', index=False)

# %% for tumor heterogenity analysis by phenograph
# spatial_interaction_with_tumor_subtypes_all
# input
interaction_data = interaction
each_row = 'tumor_het_phenograph_scaled' # tumor_het_phenograph_scaled patientid
within_barplot = 'patientid' #patientid phenotype
row_split ='tumor_het_phenograph_scaled_consolidated' # tumor_het_phenograph_scaled_consolidated patient_group

# subset the adata object based on the rows in supplied interaction df
adata = adata[adata.obs[each_row].isin(interaction_data.index)]

# interaction data
interaction_data.to_csv('for R/interaction_data.csv', index=True)

# median_cellfeatures
median_cellfeatures = pd.DataFrame(adata.obs)[['Area', 'Eccentricity', 'Extent', 'Solidity']]
median_cellfeatures = median_cellfeatures.apply(pd.to_numeric)
median_cellfeatures['phenotype'] = adata.obs[each_row].values
median_cellfeatures = median_cellfeatures.groupby('phenotype').median()
# normalize the data
median_cellfeatures = median_cellfeatures.apply(lambda x: (x - x.min()), axis=0).apply(lambda x: (x / x.max()), axis=0)
median_cellfeatures.to_csv('for R/median_cellfeatures.csv', index=True)

# phenotype distribution for ROI
#subset_yaxis = [ 'aj_IR', 'aj_BTIL', 'aj_IB', 'aj_MIS', 'aj_IT', 'aj_ET']
pheno_prop = sm.pl.stacked_barplot (adata,x_axis=each_row,y_axis=within_barplot,method='percent',return_data=True)
pheno_prop.to_csv('for R/pheno_prop.csv', index=True)

# Pheno group
pheno_group = pd.DataFrame(adata.obs[[each_row,row_split]])
pheno_group = pheno_group.drop_duplicates(subset=[each_row])
pheno_group.to_csv('for R/pheno_group.csv', index=False)


# %% stromal proportion subtyping patients
np.unique(adata.obs['phenotype_aj'])
sub_y = ['B cells', 'CD4 T cells', 'CD8 T cells',
       'Exhausted T cells', 'Follicular Dendritic cells',
       'Follicular Helper T cells', 
       'M1 Macrophages', 'M2 Macrophages', 'Myeloid Dendritic cells',
       'Other Macrophages','Regulatory T cells']
immune_composition = sm.pl.stacked_barplot (ptcl_nos,x_axis='patientid',y_axis='phenotype_aj',subset_yaxis=sub_y, method='percent', plot_tool='matplotlib',return_data=True)
immune_composition.to_csv('for R/immune_prop.csv', index=True)



# %% phenotype expression plot while normalizing for number of cells within each patient


os.chdir('/Users/aj/Dropbox (Partners HealthCare)/PTCL_Jackson/tumor_subgroups/PTCL/heatmaps') # ALCL AITL PTCL
adata = ptcl_tumor.copy()

# input
each_row = 'phenotype_tumor_het' # tumor_het_phenograph_scaled patientid tumor_het
within_barplot = 'patientid' #patientid phenotype
row_split ='myeloid_lymphoid' # tumor_het_phenograph_scaled_consolidated patient_group

# median_expression
comb_data = pd.DataFrame(adata.X, columns=adata.var.index)
comb_data['phenotype'] = adata.obs[each_row].values
median_expression = comb_data.groupby('phenotype').median()
median_expression.to_csv('for R/median_expression.csv', index=True)

# median_cellfeatures
median_cellfeatures = pd.DataFrame(adata.obs)[['Area', 'Eccentricity', 'Extent', 'Solidity']]
median_cellfeatures = median_cellfeatures.apply(pd.to_numeric)
median_cellfeatures['phenotype'] = adata.obs[each_row].values
median_cellfeatures = median_cellfeatures.groupby('phenotype').median()
# normalize the data
median_cellfeatures = median_cellfeatures.apply(lambda x: (x - x.min()), axis=0).apply(lambda x: (x / x.max()), axis=0)
median_cellfeatures.to_csv('for R/median_cellfeatures.csv', index=True)

# phenotype distribution for ROI
#subset_yaxis = [ 'aj_IR', 'aj_BTIL', 'aj_IB', 'aj_MIS', 'aj_IT', 'aj_ET']
data = pd.DataFrame(adata.obs)[[each_row,within_barplot]].astype(str)

# Normalize for total cells
total_cells = data.groupby([each_row,within_barplot]).size().unstack().fillna(0).sum(axis=0)
xx = data.groupby([each_row,within_barplot]).size().unstack().fillna(0).div(total_cells, axis=1)#.stack()
# Normalize within each cell type (all rows between 0-1)
rg = xx.div(xx.sum(axis=1), axis=0)
rg.to_csv('for R/pheno_prop.csv', index=True)

# Pheno group
pheno_group = pd.DataFrame(adata.obs[[each_row,row_split]])
pheno_group = pheno_group.drop_duplicates(subset=[each_row])
pheno_group.to_csv('for R/pheno_group.csv', index=False)



adata.var.index.get_loc('CD3D')

adata.raw[:,adata.var.index.get_loc('CD3D')].X
np.log1p(55353)





