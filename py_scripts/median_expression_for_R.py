#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 15:51:35 2020
@author: Ajit Johnson Nirmal
"""

# %% median expression of cell phenotyping for R

adata = alcl.copy()

# set wd
import os
os.chdir('/Volumes/SSD/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/re_staining_tma/5-TMA')

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


# %% Using scaled data

adata = alcl.copy()
adata = aitl.copy()
adata = adata[adata.obs['tumor_stroma'] == 'Tumor']

# input
each_row = 'tumor_het_phenograph_scaled' # tumor_het_phenograph_scaled patientid
within_barplot = 'patientid' #patientid phenotype
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
#subset_yaxis = [ 'aj_IR', 'aj_BTIL', 'aj_IB', 'aj_MIS', 'aj_IT', 'aj_ET']
pheno_prop = sm.pl.stacked_barplot (adata,x_axis=each_row,y_axis=within_barplot,method='percent',return_data=True)
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





















