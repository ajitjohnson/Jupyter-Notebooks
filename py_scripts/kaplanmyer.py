#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 21:17:12 2021

@author: aj
"""


bdata = combined_excluded_RJP.copy()

import kaplanmeier as km


# defining groups
IM = combined_excluded_RJP[combined_excluded_RJP.obs['stage_consolidated'] == 'MIS']
#IM = combined_excluded_RJP.copy()
IM = IM[IM.obs['lda_consolidated_MC'] == 'MC9']

organiz_by = 'phenotype'

# frequency calculation
mc = IM.obs[['imageid', 'patient_ID', 'lda_consolidated_MC', 'phenotype']]

# frequency calculation
x = pd.DataFrame(mc.groupby(['imageid',organiz_by]).size().fillna(0))
x.reset_index(inplace=True) 
# add total cells for each image in a new column
total_cells = combined_excluded_RJP.obs['imageid'].value_counts()
t = dict(zip(total_cells.index, total_cells.values))
# add a column to x to duplicate imageid
x['total_cells'] = x['imageid']
x = x.replace({"total_cells": t})
# normalize counts
x['Frequency'] = (x[0]/x['total_cells']) * 100
x = x[['imageid', organiz_by,'Frequency']]
# add patirnt id
p = dict(zip(mc.imageid, mc.patient_ID))
x['patient_ID'] = x['imageid']
x = x.replace({"patient_ID": p})
x = x[['patient_ID',organiz_by,'Frequency']]
x = pd.DataFrame(x.groupby(['patient_ID', organiz_by])['Frequency'].agg('sum').fillna(0))
x.reset_index(inplace=True) 
median_expression = x.pivot(index='patient_ID', columns=organiz_by, values='Frequency')



# loop
for i in list(median_expression.columns):
    median_expression[i].loc[median_expression[i] <= median_expression[i].median()] = 0
    median_expression[i].loc[median_expression[i] > median_expression[i].median()] = 1
    
    
# merge with df
df = pd.read_csv('/Users/aj/Dropbox (Partners HealthCare)/PCA_Atlas1-2020/Data/survival.csv', index_col=0)  
df.index = df.index.astype('str')
df = df.merge(median_expression, left_index=True, right_index=True, how='inner')


#list(median_expression.columns)


# import the survival data

# Compute 
for i in list(median_expression.columns): 
    out=km.fit(df['﻿time'], df['PFS'], df[i])
    km.plot(out, title= str(i))
    
    
    
# legacy

sub_y = ['MC1','MC2','MC3','MC4','MC5','MC6','MC7','MC8','MC9','MC10']
median_expression = sm.pl.stacked_barplot (IM,x_axis='patient_ID',y_axis='lda_consolidated_MC', 
                                           subset_yaxis=sub_y,order_yaxis=sub_y,
                                           method='percent', plot_tool='matplotlib', return_data=True)

median_expression = sm.pl.stacked_barplot (IM,x_axis='patient_ID',y_axis='phenotype', 
                                           method='percent', plot_tool='matplotlib', return_data=True)


