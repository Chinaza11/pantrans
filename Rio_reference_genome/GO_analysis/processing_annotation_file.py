# -*- coding: utf-8 -*-
"""
Created on Sat May 18 14:07:11 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/GO_analysis/")

df = pd.read_table("Rio_annotationInfo-reduced.txt")

subset_df = df[['locusName', 'GO']]

# remove duplicate locusName and compile the GO terms into one cell
merged_df = subset_df.groupby('locusName')['GO'].apply(lambda x: ','.join(x.astype(str))).reset_index()

# remove duplicate GO terms in each row
for index, cell in merged_df['GO'].items():
    if cell != "nan":
        lst = cell.split(',')
        if len(lst) != len(set(lst)):
            merged_df.at[index, 'GO'] = ','.join(list(set(lst)))

# convert 'nan' cells in the GO column to 'NA' for easy manipulation in R
for index, cell in merged_df['GO'].items():
    if cell == "nan":
        merged_df.at[index, 'GO'] = 'NA'
        
# merged_df.to_csv('Rio_annotationInfo-reduced_deduplicate.txt', sep='\t', index=False)

# this did not change previous GO result, seems topGO accounted for duplicates somehow