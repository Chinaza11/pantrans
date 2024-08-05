# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 09:40:01 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/GO_analysis/")

# =============================================================================
# cleaning up pantranscriptome_with_go_terms.csv "GO_annotations" column
# =============================================================================

pantranscriptome = pd.read_csv("intermediate_files/pantranscriptome_with_go_terms.csv")

print(pantranscriptome.head())

pantranscriptome['GO_annotations_cleaned'] = ""

for index, GO_term in pantranscriptome['GO_annotations'].items():
    if pd.isna(GO_term) != True:
        new = GO_term.split(',')
        unique_go_terms_lst = []
        unique_go_terms_str = ""
        for i in new:
            k = i.split('|')
            for m in k:
                m = m.strip()
                if m not in unique_go_terms_lst:
                    unique_go_terms_lst.append(m)
        unique_go_terms_str = ','.join(unique_go_terms_lst)
        pantranscriptome.at[index, 'GO_annotations_cleaned'] = unique_go_terms_str
    
    #putting NA in all cells that do not have annotations so topGO package in R reads it properly
    if pd.isna(GO_term) == True:
        pantranscriptome.at[index, 'GO_annotations_cleaned'] = "NA"
    
pantranscriptome.to_csv('pantranscriptome_plus_GO_cleaned.csv', index=False)