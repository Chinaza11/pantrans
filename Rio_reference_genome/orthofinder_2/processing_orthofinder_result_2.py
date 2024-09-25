# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 18:51:16 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/")

# =============================================================================
# keeping only the oldest node for genes with multiple duplication events
# =============================================================================

annot = pd.read_csv('orthofinder_2/annot_node_ckpt5.csv')
annot['Species.Tree.Node_2'] = annot['Species.Tree.Node']

guide = {1:'N0', 2:'N2', 3:'N4', 4:'N7', 5:'N11', 6:'N13', 7:'N14', 8:'sorghumRio'}

for index, cell in annot['Species.Tree.Node'].items():
    if pd.isna(cell) != True:
        if 'SingleCopyOrthologues' not in cell:
            if 'UnassignedGenes' not in cell:
                if ',' in cell:
                    y = cell.split(',')
                    l = []
                    for i in y:
                        for key, val in guide.items():
                            if i == val:
                                l.append(key)
                    new_entry = guide[min(l)]
                    annot.at[index, 'Species.Tree.Node_2'] = new_entry

# =============================================================================
# convert empty cells to the string value "NA" for GO analysis in R
# =============================================================================

for index, cell in annot['Species.Tree.Node_2'].items():
    if pd.isna(cell) == True:
        annot.at[index, 'Species.Tree.Node_2'] = 'NA'
        
(annot['Species.Tree.Node_2'] == 'NA').sum()


annot.to_csv('orthofinder_2/annot_node_ckpt6.csv', index=False)
