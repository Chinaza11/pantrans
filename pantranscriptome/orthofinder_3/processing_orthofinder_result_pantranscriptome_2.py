# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 10:23:10 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/orthofinder_3")

# =============================================================================
# keeping only the oldest node for genes with multiple duplication events
# =============================================================================

annot = pd.read_csv('annot_node_ckpt5.csv')
annot['Species.Tree.Node_2'] = annot['Species.Tree.Node']

guide = {1:'N0', 2:'N1', 3:'N2', 4:'N4', 5:'N6', 6:'N9', 7:'N12', 8:'N14',
         9:'N15', 10:'N16', 11:'N17'}

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


annot.to_csv('annot_node_ckpt6.csv', index=False)