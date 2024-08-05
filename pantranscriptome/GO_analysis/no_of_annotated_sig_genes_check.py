# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:25:36 2024

@author: Chinaza
"""

import os, pandas as pd, numpy as np

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# checking to see the amount of genes that are annotated among the significant genes
# =============================================================================

pantranscriptome = pd.read_csv("GO_analysis/pantranscriptome_plus_GO_cleaned.csv")
pantranscriptome.head()

# importing all the hybrid genes
header = ('Gene', 'qvals.Treatment', 'qvals.Type', 'qvals.Treatment.Type', 'GO')
hyb_genes_annot = pd.read_table("GO_analysis/hyb_genes.txt", names=header)
hyb_genes_annot.head()

#---------------------------------
# function for counting rows without annotation

def get_na_values(dataframe):

    dataframe["GO_annotations"] = np.nan
    
    na_count1 = dataframe['GO_annotations'].isna().sum()
    print(f"Number of NA values in 'GO_annotations' column (before anything is done): {na_count1}")
    
    # dealing with single genes
    for index, gene in dataframe['Gene'].items():
        for index2, gene2 in pantranscriptome['sequence_name'].items():
            if gene == gene2:
                dataframe.at[index, 'GO_annotations'] = pantranscriptome['GO_annotations_cleaned'][index2]
    
    # na count after dealing with single genes
    na_count2 = dataframe['GO_annotations'].isna().sum()
    
    print(f"Number of NA values in 'GO_annotations' column (after dealing with single genes): {na_count2}")

    # dealing with hybrid genes and getting total counts
    count = 0
    for index3, gene3 in dataframe['Gene'].items():
        for index4, gene4 in hyb_genes_annot['Gene'].items():
            if gene3 == gene4:
                count += 1
                dataframe.at[index3, 'GO_annotations'] = hyb_genes_annot['GO'][index4]
    
    print(f"Number of significant hybrid genes: {count}")
    
    # na count after dealing with hybrid genes, total counts
    na_count3 = dataframe['GO_annotations'].isna().sum()
    print(f"Number of NA values in 'GO_annotations' column (after dealing with single and hybrid genes): {na_count3}")
    
    print(f"Number of rows annotated is {na_count1-na_count3} out of {na_count1} rows; meaning, at least {na_count3} genes are not annotated")
    return dataframe

# =============================================================================
# Interaction
# =============================================================================

header2 = ('Gene', 'qvals.Treatment.Type')
sig_genes_int = pd.read_table("glmmseq/result_qval_int.txt", names=header2)
sig_genes_int = sig_genes_int.iloc[1:]
sig_genes_int.reset_index(drop=True, inplace=True)

sig_genes_int.head()

print("\n Interaction - statistics \n")
result_int = get_na_values(sig_genes_int)


# =============================================================================
# Type
# =============================================================================

header2 = ('Gene', 'qvals.Type')
sig_genes_type = pd.read_table("glmmseq/result_qval_type.txt", names=header2)
sig_genes_type = sig_genes_type.iloc[1:]
sig_genes_type.reset_index(drop=True, inplace=True)

print("\n Type - statistics \n")
result_type = get_na_values(sig_genes_type)


# =============================================================================
# Treatment
# =============================================================================

header2 = ('Gene', 'qvals.Treatment')
sig_genes_trt = pd.read_table("glmmseq/result_qval_trt.txt", names=header2)
sig_genes_trt = sig_genes_trt.iloc[1:]
sig_genes_trt.reset_index(drop=True, inplace=True)

print("\n Treatment - statistics \n")
result_trt = get_na_values(sig_genes_trt)

