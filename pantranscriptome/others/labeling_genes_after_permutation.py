# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 19:14:06 2025

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# 
# =============================================================================

type_all_annotation = pd.read_excel('final_processing_and_plotting/SUPPLEMENTARY FILE 3.xlsx', sheet_name='Pan-transcriptome Type')

type_genes = pd.read_csv("pantranscriptome/glmmseq/permutation_analysis/type_genes_df.txt", delimiter='\t')

type_genes = type_genes['type_genes'].str.split('--').explode().reset_index(drop=True).to_frame()

type_annot_subset = type_all_annotation.merge(type_genes, 
                                      left_on='sequence_name',
                                      right_on='type_genes',
                                      how='inner')

type_annot_subset = type_annot_subset.drop(columns=['type_genes'])

type_annot_subset.to_csv('pantranscriptome/others/labeling_genes/type/type_genes_description-after_permutation.csv', index=False)

# =============================================================================
# 
# =============================================================================

int_all_annotation = pd.read_excel('final_processing_and_plotting/SUPPLEMENTARY FILE 3.xlsx', sheet_name='Pan-transcriptome Interaction')

int_genes = pd.read_csv("pantranscriptome/glmmseq/permutation_analysis/int_genes_df.txt", delimiter='\t')

int_genes = int_genes['int_genes'].str.split('--').explode().reset_index(drop=True).to_frame()

int_annot_subset = int_all_annotation.merge(int_genes, 
                                      left_on='sequence_name',
                                      right_on='int_genes',
                                      how='inner')

int_annot_subset = int_annot_subset.drop(columns=['int_genes'])

int_annot_subset.to_csv('pantranscriptome/others/labeling_genes/interaction/interaction_genes_description-after_permutation.csv', index=False)