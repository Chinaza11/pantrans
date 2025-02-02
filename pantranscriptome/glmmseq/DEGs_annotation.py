# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 19:45:26 2025

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome")

# =============================================================================
# 
# =============================================================================

pantranscriptome_df = pd.read_csv("pantranscriptome_construction/pantranscriptome_plus_gene_type.csv")

core_genes = pantranscriptome_df[pantranscriptome_df['annotation'] == 'core']['sequence_name'].tolist()
shell_genes = pantranscriptome_df[pantranscriptome_df['annotation'] == 'shell']['sequence_name'].tolist()
cloud_genes = pantranscriptome_df[pantranscriptome_df['annotation'] == 'cloud']['sequence_name'].tolist()

# =============================================================================
# gene annotation function
# =============================================================================

def gene_annotation(dataframe, column_name):
    for index, gene_name in dataframe[column_name].items():
        if "--" not in gene_name:
            if gene_name in core_genes:
                dataframe.at[index, 'annotation'] = 'core'
            if gene_name in shell_genes:
                dataframe.at[index, 'annotation'] = 'shell'
            if gene_name in cloud_genes:
                dataframe.at[index, 'annotation'] = 'cloud'
        else:
            gene_types = []
            for i in range(len(gene_name.split("--"))):
                j = gene_name.split("--")[i]
                if j in core_genes:
                    gene_types.append('core')
                elif j in shell_genes:
                    gene_types.append('shell')
                elif j in cloud_genes:
                    gene_types.append('cloud')
            dataframe.at[index, 'annotation'] = ', '.join(list(set(gene_types)))            

# =============================================================================
# treatment genes
# =============================================================================

treat_genes = pd.read_csv("treat_genes_df.txt", delimiter='\t')
    
gene_annotation(treat_genes, 'treat_genes')
treat_genes['annotation'].value_counts()
treat_genes.to_csv('treat_genes_plus_gene_type.csv', index=False)

# =============================================================================
# type genes (after permutation test)
# =============================================================================

type_genes = pd.read_csv("glmmseq/permutation_analysis/type_genes_df.txt", delimiter='\t')
    
gene_annotation(type_genes, 'type_genes')
type_genes['annotation'].value_counts()
type_genes.to_csv('glmmseq/permutation_analysis/type_genes_plus_gene_type.csv', index=False)

# =============================================================================
# interaction genes (after permutation test)
# =============================================================================

int_genes = pd.read_csv("glmmseq/permutation_analysis/int_genes_df.txt", delimiter='\t')
    
gene_annotation(int_genes, 'int_genes')
int_genes['annotation'].value_counts()
int_genes.to_csv('glmmseq/permutation_analysis/int_genes_plus_gene_type.csv', index=False)
