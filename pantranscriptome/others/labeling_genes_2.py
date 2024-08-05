# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 21:42:37 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# get genes that for BLAST search on the NCBI database
# =============================================================================

def get_gene_list(df, save_file_path):
    gene_list = []
    for index, name in df['NCBI_gene_description'].items():
        if pd.isna(name):
            gene_list.append(df.at[index, 'sequence_name'])
        if 'hypothetical protein' in str(name):
            gene_list.append(df.at[index, 'sequence_name'])  
        if 'uncharacterized' in str(name):
            gene_list.append(df.at[index, 'sequence_name'])  
    
    with open(save_file_path, 'w') as fh:
        for gene in gene_list:
            fh.write(f'{gene}\n')

treat_genes = pd.read_csv("others/labeling_genes/treatment/treatment_genes_description.csv")
#get_gene_list(treat_genes, 'others/labeling_genes/treatment/treatment_genes_list.txt')

type_genes = pd.read_csv("others/labeling_genes/type/type_genes_description.csv")
#get_gene_list(type_genes, 'others/labeling_genes/type/type_genes_list.txt')

int_genes = pd.read_csv("others/labeling_genes/interaction/interaction_genes_description.csv")
#get_gene_list(int_genes, 'others/labeling_genes/interaction/interaction_genes_list.txt')
