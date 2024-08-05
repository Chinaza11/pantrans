# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:46:48 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/others/labeling_genes/")

# =============================================================================
# 
# =============================================================================

def add_blast_result_to_df(genes_df_path, fixed_effect):
    genes_df = pd.read_csv(genes_df_path, names=['genes'])
    for index, name in genes_df['genes'].items():
        print(str(index) + ':  ' + name)
        file = f'{fixed_effect}/blast_result/{name}.txt'
        column_names = ['qseqid', 'sseqid', 'pident', 'evalue', 'stitle', 'length']
        df = pd.read_table(file, index_col=None, header=None, names=column_names)
        blast_result_list = []
        for i in range(df.shape[0]):
            blast_result = str(df.iloc[i,2]) + "; " + str(df.iloc[i,3]) + "; " + str(df.iloc[i,4])
            blast_result_list.append(blast_result)
        if len(blast_result_list) > 10:
            blast_result_list = blast_result_list[0:10]
        for i in range(len(blast_result_list)):
            genes_df.at[index, f'blast_result{i+1}'] = blast_result_list[i] 
    
    genes_df = genes_df.iloc[:, :16] 
    return genes_df


treat_genes_blast = add_blast_result_to_df("treatment/treatment_genes_list.txt",
                                   "treatment")

type_genes_blast = add_blast_result_to_df("type/type_genes_list.txt",
                                   "type")

int_genes_blast = add_blast_result_to_df("interaction/interaction_genes_list.txt",
                                   "interaction")

# =============================================================================
# 
# =============================================================================

treat_genes_main = pd.read_csv("treatment/treatment_genes_description.csv")
treat_genes_main = pd.merge(treat_genes_main, treat_genes_blast, left_on='sequence_name', right_on='genes', how='left').drop('genes', axis=1)
treat_genes_main.to_csv("treatment/treatment_genes_description_2.csv", index=False, header=True)

type_genes_main = pd.read_csv("type/type_genes_description.csv")
type_genes_main = pd.merge(type_genes_main, type_genes_blast, left_on='sequence_name', right_on='genes', how='left').drop('genes', axis=1)
type_genes_main.to_csv("type/type_genes_description_2.csv", index=False, header=True)

int_genes_main = pd.read_csv("interaction/interaction_genes_description.csv")
int_genes_main = pd.merge(int_genes_main, int_genes_blast, left_on='sequence_name', right_on='genes', how='left').drop('genes', axis=1)
int_genes_main.to_csv("interaction/interaction_genes_description_2.csv", index=False, header=True)

# =============================================================================
# 
# =============================================================================