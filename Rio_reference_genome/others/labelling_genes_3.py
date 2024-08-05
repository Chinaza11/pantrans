# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 14:44:10 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/others/labeling_genes/")

# =============================================================================
# adding BLAST result to the fixed-effect_genes_description_2 dataframe
# =============================================================================

def add_blast_result_to_df(genes_df_path, fixed_effect):
    genes_df = pd.read_csv(genes_df_path)
    for index, name in genes_df['protein_genbank_id'].items():
        if pd.isna(name) != True:
            file = f'{fixed_effect}/blast_result/{name}.txt'
            column_names = ['qseqid', 'sseqid', 'pident', 'evalue', 'stitle', 'length']
            df = pd.read_table(file, index_col=None, header=None, names=column_names)
            blast_result_list = []
            for i in range(df.shape[0]):
                blast_result = str(df.iloc[i,2]) + "; " + str(df.iloc[i,3]) + "; " + str(df.iloc[i,4])
                blast_result_list.append(blast_result)
            for i in range(len(blast_result_list)):
                genes_df.at[index, f'blast_result{i+1}'] = blast_result_list[i] 
    
    genes_df = genes_df.iloc[:, :14] 
    return genes_df


treat_genes = add_blast_result_to_df("treatment/treatment_genes_description_2.csv",
                                   "treatment")
treat_genes.to_csv("treatment/treatment_genes_description_3.csv", index=False, header=True)


type_genes = add_blast_result_to_df("type/type_genes_description_2.csv",
                                   "type")
type_genes.to_csv("type/type_genes_description_3.csv", index=False, header=True)


int_genes = add_blast_result_to_df("interaction/interaction_genes_description_2.csv",
                                   "interaction")
int_genes.to_csv("interaction/interaction_genes_description_3.csv", index=False, header=True)
