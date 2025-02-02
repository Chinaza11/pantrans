# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 14:37:00 2025

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/pantranscriptome_construction")

pantranscriptome_df = pd.read_csv("pantranscriptome.csv")

# =============================================================================
# annotating cloud genes
# =============================================================================

pav_only_one_cds = pd.read_csv("presence_absence_table_only_one_cds.csv")

for index, stable_id in pantranscriptome_df['stable_id'].items():
    if stable_id in pav_only_one_cds['stable ids'].values:
        pantranscriptome_df.at[index, 'annotation'] = 'cloud'
        
# =============================================================================
# annotating core genes
# =============================================================================

pav_more_than_one_cds = pd.read_csv("presence_absence_table_more_than_one_cds.csv")
        
col_names = ['ChiAmber_cds_length', 'Grassl_cds_length', 'Leoti_cds_length', 
             'pi229841_cds_length', 'pi297155_cds_length', 'pi329311_cds_length',
             'pi506069_cds_length', 'pi510757_cds_length', 'pi655972_cds_length',
             'RioNAM_cds_length']

for index, stable_id in pav_more_than_one_cds['stable ids'].items():
    if all(not pd.isna(pav_more_than_one_cds[col][index]) for col in col_names):
        pav_more_than_one_cds.at[index, 'annotation'] = 'core'
        
core_genes_df = pav_more_than_one_cds[pav_more_than_one_cds['annotation'] == 'core']
core_genes = core_genes_df['stable ids'].tolist()

for index, stable_id in pantranscriptome_df['stable_id'].items():
    if stable_id in core_genes:
        pantranscriptome_df.at[index, 'annotation'] = 'core'
        
# =============================================================================
# annotating shell genes
# =============================================================================

for index, annotation in pav_more_than_one_cds['annotation'].items():
    if annotation == 'nan':
        pav_more_than_one_cds.at[index, 'annotation'] = 'shell'
        
shell_genes_df = pav_more_than_one_cds[pav_more_than_one_cds['annotation'] == 'shell']
shell_genes = shell_genes_df['stable ids'].tolist()

for index, stable_id in pantranscriptome_df['stable_id'].items():
    if stable_id in shell_genes:
        pantranscriptome_df.at[index, 'annotation'] = 'shell'
    
# =============================================================================
# 
# =============================================================================

pantranscriptome_df['annotation'].value_counts()

# =============================================================================
# save dataframe
# =============================================================================

pantranscriptome_df.to_csv('pantranscriptome_plus_gene_type.csv', index=False)
