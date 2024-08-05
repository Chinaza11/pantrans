# -*- coding: utf-8 -*-
"""
Created on Sun May 19 08:29:52 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result")

# =============================================================================
# 
# =============================================================================

def parse_mmquant_result(path):
    with open(path) as file:
        for line in file:
            if line.startswith("Stats") and line.startswith("Stats file") != True:
                lst = line.split(' ')
                cleaned_lst1 = [item for item in lst if item.strip()]
            if line.startswith("# hits"):
                lst = line.split(' ')
                cleaned_lst2 = [item for item in lst if item.strip()]
            if line.startswith("# uniquely mapped reads"):
                lst = line.split(' ')
                int_lst = [item for item in lst if item.strip()]
                s1 = "".join(int_lst[4:])
                cleaned_lst = s1.split(')')
                cleaned_lst3 = [item for item in cleaned_lst if item.strip()]
            if line.startswith("# ambiguous hits"):
                lst = line.split(' ')
                int_lst = [item for item in lst if item.strip()]
                s1 = "".join(int_lst[3:])
                cleaned_lst = s1.split(')')
                cleaned_lst4 = [item for item in cleaned_lst if item.strip()]
            if line.startswith("# non-uniquely mapped hits"):
                lst = line.split(' ')
                int_lst = [item for item in lst if item.strip()]
                s1 = "".join(int_lst[4:])
                cleaned_lst = s1.split(')')
                cleaned_lst5 = [item for item in cleaned_lst if item.strip()]
            if line.startswith("# unassigned hits"):
                lst = line.split(' ')
                int_lst = [item for item in lst if item.strip()]
                s1 = "".join(int_lst[3:])
                cleaned_lst = s1.split(')')
                cleaned_lst6 = [item for item in cleaned_lst if item.strip()]
    
    df = pd.DataFrame({
        'samples': cleaned_lst1[1:],
        'no_of_hits': cleaned_lst2[2:],
        'uniquely_mapped_reads': cleaned_lst3,
        'ambiguous_hits': cleaned_lst4,
        'non_uniquely_mapped_hits': cleaned_lst5,
        'unassigned_hits': cleaned_lst6,
    })
    
    df[['uniquely_mapped_reads_no', 'uniquely_mapped_reads_percent']] = df['uniquely_mapped_reads'].str.split("(", expand=True)
    df[['ambiguous_hits_no', 'ambiguous_hits_percent']] = df['ambiguous_hits'].str.split("(", expand=True)
    df[['non_uniquely_mapped_hits_no', 'non_uniquely_mapped_hits_percent']] = df['non_uniquely_mapped_hits'].str.split("(", expand=True)
    df[['unassigned_hits_no', 'unassigned_hits_percent']] = df['unassigned_hits'].str.split("(", expand=True)
    
    df = df[['samples', 'no_of_hits', 'uniquely_mapped_reads_no', 'uniquely_mapped_reads_percent', 'ambiguous_hits_no', 'ambiguous_hits_percent', 'non_uniquely_mapped_hits_no', 'non_uniquely_mapped_hits_percent', 'unassigned_hits_no', 'unassigned_hits_percent']]
    
    columns_to_convert = ['no_of_hits', 'uniquely_mapped_reads_no', 'ambiguous_hits_no', 'non_uniquely_mapped_hits_no', 'unassigned_hits_no']
    df[columns_to_convert] = df[columns_to_convert].astype(int)
    df[columns_to_convert] = round(df[columns_to_convert]/1000000, 2)
    
    return df

# =============================================================================
# Rio reference genome
# =============================================================================

rio = parse_mmquant_result("Rio_reference_genome/from_the_cluster/mmquant-slurm-8405642.out")

subset_rio_percent = rio[['samples', 'uniquely_mapped_reads_percent', 'ambiguous_hits_percent', 'non_uniquely_mapped_hits_percent', 'unassigned_hits_percent']]

melted_rio_percent = pd.melt(subset_rio_percent, id_vars=['samples'], var_name='Categories', value_name='Percentage')

melted_rio_percent['Cat'] = 'Rio'
#------------------------------------------------------------------------------

subset_rio_no = rio[['samples', 'uniquely_mapped_reads_no', 'ambiguous_hits_no', 'non_uniquely_mapped_hits_no', 'unassigned_hits_no']]

melted_rio_no = pd.melt(subset_rio_no, id_vars=['samples'], var_name='Categories', value_name='Million')

melted_rio_no['Cat'] = 'Rio'

# =============================================================================
# Pan-transcriptome
# =============================================================================

pant = parse_mmquant_result("pantranscriptome/from_the_cluster/mmquant-slurm-6734901.out")

subset_pantrans_percent = pant[['samples', 'uniquely_mapped_reads_percent', 'ambiguous_hits_percent', 'non_uniquely_mapped_hits_percent', 'unassigned_hits_percent']]

melted_pantrans_percent = pd.melt(subset_pantrans_percent, id_vars=['samples'], var_name='Categories', value_name='Percentage')

melted_pantrans_percent['Cat'] = 'PT'

#------------------------------------------------------------------------------

subset_pantrans_no = pant[['samples', 'uniquely_mapped_reads_no', 'ambiguous_hits_no', 'non_uniquely_mapped_hits_no', 'unassigned_hits_no']]

melted_pantrans_no = pd.melt(subset_pantrans_no, id_vars=['samples'], var_name='Categories', value_name='Million')

melted_pantrans_no['Cat'] = 'PT'

# =============================================================================
# Combine both result in preparation for plotting in R
# =============================================================================

combined_percent_df = pd.concat([melted_pantrans_percent, melted_rio_percent])
combined_percent_df = combined_percent_df.reset_index(drop=True)

for index, cell in combined_percent_df['Categories'].items():
    combined_percent_df.at[index,'Categories'] = cell[:-8]
    combined_percent_df.at[index,'Percentage'] = combined_percent_df.at[index,'Percentage'][:-1]

combined_percent_df['samples'] = combined_percent_df['samples'].replace({'069-n02': '069-no2'})

combined_percent_df.to_csv('mmquant_percent.csv', index=False)

#------------------------------------------------------------------------------

combined_no_df = pd.concat([melted_pantrans_no, melted_rio_no])
combined_no_df = combined_no_df.reset_index(drop=True)

for index, cell in combined_no_df['Categories'].items():
    combined_no_df.at[index,'Categories'] = cell[:-3]

combined_no_df['samples'] = combined_no_df['samples'].replace({'069-n02': '069-no2'})

combined_no_df.to_csv('mmquant_no.csv', index=False)