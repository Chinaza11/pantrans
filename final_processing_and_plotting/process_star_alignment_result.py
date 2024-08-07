# -*- coding: utf-8 -*-
"""
Created on Fri May 17 15:16:40 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result")

# =============================================================================
# Pan-transcriptome
# =============================================================================

pantrans = pd.read_csv("final_processing_and_plotting/pantranscriptome_star_alignment_plot.csv")

pantrans['total'] = pantrans[['Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped: too many mismatches', 'Unmapped: too short', 'Unmapped: other']].sum(axis=1)
pantrans['total_round'] = round(pantrans['total']/1000000)

pantrans['Uniquely mapped - Percent'] = round(pantrans['Uniquely mapped']/pantrans['total'] * 100, 1)
pantrans['Mapped to multiple loci - Percent'] = round(pantrans['Mapped to multiple loci']/pantrans['total'] * 100, 1)
pantrans['Mapped to too many loci - Percent'] = round(pantrans['Mapped to too many loci']/pantrans['total'] * 100, 1)
pantrans['Unmapped: too many mismatches - Percent'] = round(pantrans['Unmapped: too many mismatches']/pantrans['total'] * 100, 1)
pantrans['Unmapped: too short - Percent'] = round(pantrans['Unmapped: too short']/pantrans['total'] * 100, 1)
pantrans['Unmapped: other - Percent'] = round(pantrans['Unmapped: other'] / pantrans['total'] * 100, 1)

subset_pantrans = pantrans[['Category', 'Uniquely mapped - Percent', 'Mapped to multiple loci - Percent', 'Mapped to too many loci - Percent', 'Unmapped: too many mismatches - Percent', 'Unmapped: too short - Percent', 'Unmapped: other - Percent']]
subset_pantrans = subset_pantrans.rename(columns={'Category': 'Samples'})

melted_pantrans = pd.melt(subset_pantrans, id_vars=['Samples'], var_name='Categories', value_name='Percentage')

melted_pantrans['Cat'] = 'PT'

# =============================================================================
# Rio reference genome
# =============================================================================
rio = pd.read_csv("final_processing_and_plotting/rio-ref-genome_star_alignment_plot.csv")

rio['total'] = rio[['Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped: too short', 'Unmapped: other']].sum(axis=1)
rio['total_round'] = round(rio['total']/1000000)

rio['Uniquely mapped - Percent'] = round(rio['Uniquely mapped']/rio['total'] * 100, 1)
rio['Mapped to multiple loci - Percent'] = round(rio['Mapped to multiple loci']/rio['total'] * 100, 1)
rio['Mapped to too many loci - Percent'] = round(rio['Mapped to too many loci']/rio['total'] * 100, 1)
rio['Unmapped: too short - Percent'] = round(rio['Unmapped: too short']/rio['total'] * 100, 1)
rio['Unmapped: other - Percent'] = round(rio['Unmapped: other'] / rio['total'] * 100, 1)

subset_rio = rio[['Category', 'Uniquely mapped - Percent', 'Mapped to multiple loci - Percent', 'Mapped to too many loci - Percent', 'Unmapped: too short - Percent', 'Unmapped: other - Percent']]
subset_rio = subset_rio.rename(columns={'Category': 'Samples'})

melted_rio = pd.melt(subset_rio, id_vars=['Samples'], var_name='Categories', value_name='Percentage')

melted_rio['Cat'] = 'Rio'

# =============================================================================
# Combine both result in preparation for plotting in R
# =============================================================================

combined_df = pd.concat([melted_pantrans, melted_rio])
combined_df = combined_df.reset_index(drop=True)

for index, cell in combined_df['Categories'].items():
    combined_df.at[index,'Categories'] = cell[:-10]
        
replacement_dict = {
    '069-con1-exp2': '069-con1',
    '069-con2-exp2': '069-con2',
    '069-no1-exp2': '069-no1',
    '069-no2-exp2': '069-no2',
    '972-con1-exp2': '972-con1',
    '972-con3-exp2': '972-con2',
    '972-no2-exp2': '972-no1',
    '972-no3-exp2': '972-no2',
    'Gra-con1-exp2': 'Gra-con1',
    'Gra-con2-exp2': 'Gra-con2',
    'Gra-con3-exp1': 'Gra-con3',
    'Gra-no1-exp2': 'Gra-no1',
    'Gra-no2-exp2': 'Gra-no2',
    'Gra-no3-exp2': 'Gra-no3',
    'Leo-con3-exp1': 'Leo-con1',
    'Leo-con3-exp2': 'Leo-con2',
    'Leo-no2-exp1': 'Leo-no1',
    'Leo-no3-exp2': 'Leo-no2',
    'Rio-con3-exp2': 'Rio-con1',
    'Rio-no1-exp2': 'Rio-no1'
}

combined_df['Samples'] = combined_df['Samples'].replace(replacement_dict)

combined_df.to_csv('final_processing_and_plotting/star_alignment.csv', index=False)

# =============================================================================
# We decided to modify the final figure and group all unmapped categories 
# together
# =============================================================================

combined_df['Categories'] = combined_df['Categories'].str.split(':').str[0]
combined_df = combined_df.groupby(['Samples', 'Categories', 'Cat'], as_index=False)['Percentage'].sum()

combined_df.to_csv('final_processing_and_plotting/star_alignment_2.csv', index=False)
