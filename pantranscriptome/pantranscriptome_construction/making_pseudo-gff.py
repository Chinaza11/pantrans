# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 22:38:24 2024

@author: Chinaza

# =============================================================================
# making pseudo GTF and GFF files, many forms were made cause wasn't sure which 
# was gonna work
# =============================================================================
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/RNA-Seq/")

# =============================================================================
# making gtf file: had to use ' instead of " because using " had it printing " 
# twice ("") when file is opened using a notepad
# =============================================================================
pantranscriptome_df = pd.read_csv("result/pantranscriptome.csv")
pantranscriptome_df.head()

# creating the pangene column
pantranscriptome_df['origin'] = 'pseudo'
pantranscriptome_df['type'] = 'gene'
pantranscriptome_df['start_position'] = 1
pantranscriptome_df['score'] = '.'
pantranscriptome_df['strand'] = '+'
pantranscriptome_df['frame'] = '.'

pantranscriptome_df['attributes'] = "transcript_id " + pantranscriptome_df['sequence_name'] + "; "\
                                    "gene_id " + pantranscriptome_df['sequence_name'] + "; "\
                                    "gene_name " + pantranscriptome_df['sequence_name'] + ";"

columns_to_drop = ['stable_id', 'sequence_name']
pantranscriptome_df2 = pantranscriptome_df.drop(columns=columns_to_drop)

new_order = ['pangene', 'origin', 'type', 'start_position', 'sequence_length', 'score', 'strand', 'frame', 'attributes']
pantranscriptome_df2 = pantranscriptome_df2[new_order]

pantranscriptome_df2.to_csv("result/pseudo-annotation.gtf", sep='\t', index=False, header=False)

# =============================================================================
# making gtf file: 
# =============================================================================
pantranscriptome_df2['attributes'] = "transcript_id " + "'" + pantranscriptome_df2['pangene'] + "'; "\
                                    "gene_id " + "'" + pantranscriptome_df2['pangene'] + "'; "\
                                    "gene_name " + "'" + pantranscriptome_df2['pangene'] + "';"

pantranscriptome_df2.to_csv("result/pseudo-annotation2.gtf", sep='\t', index=False, header=False)

# =============================================================================
# making gff file: type=gene
# =============================================================================

pantranscriptome_df3 = pantranscriptome_df.drop(columns='stable_id')

pantranscriptome_df3['attributes'] = "ID=" + pantranscriptome_df3['sequence_name'] + ";"\
                                    "Name=" + pantranscriptome_df3['sequence_name']

pantranscriptome_df3 = pantranscriptome_df3.drop(columns='sequence_name')

new_order = ['pangene', 'origin', 'type', 'start_position', 'sequence_length', 'score', 'strand', 'frame', 'attributes']
pantranscriptome_df3 = pantranscriptome_df3[new_order]

pantranscriptome_df3.to_csv("result/pseudo-annotation.gff", sep='\t', index=False, header=False)

# =============================================================================
# making gff file: type=CDS
# =============================================================================
pantranscriptome_df4 = pantranscriptome_df.drop(columns='stable_id')
pantranscriptome_df4['type'] = 'CDS'

pantranscriptome_df4['attributes'] = "ID=" + pantranscriptome_df4['sequence_name'] + ".CDS.1;"\
                                    "Parent=" + pantranscriptome_df4['sequence_name']

pantranscriptome_df4 = pantranscriptome_df4.drop(columns='sequence_name')

new_order = ['pangene', 'origin', 'type', 'start_position', 'sequence_length', 'score', 'strand', 'frame', 'attributes']
pantranscriptome_df4 = pantranscriptome_df4[new_order]

pantranscriptome_df4.to_csv("result/pseudo-annotation2.gff", sep='\t', index=False, header=False)

# =============================================================================
# making gff file: type = gene+CDS, sorted properly
# =============================================================================

# making the idnex even numbers
pantranscriptome_df3.index = range(0, len(pantranscriptome_df3) * 2, 2)

# making the index odd numbers
pantranscriptome_df4.index = range(1, len(pantranscriptome_df4) * 2, 2)

# concatenating the two dataframes
pantranscriptome_df5 = pd.concat([pantranscriptome_df3, pantranscriptome_df4], axis=0)

# sorting the dataframe by index
pantranscriptome_df5 = pantranscriptome_df5.sort_index()

pantranscriptome_df5.to_csv("result/pseudo-annotation3.gff", sep='\t', index=False, header=False)