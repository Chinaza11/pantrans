# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 09:40:01 2024

@author: Chinaza
"""

# =============================================================================
# this task was done in the HPC environment
# =============================================================================

import os, pandas as pd

# os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# loading all CP-NAM parent interproscan result file and concantenating them
# =============================================================================

header = ["protein_accession", "md5", "seq_length", "analysis", "signature_accession", "signature_description", "start_location", "stop_location", "score", "status", "date", "interpro_annot_acce", "interpro_annot_descrip", "GO_annot"]


ChiAmber = pd.read_table("/projects/cooper_research/Ref_Genomes/ChiAmber/annotation/sorghum_chineseamber.evd.protein.interproscan.tsv", names=header)
Grassl = pd.read_table("/projects/cooper_research/Ref_Genomes/Grassl/annotation/sorghum_grassl.evd.protein.interproscan.tsv", names=header)
Leoti = pd.read_table("/projects/cooper_research/Ref_Genomes/Leoti/annotation/sorghum_leoti.evd.protein_interproscan.tsv", names=header)
pi229841 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi229841/annotation/sorghum_pi229841.evd.protein.interproscan.tsv", names=header)
pi297155 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi297155/annotation/sorghum_pi297155.evd.protein.interproscan.tsv", names=header)
pi329311 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi329311/annotation/sorghum_pi329311.evd.protein.interproscan.tsv", names=header)
pi506069 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi506069/annotation/sorghum_pi506069.evd.protein.interproscan.tsv", names=header)
pi510757 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi510757/annotation/sorghum_pi510757.evd.protein.interproscan.tsv", names=header)
pi655972 = pd.read_table("/projects/cooper_research/Ref_Genomes/pi655972/annotation/sorghum_pi655972.evd.protein.interproscan.tsv", names=header)
RioNAM = pd.read_table("/projects/cooper_research/Ref_Genomes/RioNAM/annotation/sorghum_riouncc.evd.protein_interproscan.tsv", names=header)


pop_names = ['ChiAmber', 'Grassl', 'Leoti', 'pi229841', 'pi297155', 'pi329311', 'pi506069', 'pi510757', 'pi655972', 'RioNAM']
pop_dataframes = [ChiAmber, Grassl, Leoti, pi229841, pi297155, pi329311, pi506069, pi510757, pi655972, RioNAM]

def inspection(name, dataframe):
    for i in range(len(name)):
        print(f"Dataframe for pop: {name[i]}")
        print(dataframe[i].shape)
        print(dataframe[i].info())
        print(dataframe[i].head())

# inspecting all dataframes
inspection(pop_names, pop_dataframes)

# concatenating all dataframes
total_interproscan = pd.concat([ChiAmber, Grassl, Leoti, pi229841, pi297155, pi329311, pi506069, pi510757, pi655972, RioNAM], ignore_index=True)
total_interproscan.shape
total_interproscan.info()
total_interproscan.head()

# =============================================================================
# test above code on my local PC

"""
Leoti = pd.read_table("C:/Users/nnamd/Downloads/sorghum_leoti.evd.protein_interproscan.tsv", names=header)
RioNAM = pd.read_table("C:/Users/nnamd/Downloads/sorghum_riouncc.evd.protein_interproscan.tsv", names=header)

pop_names = ["Leoti", "RioNAM"]
pop_dataframes = [Leoti, RioNAM]

inspection(pop_names, pop_dataframes)

total_interproscan = pd.concat([Leoti, RioNAM], ignore_index=True)
total_interproscan.shape
total_interproscan.info()
total_interproscan.head()
"""

# =============================================================================
# get GO terms for each of the gene/protein
# =============================================================================

pantranscriptome = pd.read_csv("data/pantranscriptome.csv")

pantranscriptome.head()

pantranscriptome['GO_annotations'] = ""

for index, sequence in pantranscriptome['sequence_name'].items():
    go_terms_list = []
    go_terms_string = ""
    
    for index2, protein in total_interproscan['protein_accession'].items():
        if sequence==protein:
            #print(total_interproscan['GO_annot'][index2])
            
            if pd.isna(total_interproscan['GO_annot'][index2]) != True and total_interproscan['GO_annot'][index2] != '-':
                #print('yes, valid GO_term here')
                go_terms_list.append(total_interproscan['GO_annot'][index2])
    
    #print(go_terms_list)
    go_terms_string = ', '.join(go_terms_list)        
    #print(go_terms_string)
    
    pantranscriptome.at[index, 'GO_annotations'] = go_terms_string

pantranscriptome.to_csv('result/intermediate_result/pantranscriptome_with_go_terms.csv', index=False)

