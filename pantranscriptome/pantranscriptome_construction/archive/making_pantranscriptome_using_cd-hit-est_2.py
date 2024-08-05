# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 14:18:20 2023

@author: Chinaza Nnamdi
"""
#-----------------------------------------
import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/RNA-Seq/")

pop_column = ['ChiAmber_cds', 'Grassl_cds', 'Leoti_cds', 'pi229841_cds', 'pi297155_cds', 'pi329311_cds', 'pi506069_cds', 'pi510757_cds', 'pi655972_cds', 'RioNAM_cds']

# =============================================================================
# presence-absence table (rows with only one cds): creating a cleaner dataframe 
# that comprise of only 3 headers (stable_id, sequence_name, sequence_length)
# =============================================================================

pav_df_one_cds = pd.read_csv("result/presence_absence_table_only_one_cds.csv")
#print(pav_df_one_cds.columns) #get column names #checkpoint
#pav_df_one_cds.head()

#-----------------------------------------
#creating an empty dataframe for storing result

num_rows = pav_df_one_cds.shape[0]

new_df_one_cds = pd.DataFrame(index=range(num_rows), columns=['stable_id', 'sequence_name', 'sequence_length'])

#-----------------------------------------
#looping through each row, getting the cds and length and putting it in the new dataframe
for index, stable_id in pav_df_one_cds['stable ids'].items():
    for pop in pop_column:
        pop_length = pop + "_length"
        if pd.isna(pav_df_one_cds[pop][index]) != True:
            #print(pav_df_one_cds[pop][index])
            #print(pav_df_one_cds[pop_length][index])
            new_df_one_cds.at[index, 'stable_id'] = stable_id
            new_df_one_cds.at[index, 'sequence_name'] = pav_df_one_cds[pop][index]
            new_df_one_cds.at[index, 'sequence_length'] = pav_df_one_cds[pop_length][index]

print(new_df_one_cds.head()) #checkpoint

# =============================================================================
# reading all cds files and putting them in a dictionary
# =============================================================================

def read_fasta(file_path):
    d = {}
    header = ""
    sequence = ""
    with open(file_path) as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if len(header) > 0:
                    d[header] = sequence
                    sequence = ""
                header = line
            else:
                sequence += line
    d[header] = sequence
    return(d)

#-----------------------------------------

"""
file_paths = ["/projects/cooper_research/Ref_Genomes/ChiAmber/annotation/sorghum_chineseamber.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/Grassl/annotation/sorghum_grassl.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/Leoti/annotation/sorghum_leoti.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi229841/annotation/sorghum_pi229841.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi297155/annotation/sorghum_pi297155.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi329311/annotation/sorghum_pi329311.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi506069/annotation/sorghum_pi506069.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi510757/annotation/sorghum_pi510757.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/pi655972/annotation/sorghum_pi655972.evd.cds.fasta",
              "/projects/cooper_research/Ref_Genomes/RioNAM/annotation/sorghum_riouncc.evd.cds.fasta"]
"""

file_paths = ["G:/My Drive/PhD/project/RNA-Seq/data/testing_chineseamber.fasta",
              "G:/My Drive/PhD/project/RNA-Seq/data/testing_grassl.fasta",
              "G:/My Drive/PhD/project/RNA-Seq/data/testing_riouncc.fasta",
              "G:/My Drive/PhD/project/RNA-Seq/data/testing_leoti.fasta"]

all_genome_cds = {} #dictionary for storing all genome

#getting all genomes, storing it in all_genome_cds and naming it based on its file path
for file_path in file_paths:
    cds_data = read_fasta(file_path)
    
    split_file_path = file_path.split('/') #splitting file_path
    genome_name = split_file_path[6] #naming genome based on its file_path
    all_genome_cds[genome_name] = cds_data

# =============================================================================
# getting the sequence of genes that have been selected to have the maximum 
# length for each locus, this was then put in a new file 
# =============================================================================

# iterating through each row of the table
for index, sequence in new_df_one_cds['sequence_name'].items():

    # iterating through each population
    for key2,value2 in all_genome_cds.items():
        
        # iterating through each gene in each population
        for key3,value3, in all_genome_cds[key2].items():
            key3_split = key3.split('>')
            
            # checking if the gene in the population correspond to the sequence in the row of new_df_one_cds and a file
            if new_df_one_cds['sequence_name'][index] == key3_split[1]:
                to_be_saved = f"\n{key3}\n{value3}"
                print(to_be_saved)
                
                with open('final_result/rows_with_only_one_cds.fasta', 'a') as fh:
                    fh.write(to_be_saved)
                    