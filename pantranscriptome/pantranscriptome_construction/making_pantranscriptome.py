# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 12:08:11 2023

@author: Chinaza Nnamdi
"""

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
            
            # getting each pop gene name(s) and splitting cells with more than one gene name 
            gene_names = pav_df_one_cds[pop][index]
            b = gene_names.split(',')

            # dealing with cells with only one gene name: this was added directly to the table
            if len(b) == 1:
                new_df_one_cds.at[index, 'stable_id'] = stable_id
                new_df_one_cds.at[index, 'sequence_name'] = pav_df_one_cds[pop][index]
                new_df_one_cds.at[index, 'sequence_length'] = pav_df_one_cds[pop_length][index]
                
            # dealing with cells that has more than one gene name: this was split up and added to the dictionary
            elif len(b) > 1:
                temp_dict = {} 
                
                #removing extra spaces between gene names
                gene_names_list = []
                for i in b:
                    z = i.strip()
                    gene_names_list.append(z)
                    
                length_multi = pav_df_one_cds[pop_length][index] #getting length
                length_multi_split = length_multi.split(',') #splitting up length 
                
                ##removing extra spaces between length
                length_list = []
                for j in length_multi_split:
                    k = j.strip()
                    length_list.append(k)
                                
                #adding them to dict
                for i in range(len(b)): 
                    key2 = gene_names_list[i]
                    length2 = length_list[i]
                    temp_dict[key2] = length2
                
                max_value = float('-inf')  # Initialize with negative infinity
                max_key = None
                
                for key, value in temp_dict.items():
                    value = int(value)
                    if value > max_value:
                        max_value = value
                        max_key = key
                
                new_df_one_cds.at[index, 'stable_id'] = stable_id
                new_df_one_cds.at[index, 'sequence_name'] = max_key
                new_df_one_cds.at[index, 'sequence_length'] = max_value
                        
#print(new_df_one_cds)
print(new_df_one_cds.head()) #checkpoint

# =============================================================================
# presence-absence table (rows with more than one cds): getting the cds with the 
# maximum length in each row and putting it in a new dataframe that comprise of 
# only 3 headers (stable_id, sequence_name, sequence_length). 
# an additional block of code at the end makes Rio the default sequence to pick 
# if all the population sequences have the same length 
# =============================================================================
pav_df_more_cds = pd.read_csv("result/presence_absence_table_more_than_one_cds.csv")
#print(pav_df_more_cds.columns) #get column names #checkpoint

#-----------------------------------------
#creating an empty dataframe for storing result

num_rows = pav_df_more_cds.shape[0]

new_df_more_cds = pd.DataFrame(index=range(num_rows), columns=['stable_id', 'sequence_name', 'sequence_length'])
#-----------------------------------------
for index, stable_id in pav_df_more_cds['stable ids'].items():
    main_dict = {} # dict for storing sequences name and length for each row
    
    # looping throw each column
    for pop in pop_column:
        pop_length = pop + "_length"
        
        if pd.isna(pav_df_more_cds[pop][index]) != True: # ignoring empty cells
            #print(pav_df_more_cds[pop][index])
            #print(pav_df_more_cds[pop_length][index])
            
            # getting each pop gene name(s) and splitting cells with more than one gene name 
            gene_names = pav_df_more_cds[pop][index]
            b = gene_names.split(',')
            
            # dealing with cells with only one gene name: this was added directly to the dictionary
            if len(b) == 1:
                key = pav_df_more_cds[pop][index]
                length = pav_df_more_cds[pop_length][index]
                main_dict[key] = length
            
            # dealing with cells that has more than one gene name: this was split up and added to the dictionary
            elif len(b) > 1:
                gene_names_list = []
                for i in b:
                    z = i.strip()
                    gene_names_list.append(z)
                    
                length_multi = pav_df_more_cds[pop_length][index] #getting length
                length_multi_split = length_multi.split(',') #splitting up length 
                length_list = []
                for j in length_multi_split:
                    k = j.strip()
                    length_list.append(k)
                
                #adding them to dict
                for i in range(len(b)): 
                    key2 = gene_names_list[i]
                    length2 = length_list[i]
                    main_dict[key2] = length2


    max_value = float('-inf')  # Initialize with negative infinity
    max_key = None                                

    # getting the sequence with max length for Rio pop
    for key, value in main_dict.items():
        if key.startswith("SbiPI563295"):
            value = int(value)
            if value > max_value:
                max_value = value
                max_key = key
    
    # if another pop has the maximum length, this will replace the sequence from Rio selected as max length
    for key, value in main_dict.items():
        value = int(value)
        if value > max_value:
            max_value = value
            max_key = key
            
    new_df_more_cds.at[index, 'stable_id'] = stable_id
    new_df_more_cds.at[index, 'sequence_name'] = max_key
    new_df_more_cds.at[index, 'sequence_length'] = max_value
    
print(new_df_more_cds)
#print(new_df_more_cds.head()) #checkpoint


# =============================================================================
# joining new_df_one_cds and new_df_more_cds dataframes together, 
# sorting dataframe by 'stable_id' column and resetting the index
# creating the pangene column
# =============================================================================

joined_df = pd.concat([new_df_one_cds, new_df_more_cds], ignore_index=True)
print(joined_df)

# sorting dataframe and reseting the index
joined_df = joined_df.sort_values(by='stable_id')
joined_df = joined_df.reset_index(drop=True)

# creating the pangene column
joined_df['StringColumn'] = 'pangene'
joined_df['CountColumn'] = range(1, len(joined_df) + 1)
joined_df['pangene'] = joined_df['StringColumn'] + joined_df['CountColumn'].astype(str)
print(joined_df)

# dropping irrelevant columns
columns_to_drop = ['StringColumn', 'CountColumn']
joined_df.drop(columns=columns_to_drop, inplace=True)
print(joined_df)

# write table out
joined_df.to_csv('result/pantranscriptome.csv', index=False)

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
for index, sequence in joined_df['sequence_name'].items():

    # iterating through each population
    for key2,value2 in all_genome_cds.items():
        
        # iterating through each gene in each population
        for key3,value3, in all_genome_cds[key2].items():
            key3_split = key3.split('>')
            
            # checking if the gene in the population correspond to the sequence in the row of joined_df and adding it to consensus_sequence.fasta file if it does
            if joined_df['sequence_name'][index] == key3_split[1]:
                header = '>' + joined_df['pangene'][index] + ' Pop_specific_name:' + key3_split[1] + ' Stable_ID:' + joined_df['stable_id'][index] + ' Sequence_length:' + str(joined_df['sequence_length'][index])
                to_be_saved = f"\n{header}\n{value3}"
                print(to_be_saved)
                
                with open('result/pantranscriptome.fasta', 'a') as fh:
                    fh.write(to_be_saved)