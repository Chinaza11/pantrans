# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 12:10:28 2023

@author: Chinaza Nnamdi
"""

import os, pandas as pd, subprocess

os.chdir("G:/My Drive/PhD/project/RNA-Seq/")

#-----------------------------------------
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
              "G:/My Drive/PhD/project/RNA-Seq/data/testing_riouncc.fasta"]

all_genome_cds = {} #dictionary for storing all genome

#getting all genomes, storing it in all_genome_cds and naming it based on its file path
for file_path in file_paths:
    cds_data = read_fasta(file_path)
    
    split_file_path = file_path.split('/') #splitting file_path
    genome_name = split_file_path[6] #naming genome based on its file_path
    all_genome_cds[genome_name] = cds_data

#-----------------------------------------

new_pav_df = pd.read_csv("result/presence_absence_table_more_than_one_cds.csv")
new_pav_df.head()

pop_column = ['ChiAmber_cds', 'Grassl_cds', 'Leoti_cds', 'pi229841_cds', 'pi297155_cds', 'pi329311_cds', 'pi506069_cds', 'pi510757_cds', 'pi655972_cds', 'RioNAM_cds']

#-----------------------------------------
#extracting all coding sequences associated with each stable ID
def extract_name_and_sequence(current_index,pop_column):
    
    #####----for trouble_shooting-----#####
    #current_index = 0 #3 #243 #20628
    #pop_column = "ChiAmber_cds"

    storage_dict = {}
    if pd.isna(new_pav_df[pop_column][current_index]) != True: #skipping all empty cells
        x = new_pav_df[pop_column][current_index] #getting the content of the cell
        x = x.split(", ") #splitting up cell contents if it has more than one sequence name
        
        #adding ">" before sequence name in preparation for search
        for i in x:
            cds_name = ">" + i
            
            #looping through all populations sequences in search of cds_name
            for value in all_genome_cds.values():
                if cds_name in value:
                    sequence = value.get(cds_name)
                    storage_dict[cds_name] = sequence
    
    #print(storage_dict)
    return(storage_dict)


#looping through all stable IDs, applying extract_name_and_sequence function, putting result in fasta file (input for cdhit), calling cdhit using default settings, and saving all cdhit output in one file
for index, gene in new_pav_df['stable ids'].iteritems():
    data = {}
    for pop in pop_column:
        function_output = extract_name_and_sequence(index, pop)
        data.update(function_output)
    
    cdhit_input_file = f"cdhit_result/{gene}_input.fasta"
    with open(cdhit_input_file, "w") as temp_file:
        for key,value in data.items():
            temp_file.write(f"{key}\n{value}\n")
            
    #calling muscle for sequence alignment
    cdhit_output_file = f"cdhit_result/{gene}_output.fasta"
    command = ["cd-hit-est", "-i", cdhit_input_file, "-o", cdhit_output_file, "-T", "0"]
    subprocess.call(command)
    
    #adding cdhit output to pantranscriptome.fasta
    with open(cdhit_output_file, 'r') as cdhit_output:
        with open("final_result/pantranscriptome.fasta", 'a') as pantranscriptome:
            pantranscriptome.write(f"\n{cdhit_output.read()}")
    
    