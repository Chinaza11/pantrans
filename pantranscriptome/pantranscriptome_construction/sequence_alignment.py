import os, pandas as pd, subprocess

os.chdir("G:/My Drive/PhD/project/RNA-Seq/")

pav_df = pd.read_csv("result/from_the_cluster/presence_absence_table.csv")
#pav_df = pd.read_csv("C:/Users/nnamd/Downloads/testing.csv")

pav_df.columns
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

all_genome_cds = {} #dictionary for storing all genomes

#getting all genomes, storing it in all_genome_cds and naming it based on its file path
for file_path in file_paths:
    cds_data = read_fasta(file_path)
    
    split_file_path = file_path.split('/')
    genome_name = split_file_path[6] #naming genome based on its file_path
    all_genome_cds[genome_name] = cds_data

#-----------------------------------------
pop_column = ['ChiAmber_cds', 'Grassl_cds', 'Leoti_cds', 'pi229841_cds', 'pi297155_cds', 'pi329311_cds', 'pi506069_cds', 'pi510757_cds', 'pi655972_cds', 'RioNAM_cds']


# =============================================================================
# cleaning up presence-absence table: removing rows without cds in any of the population
# =============================================================================

del_indices = [] 

for index, stable_id in pav_df['stable ids'].iteritems():

    #getting indices of rows without any cds
    if all(pd.isna(pav_df.at[index,pop]) for pop in pop_column):
        del_indices.append(index)
        
print(del_indices) #checkpoint
pav_df.drop(del_indices, inplace=True) #dropping these rows from the dataframe
pav_df = pav_df.reset_index(drop=True) #resetting the index

pav_df.to_csv('result/presence_absence_table_without_empty_cds_rows.csv', index=False)

# =============================================================================
# cleaning up presence-absence table: separating rows with only one cds from the dataframe
# =============================================================================

one_cds_indices = []

for index, stable_id in pav_df['stable ids'].iteritems():
   
    #putting the length of all pop cds in each row of the dataframe in a list
    one_cds_lst = []
    for pop in pop_column:
        pop_cds_length = pop + '_length'
        if pd.isna(pav_df[pop_cds_length][index]) != True:
            count = pav_df[pop_cds_length][index]
            one_cds_lst.append(count)
    #print(one_cds_lst)
    
    #getting the indices of rows where there is only one cds
    if len(one_cds_lst) == 1:
        one_cds_indices.append(index)
print(one_cds_indices) #checkpoint

one_cds_df = pav_df.loc[one_cds_indices] #making a dataframe that comprise of rows with only one cds
one_cds_df.to_csv('result/presence_absence_table_only_one_cds.csv', index=False)


new_pav_df = pav_df.drop(one_cds_indices) #dropping rows with only one cds from the dataframe
new_pav_df = new_pav_df.reset_index(drop=True) #resetting the index
new_pav_df.to_csv('result/presence_absence_table_more_than_one_cds.csv', index=False)

# =============================================================================
# #extracting all coding sequences associated with each stable ID
# =============================================================================
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
            #print(cds_name)
            
            #loopig through all populations sequences in search of cds_name
            for value in all_genome_cds.values():
                #print(value)
                if cds_name in value:
                    #print('yes')
                    sequence = value.get(cds_name)
                    storage_dict[cds_name] = sequence
    
    #print(storage_dict)
    return(storage_dict)
                
# =============================================================================
# looping through all stable IDs, applying extract_name_and_sequence function, 
# putting result in fasta file (input for muscle), conducting alignment with 
# muscle and generating percent identity matrix identity with clustal omega
# =============================================================================
for index, gene in new_pav_df['stable ids'].iteritems():
    data = {}
    for pop in pop_column:
        function_output = extract_name_and_sequence(index, pop)
        data.update(function_output)
    
    muscle_input_file = f"muscle_result/{gene}_input.fasta"
    with open(muscle_input_file, "w") as temp_file:
        for key,value in data.items():
            temp_file.write(f"{key}\n{value}\n")
            
    #calling muscle for sequence alignment
    muscle_output_file = f"muscle_result/{gene}_output.fasta"
    command = ["muscle", "-in", muscle_input_file, "-out", muscle_output_file]
    subprocess.call(command)
    
    #calling clustal omega and generating percent identity matrix
    clustal_output_file = f"muscle_result/{gene}_percent_identity.txt"
    command2 = ["clustalo", "--infile", muscle_output_file, "--full", "--distmat-out", clustal_output_file, "--percent-id"]
    subprocess.call(command2)
