import os, pandas as pd

print(os.getcwd())

#setting header before importing files
header = ["index", "stable_id", "pop_specific_id", "chromosome", "start", "stop"]

#importing stable id files of all cp-nams populations
ChiAmber = pd.read_table("stable_ids/ChineseAmber.assign_stable_IDs", names=header)
Grassl = pd.read_table("stable_ids/Grassel.assign_stable_IDs", names=header)
Leoti = pd.read_table("stable_ids/Leoti.assign_stable_IDs", names=header)
pi229841 = pd.read_table("stable_ids/pi229841.assign_stable_IDs", names=header)
pi297155 = pd.read_table("stable_ids/pi297155.assign_stable_IDs", names=header)
pi329311 = pd.read_table("stable_ids/pi329311.assign_stable_IDs", names=header)
pi506069 = pd.read_table("stable_ids/pi506069.assign_stable_IDs", names=header)
pi510757 = pd.read_table("stable_ids/pi510757.assign_stable_IDs", names=header)
pi655972 = pd.read_table("stable_ids/pi655972.assign_stable_IDs", names=header)
RioNAM = pd.read_table("stable_ids/Rio.assign_stable_IDs", names=header)

#print(RioNAM)
print(RioNAM.head(3)) #checkpoint

all_stable_ids = []

#converting all stable_id columns in each of the dataframes to a list
ChiAmber_stable_ids = ChiAmber['stable_id'].tolist()
Grassl_stable_ids = Grassl['stable_id'].tolist()
Leoti_stable_ids = Leoti['stable_id'].tolist()
pi229841_stable_ids = pi229841['stable_id'].tolist()
pi297155_stable_ids = pi297155['stable_id'].tolist()
pi329311_stable_ids = pi329311['stable_id'].tolist()
pi506069_stable_ids = pi506069['stable_id'].tolist()
pi510757_stable_ids = pi510757['stable_id'].tolist()
pi655972_stable_ids = pi655972['stable_id'].tolist()
RioNAM_stable_ids = RioNAM['stable_id'].tolist()

#storing all the stable ids in one list
all_stable_ids = ChiAmber_stable_ids + Grassl_stable_ids + Leoti_stable_ids + pi229841_stable_ids + pi297155_stable_ids + pi329311_stable_ids + pi506069_stable_ids + pi510757_stable_ids + pi655972_stable_ids + RioNAM_stable_ids

print(len(all_stable_ids)) #checkpoint

#keeping only the unique stable ids name
unique_stable_ids = []
for i in all_stable_ids:
    if i not in unique_stable_ids:
        unique_stable_ids.append(i)

print(len(unique_stable_ids)) #checkpoint

unique_stable_ids.sort() #ordering the unique gene

#putting the unique genes in a dataframe
pav_df = pd.DataFrame(unique_stable_ids, columns=["stable ids"])

#list of column names for pav_df
column_names = ["ChiAmber_matching_id", "ChiAmber_cds", "ChiAmber_cds_length", 
              "Grassl_matching_id", "Grassl_cds", "Grassl_cds_length", 
              "Leoti_matching_id", "Leoti_cds", "Leoti_cds_length", 
              "pi229841_matching_id", "pi229841_cds", "pi229841_cds_length", 
              "pi297155_matching_id", "pi297155_cds", "pi297155_cds_length", 
              "pi329311_matching_id", "pi329311_cds", "pi329311_cds_length", 
              "pi506069_matching_id", "pi506069_cds", "pi506069_cds_length", 
              "pi510757_matching_id", "pi510757_cds", "pi510757_cds_length", 
              "pi655972_matching_id", "pi655972_cds", "pi655972_cds_length", 
              "RioNAM_matching_id", "RioNAM_cds", "RioNAM_cds_length"]

pav_df[column_names] = ""

#-----------------------------------------
#matching stable_ids with population specific id and extracting this into pav_df
def extracting_pop_id(pop_df, pav_df_column):
    for index, gene in pav_df['stable ids'].iteritems():
        l = []
        s = ""
        #print(f"{index}, {gene}")
        x = pop_df[pop_df['stable_id'] == gene]
        
        #for-loop helps account for when stable_id match to more than one pop specific id         
        for j in x['pop_specific_id']:
            l.append(j)
            #print(l)
            s = ', '.join(l)
            #print(s)
        if len(l) > 0:
            pav_df.at[index, pav_df_column] = s

pop_df_list = [ChiAmber, Grassl, Leoti, pi229841, pi297155, pi329311, pi506069, pi510757, pi655972, RioNAM]
pav_df_column_names = ['ChiAmber_matching_id', 'Grassl_matching_id', 'Leoti_matching_id', 'pi229841_matching_id', 'pi297155_matching_id', 'pi329311_matching_id', 'pi506069_matching_id', 'pi510757_matching_id', 'pi655972_matching_id', 'RioNAM_matching_id']

for i in range(len(pav_df_column_names)):
    extracting_pop_id(pop_df_list[i], pav_df_column_names[i])
    
print(pav_df.head(3)) #checkpoint
    
#-----------------------------------------
#to read fasta files
def read_fasta(file_path):
    d = {}
    header = ""
    sequence = ""
    with open(file_path) as file:
        for line in file:
            #print(line)
            line = line.strip()
            if line.startswith(">"):
                if len(header) > 0:
                    d[header] = sequence
                    sequence = ""
                header = line
                #print(header)
            else:
                sequence += line
    d[header] = sequence
    #print(d)
    return(d)

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

all_genome_cds = {} #dictionary for storing all genome

#getting all genomes, storing it in all_genome_cds and naming it based on its file path
for file_path in file_paths:
    cds_data = read_fasta(file_path)
    
    split_file_path = file_path.split('/') #splitting file_path
    genome_name = split_file_path[4] #naming genome based on its file_path
    all_genome_cds[genome_name] = cds_data

print(all_genome_cds.keys()) #checkpoint
#all_genome_cds.values()

#-----------------------------------------
#function for extracting cds name and length from fasta file that has been turned to a dictionary
def extracting_cds_name_and_length(pop_name):   
    for index, row in pav_df.iterrows():
        matching_id_column = pop_name + '_matching_id'
        cds_column = pop_name + '_cds'
        cds_length = pop_name + '_cds_length'
    
        cds_name_list = []
        cds_name_str = ""
        seq_len_list = []
        seq_len_str = ""
        
        #accessing genes in pop matching id column
        cds = pav_df[matching_id_column][index]
    
        #if statement and for-loop that follows splits rows that have multiple genes so that each is handled one at a time 
        if len(cds) > 0:
            cds = cds.split(',')
    
            for i in cds:
                i = i.strip()
                
                #for-loop loops through dict of populations and get matching cds name and length
                for cds_name, sequence in all_genome_cds[pop_name].items(): 
                    if cds_name.find(i) == True:
                        cds_name = cds_name.strip(">")
                        cds_name_list.append(cds_name)
                        cds_name_str = ', '.join(cds_name_list)
                        pav_df.at[index,cds_column] = cds_name_str
                        
                        seq_len = str(len(sequence))
                        seq_len_list.append(seq_len)
                        seq_len_str = ', '.join(seq_len_list)
                        pav_df.at[index,cds_length] = seq_len_str

pop_names = ['ChiAmber', 'Grassl', 'Leoti', 'pi229841', 'pi297155', 'pi329311', 'pi506069', 'pi510757', 'pi655972', 'RioNAM']

for i in range(len(pop_names)):
    extracting_cds_name_and_length(pop_names[i])

pav_df.to_csv('result/presence_absence_table.csv', index=False)


