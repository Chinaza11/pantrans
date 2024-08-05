# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:52:00 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/others/labeling_genes/")

# =============================================================================
# identifying and getting genes that generated error and provided no result in 
# previous run on the HPC environment. They will be run again.
# =============================================================================

type_genes_blast = pd.read_table("type/type_genes_list.txt", names=['genes'])

l = []
for index, name in type_genes_blast['genes'].items():
    print(str(index) + ':  ' + name)
    with open('type/missing.txt', 'w') as fh:
        if os.path.isfile(f'type/blast_result/{name}.txt') != True:
            l.append(name)
            fh.write(f'{name}\n')
    
with open('type/missing.txt', 'w') as fh:    
    for i in l:
        fh.write(f'{i}\n')

# =============================================================================
# splitting the gene list file into several small files
# =============================================================================

def split_file(input_file, lines_per_file, output_file_name):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Iterate over lines in chunks
    for i in range(0, len(lines), lines_per_file):
        #print(i)
        chunk = lines[i: i+lines_per_file]
        output_file = f'{output_file_name}_{i // lines_per_file + 1}.txt'
        with open(output_file, 'w') as file:
            file.writelines(chunk)
        print(f'Created {output_file}')

# Usage
input_file = 'treatment/uncharacterized_protein.txt'
lines_per_file = 5
output_file_name = 'treatment/treatment_uncharacterized_protein'

split_file(input_file, lines_per_file, output_file_name)

# =============================================================================
# code for running several scripts in parallel in the HPC environment
# =============================================================================

'''
for i in {1..16}; do cp multi_blast.slurm multi_blast_$i.slurm; done
    
for i in {1..16}; do sed -i "s/python blast/python multi_blast_$i/g" multi_blast_$i.slurm; done
    
module load parallel

parallel sbatch ::: *.slurm
'''