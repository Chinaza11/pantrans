# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 16:20:58 2024

@author: Chinaza
"""

import os

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/others/labeling_genes/")

# =============================================================================
# splitting the gene list file into several small files
# =============================================================================

def split_file(input_file, lines_per_file, output_path):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Iterate over lines in chunks
    for i in range(0, len(lines), lines_per_file):
        #print(i)
        chunk = lines[i: i+lines_per_file]
        output_file = f'{output_path}_{i // lines_per_file + 1}.txt'
        with open(output_file, 'w') as file:
            file.writelines(chunk)
        print(f'Created {output_file}')

# Usage
input_file = 'treatment/treatment_protein_id_for_hpc.csv'
lines_per_file = 11
output_path = 'treatment/treatment_gene_list'
split_file(input_file, lines_per_file, output_path)

# =============================================================================
# code for running several scripts in parallel in the HPC environment
# =============================================================================

'''
for i in {1..16}; do cp multi_blast.slurm multi_blast_$i.slurm; done
    
for i in {1..16}; do sed -i "s/python blast/python multi_blast_$i/g" multi_blast_$i.slurm; done
    
module load parallel

parallel sbatch ::: *.slurm
'''