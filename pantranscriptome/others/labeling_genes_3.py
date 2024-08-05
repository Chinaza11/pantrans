# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 07:35:55 2024

@author: Chinaza
"""

import os, subprocess

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# running BLAST
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

pantranscriptome = read_fasta("pantranscriptome_construction/pantranscriptome_2.fasta")


def blast(gene_list, input_file, output_folder):
    with open(gene_list) as fh:
        for line in fh:
            for key,val in pantranscriptome.items():
                if line.strip() in key:
                    with open(input_file, 'w') as fw:
                        print(f'Processing: {line}\n{key}\n-------')
                        fw.write(f"{key}\n{val}")
                    blast_cmd = f"blastn -query {input_file} -db /scratch/cnnamdi/blast/nt_db/nt -out blast_result/{output_folder}/{line.strip()}.txt -max_target_seqs 10 -num_threads 40 -outfmt '6 qseqid sseqid pident evalue stitle length'"
                    subprocess.run(blast_cmd, shell=True, check=True)


'''
blast('treatment_genes_list.txt', 'trt_blast_input.fasta', 'treatment')

blast('type_genes_list.txt', 'type_blast_input.fasta', 'type')

blast('interaction_genes_list.txt', 'int_blast_input.fasta', 'interaction')
'''