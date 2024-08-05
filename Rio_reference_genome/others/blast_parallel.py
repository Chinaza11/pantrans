# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:27:59 2024

@author: Chinaza
"""

import subprocess
from concurrent.futures import ThreadPoolExecutor

def process_id(protein_id):
    try:
        print(f"Processing: {protein_id}")

        # Get protein sequence in FASTA format
        efetch_cmd = f"efetch -db protein -id {protein_id} -format fasta > input_{protein_id}.fasta"
        subprocess.run(efetch_cmd, shell=True, check=True)

        # blastp search of protein sequence
        blast_cmd = f"blastp -query input_{protein_id}.fasta -db nr -out blast_result/{protein_id}.txt -max_target_seqs 10 -outfmt '6 qseqid sseqid pident evalue stitle length' -remote"
        subprocess.run(blast_cmd, shell=True, check=True)
        
        print(f"Completed: {protein_id}")

    except subprocess.CalledProcessError as e:
        print(f"Error processing {protein_id}: {e}")

def main():
    # Read the list of IDs
    with open("id_list.txt") as f:
        protein_ids = [line.strip() for line in f if line.strip()]

    # Use ThreadPoolExecutor to parallelize the task
    with ThreadPoolExecutor(max_workers=32) as executor:
        executor.map(process_id, protein_ids)

if __name__ == "__main__":
    main()


"""
Didn't use it because 'efetch' generates error when used concurrently several times

"""