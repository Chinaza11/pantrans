# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 06:28:24 2024

@author: Chinaza
"""

# =============================================================================
# making pantranscriptome that works with STAR
# =============================================================================

import os

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

#-----------------------------------------
pantranscriptome = read_fasta("result/pantranscriptome_2.fasta")

for key, value in pantranscriptome.items():
    new_header = key.split(" ")[0]
    to_be_saved = f"{new_header}\n{value}\n"
    #print(to_be_saved)
                
    with open('result/pantranscriptome_for_star.fasta', 'a') as fh:
        fh.write(to_be_saved)

