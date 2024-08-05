# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 03:50:51 2023

@author: Chinaza Nnamdi
"""

lst = ["pantranscriptome", "rows_with_only_one_cds"]

for i in lst:
    input_files = f"final_result/{i}.fasta"
    output_file = "final_result/final_pantranscriptome.fasta"
    #print(f"\n{i}")
    with open(input_files, 'r') as input_f:
        with open(output_file, 'a') as output_f:
            output_f.write(f"\n{input_f.read()}")
            
            
