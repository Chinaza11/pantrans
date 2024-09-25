# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 10:31:07 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/")

# =============================================================================
# annotate duplicate gene copies with their species tree node/age of duplication
# =============================================================================

annot = pd.read_table("GO_analysis/Rio_annotationInfo-reduced_deduplicate.txt")
annot.rename(columns={'locusName': 'sequence_name'}, inplace=True)

duplications = pd.read_csv("orthofinder_3/Duplications.tsv", sep='\t')

def duplication_info(gene_column_duplication_file, STN_column_annot_file, orthogrp_column_annot_file):
    for index, gene in annot['sequence_name'].items():
        matching_rows = duplications[duplications[gene_column_duplication_file].str.contains(gene, case=False, na=False, regex=False)]
        matching_indices = matching_rows.index
        matching_STN = duplications.loc[matching_indices, 'Species Tree Node'].tolist()
        annot.at[index, STN_column_annot_file] = ','.join(matching_STN)
        
        matching_O = duplications.loc[matching_indices, 'Orthogroup'].tolist()
        annot.at[index, orthogrp_column_annot_file] = ','.join(matching_O)


duplication_info('Genes 1', 'Genes1 Species Tree Node', 'Genes1 Orthogroup')

duplication_info('Genes 2', 'Genes2 Species Tree Node', 'Genes2 Orthogroup')
    
annot = annot[['sequence_name','Genes1 Species Tree Node','Genes2 Species Tree Node','Genes1 Orthogroup','Genes2 Orthogroup']]

annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt1.csv', index=False)

# =============================================================================
# annotate all genes with their orthogroups gotten from Orthogroups.tsv file
# =============================================================================

orthogroups = pd.read_csv("orthofinder_3/Orthogroups.tsv", sep='\t')

for index, gene in annot['sequence_name'].items():
    matching_rows = orthogroups[orthogroups['sorghumRio'].str.contains(gene, case=False, na=False, regex=False)]
    matching_O = matching_rows['Orthogroup'].tolist()
    annot.at[index, 'orthogroup_from_orthogroups.tsv'] = ','.join(matching_O)

annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt2.csv', index=False)

# =============================================================================
# annotate genes that are SingleCopyOrthologues
# =============================================================================

header = ["group_name"]
Orthogroups_SingleCopyOrthologues = pd.read_table("orthofinder_3/Orthogroups_SingleCopyOrthologues.txt", names=header)

for index, orthogroup in annot['orthogroup_from_orthogroups.tsv'].items():
    if pd.isna(orthogroup) != True:
        for index2, orthogroup2 in Orthogroups_SingleCopyOrthologues['group_name'].items():
            if orthogroup == orthogroup2:
                annot.at[index, "SingleCopyOrthologues"] = "SingleCopyOrthologues"

annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt3.csv', index=False)

# =============================================================================
# annotate genes that are unassigned to any orthogroup)
# =============================================================================

Orthogroups_UnassignedGenes = pd.read_csv("orthofinder_3/Orthogroups_UnassignedGenes.tsv", sep='\t', low_memory=False)

for index, gene in annot['sequence_name'].items():
    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumRio'].items():
        if pd.isna(gene2) != True:
            if gene in gene2:
                annot.at[index, "UnassignedGenes"] = "UnassignedGenes"
            
annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt4.csv', index=False)

# =============================================================================
# descriptive statistics
# =============================================================================

annot = pd.read_csv("orthofinder_3/intermediate_files/annot_node_ckpt4.csv")

pd.notna(annot['Genes1 Species Tree Node']).sum()
pd.notna(annot['Genes2 Species Tree Node']).sum()
pd.notna(annot['SingleCopyOrthologues']).sum()
pd.notna(annot['orthogroup_from_orthogroups.tsv']).sum()
pd.notna(annot['UnassignedGenes']).sum()

# 15950 rows were annotated from Genes1 column of the Duplications.tsv file
# 17145 rows were annotated from Genes2 column of the Duplications.tsv file
# 20 genes were SingleCopyOrthologues
# 1295 genes were UnassignedGenes
# 34195 rows from orthogroup_from_orthogroups.tsv column + 1295 rows from UnassignedGenes column = 35490

# =============================================================================
# join Genes1 Species Tree Node, Genes2 Species Tree Node, SingleCopyOrthologues
# UnassignedGenes and remove duplicate entries
# =============================================================================

for index, cell in annot['Genes1 Species Tree Node'].items():       
    l = str(cell) + ',' + str(annot.at[index,'Genes2 Species Tree Node']) + ',' + str(annot.at[index,'SingleCopyOrthologues']) + ',' + str(annot.at[index,'UnassignedGenes'])
    cleaned_list = [item for item in l.split(',') if item != 'nan']
    unique_list = list(set(cleaned_list))
    annot.at[index, 'Species.Tree.Node'] = ','.join(unique_list)
    
# =============================================================================
# convert empty cells to the string value "NA" for GO analysis in R
# =============================================================================

for index, cell in annot['Species.Tree.Node'].items():
    if len(annot.at[index, 'Species.Tree.Node']) == 0:
        annot.at[index, 'Species.Tree.Node'] = 'NA'
        
(annot['Species.Tree.Node'] == 'NA').sum()

# 23301 genes annotated (leaving 12189 not yet annotated)

annot.rename(columns={'sequence_name': 'locusName'}, inplace=True)
    
annot.to_csv('orthofinder_3/annot_node_ckpt5.csv', index=False)
