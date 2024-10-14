# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 04:34:17 2024

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# annotate duplicate gene copies with their species tree node/age of duplication
# =============================================================================

annot = pd.read_csv("GO_analysis/pantranscriptome_plus_GO_cleaned.csv")

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
    
annot = annot[['stable_id','sequence_name','pangene','Genes1 Species Tree Node','Genes2 Species Tree Node','Genes1 Orthogroup','Genes2 Orthogroup']]

# annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt1.csv', index=False)

# =============================================================================
# annotate all genes with their orthogroups gotten from Orthogroups.tsv file
# =============================================================================

orthogroups = pd.read_csv("orthofinder_3/Orthogroups.tsv", sep='\t')

def add_orthogroup_info(gene, index, column_name):
    matching_rows = orthogroups[orthogroups[column_name].str.contains(gene, case=False, na=False, regex=False)]
    matching_O = matching_rows['Orthogroup'].tolist()
    annot.at[index, 'orthogroup_from_orthogroups.tsv'] = ','.join(matching_O)

prefix_to_column = {
    'SbiCamber': 'sorghumChineseamber',
    'SbiGrassl': 'sorghumGrassl',
    'SbiLeoti': 'sorghumLeoti',
    'SbiPI229841': 'sorghumPI229841',
    'SbiPI297155': 'sorghumPI297155',
    'SbiPI329311': 'sorghumPI329311',
    'SbiPI506069': 'sorghumPI506069',
    'SbiPI510757': 'sorghumPI510757',
    'SbiPI655972': 'sorghumPI655972',
    'SbiPI563295': 'sorghumRioNAM'
}

for index, gene in annot['sequence_name'].items():
    for prefix, column_name in prefix_to_column.items():
        if gene.startswith(prefix):
            add_orthogroup_info(gene, index, column_name)
            break

# annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt2.csv', index=False)

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

# annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt3.csv', index=False)

# =============================================================================
# annotate genes that are unassigned to any orthogroup
# =============================================================================

Orthogroups_UnassignedGenes = pd.read_csv("orthofinder_3/Orthogroups_UnassignedGenes.tsv", sep='\t', low_memory=False)

def unassigned_genes(gene, index, column_name):
    for index2, gene2 in Orthogroups_UnassignedGenes[column_name].items():
        if gene == gene2:
            annot.at[index, "UnassignedGenes"] = "UnassignedGenes"

for index, gene in annot['sequence_name'].items():
    for prefix, column_name in prefix_to_column.items():
        if gene.startswith(prefix):
            unassigned_genes(gene, index, column_name)
            break

# annot.to_csv('orthofinder_3/intermediate_files/annot_node_ckpt4.csv', index=False)

# =============================================================================
# descriptive statistics
# =============================================================================

annot = pd.read_csv("orthofinder_3/intermediate_files/annot_node_ckpt4.csv")

pd.notna(annot['Genes1 Species Tree Node']).sum()
pd.notna(annot['Genes2 Species Tree Node']).sum()
pd.notna(annot['SingleCopyOrthologues']).sum()
pd.notna(annot['orthogroup_from_orthogroups.tsv']).sum()
pd.notna(annot['UnassignedGenes']).sum()

# 25314 rows were annotated from Genes1 column of the Duplications.tsv file
# 27307 rows were annotated from Genes2 column of the Duplications.tsv file
# 4 genes were SingleCopyOrthologues
# 2265 genes were UnassignedGenes
# 59779 rows from orthogroup_from_orthogroups.tsv column + 2265 rows from UnassignedGenes column = 62044

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

# 35992 genes annotated (leaving 26052 not yet annotated)

# =============================================================================
# changing all nodes in CP-NAM to N17
# =============================================================================

replacement_dict = {
    'N19': 'N17',
    'N21': 'N17',
    'N23': 'N17',
    'N24': 'N17',
    'N26': 'N17',
    'N27': 'N17',
    'N28': 'N17',
    'N29': 'N17',
    'sorghumChineseamber': 'N17',
    'sorghumGrassl': 'N17',
    'sorghumLeoti': 'N17',
    'sorghumPI229841': 'N17',
    'sorghumPI297155': 'N17',
    'sorghumPI329311': 'N17',
    'sorghumPI510757': 'N17',
    'sorghumPI506069': 'N17',
    'sorghumPI655972': 'N17',
    'sorghumRioNAM': 'N17'}

for key, value in replacement_dict.items():
    annot['Species.Tree.Node'] = annot['Species.Tree.Node'].str.replace(key,value)
    
annot.rename(columns={'sequence_name': 'locusName'}, inplace=True)
    
annot.to_csv('orthofinder_3/annot_node_ckpt5.csv', index=False)

