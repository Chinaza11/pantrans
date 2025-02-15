# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 07:43:55 2024

@author: Chinaza
"""

import os, pandas as pd, matplotlib.pyplot as plt
from matplotlib_venn import venn2

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# importing all dataframes
# =============================================================================

# importing dataframe with all orthogroups
Orthogroups = pd.read_csv("pantranscriptome/orthofinder_2/Orthogroups.tsv", sep='\t')
Orthogroups.head()
Orthogroups.info()

# -----------------------------------------------------------------------------
# importing dataframe with UnassignedGenes (genes unassigned to any orthogroup)
Orthogroups_UnassignedGenes = pd.read_csv("pantranscriptome/orthofinder_2/Orthogroups_UnassignedGenes.tsv", sep='\t', low_memory=False)
Orthogroups_UnassignedGenes.head()
Orthogroups_UnassignedGenes.info()

# -----------------------------------------------------------------------------
# DEGs from analysis that used the pantranscriptome
pan_type = pd.read_table("pantranscriptome/type_genes_df.txt")
pan_type.head()
pan_type['orthogroup'] = ""

pan_treat = pd.read_table("pantranscriptome/treat_genes_df.txt")
pan_treat.head()
pan_treat['orthogroup'] = ""

pan_int = pd.read_table("pantranscriptome/int_genes_df.txt")
pan_int.head()
pan_int['orthogroup'] = ""

# -----------------------------------------------------------------------------
# DEGs from analysis that used the Rio reference genome
rio_type = pd.read_table("Rio_reference_genome/type_genes_df.txt")
rio_type.head()
rio_type['orthogroup'] = ""

rio_treat = pd.read_table("Rio_reference_genome/treat_genes_df.txt")
rio_treat.head()
rio_treat['orthogroup'] = ""

rio_int = pd.read_table("Rio_reference_genome/int_genes_df.txt")
rio_int.head()
rio_int['orthogroup'] = ""

# =============================================================================
# annotating DEGs gotten using Rio reference genome
# =============================================================================

# -----------------------------------------------------------------------------
# annotating from orthogroups file (Orthogroups.tsv)
def annot_assigned_rio_genes(df, column):
    for index, gene in df[column].items():
        
        if gene.find("--") != -1:
            genes = gene.split("--")
            orthogroup_list = []
            for value in genes:
                for index2, gene2 in Orthogroups['sorghumRio'].items():
                    if value in str(gene2):
                        orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                        
            deduplicate_orthogroup_list = []
            for val in orthogroup_list:
                if val not in deduplicate_orthogroup_list:
                    deduplicate_orthogroup_list.append(val)
                    
            df.at[index,"orthogroup"] = ", ".join(deduplicate_orthogroup_list)
        
        else:
            for index2, gene2 in Orthogroups['sorghumRio'].items():
                if gene in str(gene2):
                    df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
    
annot_assigned_rio_genes(rio_type,"type_genes")
annot_assigned_rio_genes(rio_treat,"treat_genes")
annot_assigned_rio_genes(rio_int,"int_genes")

# -----------------------------------------------------------------------------
# annotating single genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_single_rio_genes(df, column):
    for index, group in df["orthogroup"].items():
        if len(group) == 0 :
            gene = df.at[index, column]
            
            for index2, gene2 in Orthogroups_UnassignedGenes['sorghumRio'].items():
                if pd.isna(gene2) != True:
                    if gene in str(gene2):
                        df.at[index, "orthogroup"] = Orthogroups_UnassignedGenes.at[index2, "Orthogroup"]
                                            
annot_unassigned_single_rio_genes(rio_type,"type_genes")
annot_unassigned_single_rio_genes(rio_treat,"treat_genes")
annot_unassigned_single_rio_genes(rio_int,"int_genes")  

# -----------------------------------------------------------------------------
# annotating hybrid genes unassigned to any orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_hybrid_rio_genes(df, column):
    for index, gene in df[column].items():        
        if gene.find("--") != -1:
            genes = gene.split("--")
            additional_orthogroup_list = []
            for value in genes:
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumRio'].items():
                    if value in str(gene2):
                        additional_orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                 
            if len(additional_orthogroup_list) > 0:
                deduplicate_orthogroup_list = []
                for val in additional_orthogroup_list:
                    if val not in deduplicate_orthogroup_list:
                        deduplicate_orthogroup_list.append(val)
                    
                deduplicate_orthogroup_str = ', ' + ", ".join(deduplicate_orthogroup_list)
                df.at[index,"orthogroup"] += deduplicate_orthogroup_str
            
annot_unassigned_hybrid_rio_genes(rio_type,"type_genes")
annot_unassigned_hybrid_rio_genes(rio_treat,"treat_genes")
annot_unassigned_hybrid_rio_genes(rio_int,"int_genes")


# =============================================================================
# annotating DEGs gotten using pantranscriptome
# =============================================================================

group_map = {
    "SbiCamber": "sorghumChineseamber",
    "SbiGrassl": "sorghumGrassl",
    "SbiLeoti": "sorghumLeoti",
    "SbiPI229841": "sorghumPI229841",
    "SbiPI297155": "sorghumPI297155",
    "SbiPI329311": "sorghumPI329311",
    "SbiPI506069": "sorghumPI506069",
    "SbiPI510757": "sorghumPI510757",
    "SbiPI655972": "sorghumPI655972",
    "SbiPI563295": "sorghumRioNAM"
}

# -----------------------------------------------------------------------------
# annotating from orthogroups file (Orthogroups.tsv)
def annot_assigned_pan_genes(df, column):
    
    # Helper function to find orthogroups for a gene based on its prefix
    def find_orthogroup(value, group_key):
        orthogroup_list = []
        for index2, gene2 in Orthogroups[group_key].items():
            if value in str(gene2):
                orthogroup_list.append(Orthogroups.at[index2, "Orthogroup"])
        return orthogroup_list

    for index, gene in df[column].items():
        orthogroup_list = []

        if "--" in gene:
            # Split gene string if it contains multiple genes separated by '--'
            genes = gene.split("--")
            for value in genes:
                for prefix, group_key in group_map.items():
                    if value.startswith(prefix):
                        orthogroup_list.extend(find_orthogroup(value, group_key))
        else:
            # Single gene case
            for prefix, group_key in group_map.items():
                if gene.startswith(prefix):
                    orthogroup_list = find_orthogroup(gene, group_key)
                    break

        # Deduplicate and assign orthogroups to the dataframe
        if orthogroup_list:
            df.at[index, "orthogroup"] = ", ".join(sorted(set(orthogroup_list)))

annot_assigned_pan_genes(pan_type,"type_genes")
annot_assigned_pan_genes(pan_treat,"treat_genes")
annot_assigned_pan_genes(pan_int,"int_genes")

# -----------------------------------------------------------------------------
# annotating single genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_single_pan_genes(df, column):
    
    # Helper function to find and assign orthogroup for unassigned genes
    def assign_orthogroup(gene, group_key, index):
        for index2, gene2 in Orthogroups_UnassignedGenes[group_key].items():
            if pd.notna(gene2) and gene in str(gene2):
                df.at[index, "orthogroup"] = Orthogroups_UnassignedGenes.at[index2, "Orthogroup"]
                break

    for index, group in df["orthogroup"].items():
        if not group:
            gene = df.at[index, column]
            for prefix, group_key in group_map.items():
                if gene.startswith(prefix):
                    assign_orthogroup(gene, group_key, index)
                    break
                                            
annot_unassigned_single_pan_genes(pan_type,"type_genes")
annot_unassigned_single_pan_genes(pan_treat,"treat_genes")
annot_unassigned_single_pan_genes(pan_int,"int_genes")  


# -----------------------------------------------------------------------------
# annotating hybrid genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_hybrid_pan_genes(df, column):
    # Helper function to find and append orthogroups for each gene
    def find_orthogroup(value, group_key):
        orthogroup_list = []
        for index2, gene2 in Orthogroups_UnassignedGenes[group_key].items():
            if value in str(gene2):
                orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2, "Orthogroup"])
        return orthogroup_list

    for index, gene in df[column].items():
        if "--" in gene:
            genes = gene.split("--")
            orthogroup_list = []

            # Loop through each split gene and assign orthogroups based on prefixes
            for value in genes:
                for prefix, group_key in group_map.items():
                    if value.startswith(prefix):
                        orthogroup_list.extend(find_orthogroup(value, group_key))
                        break  # Exit once a matching prefix is found

            # Deduplicate the orthogroup list and update the dataframe
            if orthogroup_list:
                deduplicated_orthogroups = list(set(orthogroup_list))  # Removes duplicates
                deduplicate_orthogroup_str = ', ' + ", ".join(deduplicated_orthogroups)
                df.at[index, "orthogroup"] += deduplicate_orthogroup_str      

annot_unassigned_hybrid_pan_genes(pan_type,"type_genes")
annot_unassigned_hybrid_pan_genes(pan_treat,"treat_genes")
annot_unassigned_hybrid_pan_genes(pan_int,"int_genes")  

# =============================================================================
# comparing DEGs orthogroups from rio-ref-genome and pantranscriptome and plotting it
# =============================================================================

def comparison(df1, df1_name, df2, df2_name, name):
    df1_lst = []
    for index1, group1 in df1['orthogroup'].items():
        split_group1 = group1.split(',')
        for i in split_group1:
            j = i.strip()
            if j not in df1_lst:
                df1_lst.append(j)
    print(f"Unique orthogroups in {df1_name} is: {len(df1_lst)}")
    
    df2_lst = []
    for index2, group2 in df2['orthogroup'].items():
        split_group2 = group2.split(',')
        for i in split_group2:
            j = i.strip()
            if j not in df2_lst:
                df2_lst.append(j)            
    print(f"Unique orthogroups in {df2_name} is: {len(df2_lst)}")
    
    #plot venn diagram
    venn2([set(df1_lst), set(df2_lst)], set_labels=(f"{df1_name}", f"{df2_name}"))
    plt.title(f'DEGs Orthogroups ({name})')
    plt.savefig(f'final_processing_and_plotting/DEGs_Orthogroups_({name})', dpi=300) 
    plt.show()
    
comparison(rio_type, 'Rio-Ref-Genome', pan_type, 'Pan-transcriptome', 'Type')
comparison(rio_treat, 'Rio-Ref-Genome', pan_treat, 'Pan-transcriptome', 'Treatment')
comparison(rio_int, 'Rio-Ref-Genome', pan_int, 'Pan-transcriptome', 'Interaction')

# =============================================================================
# Write out dataframes for plotting in R
# =============================================================================

rio_treat.to_csv("Rio_reference_genome/treat_genes_orthogroups.csv", index=False)
rio_type.to_csv("Rio_reference_genome/type_genes_orthogroups.csv", index=False)
rio_int.to_csv("Rio_reference_genome/int_genes_orthogroups.csv", index=False)

pan_treat.to_csv("pantranscriptome/treat_genes_orthogroups.csv", index=False)
pan_type.to_csv("pantranscriptome/type_genes_orthogroups.csv", index=False)
pan_int.to_csv("pantranscriptome/int_genes_orthogroups.csv", index=False)