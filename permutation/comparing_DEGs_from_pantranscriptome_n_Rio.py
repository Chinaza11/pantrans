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
Orthogroups = pd.read_csv("pantranscriptome/orthofinder/Orthogroups.tsv", sep='\t')

# -----------------------------------------------------------------------------
# importing dataframe with UnassignedGenes (genes unassigned to any orthogroup)
Orthogroups_UnassignedGenes = pd.read_csv("pantranscriptome/orthofinder/Orthogroups_UnassignedGenes.tsv", sep='\t', low_memory=False)

# -----------------------------------------------------------------------------
# DEGs: using pantranscriptome
pan_type = pd.read_table("pantranscriptome/glmmseq/permutation_analysis/type_genes_df.txt")
pan_type.head()
pan_type['orthogroup'] = ""

pan_int = pd.read_table("pantranscriptome/glmmseq/permutation_analysis/int_genes_df.txt")
pan_int.head()
pan_int['orthogroup'] = ""

# -----------------------------------------------------------------------------
# DEGs: using Rio reference genome
rio_type = pd.read_table("Rio_reference_genome/glmmseq/permutation_analysis/type_genes_df.txt")
rio_type.head()
rio_type['orthogroup'] = ""

rio_int = pd.read_table("Rio_reference_genome/glmmseq/permutation_analysis/int_genes_df.txt")
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
annot_unassigned_single_rio_genes(rio_int,"int_genes")  

# -----------------------------------------------------------------------------
# annotating hybrid genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
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
annot_unassigned_hybrid_rio_genes(rio_int,"int_genes")


# =============================================================================
# annotating DEGs gotten using pantranscriptome
# =============================================================================

# -----------------------------------------------------------------------------
# annotating from orthogroups file (Orthogroups.tsv)
def annot_assigned_pan_genes(df, column):
    for index, gene in df[column].items():
                
        if gene.find("--") != -1:
            genes = gene.split("--")
            orthogroup_list = []
            for value in genes:
                if value.startswith("SbiCamber"):
                    for index2, gene2 in Orthogroups['sorghumChineseamber'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiGrassl"):
                    for index2, gene2 in Orthogroups['sorghumGrassl'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiLeoti"):
                    for index2, gene2 in Orthogroups['sorghumLeoti'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI229841"):
                    for index2, gene2 in Orthogroups['sorghumPI229841'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])                      
                elif value.startswith("SbiPI297155"):
                    for index2, gene2 in Orthogroups['sorghumPI297155'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI329311"):
                    for index2, gene2 in Orthogroups['sorghumPI329311'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI506069"):
                    for index2, gene2 in Orthogroups['sorghumPI506069'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI510757"):
                    for index2, gene2 in Orthogroups['sorghumPI510757'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])                       
                elif value.startswith("SbiPI655972"):
                    for index2, gene2 in Orthogroups['sorghumPI655972'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI563295"):
                    for index2, gene2 in Orthogroups['sorghumRioNAM'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups.at[index2,"Orthogroup"])
            
            if len(orthogroup_list) > 0:
                deduplicate_orthogroup_list = []
                for val in orthogroup_list:
                    if val not in deduplicate_orthogroup_list:
                        deduplicate_orthogroup_list.append(val)
                        
                df.at[index,"orthogroup"] = ", ".join(deduplicate_orthogroup_list)
            
        else:
            if gene.startswith("SbiCamber"):
                for index2, gene2 in Orthogroups['sorghumChineseamber'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiGrassl"):
                for index2, gene2 in Orthogroups['sorghumGrassl'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiLeoti"):
                for index2, gene2 in Orthogroups['sorghumLeoti'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI229841"):
                for index2, gene2 in Orthogroups['sorghumPI229841'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]                        
            elif gene.startswith("SbiPI297155"):
                for index2, gene2 in Orthogroups['sorghumPI297155'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI329311"):
                for index2, gene2 in Orthogroups['sorghumPI329311'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI506069"):
                for index2, gene2 in Orthogroups['sorghumPI506069'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI510757"):
                for index2, gene2 in Orthogroups['sorghumPI510757'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]                         
            elif gene.startswith("SbiPI655972"):
                for index2, gene2 in Orthogroups['sorghumPI655972'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI563295"):
                for index2, gene2 in Orthogroups['sorghumRioNAM'].items():
                    if gene in str(gene2):
                        df.at[index,"orthogroup"] = Orthogroups.at[index2,"Orthogroup"]                          


annot_assigned_pan_genes(pan_type,"type_genes")
annot_assigned_pan_genes(pan_int,"int_genes")


# -----------------------------------------------------------------------------
# annotating single genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_single_pan_genes(df, column):
    for index, group in df["orthogroup"].items():
        if len(group) == 0 :
            gene = df.at[index, column]
            
            if gene.startswith("SbiCamber"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumChineseamber'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]
            elif gene.startswith("SbiGrassl"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumGrassl'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]            
            elif gene.startswith("SbiLeoti"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumLeoti'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI229841"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI229841'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]  
            elif gene.startswith("SbiPI297155"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI297155'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI329311"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI329311'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]            
            elif gene.startswith("SbiPI506069"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI506069'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI510757"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI510757'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"] 
            elif gene.startswith("SbiPI655972"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI655972'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"]
            elif gene.startswith("SbiPI563295"):
                for index2, gene2 in Orthogroups_UnassignedGenes['sorghumRioNAM'].items():
                    if pd.isna(gene2) != True:
                        if gene in str(gene2):
                            df.at[index,"orthogroup"] = Orthogroups_UnassignedGenes.at[index2,"Orthogroup"] 
                                            
annot_unassigned_single_pan_genes(pan_type,"type_genes")
annot_unassigned_single_pan_genes(pan_int,"int_genes")  


# -----------------------------------------------------------------------------
# annotating hybrid genes unassigned to no orthogroups using this file (Orthogroups_UnassignedGenes.tsv)
def annot_unassigned_hybrid_pan_genes(df, column):
    for index, gene in df[column].items():    
        if gene.find("--") != -1:
            genes = gene.split("--")
            orthogroup_list = []
            for value in genes:
                if value.startswith("SbiCamber"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumChineseamber'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiGrassl"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumGrassl'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiLeoti"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumLeoti'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI229841"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI229841'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])                      
                elif value.startswith("SbiPI297155"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI297155'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI329311"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI329311'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI506069"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI506069'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI510757"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI510757'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])                       
                elif value.startswith("SbiPI655972"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumPI655972'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
                elif value.startswith("SbiPI563295"):
                    for index2, gene2 in Orthogroups_UnassignedGenes['sorghumRioNAM'].items():
                        if value in str(gene2):
                            orthogroup_list.append(Orthogroups_UnassignedGenes.at[index2,"Orthogroup"])
            if len(orthogroup_list) > 0:
                deduplicate_orthogroup_list = []
                for val in orthogroup_list:
                    if val not in deduplicate_orthogroup_list:
                        deduplicate_orthogroup_list.append(val)
                    
                deduplicate_orthogroup_str = ', ' + ", ".join(deduplicate_orthogroup_list)
                df.at[index,"orthogroup"] += deduplicate_orthogroup_str            

    
annot_unassigned_hybrid_pan_genes(pan_type,"type_genes")
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
    #plt.savefig(f'DEGs_Orthogroups_({name})', dpi=300) 
    plt.show()
    
comparison(rio_type, 'Rio-Ref-Genome', pan_type, 'Pan-transcriptome', 'Type')
comparison(rio_int, 'Rio-Ref-Genome', pan_int, 'Pan-transcriptome', 'Interaction')


# =============================================================================
# Write out dataframes for plotting in R
# =============================================================================

rio_type.to_csv("Rio_reference_genome/permutation/type_genes_orthogroups.csv", index=False)
rio_int.to_csv("Rio_reference_genome/permutation/int_genes_orthogroups.csv", index=False)

pan_type.to_csv("pantranscriptome/permutation/type_genes_orthogroups.csv", index=False)
pan_int.to_csv("pantranscriptome/permutation/int_genes_orthogroups.csv", index=False)
