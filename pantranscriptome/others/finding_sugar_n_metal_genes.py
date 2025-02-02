# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:53:09 2025

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# import sugar genes and metal genes
# =============================================================================

sugar_genes = pd.read_excel("others/labeling_genes/Sorghum Gene Families of Interest.xlsx", sheet_name="Sugar Transport & Synthesis")

metal_genes = pd.read_excel("others/labeling_genes/Sorghum Gene Families of Interest.xlsx", sheet_name="Metal Transport Genes")

# =============================================================================
# import treatment, type and int genes
# =============================================================================

treatment_genes = pd.read_excel("others/labeling_genes/treatment/treatment_genes_description_final.xlsx", sheet_name="treatment_genes_description_2")

type_genes = pd.read_csv("others/labeling_genes/type/type_genes_description-after_permutation.csv")

int_genes = pd.read_csv("others/labeling_genes/interaction/interaction_genes_description-after_permutation.csv")

# =============================================================================
# function to find sugar and metal genes present in the fixed effects
# =============================================================================

def find_matching_rows(sugar_or_metal_genes_df, fixed_effect_df, fixed_effect):
    matches_dict = {}
    for i in sugar_or_metal_genes_df["BTx623 ID"]:
        x = i.split('.')[1]
        match = fixed_effect_df[fixed_effect_df["stable_id"].str.contains(x, na=False)]
        if len(match) > 1:
            print("inspect dataframe for multiple matches")
        elif len(match) == 1:
            matches_dict[i] = match["stable_id"].tolist()
            
    table_df = pd.DataFrame([
        {"BTx623_ID": key, "Stable_ID": value} 
        for key, values in matches_dict.items() 
        for value in values
    ])

    table_df["Fixed_Effect"] = fixed_effect
    return table_df

# =============================================================================
# find sugar genes and merge in one df
# =============================================================================

sugar_in_trt = find_matching_rows(sugar_genes, treatment_genes, "Treatment")
sugar_in_type = find_matching_rows(sugar_genes, type_genes, "Type")
sugar_in_int = find_matching_rows(sugar_genes, int_genes, "Interaction")

combined_sugar = pd.concat([sugar_in_trt, 
                            sugar_in_type, 
                            sugar_in_int], ignore_index=True)

combined_sugar["Gene_Type"] = "sugar"

# =============================================================================
# find metal genes and merge in one df
# =============================================================================

metal_in_trt = find_matching_rows(metal_genes, treatment_genes, "Treatment")
metal_in_type = find_matching_rows(metal_genes, type_genes, "Type")
metal_in_int = find_matching_rows(metal_genes, int_genes, "Interaction")

combined_metal = pd.concat([metal_in_trt, 
                            metal_in_type, 
                            metal_in_int], ignore_index=True)

combined_metal["Gene_Type"] = "metal"

# =============================================================================
# combine sugar and metal df; merge with relevant info in other df
# =============================================================================

combined = pd.concat([combined_sugar, 
                      combined_metal], ignore_index=True)

combined = pd.merge(combined, sugar_genes, left_on="BTx623_ID", right_on="BTx623 ID", how="left", suffixes=("", ""))
combined.drop(columns=["BTx623 ID"], inplace=True)

combined = pd.merge(combined, metal_genes, left_on="BTx623_ID", right_on="BTx623 ID", how="left", suffixes=("", "_"))

combined.drop(columns=["BTx623 ID", "Unnamed: 3", 
                       "Unnamed: 4", "Unnamed: 5",
                       "Rio Ortholog"], inplace=True)

combined["Gene_Name"] = combined["Gene Name"].fillna(combined["Gene Name_"])
combined.drop(columns=["Gene Name", "Gene Name_"], inplace=True)

combined = pd.merge(combined, 
                    treatment_genes[["stable_id", 
                                     "NCBI_gene_description"]], 
                    left_on="Stable_ID", 
                    right_on="stable_id",
                    how="left")

combined = pd.merge(combined, 
                    type_genes[["stable_id", 
                                     "NCBI_gene_description"]], 
                    left_on="Stable_ID", 
                    right_on="stable_id",
                    how="left",
                    suffixes=("", "_"))

combined["NCBI_gene_description"] = combined["NCBI_gene_description"].fillna(combined["NCBI_gene_description_"])
combined.drop(columns=["stable_id", "stable_id_", 
                       "NCBI_gene_description_"], inplace=True)


combined.rename(columns={"Unnamed: 2": "Gene_Description_From_File_Dr_Cooper_Shared"}, inplace=True)
combined.rename(columns={"NCBI_gene_description": "NCBI_Gene_Description"}, inplace=True)


agg_dict = {                       
    'Stable_ID': 'first',
    'Fixed_Effect': lambda x: ', '.join(map(str, x)),
    'Gene_Type': 'first',                         
    'Gene_Description_From_File_Dr_Cooper_Shared': 'first',
    'Gene_Name': 'first',                         
    'NCBI_Gene_Description': 'first'
}

combined = combined.groupby('BTx623_ID').agg(agg_dict).reset_index()

combined = combined.sort_values(by='Gene_Type', ascending=False).reset_index(drop=True)

combined.to_csv('others/labeling_genes/sugar_and_metal_genes.csv', index=False)