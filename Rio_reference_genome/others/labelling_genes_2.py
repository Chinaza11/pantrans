# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:17:06 2024

@author: Chinaza
"""

import requests, re, os, pandas as pd
from bs4 import BeautifulSoup

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/")

# =============================================================================
# 
# =============================================================================

def get_protein_id(locus_tag):
    url = f"https://www.ncbi.nlm.nih.gov/search/all/?term={locus_tag}"
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")
        soup_str = str(soup)
        i = re.search('Sequence length:', soup_str)
        if i:
            r = soup_str.find('<li>', i.end()) + 4
            s = soup_str.find('</li>', i.end())
            genbank_id = soup_str[r:s]
            if genbank_id:
                return genbank_id
            else:
                print("Failed to retrieve genbank id")
        else:
            print("Genebank ID not found on webpage")
    else:
        print("Failed to retrieve data from NCBI website")


def annot_df(df):
    for index, name in df['NCBI_locus_tag'].items():
        res = get_protein_id(name)
        if res:
            df.at[index, 'protein_genbank_id'] = res
    return df

# =============================================================================
# 
# =============================================================================

treat_df = pd.read_csv("others/labeling_genes/treatment/treatment_genes_description.csv")
treat_df = annot_df(treat_df)
treat_df.to_csv("others/labeling_genes/treatment/treatment_genes_description_2.csv", index=False)
treat_df_2 = treat_df.iloc[:, 3:]
treat_df_2 = treat_df_2[treat_df_2['protein_genbank_id'] != 'nan']
treat_df_2.to_csv("others/labeling_genes/treatment/treatment_protein_id_for_hpc.csv", index=False, header=None)

type_df = pd.read_csv("others/labeling_genes/type/type_genes_description.csv")
type_df = annot_df(type_df)
type_df.to_csv("others/labeling_genes/type/type_genes_description_2.csv", index=False)
type_df_2 = type_df.iloc[:, 3:]
type_df_2 = type_df_2[type_df_2['protein_genbank_id'] != 'nan']
type_df_2.to_csv("others/labeling_genes/type/type_protein_id_for_hpc.csv", index=False, header=None)

int_genes = pd.read_csv("others/labeling_genes/interaction/interaction_genes_description.csv")
int_genes = annot_df(int_genes)
int_genes.to_csv("others/labeling_genes/interaction/interaction_genes_description_2.csv", index=False)
int_genes_2 = int_genes.iloc[:, 3:]
int_genes_2.to_csv("others/labeling_genes/interaction/int_protein_id_for_hpc.csv", index=False, header=None)
