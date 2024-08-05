# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:46:14 2024

@author: Chinaza
"""

import requests, re, os, pandas as pd
from bs4 import BeautifulSoup

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# functions for getting gene description
# =============================================================================

def get_gene_id(gene_name):
    url = f"https://www.ncbi.nlm.nih.gov/search/all/?term={gene_name}" # Construct the URL    
    response = requests.get(url) # Send a GET request to the URL
    if response.status_code == 200: # Check if the request was successful (status code 200)
        # Parse the HTML content of the page
        soup = BeautifulSoup(response.content, "html.parser")
        soup_str = str(soup) # convert HTML content to string
        i = re.search('Gene ID: ', soup_str) # look for position of exact match for substring 
        if i:
            r = soup_str.find('</li>', i.start()) # find end position of info to extract
            gene_id = soup_str[i.end():r] #extract gene id
            if gene_id:
                return gene_id
            else:
                print("Failed to retrieve Gene ID")
        else:
            print(f"Gene ID info for {gene_name} not found on webpage")
    else:
        print("Failed to retrieve data from NCBI website: gene ID function")
        

def get_gene_description(gene_id):
    url = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")
        soup_str = str(soup)
        i = re.search('Gene description', soup_str)
        if i:
            r = soup_str.find('<dd>', i.end()) + 4
            s = soup_str.find('</dd>', i.end())
            gene_description = soup_str[r:s]
            if gene_description:
                return gene_description
            else:
                print("Failed to retrieve gene description")
        else:
            print("Gene description info not found on webpage")
    else:
        print("Failed to retrieve data from NCBI website: gene description function")


def main_function(gene_name):
    gene_id = get_gene_id(gene_name)
    if gene_id:
        gene_description = get_gene_description(gene_id)
        return gene_description

# -----------------------------------------------------------------------------
def get_gene_description2(gene_name):
    """annotating hypothetical proteins"""
    url = f"https://www.ncbi.nlm.nih.gov/search/all/?term={gene_name}"
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")        
        soup_str = str(soup)
        i = re.search('hypothetical protein', soup_str)
        if i:
            return "hypothetical protein"
        else:
            print(f"Gene description info for {gene_name} not found on webpage")
    else:
        print(f"Failed to retrieve data from NCBI website for {gene_name}")


# =============================================================================
# testing
# =============================================================================

gene_name = 'SORBI_3002G226000'
result = main_function(gene_name)
print(result)

gene_name = 'SORBI_3010G278500'
result = main_function(gene_name)
print(result)

# =============================================================================
# import pantranscriptome
# =============================================================================

pantranscriptome = pd.read_csv("pantranscriptome_construction/pantranscriptome.csv")

# =============================================================================
# import gene list and left join it to the pantranscriptome
# =============================================================================

def process(file_path):
    df = pd.read_table(file_path)    
    colname = df.columns[0]
    
    # make hybrid genes be split into new rows
    df = df[colname].str.split('--').explode().reset_index(drop=True)
    
    # rename header in preparation for merger and drop duplicates
    df = df.rename('sequence_name')
    df = df.drop_duplicates()
    
    # merge gene list and pantranscritome file
    df = pd.merge(df, pantranscriptome, on='sequence_name', how='left')
    return df

treat_genes_df = process("treat_genes_df.txt")
type_genes_df = process("type_genes_df.txt")
int_genes_df = process("int_genes_df.txt")

# =============================================================================
# annotate genes with their description
# =============================================================================

def annot_genes(df):
    for index, name in df['stable_id'].items():
        res = main_function(name)
        if res:
            df.at[index, 'NCBI_gene_description'] = res
    return df        
    

treat_genes_df = annot_genes(treat_genes_df)
nan = (treat_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 545 out of 1324 not yet annotated

type_genes_df = annot_genes(type_genes_df)
nan = (type_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 1109 out of 1891 not yet annotated

int_genes_df = annot_genes(int_genes_df)
nan = (int_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 175 out of 273 not yet annotated

# =============================================================================
# annotate hypothetical proteins
# =============================================================================

def annot_genes2(df):
    for index, name in df['stable_id'].items():
        if df.at[index, 'NCBI_gene_description'] == 'nan':
            res = get_gene_description2(name)
            if res:
                df.at[index, 'NCBI_gene_description'] = res
    return df

treat_genes_df = annot_genes2(treat_genes_df)        
nan = (treat_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 447 out of 1324 not yet annotated

type_genes_df = annot_genes2(type_genes_df)        
nan = (type_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 928 out of 1891 not yet annotated

int_genes_df = annot_genes2(int_genes_df)        
nan = (int_genes_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan}') # 156 out of 273 not yet annotated

# treat_genes_df.to_csv("others/labeling_genes/treatment/treatment_genes_description.csv", index=False)
# type_genes_df.to_csv("others/labeling_genes/type/type_genes_description.csv", index=False)
# int_genes_df.to_csv("others/labeling_genes/interaction/interaction_genes_description.csv", index=False)
