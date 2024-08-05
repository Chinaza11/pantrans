# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 20:38:28 2024

@author: Chinaza
"""

import requests, re, os, pandas as pd
from bs4 import BeautifulSoup

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/")

# =============================================================================
# import GTF file
# =============================================================================

column_names = ['seqname','source','feature','start','end','score','strand','frame','attribute']

genomic_gtf = pd.read_csv("others/genomic.gtf", sep='\t', comment='#', names=column_names, header=None, low_memory=False)

# =============================================================================
# process Treatment, Type and Interaction files
# =============================================================================

def process(file_path):    
    genes_df = pd.read_table(file_path)
    colname = genes_df.columns[0]
    
    # make hybrid genes be split into new rows
    expanded_genes = genes_df[colname].str.split('--').explode().reset_index(drop=True)
    
    # remove duplicate rows
    expanded_genes = expanded_genes.drop_duplicates().reset_index(drop=True)
    
    # convert panda series to a dataframe
    expanded_genes = expanded_genes.to_frame()
    
    return expanded_genes

treat_df = process("treat_genes_df.txt")
type_df = process("type_genes_df.txt")
int_df = process("int_genes_df.txt")

# =============================================================================
# get NCBI locus tag
# =============================================================================

def ncbi_locus_tag(df):
    colname = df.columns[0]
    for index, name in df[colname].items():
        matching_row = genomic_gtf[genomic_gtf['attribute'].str.contains(name, case=False, na=False, regex=False)]
        gene_row = matching_row[matching_row['feature'] == 'gene']
        i = gene_row.at[gene_row.index[0], 'attribute']
        j = i.split('locus_tag')
        k = j[1].strip()
        l = k[1:-2]
        df.at[index, 'NCBI_locus_tag'] = l
    for index, name in df['NCBI_locus_tag'].items():
        if 'partial' in name:
            a = name.split('"')[0]
            df.at[index, 'NCBI_locus_tag'] = a
    return df
            
treat_df = ncbi_locus_tag(treat_df)
type_df = ncbi_locus_tag(type_df)
int_df = ncbi_locus_tag(int_df)
        
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

gene_name = 'BDA96_01G006200'
result = main_function(gene_name)
print(result)

gene_name = 'BDA96_01G006200'
result = get_gene_description2(gene_name)
print(result)


# =============================================================================
# annotate genes with their description
# =============================================================================

def annot_gene_description(df):
    for index, name in df['NCBI_locus_tag'].items():
        res = main_function(name)
        if res:
            df.at[index, 'NCBI_gene_description'] = res
    return df

treat_df = annot_gene_description(treat_df)
nan = (treat_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {treat_df.shape[0]}')

type_df = annot_gene_description(type_df)
nan = (type_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {type_df.shape[0]}')

int_df = annot_gene_description(int_df)
nan = (int_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {int_df.shape[0]}')

""" COMMENT: All genes had no description which is weird """

# =============================================================================
# annotate hypothetical proteins
# =============================================================================

def annot_hypothetical_protein(df):
    for index, name in df['NCBI_locus_tag'].items():
        res = get_gene_description2(name)
        if res:
            df.at[index, 'NCBI_gene_description'] = res
    return df

treat_df = annot_hypothetical_protein(treat_df)
nan = (treat_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {treat_df.shape[0]}')
            
type_df = annot_hypothetical_protein(type_df)
nan = (type_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {type_df.shape[0]}')

int_df = annot_hypothetical_protein(int_df)
nan = (int_df['NCBI_gene_description'] == 'nan').sum()
print(f'Genes not yet annotated are {nan} out of {int_df.shape[0]}')

"""
COMMENT: All genes were labeled as hypothetical protein which does not make sense.
NEW PLAN: Each sequence will be blasted on NCBI and hits will be recorded. And I need the gene bank ID of each of the protein to achieve this. Script labelling_genes_2.py will extract the gene bank ID.
"""

# treat_df.to_csv("others/labeling_genes/treatment/treatment_genes_description.csv", index=False)
# type_df.to_csv("others/labeling_genes/type/type_genes_description.csv", index=False)
# int_df.to_csv("others/labeling_genes/interaction/interaction_genes_description.csv", index=False)


