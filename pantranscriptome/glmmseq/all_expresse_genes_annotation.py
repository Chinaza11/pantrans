# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 18:54:10 2025

@author: Chinaza
"""

import os, pandas as pd

os.chdir("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome")

# =============================================================================
# get list of genes expressed
# =============================================================================

raw_counts = pd.read_csv("glmmseq/rawCounts20_pantranscriptome.txt", delimiter='\t')

# mmquant only outputs expressed genes so raw_count shouldn't have any row with
# zero counts. However, below code will double-check this.
raw_counts[raw_counts[list(raw_counts.columns)].eq(0).all(axis=1)] # no row with zero counts

raw_counts["Gene"] = raw_counts["Gene"].str.split("--")
raw_counts = raw_counts.explode("Gene", ignore_index=True)

raw_counts_genes = list(set(raw_counts["Gene"].tolist()))

# =============================================================================
# get list of pantranscriptome core, shell and cloud genes
# =============================================================================

pantranscriptome_df = pd.read_csv("pantranscriptome_construction/pantranscriptome_plus_gene_type.csv")

core_genes = pantranscriptome_df[pantranscriptome_df["annotation"] == "core"]["sequence_name"].tolist()
shell_genes = pantranscriptome_df[pantranscriptome_df["annotation"] == "shell"]["sequence_name"].tolist()
cloud_genes = pantranscriptome_df[pantranscriptome_df["annotation"] == "cloud"]["sequence_name"].tolist()

# =============================================================================
# get no of genes expressed that are core, shell and cloud
# =============================================================================

core_genes_expressed = list(set(raw_counts_genes) & set(core_genes))
shell_genes_expressed = list(set(raw_counts_genes) & set(shell_genes))
cloud_genes_expressed = list(set(raw_counts_genes) & set(cloud_genes))
