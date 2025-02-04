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
raw_counts_explode = raw_counts.explode("Gene", ignore_index=True)

raw_counts_genes = list(set(raw_counts_explode["Gene"].tolist()))

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

# =============================================================================
# dataframe for piechart in R
# =============================================================================

expressed_genes_annotation = pd.DataFrame(
    [(gene, "core") for gene in core_genes_expressed
    ] + [(gene, "shell") for gene in shell_genes_expressed
    ] + [(gene, "cloud") for gene in cloud_genes_expressed
    ], columns=["gene", "annotation"])

expressed_genes_annotation.to_csv("glmmseq/expressed_genes_annotation.csv", index=False)

# =============================================================================
# dataframe for violin plot in R
# =============================================================================

for index, gene in raw_counts['Gene'].items():
    if "--" not in gene:
        if gene in core_genes:
            raw_counts.at[index, 'annotation'] = 'core'
        elif gene in shell_genes:
            raw_counts.at[index, 'annotation'] = 'shell'
        elif gene in cloud_genes:
            raw_counts.at[index, 'annotation'] = 'cloud'
    elif "--" in gene:
        hybrid_genes = []
        for i in gene.split("--"):
            if i in core_genes:
                hybrid_genes.append('core')
            elif i in shell_genes:
                hybrid_genes.append('shell')
            elif i in cloud_genes:
                hybrid_genes.append('cloud')
        raw_counts.at[index, 'annotation'] = ", ".join(list(set(hybrid_genes)))

raw_counts["annotation"] = raw_counts["annotation"].str.split(", ")
raw_counts = raw_counts.explode("annotation", ignore_index=True)

raw_counts["Average"] = raw_counts.iloc[:, 1:-1].mean(axis=1)

raw_counts.to_csv('glmmseq/raw_count_for_violin_plot.csv', index=False)
