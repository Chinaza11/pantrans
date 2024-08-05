rm(list=ls())
ls()

library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(ggupset)
library(ggplot2)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# read in the annotation file and subset it
# =============================================================================

annot = read.csv("orthofinder_2/annot_node_ckpt5.csv")

subset_annot = annot[c("locusName", "Species.Tree.Node")]

# =============================================================================
# import glmmseq result and prepare it
# =============================================================================

stats = read.table("glmmseq/glmm-allGenes-results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# creating the gene universe
gene_universe = rownames(stats)

#------------------------------------------------------------------------------
treatList = stats$qvals.Treatment
names(treatList) = gene_universe
subset_treatList = treatList[treatList <= 0.05]
treat_genes = names(subset_treatList)

#------------------------------------------------------------------------------
result_qval_type = read.table("glmmseq/permutation_analysis/result_qval_type.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

qval_type_gene_names = result_qval_type$Gene
qval_type = result_qval_type$perm_pval_type
names(qval_type) = qval_type_gene_names
sig_qval_type = qval_type[qval_type <= 0.05]
type_genes = names(sig_qval_type)

#------------------------------------------------------------------------------
result_qval_trt_type = read.table("glmmseq/permutation_analysis/result_qval_trt_type.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

qval_trt_type_gene_names = result_qval_trt_type$Gene
qval_trt_type = result_qval_trt_type$perm_pval_trt_type
names(qval_trt_type) = qval_trt_type_gene_names
sig_qval_trt_type = qval_trt_type[qval_trt_type <= 0.05]

int_genes = names(sig_qval_trt_type)

#------------------------------------------------------------------------------
# genes were written out to a .txt file so I could compare it with genes from DEGs analysis using the Rio reference genome (comparison of their orthogroups)

# type_genes_df = data.frame(type_genes)
# treat_genes_df = data.frame(treat_genes)
# int_genes_df = data.frame(int_genes)
# 
# write.table(type_genes_df, "type_genes_df.txt", quote=FALSE, sep="\t")
# write.table(treat_genes_df, "treat_genes_df.txt", quote=FALSE, sep="\t")
# write.table(int_genes_df, "int_genes_df.txt", quote=FALSE, sep="\t")
#------------------------------------------------------------------------------
# annotating all hybrid genes
# hybrid genes here means reads that mapped to more than one locus, these genes are joined by "--" e.g. SbRio.02G124800--SbRio.07G193200

# finding all rownames that has "--" and putting this in a list
hyb_index = grep("--", rownames(stats))

# create an empty dataframe to store the hybrid genes
hyb_genes = data.frame() 

# finding all hybrid genes and adding them to the new dataframe created
for (i in hyb_index) {
  selected_row = stats[i,] # select the row based on the index
  hyb_genes = rbind(hyb_genes, selected_row) # append the selected row to the new dataframe
}

# subseting dataframe to include only the q-vals for type, treatment and interactions
hyb_genes = hyb_genes[,c(22,23,24)]

# annotating all hybrid genes
for (i in 1:length(rownames(hyb_genes))) {
  temp = rownames(hyb_genes)[i]
  genes = unlist(strsplit(temp, split="--"))
  m = match(genes, subset_annot$locusName)
  hyb_genes$Species.Tree.Node[i] = paste0(unique(subset_annot$Species.Tree.Node[m]), collapse=",")
}

# sub-setting annotated hybrid genes file and selecting only the gene names and Species.Tree.Node annotation
subset_hyb_genes = data.frame(cbind(rownames(hyb_genes), hyb_genes$Species.Tree.Node))

# renaming columns in preparation for merger with main annotation file
colnames(subset_hyb_genes)[1] = "locusName"
colnames(subset_hyb_genes)[2] = "Species.Tree.Node"

# removing duplicate entries in each cell of subset_hyb_genes
for (i in 1:nrow(subset_hyb_genes)) {
  cell_content = subset_hyb_genes$Species.Tree.Node[i]
  unique_entries <- unique(unlist(strsplit(cell_content, ",")))
  subset_hyb_genes$Species.Tree.Node[i] = paste(unique_entries, collapse = ",")
}
#------------------------------------------------------------------------------

# merging hybrid genes annotation with other genes annotation
subset_annot = rbind(subset_annot, subset_hyb_genes)

# making each annotation separated by "," to be a row of its own
subset_annot = separate_rows(subset_annot, Species.Tree.Node, sep = ",")

# switching the columns placement in preparation for enrichment analysis. enricher() and GSEA() functions requires the file to be in TERM2GENE format
subset_annot = subset_annot[, c("Species.Tree.Node", "locusName")]

na_count = sum(is.na(subset_annot$Species.Tree.Node))

# making sure all NA annotations (from hybrid genes) reads as NA and not as another annotation category called "NA" 
subset_annot$Species.Tree.Node[subset_annot$Species.Tree.Node == "NA"] = NA

na_count = sum(is.na(subset_annot$Species.Tree.Node))

# =============================================================================
# over representation analysis using hypergeometric test. maxGSSize was set as 20,000 
# if not the default of 500 annotation will exclude several genes and lead to wrong
# result. 
# Max annotation for any single node is less than 20,000. Use this code to confirm:
# table(subset_annot$Species.Tree.Node)
# =============================================================================

treat_enricher = enricher(treat_genes, minGSSize=1, maxGSSize=20000, universe=gene_universe, TERM2GENE=subset_annot)
treat_enricher_df = data.frame(treat_enricher)

barplot(treat_enricher, showCategory=20) 

dotplot(treat_enricher, showCategory=30) + 
  ggtitle("Treatment") + 
  theme(
    plot.title = element_text(size = 24, face ="bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face ="bold"), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 18, face ="bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
    )

#------------------------------------------------------------------------------

type_enricher = enricher(type_genes, minGSSize=1, maxGSSize=20000, universe=gene_universe, TERM2GENE=subset_annot)
type_enricher_df = data.frame(type_enricher)
write.table(type_enricher_df, "orthofinder_2/permutation/type_enricher_df.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)
barplot(type_enricher, showCategory=20) 

dotplot(type_enricher, showCategory=30) + 
  ggtitle("Type") + 
  theme(
    plot.title = element_text(size = 24, face ="bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face ="bold"), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 18, face ="bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

#------------------------------------------------------------------------------

int_enricher = enricher(int_genes, minGSSize=1, maxGSSize=20000, universe=gene_universe, TERM2GENE=subset_annot)
int_enricher_df = data.frame(int_enricher)

## no result so code below was commented out

# write.table(int_enricher_df, "orthofinder_2/permutation/int_enricher_df.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)
# barplot(int_enricher, showCategory=20) 
# 
# dotplot(int_enricher, showCategory=30) + 
#   ggtitle("Interaction") + 
#   theme(
#     plot.title = element_text(size = 24, face ="bold", hjust = 0.5),
#     axis.title.x = element_text(size = 14, face ="bold"), 
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 18, face ="bold"),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 14)
#   )