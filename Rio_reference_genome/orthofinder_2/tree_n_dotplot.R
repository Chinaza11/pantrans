rm(list=ls())
ls()

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtree)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/orthofinder_2")

# =============================================================================
# Tree with node labels
# =============================================================================

tree = read.tree("SpeciesTree_rooted_node_labels.txt")
tree

png(file="tree.png", width=9, height=8, units="in", res=500)
ggtree(tree, branch.length = "none") + 
  geom_tiplab() +
  geom_nodelab(geom='label') + hexpand(.05)
dev.off()

# =============================================================================
# Tree with no node labels
# =============================================================================

png(file="tree_2.png", width=9, height=8, units="in", res=500)
ggtree(tree, branch.length = "none") + 
  geom_tiplab() + hexpand(.05)
dev.off()

# =============================================================================
# Tree with customized node labels
# =============================================================================

tree$node.label[tree$node.label == "N0"] = "Angiosperms"
tree$node.label[tree$node.label == "N2"] = "Monocots"
tree$node.label[tree$node.label == "N4"] = "Poaceae (Grasses)"
tree$node.label[tree$node.label == "N7"] = "Panicoideae"
tree$node.label[tree$node.label == "N11"] = "Andropogoneae"
tree$node.label[tree$node.label == "N13"] = "Sorghum"
tree$node.label[tree$node.label == "N14"] = "Sorghum (subgroup)"

nodes_to_replace_with_na = c("N1","N3","N5","N6","N8","N9","N10","N12","N15","N16")
tree$node.label[tree$node.label %in% nodes_to_replace_with_na] = NA

png(file="tree_3.png", width=9, height=9, units="in", res=500)
tree_3 = ggtree(tree, branch.length = "none") + 
  geom_tiplab() +
  geom_nodelab(geom='label', size=3.4) + hexpand(.07) + xlim(-0.2, NA)

tree_3
dev.off()

# =============================================================================
# Importing the annotation file and getting all the unique nodes
# =============================================================================

annot = read.csv("annot_node_ckpt5.csv")
subset_annot = annot[c("locusName", "Species.Tree.Node")]

subset_annot = separate_rows(subset_annot, Species.Tree.Node, sep=',')
unique_nodes = subset_annot$Species.Tree.Node %>% unique()

unique_nodes = sort(unique_nodes)
unique_nodes = unique_nodes[!is.na(unique_nodes)]

# =============================================================================
# Importing the enrichment result and combining them into one dataframe
# =============================================================================

trt_result = read.table("treat_enricher_df.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
subset_trt_result = trt_result[c("ID", "p.adjust", "Count")]
subset_trt_result$factor = 'Treatment'
row.names(subset_trt_result) = NULL

type_result = read.table("type_enricher_df.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
subset_type_result = type_result[c("ID", "p.adjust", "Count")]
subset_type_result$factor = 'Type'
row.names(subset_type_result) = NULL

int_result = read.table("int_enricher_df.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
subset_int_result = int_result[c("ID", "p.adjust", "Count")]
subset_int_result$factor = 'Interaction'
row.names(subset_int_result) = NULL

df = rbind(subset_trt_result,subset_type_result,subset_int_result)

for (i in unique_nodes){
  if(!(i %in% df$ID)){
    df = rbind(df, c(ID = i, p.adjust=NA, Count=NA, factor='Treatment,Type,Interaction'))
  }
}

df = separate_rows(df, factor, sep=',')

df$Count = as.numeric(df$Count)
df$p.adjust = as.numeric(df$p.adjust)

df = df %>% rename(node = ID)
df = df %>% rename(GeneCount = Count)

df$node[df$node == "N0"] = "Angiosperms"
df$node[df$node == "N2"] = "Monocots"
df$node[df$node == "N4"] = "Poaceae (Grasses)"
df$node[df$node == "N7"] = "Panicoideae"
df$node[df$node == "N11"] = "Andropogoneae"
df$node[df$node == "N13"] = "Sorghum"
df$node[df$node == "N14"] = "Sorghum (subgroup)"
df$node[df$node == "sorghumRio"] = "SorghumRio"
df$node[df$node == "SingleCopyOrthologues"] = "Single Copy Orthologues"
df$node[df$node == "UnassignedGenes"] = "No Orthorgroup"

# =============================================================================
# Plotting enrichment result
# =============================================================================

df$factor = factor(df$factor, levels = c("Treatment", "Type", "Interaction"))
df$node = factor(df$node, levels = rev(c("Angiosperms", "Monocots", "Poaceae (Grasses)", "Panicoideae", "Andropogoneae", "Sorghum", "Sorghum (subgroup)", "SorghumRio", "Single Copy Orthologues", "No Orthorgroup")))

png(file="enrichment_dotplot.png", width=9, height=8, units="in", res=500)
enrichment_dotplot=ggplot(df, aes(x=factor, y=node, color=p.adjust, size=GeneCount)) +
                    geom_point() +
                    xlab('') +
                    ylab('') + 
                    scale_color_gradient(low = "darkred", high = "blue") +
                    cowplot::theme_cowplot() + 
                    theme(
                      panel.grid.major = element_line(color = "grey90"),
                      panel.grid.minor = element_line(color = "grey98")) +
                    theme(axis.line  = element_blank()) +
                    theme(axis.ticks = element_blank())
enrichment_dotplot
dev.off()

# =============================================================================
# Combining tree and enrichment result in one layout
# =============================================================================

png(file="tree_n_dotplot.png", width=12, height=12, units="in", res=300)

ggpubr::ggarrange(tree_3, enrichment_dotplot, 
                  nrow=2,
                  labels = "AUTO")

dev.off()
