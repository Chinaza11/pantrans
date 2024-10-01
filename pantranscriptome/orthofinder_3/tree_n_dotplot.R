rm(list=ls())
ls()

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtree)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/orthofinder_3")

# =============================================================================
# Tree with node label
# =============================================================================

tree = read.tree("SpeciesTree_rooted_node_labels.txt")
tree

png(file="tree.png", width=9, height=8, units="in", res=500)
ggtree(tree, branch.length = "none") + 
  geom_tiplab() +
  geom_nodelab(geom='label') + hexpand(.05)
dev.off()

# =============================================================================
# Tree with no node label
# =============================================================================

png(file="tree_2.png", width=9, height=8, units="in", res=500)
ggtree(tree, branch.length = "none") + 
  geom_tiplab() + hexpand(.15)
dev.off()

# =============================================================================
# Tree with customized node labels
# =============================================================================

tree$node.label[tree$node.label == "N0"] = "Plantae"
tree$node.label[tree$node.label == "N1"] = "Angiosperms"
tree$node.label[tree$node.label == "N2"] = "Eudicots"
tree$node.label[tree$node.label == "N4"] = "Monocots"
tree$node.label[tree$node.label == "N6"] = "Poales"
tree$node.label[tree$node.label == "N9"] = "Poaceae (Grasses)"
tree$node.label[tree$node.label == "N12"] = "Panicoideae"
tree$node.label[tree$node.label == "N14"] = "Andropogoneae"
tree$node.label[tree$node.label == "N16"] = "Sorghum"
tree$node.label[tree$node.label == "N17"] = "CP-NAM Pop. Parents"

nodes_to_replace_with_na = c("N3","N5","N7","N8","N10","N11","N13","N18","N19",
                             "N20","N21","N22","N23","N24","N25","N26","N27",
                             "N28","N29")

tree$node.label[tree$node.label %in% nodes_to_replace_with_na] = NA

png(file="tree_3.png", width=13, height=8, units="in", res=500)
tree_3 = ggtree(tree, branch.length = "none") + 
          geom_tiplab() +
          geom_nodelab(geom='label', size=3.4) + hexpand(.07)
tree_3
dev.off()

# =============================================================================
# get all unique species tree nodes/age of duplication event
# =============================================================================

annot = read.csv("annot_node_ckpt5.csv")
subset_annot = annot[c("locusName", "Species.Tree.Node")]

subset_annot = separate_rows(subset_annot, Species.Tree.Node, sep=',')
unique_nodes = subset_annot$Species.Tree.Node %>% unique()

unique_nodes = sort(unique_nodes)

# =============================================================================
# prepare enrichment result for plotting
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

df$node[df$node == "N0"] = "Plantae"
df$node[df$node == "N1"] = "Angiosperms"
df$node[df$node == "N2"] = "Eudicots"
df$node[df$node == "N4"] = "Monocots"
df$node[df$node == "N6"] = "Poales"
df$node[df$node == "N9"] = "Poaceae (Grasses)"
df$node[df$node == "N12"] = "Panicoideae"
df$node[df$node == "N14"] = "Andropogoneae"
df$node[df$node == "N16"] = "Sorghum"
df$node[df$node == "N17"] = "CP-NAM Pop. Parents"
df$node[df$node == "SingleCopyOrthologues"] = "Single Copy Orthologues"
df$node[df$node == "UnassignedGenes"] = "No Orthorgroup"

# =============================================================================
# plotting enrichment result
# =============================================================================

df$factor = factor(df$factor, levels = c("Treatment", "Type", "Interaction"))
df$node = factor(df$node, levels = rev(c("Plantae", "Angiosperms", "Eudicots", 
                                         "Monocots", "Poales", "Poaceae (Grasses)",
                                         "Panicoideae", "Andropogoneae", "N15",
                                         "Sorghum", "CP-NAM Pop. Parents", 
                                         "Single Copy Orthologues", "No Orthorgroup")))

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
# combining tree and enrichment result in one layout
# =============================================================================

png(file="tree_n_dotplot.png", width=12, height=12, units="in", res=300)

ggpubr::ggarrange(tree_3, enrichment_dotplot, 
                  nrow=2,
                  labels = "AUTO")

dev.off()

# =============================================================================
# Tree for Adobe Indesign
# =============================================================================

png(file="tree_4.png", width=11.5, height=11.5, units="in", res=500)
tree_4 = ggtree(tree, branch.length = "none") + 
  geom_tiplab(fontface='bold', size=3.75) +
  geom_nodelab(geom='label', size=4, angle=90, fontface='bold') + hexpand(.07)

tree_4
dev.off()


# =============================================================================
# Calculate diameter of enrichment circles for Adobe Indesign plot
# =============================================================================

get_circle_diameter = function(value){
  new_value = -log10(value)
  lower_bound = -log10(0.05)
  upper_bound = -log10(1*10^-10)
  
  value_percentile = 100 * ((new_value - lower_bound)/(upper_bound - lower_bound))
  
  adobe_scale_lower_bound = 0.07
  adobe_scale_upper_bound = 0.25
  
  diameter_of_circle = adobe_scale_lower_bound + 
    ((value_percentile/100) * (adobe_scale_upper_bound - adobe_scale_lower_bound))
  return(round(diameter_of_circle,3))
}

df_new = na.omit(df)
df_new$adobe_circle_diameter = get_circle_diameter(df_new$p.adjust)
