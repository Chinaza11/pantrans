rm(list=ls())
ls()

library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# 
# =============================================================================

pantranscriptome <- read.csv("pantranscriptome/pantranscriptome_construction/pantranscriptome_plus_gene_type.csv")
annotation_count <- pantranscriptome %>% count(annotation, name="count")

annotation_count$percent <- round(annotation_count$count / sum(annotation_count$count) * 100, 1)
annotation_count$label <- paste0(annotation_count$count, " (", annotation_count$percent, "%)")

ggplot(annotation_count, aes(x="", y=count, fill=annotation)) +
  geom_bar(stat="identity", width=1) +  # Create bar chart as pie base
  coord_polar("y", start=0) +  # Convert to pie chart
  geom_text(aes(label=label), position=position_stack(vjust=0.5), size=5) +  # Add labels
  theme_void() + 
  theme(legend.title = element_blank()) +
  scale_fill_brewer(type="qual", palette=1)

# =============================================================================
# 
# =============================================================================

