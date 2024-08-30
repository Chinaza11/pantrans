rm(list=ls())
ls()

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# UpSetR plot
# =============================================================================

library(UpSetR)

treat = read.table("treat_genes_df.txt", header=TRUE)
type = read.table("type_genes_df.txt", header=TRUE)
int = read.table("int_genes_df.txt", header=TRUE)

listInput = list(Treatment=treat$treat_genes, Type=type$type_genes, Interaction=int$int_genes)

png("glmmseq/UpSetR_plot.png", width=9, height=9, units="in", res=500)

upset(fromList(listInput), keep.order = T, sets = c("Interaction", "Type", "Treatment"), order.by="freq", mainbar.y.label="DEGs Intersections", sets.x.label="DEGs Per Effect", text.scale=c(1.5, 1.5, 1.25, 1.1, 2, 2), sets.bar.color=c('steelblue','purple','orange'), set_size.scale_max = 2000, shade.color = "gray50")

dev.off()
