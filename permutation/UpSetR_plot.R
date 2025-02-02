rm(list=ls())
ls()

library(UpSetR)
library(cowplot)
library(gridExtra)
library(ggplot2)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# Rio reference genome
# =============================================================================

treat = read.table("Rio_reference_genome/treat_genes_df.txt", header=TRUE)
type = read.table("Rio_reference_genome/glmmseq/permutation_analysis/type_genes_df.txt", header=TRUE)
int = read.table("Rio_reference_genome/glmmseq/permutation_analysis/int_genes_df.txt", header=TRUE)

rio_listInput = list(Treatment=treat$treat_genes, Type=type$type_genes, Interaction=int$int_genes)

rio = upset(fromList(rio_listInput), 
            keep.order = T, 
            sets = c("Interaction", "Type", "Treatment"),
            mainbar.y.max = 1200,
            order.by="freq", 
            mainbar.y.label="DEGs Intersections", 
            sets.x.label="DEGs Per Effect", 
            text.scale=c(1.25, 1.25, 1, 0.85, 1.5, 1.75), 
            sets.bar.color=c('steelblue','purple','orange'),
            set_size.scale_max = 2000,
            shade.color = "gray50")
rio

# =============================================================================
# Pan-transcriptome
# =============================================================================

treat = read.table("pantranscriptome/treat_genes_df.txt", header=TRUE)
type = read.table("pantranscriptome/glmmseq/permutation_analysis/type_genes_df.txt", header=TRUE)
int = read.table("pantranscriptome/glmmseq/permutation_analysis/int_genes_df.txt", header=TRUE)

pant_listInput = list(Treatment=treat$treat_genes, Type=type$type_genes, Interaction=int$int_genes)

pan_t = upset(fromList(pant_listInput), 
              keep.order = T, sets = c("Interaction", "Type", "Treatment"), 
              mainbar.y.max = 1200,
              order.by="freq", 
              mainbar.y.label="DEGs Intersections", 
              sets.x.label="DEGs Per Effect", 
              text.scale=c(1.25, 1.25, 1, 0.85, 1.5, 1.75), 
              sets.bar.color=c('steelblue','purple','orange'),
              set_size.scale_max = 2000,
              shade.color = "gray50")
pan_t

png("permutation/pantranscriptome_UpSetR_plot.jpg", width=20, height=20, units="cm", res=300)
pan_t
dev.off()

# =============================================================================
# Both plots in one layout
# =============================================================================

png("permutation/UpSetR_plot.jpg", width=20, height=20, units="cm", res=300)

rio_x = plot_grid(NULL, rio$Main_bar, rio$Sizes, rio$Matrix, 
                  nrow=2, align='hv', rel_heights = c(3.5,1.5), 
                  rel_widths = c(1.5,3.5), labels="A", label_x=0.3)
pan_t_x = plot_grid(NULL, pan_t$Main_bar, pan_t$Sizes, pan_t$Matrix, 
                    nrow=2, align='hv', rel_heights = c(3.5,1.5), 
                    rel_widths = c(1.5,3.5),  labels="B", label_x=0.3)

grid.arrange(rio_x + theme(plot.margin = unit(c(0, 0, 0, -2), "cm")),
             pan_t_x + theme(plot.margin = unit(c(0, 0, 0, -2), "cm")), 
             nrow=2, ncol=1)
dev.off()
