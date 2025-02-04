rm(list=ls())
ls()

library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)


setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# 
# =============================================================================

pie_chart <- function(df){
  annotation_count <- df %>% count(annotation, name="count")
  
  annotation_count$percent <- round(annotation_count$count / sum(annotation_count$count) * 100, 1)
  annotation_count$label <- paste0(annotation_count$count, " (", annotation_count$percent, "%)")
  
  ggplot(annotation_count, aes(x="", y=count, fill=annotation)) +
    geom_bar(stat="identity", width=1) +  # Create bar chart as pie base
    coord_polar("y", start=0) +  # Convert to pie chart
    geom_text(aes(label=label), position=position_stack(vjust=0.5), size=2) +  # Add labels
    theme_void() + 
    theme(legend.title = element_blank()) +
    scale_fill_brewer(type="qual", palette=1)
}

# =============================================================================
# 
# =============================================================================

pantranscriptome <- read.csv("pantranscriptome/pantranscriptome_construction/pantranscriptome_plus_gene_type.csv")

pantranscriptome_pie_chart <- pie_chart(pantranscriptome)
pantranscriptome_pie_chart

# =============================================================================
# 
# =============================================================================

expressed_genes <- read.csv("pantranscriptome/glmmseq/expressed_genes_annotation.csv")

expressed_genes_pie_chart <- pie_chart(expressed_genes)
expressed_genes_pie_chart

# =============================================================================
# 
# =============================================================================

ggarrange(pantranscriptome_pie_chart, expressed_genes_pie_chart, 
          nrow=1, common.legend=T, legend="right",
          labels = "AUTO")

# =============================================================================
# 
# =============================================================================

df_for_violin_plot <- read.csv("pantranscriptome/glmmseq/raw_count_for_violin_plot.csv")

violin_plot <- ggplot(df_for_violin_plot, aes(x=annotation, y=log10(Average), fill=annotation)) +
                  geom_violin(trim=F) +
                  geom_boxplot(width=0.1, fill="white") +
                  scale_fill_brewer(type="qual", palette=1) +
                  theme_classic2() + 
                  theme(axis.text.x = element_blank(),  
                        axis.ticks.x = element_blank(),
                        legend.title = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size = 7.5)) +
  labs(y = expression(Log[10] ~ "(average read count across the samples)"))

violin_plot

# =============================================================================
# 
# =============================================================================

# Extract the legend from pantranscriptome_pie_chart
legend <- cowplot::get_legend(pantranscriptome_pie_chart + theme(legend.position="right"))

# Create the main plot layout
main_plot <- plot_grid(
  # Top row: Two pie charts
  plot_grid(pantranscriptome_pie_chart + theme(legend.position="none"), 
            expressed_genes_pie_chart + theme(legend.position="none"), 
            ncol = 2, labels = c('A', 'B'), label_x=0.1),
  
  # Bottom row: Centering the violin plot using blank space
  plot_grid(NULL, violin_plot + theme(legend.position="none"), NULL, 
            ncol = 3, rel_widths = c(1, 2, 1), labels = c('', 'C', ''), label_x=-0.15),
  nrow = 2
)

# Combine main plot with the legend
final_plot <- plot_grid(main_plot, legend, ncol = 2, rel_widths = c(10, 1))


png(file="final_processing_and_plotting/pantranscrptome_stats.png", width=20, height=15, units="cm", res=300)

final_plot

dev.off()
