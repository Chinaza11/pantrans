rm(list=ls())
ls() 

library(ggvenn)
library(tidyr)
library(ggpubr)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# Type
# =============================================================================

rio_type = read.csv("Rio_reference_genome/permutation/type_genes_orthogroups.csv")
rio_type = rio_type[c("orthogroup")]
rio_type = separate_rows(rio_type, orthogroup, sep=', ')
rio_type_lst = rio_type$orthogroup %>% unique()

pan_type = read.csv("pantranscriptome/permutation/type_genes_orthogroups.csv")
pan_type = pan_type[c("orthogroup")]
pan_type = separate_rows(pan_type, orthogroup, sep=', ')
pan_type_lst = pan_type$orthogroup %>% unique()

type = list(rio = rio_type_lst, pan_t = pan_type_lst)

type_venn = ggvenn(
              type, columns = c("rio", "pan_t"),
              stroke_linetype = 'blank',
              show_percentage = F,
              set_name_color = 'white',
              fill_color = c("#CC79A7", "#0072B2"),
              text_size = 4,
              auto_scale = T)
type_venn

# =============================================================================
# Interaction
# =============================================================================

rio_int = read.csv("Rio_reference_genome/permutation/int_genes_orthogroups.csv")
rio_int = rio_int[c("orthogroup")]
rio_int = separate_rows(rio_int, orthogroup, sep=', ')
rio_int_lst = rio_int$orthogroup %>% unique()

pan_int = read.csv("pantranscriptome/permutation/int_genes_orthogroups.csv")
pan_int = pan_int[c("orthogroup")]
pan_int = separate_rows(pan_int, orthogroup, sep=', ')
pan_int_lst = pan_int$orthogroup %>% unique()

int = list(rio = rio_int_lst, pan_t = pan_int_lst)

int_venn = ggvenn(
            int, columns = c("rio", "pan_t"),
            stroke_linetype = 'blank',
            show_percentage = F,
            set_name_color = 'white',
            fill_color = c("#CC79A7", "#0072B2"),
            text_size = 4,
            auto_scale = T)
int_venn

# =============================================================================
# All plots in one
# =============================================================================

png(file="permutation/DEGs_orthogroup_venn.jpg", width=20, height=20, units="cm", res=300)

ggarrange(type_venn, int_venn, nrow = 2, labels = 'AUTO')

dev.off()
