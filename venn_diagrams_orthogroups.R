rm(list=ls())
ls() 

library(ggvenn)
library(tidyr)
library(ggpubr)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# Treatment
# =============================================================================

rio_treat = read.csv("Rio_reference_genome/treat_genes_orthogroups.csv")
rio_treat = rio_treat[c("orthogroup")]
rio_treat = separate_rows(rio_treat, orthogroup, sep=', ')
rio_treat_lst = rio_treat$orthogroup %>% unique()

pan_treat = read.csv("pantranscriptome/treat_genes_orthogroups.csv")
pan_treat = pan_treat[c("orthogroup")]
pan_treat = separate_rows(pan_treat, orthogroup, sep=', ')
pan_treat_lst = pan_treat$orthogroup %>% unique()

treat = list(rio = rio_treat_lst, pan_t = pan_treat_lst)

treat_venn = ggvenn(
              treat, columns = c("rio", "pan_t"),
              stroke_linetype = 'blank',
              show_percentage = F,
              set_name_color = 'white',
              fill_color = c("#CC79A7", "#0072B2"),
              text_size = 6,
              auto_scale = T)
treat_venn

# =============================================================================
# Type
# =============================================================================

rio_type = read.csv("Rio_reference_genome/type_genes_orthogroups.csv")
rio_type = rio_type[c("orthogroup")]
rio_type = separate_rows(rio_type, orthogroup, sep=', ')
rio_type_lst = rio_type$orthogroup %>% unique()

pan_type = read.csv("pantranscriptome/type_genes_orthogroups.csv")
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
              text_size = 6,
              auto_scale = T)
type_venn

# =============================================================================
# Interaction
# =============================================================================

rio_int = read.csv("Rio_reference_genome/int_genes_orthogroups.csv")
rio_int = rio_int[c("orthogroup")]
rio_int = separate_rows(rio_int, orthogroup, sep=', ')
rio_int_lst = rio_int$orthogroup %>% unique()

pan_int = read.csv("pantranscriptome/int_genes_orthogroups.csv")
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
            text_size = 6,
            auto_scale = T)
int_venn

# =============================================================================
# All plots in one
# =============================================================================

png(file="DEGs_orthogroup_venn.png", width=12, height=12, units="in", res=300)

blank <- ggplot() + theme_void()

ggarrange(
  ggarrange(treat_venn, type_venn, blank, nrow = 1,
            widths = c(1, 1, 0), labels = c('A','B')),
  ggarrange(blank, int_venn, blank, nrow = 1,
            widths = c(0.5, 1, 0.5), labels = c('C'), label.x=1),
  nrow = 2)

dev.off()
