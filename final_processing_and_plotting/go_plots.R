rm(list=ls())
ls()

library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/")

# =============================================================================
# read in data
# =============================================================================

rio_treat_df = read.table("Rio_reference_genome/GO_analysis/go_figure/treat/molecular_function_result_treat.tsv", header=TRUE, sep="\t")

rio_type_df = read.table("Rio_reference_genome/GO_analysis/go_figure/type/molecular_function_result_type.tsv", header=TRUE, sep="\t")

rio_int_df = read.table("Rio_reference_genome/GO_analysis/go_figure/int/molecular_function_result_int.tsv", header=TRUE, sep="\t")

pant_treat_df = read.table("pantranscriptome/GO_analysis/go_figure/treat/molecular_function_result_treat.tsv", header=TRUE, sep="\t")

pant_type_df = read.table("pantranscriptome/GO_analysis/go_figure/type/molecular_function_result_type.tsv", header=TRUE, sep="\t")

pant_int_df = read.table("pantranscriptome/GO_analysis/go_figure/int/molecular_function_result_int.tsv", header=TRUE, sep="\t")

# =============================================================================
# set minimum and maximum values for x and y-axis
# =============================================================================

df_list = list(rio_treat_df, rio_type_df, rio_int_df,
            pant_treat_df, pant_type_df, pant_int_df)

xmax_limit = max(sapply(df_list, function(i) max(i$x)))+.5
xmin_limit = min(sapply(df_list, function(i) min(i$x)))-.5
ymax_limit = max(sapply(df_list, function(i) max(i$y)))+.5
ymin_limit = min(sapply(df_list, function(i) min(i$y)))-.5

# =============================================================================
# set size limits for the points
# =============================================================================

all_sizes <- c(rio_treat_df$colour, rio_type_df$colour, rio_int_df$colour, pant_treat_df$colour, pant_type_df$colour, pant_int_df$colour)

common_max = -(min(all_sizes)) 
common_min = -(max(all_sizes))

# =============================================================================
# compare Rio reference genome and the pantranscriptome GO terms and create a
# column from this for each of the dataframes
# =============================================================================

compare_go_terms = function(rio_df, pant_df){
  for(i in 1:nrow(rio_df)){
    rio_term = rio_df[i, 'representative']
    if(rio_term %in% pant_df$representative){
      rio_df[i, 'comparison'] = 'not_unique'
    } else {
      rio_df[i, 'comparison'] = 'unique'
    }
  }
  for(i in 1:nrow(pant_df)){
    pant_term = pant_df[i, 'representative']
    if(pant_term %in% rio_df$representative){
      pant_df[i, 'comparison'] = 'not_unique'
    } else {
      pant_df[i, 'comparison'] = 'unique'
    }
  }
  return(list(rio_df = rio_df, pant_df = pant_df))
}

rio_int_df = compare_go_terms(rio_int_df, pant_int_df)$rio_df
pant_int_df = compare_go_terms(rio_int_df, pant_int_df)$pant_df

rio_type_df = compare_go_terms(rio_type_df, pant_type_df)$rio_df
pant_type_df = compare_go_terms(rio_type_df, pant_type_df)$pant_df

rio_treat_df = compare_go_terms(rio_treat_df, pant_treat_df)$rio_df
pant_treat_df = compare_go_terms(rio_treat_df, pant_treat_df)$pant_df

# =============================================================================
# function for plotting
# =============================================================================

get_plot = function(df, source){
  names(df)[names(df) == "colour"] = "log_pval"
  df = df %>% select(-size)
  
  df = df %>%
    mutate(members_count = str_count(members, "GO:"))
  
  if(source=="rio"){
    color="#B9638A"
    } else if(source=="pantranscriptome"){
    color="#0072B2"
  }

  plot = ggplot(df, aes(x, y, color=comparison)) +
    geom_point(aes(size=-(log_pval)), alpha=0.5) +
    geom_text_repel(data=head(df,15), aes(label=description), size=3.5, color='black') +
    scale_color_manual(values = c('unique'=color, 'not_unique'='#D3D3D3')) + 
    xlab('') +
    ylab('') +
    xlim(xmin_limit,xmax_limit) +
    ylim(ymin_limit,ymax_limit) +  
    scale_size(range=c(1,30),
               limits=c(common_min,common_max),
               guide="none") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill="white",color=NA),
      panel.border = element_rect(color="black",fill=NA))
  return(plot)
}

# =============================================================================
# plotting
# =============================================================================

rio_treat_plot = get_plot(rio_treat_df, 'rio')
rio_type_plot = get_plot(rio_type_df, 'rio')
rio_int_plot = get_plot(rio_int_df, 'rio')
pant_treat_plot = get_plot(pant_treat_df, 'pantranscriptome')
pant_type_plot = get_plot(pant_type_df, 'pantranscriptome')
pant_int_plot = get_plot(pant_int_df, 'pantranscriptome')

# =============================================================================
# saving plots
# =============================================================================

png(file="final_processing_and_plotting/terms_treatment.png", width=16, height=10, units="in", res=300)
ggpubr::ggarrange(rio_treat_plot, pant_treat_plot,
                  ncol=2, nrow=1,
                  labels = "AUTO")
dev.off()

png(file="final_processing_and_plotting/terms_type.png", width=16, height=10, units="in", res=300)
ggpubr::ggarrange(rio_type_plot, pant_type_plot,
                  ncol=2, nrow=1,
                  labels = "AUTO")
dev.off()

png(file="final_processing_and_plotting/terms_interaction.png", width=16, height=10, units="in", res=300)
ggpubr::ggarrange(rio_int_plot, pant_int_plot,
                  ncol=2, nrow=1,
                  labels = "AUTO")
dev.off()

png(file="final_processing_and_plotting/terms_all.png", width=16, height=10, units="in", res=300)
ggpubr::ggarrange(rio_treat_plot, pant_treat_plot,
                  rio_type_plot, pant_type_plot,
                  rio_int_plot, pant_int_plot,
                  ncol=2, nrow=3,
                  labels = "AUTO")
dev.off()

png(file="final_processing_and_plotting/terms_treatment_n_type.png", width=16, height=16, units="in", res=300)
ggpubr::ggarrange(rio_treat_plot, pant_treat_plot,
                  rio_type_plot, pant_type_plot,
                  ncol=2, nrow=2,
                  labels = "AUTO")
dev.off()
