rm(list=ls())
ls()

# =============================================================================
# Function
# =============================================================================

merge_and_save = function(top_go_path, go_figure_path, merge_path){
  top_go = read.table(top_go_path, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  colnames(top_go)[colnames(top_go) == "classicFisher"] = "Classic.Fisher.Pvalue"
  top_go = top_go[, !(names(top_go) %in% "Term")]
  
  go_figure = read.table(go_figure_path, header=TRUE, sep="\t")
  go_figure = go_figure[, !(names(go_figure) %in% c("Member.P.value", "Member.user.defined.value"))]
  
  merge = merge(top_go, go_figure, by.x='GO.ID', by.y='Cluster.member')
  merge = merge[order(merge$Cluster.representative), ]
  colnames(merge)[colnames(merge) == "Cluster.representative"] = "Cluster.Representative"
  colnames(merge)[colnames(merge) == "Member.IC"] = "Information.Content"
  colnames(merge)[colnames(merge) == "Member.frequency"] = "Frequency"
  colnames(merge)[colnames(merge) == "Cluster.member.description"] = "Term.Description"
  
  merge = merge[, c("GO.ID", "Term.Description", "Annotated", "Significant", "Expected", "Classic.Fisher.Pvalue", "Cluster.Representative", "Information.Content", "Frequency")]
  
  write.table(merge, merge_path, quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)
}

# =============================================================================
# Rio ref genome
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/GO_analysis/")

# Treatment
merge_and_save("res.treat.mf.txt",
               "go_figure/treat/molecular_function_full_table_result_treat.tsv",
               "topGO_n_GOFigure_treat.txt")

# Type
merge_and_save("res.type.mf.txt",
               "go_figure/type/molecular_function_full_table_result_type.tsv",
               "topGO_n_GOFigure_type.txt")

# Interaction
merge_and_save("res.int.mf.txt",
               "go_figure/int/molecular_function_full_table_result_int.tsv",
               "topGO_n_GOFigure_int.txt")

# =============================================================================
# Pan-transcriptome
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/GO_analysis/")

# Treatment
merge_and_save("res.treat.mf.txt",
               "go_figure/treat/molecular_function_full_table_result_treat.tsv",
               "topGO_n_GOFigure_treat.txt")

# Type
merge_and_save("res.type.mf.txt",
               "go_figure/type/molecular_function_full_table_result_type.tsv",
               "topGO_n_GOFigure_type.txt")

# Interaction
merge_and_save("res.int.mf.txt",
               "go_figure/int/molecular_function_full_table_result_int.tsv",
               "topGO_n_GOFigure_int.txt")

# =============================================================================
# Rio ref genome - permutation
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/GO_analysis/permutation/")

# Type
merge_and_save("res.type.mf.txt",
               "go_figure/permutation_type/molecular_function_full_table_result_permutation_type.tsv",
               "topGO_n_GOFigure_type.txt")

# Interaction
merge_and_save("res.int.mf.txt",
               "go_figure/permutation_int/molecular_function_full_table_result_permutation_int.tsv",
               "topGO_n_GOFigure_int.txt")

# =============================================================================
# Pan-transcriptome - permutation
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/GO_analysis/permutation/")

# Type
merge_and_save("res.type.mf.txt",
               "go_figure/permutation_type/molecular_function_full_table_result_permutation_type.tsv",
               "topGO_n_GOFigure_type.txt")

# Interaction
merge_and_save("res.int.mf.txt",
               "go_figure/permutation_int/molecular_function_full_table_result_permutation_int.tsv",
               "topGO_n_GOFigure_int.txt")
