#getwd() #------ in HPC environ: uncomment

#.libPaths("~/rlibs") #------ in HPC environ: uncomment

library(dplyr)
library(stringr)
library(tibble)

rm(list=ls()) #------ in HPC environ: comment out
ls() #------ in HPC environ: comment out

# =============================================================================
# Set up the working directory, and read in raw counts from mmquant
# The format is rows=genes; columns=samples
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/glmmseq") #------ in HPC environ: comment out
cts = as.matrix(read.table("rawCounts20.txt", header=TRUE, row.names="Gene", sep="\t"))
colnames(cts) = gsub("X", "pi", colnames(cts))
colnames(cts) = gsub("pi069.n02", "pi069.no2", colnames(cts))

# =============================================================================
# Set up the sample information table 
# (pull genotype, treatment, and type info from the sample names)
# =============================================================================

SampleID = colnames(cts)
Genotype = unlist(lapply(strsplit(SampleID, split="\\."), `[[`, 1)) 
Condition = unlist(lapply(strsplit(SampleID, split="\\."), `[[`, 2)) 
Rep = substr(Condition,nchar(Condition),nchar(Condition))
Condition = gsub("0","o", Condition)
Condition = substr(Condition, 1, nchar(Condition)-1) 
Condition = as.factor(gsub("con", "C", gsub("no", "N", Condition)))
GenID = paste0(Genotype, Rep)
Genotype[-c(grep("[0-9]",Genotype))] = "Sweet"
Genotype[grep("[0-9]",Genotype)] = "Biomass" 

# =============================================================================
# Now put the counts and sample info into a DESeq object
# This requires making a metadata table (colData) and specifying the model 
# (in this case counts as a function of iron treatment and sweet/non-sweet type)
# =============================================================================

library(DESeq2)

colData = data.frame(SAMID = as.factor(SampleID),
                     GENID = as.factor(GenID),
                     Treatment = Condition,
                     Type = as.factor(Genotype))

rownames(colData) = SampleID

print(colData)

# =============================================================================
# This section is for preparing the permutation groups
# =============================================================================

library(gtools)
library(tidyr)

set.seed(1234)

Type = c("pi069", "pi972", "Gra", "Leo", "Rio")
combinations1 = data.frame(combinations(5,3,Type))
# combinations2 = data.frame(combinations(5,2,Type)) #this can't be merged with combination1

new_column_names = c("Biomass1", "Biomass2", "Biomass3")
colnames(combinations1) = new_column_names

combinations1$Sweet1 = NA
combinations1$Sweet2 = NA

# function to add individuals not already in the permutation groups
complete_df = function(df,column){
  for (i in 1:nrow(df)){
    for (str in Type){
      if (str %in% df[i,]) {
      } else {
        df[i, column] = str
      }}}
  return(df)
}

combinations1 = complete_df(combinations1,'Sweet1')
combinations1 = complete_df(combinations1,'Sweet2')

# duplicating the genotypes in the cell
for (col_name in colnames(combinations1)) {
  combinations1[[col_name]] <- ifelse(combinations1[[col_name]] %in% "Gra", 
      paste(combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], sep = ","), 
      paste(combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], combinations1[[col_name]], sep = ","))
}

# function to process each cell in a column: labeling each duplicated individual
process_column <- function(column) {
  modified <- sapply(column, function(cell) {
    parts <- unlist(strsplit(cell, ","))
    modified_parts <- sapply(seq_along(parts), function(i) {
      part <- parts[i]
      if (length(parts) == 6) {
        if (i %% 6 == 1) {
          return(paste0(part, "1"))
        } else if (i %% 6 == 2) {
          return(paste0(part, "2"))
        } else if (i %% 6 == 3) {
          return(paste0(part, "3"))
        } else if (i %% 6 == 4) {
          return(paste0(part, "4"))
        } else if (i %% 6 == 5) {
          return(paste0(part, "5"))
        } else {
          return(paste0(part, "6"))
        }
      } else if (length(parts) == 4) {
        if (i %% 4 == 1) {
          return(paste0(part, "1"))
        } else if (i %% 4 == 2) {
          return(paste0(part, "2"))
        } else if (i %% 4 == 3) {
          return(paste0(part, "3"))
        } else {
          return(paste0(part, "4"))
        }
      } else {
        # Handle other cases as needed
        return(part)
      }
    })
    return(paste(modified_parts, collapse = ", "))  # rejoin the modified parts
  })
  return(modified)
}

# applying the process_column function to each column in the dataframe
combinations1 <- combinations1 %>%
  mutate_all(process_column)

# combining the biomass columns into 1 and the sweet column into 1
combined_column1 = do.call(paste, c(combinations1[1:3], sep = ","))
combined_column2 = do.call(paste, c(combinations1[4:5], sep = ","))
combined_column = data.frame(Biomass=combined_column1, Sweet = combined_column2)
combined_column

# transposing dataframe
combn_t = t(combinations1)
combn_t

# adding column names
colnames_list = paste0("permutation", 1:10)
colnames(combn_t) = colnames_list

# making a dataframe from each of the columns from combn_t and putting this in a list for further manipulation
new_dfs <- list() # Create an empty list to store the new dataframes
for (col_name in colnames(combn_t)) {
  new_df <- data.frame(RowName = row.names(combn_t), Value = combn_t[, col_name])
  new_dfs[[length(new_dfs) + 1]] <- new_df
}
new_dfs

# splitting the strings in each cell into rows
new_dfs2 <- list()
for (i in 1:(length(new_dfs))){
  split_df = new_dfs[[i]] %>%
    separate_rows(Value, sep = ",")
  new_dfs2[[i]] <- split_df
}
print(new_dfs2[[3]], n=22)

# dropping the "Rio3" and "Rio4" rows
name_to_drop <- c(" Rio3", " Rio4")
filtered_dataframes <- list() # create a list to store the filtered dataframes

# loop through each dataframe in the list
for (i in 1:length(new_dfs2)) {
  # filter out rows with the specified name
  filtered_df <- new_dfs2[[i]] %>%
    filter(Value != name_to_drop)
  
  # append the filtered dataframe to the list
  filtered_dataframes[[i]] <- filtered_df
}
print(filtered_dataframes[[3]], n=22)

# removing extra spaces from the "Value" column
column_to_clean <- "Value"
cleaned_dataframes <- list() # create a list to store the cleaned dataframes
# loop through each dataframe in the list
for (i in 1:length(filtered_dataframes)) {
  # remove extra spaces in the specified column using gsub
  cleaned_df <- filtered_dataframes[[i]]
  cleaned_df[[column_to_clean]] <- gsub("\\s+", "", cleaned_df[[column_to_clean]])
  
  # append the cleaned dataframe to the list
  cleaned_dataframes[[i]] <- cleaned_df
}
print(cleaned_dataframes[[3]], n=22)


# renaming Biomass1, Biomass2 and Biomass3 to Biomass. Renaming Sweet1 and Sweet2 to Sweet
column_to_edit <- "RowName"
cleaned_dataframes2 <- list()
for (i in 1:length(cleaned_dataframes)) {
  cleaned_df2 <- cleaned_dataframes[[i]]
  cleaned_df2[[column_to_edit]] <- gsub("Biomass[0-9]", "Biomass", cleaned_df2[[column_to_edit]])
  cleaned_df2[[column_to_edit]] <- gsub("Sweet[0-9]", "Sweet", cleaned_df2[[column_to_edit]])
  
  # append the cleaned dataframe to the list
  cleaned_dataframes2[[i]] <- cleaned_df2
}
print(cleaned_dataframes2[[3]], n=22)


# renaming genotype name in preparation for merger
# define a mapping of strings to find and their replacements
find_replace_map <- list(
  list(find = "Gra1", replace = "Gra.con1"),
  list(find = "Gra2", replace = "Gra.con2"),
  list(find = "Gra3", replace = "Gra.con3"),
  list(find = "Gra4", replace = "Gra.no1"),
  list(find = "Gra5", replace = "Gra.no2"),
  list(find = "Gra6", replace = "Gra.no3"),
  list(find = "Leo1", replace = "Leo.con1"),
  list(find = "Leo2", replace = "Leo.con2"),
  list(find = "Leo3", replace = "Leo.no1"),
  list(find = "Leo4", replace = "Leo.no2"),
  list(find = "Rio1", replace = "Rio.con1"),
  list(find = "Rio2", replace = "Rio.no1"),
  list(find = "pi9721", replace = "pi972.con1"),
  list(find = "pi9722", replace = "pi972.con2"),
  list(find = "pi9723", replace = "pi972.no1"),
  list(find = "pi9724", replace = "pi972.no2"),
  list(find = "pi0691", replace = "pi069.con1"),
  list(find = "pi0692", replace = "pi069.con2"),
  list(find = "pi0693", replace = "pi069.no1"),
  list(find = "pi0694", replace = "pi069.no2")
)

# loop through each dataframe in the list
for (i in 1:length(cleaned_dataframes2)) {
  # loop through the find-replace mapping
  for (mapping in find_replace_map) {
    find_str <- mapping$find
    replace_str <- mapping$replace
    
    # find and replace the string in the "Name" column
    cleaned_dataframes2[[i]]$Value <- str_replace(cleaned_dataframes2[[i]]$Value, find_str, replace_str)
  }
}
print(cleaned_dataframes2[[3]], n=22)


# renaming column names
new_column_names <- c("Type", "SAMID") # Define new column names
dataframes_with_new_names <- list() # Create a list to store the dataframes with renamed columns

# Loop through each dataframe in the list
for (i in 1:length(cleaned_dataframes2)) {
  # Rename columns using the colnames() function
  colnames(cleaned_dataframes2[[i]]) <- new_column_names
  
  # Append the dataframe with renamed columns to the list
  dataframes_with_new_names[[i]] <- cleaned_dataframes2[[i]]
}
print(dataframes_with_new_names[[3]], n=22)

#duplicating colData(original meta_data) and renaming the "Type" column so it can be easily identified after merger 
colData_for_merge = colData
colnames(colData_for_merge)[4] = "Type_original_grouping"
colData_for_merge
str(colData_for_merge)


#merging colData_for_merge with individual permutation groups in the dataframes in the list 
# Loop through each dataframe in the list and merge with the additional dataframe
for (i in 1:length(dataframes_with_new_names)) {
  merged_df <- merge(x=dataframes_with_new_names[[i]], y=colData_for_merge, by="SAMID")
  # Replace the dataframe in the list with the merged dataframe
  dataframes_with_new_names[[i]] <- merged_df
}
print(dataframes_with_new_names[[3]])


rownames_column <- "SAMID" # Define the column to be used as rownames
# Loop through each dataframe in the list and set the rownames
for (i in 1:length(dataframes_with_new_names)) {
  df <- dataframes_with_new_names[[i]]
  rownames(df) <- df[, rownames_column]
  dataframes_with_new_names[[i]] <- df
}
print(dataframes_with_new_names[[3]])
str(dataframes_with_new_names[[3]])


columns_to_factor <- c("SAMID", "Type") # Define the column to be converted to a factor
# Loop through each dataframe in the list and convert the column to a factor
for (i in 1:length(dataframes_with_new_names)) {
  df <- dataframes_with_new_names[[i]]
  # Convert the specified columns to factors
  for (col_name in columns_to_factor) {
  df[[col_name]] <- as.factor(df[[col_name]])
  }
  dataframes_with_new_names[[i]] <- df
}
print(dataframes_with_new_names[[3]])
str(dataframes_with_new_names[[3]])


column_to_drop <- "Type_original_grouping" # Define the column to be dropped
final_df <- list()
# Loop through each dataframe in the list and drop the column
for (i in 1:length(dataframes_with_new_names)) {
  df <- dataframes_with_new_names[[i]]
  
  # Drop the specified column
  df <- df[, !names(df) %in% column_to_drop]
  
  final_df[[i]] <- df
}

print(final_df[[3]])
str(final_df[[3]])

# =============================================================================
# creating the other set of permutation groups by duplicating previous group 
# and swiping Type category
# =============================================================================

final_df2 = final_df 
print(final_df2[[3]])

str(final_df2[[3]])

# Define a mapping of strings to find and their replacements
find_replace <- list(
  list(find = "Biomass", replace = "Biomass_to_be_replaced"),
  list(find = "Sweet", replace = "Sweet_to_be_replaced"))

# Loop through each dataframe in the list
for (i in 1:length(final_df2)) {
  # Loop through the find-replace mapping
  for (mapping in find_replace) {
    find_str <- mapping$find
    replace_str <- mapping$replace
    
    # Find and replace the string in the "Type" column
    final_df2[[i]]$Type <- str_replace(final_df2[[i]]$Type, find_str, replace_str)
  }
}
print(final_df2[[3]])


# Define a mapping of strings to find and their replacements
find_replace <- list(
  list(find = "Biomass_to_be_replaced", replace = "Sweet"),
  list(find = "Sweet_to_be_replaced", replace = "Biomass"))

# Loop through each dataframe in the list
for (i in 1:length(final_df2)) {
  # Loop through the find-replace mapping
  for (mapping in find_replace) {
    find_str <- mapping$find
    replace_str <- mapping$replace
    
    # Find and replace the string in the "Type" column
    final_df2[[i]]$Type <- str_replace(final_df2[[i]]$Type, find_str, replace_str)
  }
}
print(final_df2[[3]])
str(final_df2[[3]])


columns_to_factor <- c("Type") # Define the column to be converted to a factor
# Loop through each dataframe in the list and convert the column to a factor
for (i in 1:length(final_df2)) {
  df <- final_df2[[i]]
  
  # Convert the specified columns to factors
  for (col_name in columns_to_factor) {
    df[[col_name]] <- as.factor(df[[col_name]])
  }
  final_df2[[i]] <- df
}
str(final_df2[[3]])

final_df <- c(final_df, final_df2)

str(final_df[[3]])
str(final_df[[12]])

# =============================================================================
# looping through each of the dataframes in the list and running deseq2 analysis and extracting the normcts and dispersion data
# =============================================================================

# normalized_counts_list <- list()
# dispersions_list <- list()
# 
# 
# for (i in 1:length(final_df)) {
#   df = final_df[[i]]
#   cts = cts[, rownames(final_df[[i]])] #sorting the colnames of cts in the same order as the rownames of coldata. This needs to be done or later functions will return an error message
# 
#   dds_name <- paste0("dds", i)  # Define the variable name for the DESeqDataSet
#   dds <- DESeqDataSetFromMatrix(countData = cts, colData = df, design = ~ Treatment + Type)
# 
#   # Estimate size factors for the DESeqDataSet
#   dds <- estimateSizeFactors(dds)
# 
#   dds <- dds[rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3,]
#   norm_counts <- counts(dds, normalized=TRUE)
#   dds <- DESeq(dds)
#   dispersions <- setNames(dispersions(dds), rownames(norm_counts))
# 
#   # Save the DESeqDataSet as a variable with a unique name (e.g., dds1, dds2, dds3)
#   assign(dds_name, dds, envir = .GlobalEnv)
# 
#   # Store the normalized counts and dispersions in their respective lists
#   normalized_counts_list[[i]] <- norm_counts
#   dispersions_list[[i]] <- setNames(dispersions, rownames(norm_counts))
# }
# 
# save(normalized_counts_list, file = "permutation_analysis/normalized_counts_list.Rdata")
# save(dispersions_list, file = "permutation_analysis/dispersions_list.Rdata")

load(file = "permutation_analysis/normalized_counts_list.Rdata")
load(file = "permutation_analysis/dispersions_list.Rdata")

# =============================================================================
# validating for-loop above
# =============================================================================

# #------------------------------------------------------------------------------
# #permutation #1 (final_df[[1]])

# cts <- cts[, rownames(final_df[[1]])]
# all(rownames(final_df[[1]]) == colnames(cts))
# 
# dds1 <- DESeqDataSetFromMatrix(countData = cts, colData = final_df[[1]], design = ~ Treatment + Type)
# 
# dds1 <- estimateSizeFactors(dds1)
# dds1 <- dds1[rowSums(counts(dds1, normalized=TRUE) >= 5 ) >= 3,]
# normcts1 = counts(dds1, normalized=TRUE)
# dds1 <- DESeq(dds1)
# dispersions1 <- setNames(dispersions(dds1), rownames(normcts1))
# head(dispersions1, n=10)
# 
# head(dispersions_list[[1]], n=10)
# 
# identical(normcts1, normalized_counts_list[[1]]) #they are same
# identical(dispersions1, dispersions_list[[1]]) #they seem to be different but visual examination of head(dispersions1) and head(dispersions_list[[1]]) shows they might be similar
# 
# 
# #------------------------------------------------------------------------------
# #permutation #6 (final_df[[6]])

# cts <- cts[, rownames(final_df[[6]])]
# all(rownames(final_df[[6]]) == colnames(cts))
# 
# dds6 <- DESeqDataSetFromMatrix(countData = cts, colData = final_df[[6]], design = ~ Treatment + Type)
# 
# dds6 <- estimateSizeFactors(dds6)
# dds6 <- dds6[rowSums(counts(dds6, normalized=TRUE) >= 5 ) >= 3,]
# normcts6 = counts(dds6, normalized=TRUE)
# dds6 <- DESeq(dds6)
# dispersions6 <- setNames(dispersions(dds6), rownames(normcts6))
# 
# identical(normcts6, normalized_counts_list[[6]]) #they are same
# identical(dispersions6, dispersions_list[[6]]) #they seem to be different but visual examination of head(dispersions6) and head(dispersions_list[[6]]) shows they might be similar
# 
# 
# identical(normcts1, normcts6) #the normcts for all the different combinations is same
# identical(dispersions1, dispersions6) #the dispersions for all the different combinations is different
# 
# #------------------------------------------------------------------------------
# #permutation #20 (final_df[[20]])

# cts <- cts[, rownames(final_df[[20]])]
# all(rownames(final_df[[20]]) == colnames(cts))
# 
# dds20 <- DESeqDataSetFromMatrix(countData = cts, colData = final_df[[20]], design = ~ Treatment + Type)
# 
# dds20 <- estimateSizeFactors(dds20)
# dds20 <- dds20[rowSums(counts(dds20, normalized=TRUE) >= 5 ) >= 3,]
# normcts20 = counts(dds20, normalized=TRUE)
# dds20 <- DESeq(dds20)
# dispersions20 <- setNames(dispersions(dds20), rownames(normcts20))
# 
# identical(normcts20, normalized_counts_list[[20]]) #they are same
# identical(dispersions20, dispersions_list[[20]]) #they seem to be different but visual examination of head(dispersions20) and head(dispersions_list[[20]]) shows they might be similar
# 

# =============================================================================
# deseq2 on the original dataset, this will be compared to the default dispersion 
# numbers coming from code below and the for-loop
# =============================================================================

# cts <- cts[, rownames(colData)] 
# all(rownames(colData) == colnames(cts))
# 
# dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ Treatment + Type)
# 
# dds <- estimateSizeFactors(dds)  
# idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
# dds <- dds[idx,] #this filter ensures at least 3 samples with a count of 5 or more
# normcts = counts(dds, normalized=TRUE)
# dds <- DESeq(dds)
# dispersions <- setNames(dispersions(dds), rownames(normcts))
# 
# #------------------------------------------------------------------------------
# #validating for-loop above: permutation #13 (the default group) (final_df[[13]])

# cts <- cts[, rownames(final_df[[13]])]
# all(rownames(final_df[[13]]) == colnames(cts))
# 
# dds13 <- DESeqDataSetFromMatrix(countData = cts, colData = final_df[[13]], design = ~ Treatment + Type)
# 
# dds13 <- estimateSizeFactors(dds13)
# dds13 <- dds13[rowSums(counts(dds13, normalized=TRUE) >= 5 ) >= 3,]
# normcts13 = counts(dds13, normalized=TRUE)
# dds13 <- DESeq(dds13)
# dispersions13 <- setNames(dispersions(dds13), rownames(normcts13))
# 
# identical(dispersions, dispersions13) #they seem to be different but visual examination of head(dispersions) and head(dispersions_list[[13]]) shows they might be similar
# head(dispersions, n=10)
# head(dispersions13, n=10)
# 
# identical(normcts13, normalized_counts_list[[13]]) #they are same
# identical(dispersions13, dispersions_list[[13]]) #they seem to be different but visual examination of head(dispersions13) and head(dispersions_list[[13]]) shows they might be similar

# =============================================================================
# running glmmseq on the different permutation groups
# =============================================================================

library(glmmSeq)
set.seed(1234)

# results_list <- list()
# 
# # Loop through the list of dataframes
# for (i in 1:length(dispersions_list)) {
#   # Extract the current dataframe from the list
#   metadata = final_df[[i]]
#   countdata = normalized_counts_list[[i]]
#   dispersions <- dispersions_list[[i]]
#   
#   # Perform glmmSeq for the current dataframe
#   current_result <- glmmSeq(~ Treatment * Type + (1 | GENID),
#                             id = "GENID",
#                             countdata = countdata,
#                             metadata = metadata,
#                             dispersion = dispersions,
#                             removeSingles = FALSE,
#                             cores = 4, #change to 20 cores for the cluster
#                             progress = TRUE)
#   
#   # Store the result in the results_list
#   results_list[[i]] <- current_result
# }
# 
# save(results_list, file = "permutation_analysis/results_list.Rdata")
load(file = "permutation_analysis/results_list.Rdata")

results_list[[1]]@modelData

results13 = glmmQvals(results_list[[13]])
results20 = glmmQvals(results_list[[20]])

# getting qvalue for all the permutation groups and putting it in a list
result_qval_summary = list()
for (i in 1:length(results_list)) {
  df = results_list[[i]]
  
  result_qval = glmmQvals(df)
  result_qval_summary[[i]] = result_qval
}

# getting the stats for all the permutation groups and putting it in a list
stats_list = list()
for (i in 1:length(result_qval_summary)) {
  df = result_qval_summary[[i]]
  
  stats = data.frame(df@stats)
  stats_list[[i]] = stats
}

colnames(stats_list[[2]])

# =============================================================================
# loading the result from the original analysis (default grouping), the q_vals 
# of type and treatment:type interaction will be extracted for comparison with 
# the result from the permutation analysis
# =============================================================================

load(file = "results.Rdata")

results@modelData

# Get Q values and a significance summary table with correct P-values
results <- glmmQvals(results)
stats = data.frame(results@stats)

stats_default = stats

#------------------------------------------------------------------------------
# extracting the type q_vals of the default grouping for comparison with the 
# permutation analysis result
#------------------------------------------------------------------------------

q_vals_type = data.frame(stats_default$qvals.Type)
colnames(q_vals_type) = "default_qvals.type"
rownames(q_vals_type) = rownames(stats_default)

#------------------------------------------------------------------------------
# renaming qvals.Type in the permutation result list of dataframes to qvals.Type_permutation(i)

colnames(stats_list[[2]])

for (i in 1:length(stats_list)) {
  # Get the current dataframe
  current_df <- stats_list[[i]]
  
  # Change the column name
  colnames(current_df)[colnames(current_df) == "qvals.Type"] <- paste0("qvals.Type_permutation_", i)
  
  # Update the dataframe in the list
  stats_list[[i]] <- current_df
}

colnames(stats_list[[2]])

#------------------------------------------------------------------------------
# merging q_vals.Type result from the permutations to q_vals_type dataframe
q_vals_type <- q_vals_type %>% rownames_to_column(var = "Gene") # Convert row names to a column in 'q_vals_type'

# Loop through the list of permutation(stats_list) and convert row names to a column
for (i in 1:length(stats_list)) {
  stats_list[[i]] <- stats_list[[i]] %>% rownames_to_column(var = "Gene")
}  
colnames(stats_list[[2]])

# Loop through the list of dataframes and merge the "qvals.Type_permutation_(i)" column to 'q_vals_type' by row names
for (i in 1:length(stats_list)) {
  q_vals_type <- merge(q_vals_type, stats_list[[i]][, c("Gene", paste0("qvals.Type_permutation_", i))], by.x = "Gene", by.y = "Gene", all.x = TRUE)
}

#------------------------------------------------------------------------------
#comparing result from default grouping ran normally to result from default grouping ran alongside other permutation tests

colnames(q_vals_type)

subset_q_vals_type <- q_vals_type[, c("Gene", "default_qvals.type", "qvals.Type_permutation_13")]

subset_q_vals_type$comparison = ifelse(subset_q_vals_type$default_qvals.type <= 0.05 & subset_q_vals_type$qvals.Type_permutation_13 <= 0.05, "both siginificant", ifelse(subset_q_vals_type$default_qvals.type > 0.05 & subset_q_vals_type$qvals.Type_permutation_13 > 0.05, "both not significant", "result differ"))

table(subset_q_vals_type$comparison) #they only differ by 21 genes which means the permutation analysis is pretty valid. 507 of the 514 genes in the original grouping is significant in both.

#------------------------------------------------------------------------------
#calculating the permutation p-value

colnames(q_vals_type)

#dropping permutation result from the default grouping
q_vals_type = subset(q_vals_type, select = -qvals.Type_permutation_13)
colnames(q_vals_type)

#filtered q_vals: focusing on only the significant genes from the default grouping
q_vals_type_filtered <- q_vals_type %>% filter(default_qvals.type <= 0.05)

# instances where: # of time(permutation group <= observed default group)/# of permutation
df_qval_type <- q_vals_type_filtered %>%
  select(-Gene, -default_qvals.type) %>%
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  mutate(perm_pval_type = rowMeans(. <= q_vals_type_filtered$default_qvals.type, na.rm = TRUE))

result_qval_type <- q_vals_type_filtered %>%
  mutate(perm_pval_type = df_qval_type$perm_pval_type)

count(result_qval_type$perm_pval_type <= 0.05)

#write.table(result_qval_type, "permutation_analysis/result_qval_type.txt", quote=FALSE, col.names = TRUE, row.names=TRUE, sep="\t")

# save gene names
result_qval_type <- subset(result_qval_type, perm_pval_type <= 0.05)
result_qval_type = data.frame(result_qval_type$Gene)
colnames(result_qval_type) = 'type_genes'
result_qval_type$type_genes = gsub("\\.v2\\.1","", result_qval_type$type_genes)
# write.table(result_qval_type, "permutation_analysis/type_genes_df.txt", quote=FALSE, sep="\t")

#------------------------------------------------------------------------------
# extracting the treatment:type q_vals of the default grouping for comparison 
# with the permutation analysis result
#------------------------------------------------------------------------------

colnames(stats_default)

q_vals_trt_type = data.frame(stats_default$qvals.Treatment.Type)
colnames(q_vals_trt_type) = "default_qvals.trt.type"
rownames(q_vals_trt_type) = rownames(stats_default)

table(rownames(q_vals_trt_type) %in% rownames(stats_default))
table(q_vals_trt_type$default_qvals.trt.type %in% stats_default$qvals.Treatment.Type)

#------------------------------------------------------------------------------
#renaming qvals.Treatment.Type in the permutation result list of dataframes to qvals.Treatment.Type_permutation(i)

colnames(stats_list[[2]])

for (i in 1:length(stats_list)) {
  # Get the current dataframe
  current_df <- stats_list[[i]]
  
  # Change the column name
  colnames(current_df)[colnames(current_df) == "qvals.Treatment.Type"] <-paste0("qvals.Treatment.Type_permutation_", i)
  
  # Update the dataframe in the list
  stats_list[[i]] <- current_df
}

colnames(stats_list[[2]])

#------------------------------------------------------------------------------
#merging q_vals_trt_type result from the permutations to q_vals_trt_type dataframe
q_vals_trt_type <- q_vals_trt_type %>% rownames_to_column(var = "Gene") #Convert row names to a column in 'q_vals_trt_type'

# Loop through the list of dataframes and merge the "qvals.Treatment.Type_permutation_(i)" column to 'q_vals_trt_type' by row names
for (i in 1:length(stats_list)) {
  q_vals_trt_type <- merge(q_vals_trt_type, stats_list[[i]][, c("Gene", paste0("qvals.Treatment.Type_permutation_", i))], by.x = "Gene", by.y = "Gene", all.x = TRUE)
}

#------------------------------------------------------------------------------
#comparing result from default grouping ran normally to result from default grouping ran alongside other permutation tests

colnames(q_vals_trt_type)

subset_q_vals_trt_type <- q_vals_trt_type[, c("Gene", "default_qvals.trt.type", "qvals.Treatment.Type_permutation_13")]

subset_q_vals_trt_type$comparison = ifelse(subset_q_vals_trt_type$default_qvals.trt.type <= 0.05 & subset_q_vals_trt_type$qvals.Treatment.Type_permutation_13 <= 0.05, "both siginificant", ifelse(subset_q_vals_trt_type$default_qvals.trt.type > 0.05 & subset_q_vals_trt_type$qvals.Treatment.Type_permutation_13 > 0.05, "both not significant", "result differ"))

table(subset_q_vals_trt_type$comparison) #they only differ by 20 genes which means the permutation analysis is pretty valid. 67 of the 76 genes in the original grouping is significant in both.

#------------------------------------------------------------------------------
#calculating the permutation p-value

colnames(q_vals_trt_type)

#dropping permutation result from the default grouping
q_vals_trt_type = subset(q_vals_trt_type, select = -qvals.Treatment.Type_permutation_13)
colnames(q_vals_trt_type)

#filtered q_vals: focusing on only the significant genes from the default grouping
q_vals_trt_type_filtered <- q_vals_trt_type %>% filter(default_qvals.trt.type <= 0.05)

# instances where: # of time(permutation group <= observed default group)/# of permutation
df_qval_trt_type <- q_vals_trt_type_filtered %>%
  select(-Gene, -default_qvals.trt.type) %>%
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  mutate(perm_pval_trt_type = rowMeans(. <= q_vals_trt_type_filtered$default_qvals.trt.type, na.rm = TRUE))

result_qval_trt_type <- q_vals_trt_type_filtered %>%
  mutate(perm_pval_trt_type = df_qval_trt_type$perm_pval_trt_type) 

count(result_qval_trt_type$perm_pval_trt_type <= 0.05)

#write.table(result_qval_trt_type, "permutation_analysis/result_qval_trt_type.txt", quote=FALSE, col.names = TRUE, row.names=TRUE, sep="\t")

# save gene names
result_qval_trt_type <- subset(result_qval_trt_type, perm_pval_trt_type <= 0.05)
result_qval_trt_type = data.frame(result_qval_trt_type$Gene)
colnames(result_qval_trt_type) = 'int_genes'
result_qval_trt_type$int_genes = gsub("\\.v2\\.1","", result_qval_trt_type$int_genes)
# write.table(result_qval_trt_type, "permutation_analysis/int_genes_df.txt", quote=FALSE, sep="\t")

