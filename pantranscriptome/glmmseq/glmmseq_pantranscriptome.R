rm(list=ls())
ls() 

library(dplyr)
library(DESeq2)

# =============================================================================
# Set up the working directory, and read in raw counts from mmquant
# The format is rows=genes; columns=samples
# =============================================================================

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/glmmseq/")
cts = as.matrix(read.table("rawCounts20_pantranscriptome.txt", header=TRUE, row.names="Gene", sep="\t"))
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

colData = data.frame(SAMID = as.factor(SampleID),
                     GENID = as.factor(GenID),
                     Treatment = Condition,
                     Type = as.factor(Genotype))

rownames(colData) = SampleID

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ Treatment + Type)


# =============================================================================
# Use DESeq to 1) remove genes with low counts
# 2) get normalized count values
# and 3) estimate dispersion factors
# =============================================================================

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,] #this filter ensures at least 3 samples with a count of 5 or more
normcts = counts(dds, normalized=TRUE)
dds <- DESeq(dds)
dispersions <- setNames(dispersions(dds), rownames(normcts))

# =============================================================================
# Now we have everything we need for glmmSeq
# Load glmmSeq package and set a random number seed
# =============================================================================

library(glmmSeq)
set.seed(1234)

# =============================================================================
# Fit the model with Fixed effects plus Random Effects. Treatment and Type are 
# the fixed effects, while the GENID (genotype) is the random effect
# This model fitting step may take 15-20 minutes to run depending on the no of 
# cores available.
# =============================================================================

# results <- glmmSeq(~ Treatment * Type + (1 | GENID),
#                     id = "GENID",
#                     countdata = normcts,
#                     metadata = colData,
#                     dispersion = dispersions,
#                     removeSingles=FALSE,
#                     cores = 4,
#                     progress = TRUE)

# The variables used by the model are in the @modeldata:

# save(results, file = "results.Rdata")
load(file = "results.Rdata")

results@modelData

# Get Q values and a significance summary table with correct P-values
results <- glmmQvals(results)
stats = data.frame(results@stats)

# Write this table to output for future loading
# write.table(stats, "glmm-allGenes-results.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# Getting significant genes in Treatment, Type and Interaction
# Venn diagram of the three fixed effects.
# =============================================================================

type.sig = stats[which(stats$qvals.Type < 0.05),]
treat.sig = stats[which(stats$qvals.Treatment < 0.05),]
inter.sig = stats[which(stats$qvals.Treatment.Type < 0.05),]

# Re-order each set of results to be in order of significance
type.sig = type.sig[order(type.sig$qvals.Type),]
treat.sig = treat.sig[order(treat.sig$qvals.Treatment),]
inter.sig = inter.sig[order(inter.sig$qvals.Treatment.Type),]

# Create Venn Diagram of overlapping genes in each set
library(VennDiagram)
x <- list(
  Type = rownames(type.sig), 
  Treatment = rownames(treat.sig), 
  Interaction = rownames(inter.sig)
)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

png(file="Venn.png", width=10, height=10, units="in", res=300)

display_venn(
  x,
  category.names = c("Type" , "Treatment" , "Interaction (Type:Treatment)"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("purple", "orange", "steelblue"),
  # Numbers
  cex = 2,
  fontface = "italic",
  # Set names
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05, 0.05)
)

dev.off()

# Get genes that overlap among the 3 sets
x24 = intersect(rownames(type.sig), rownames(treat.sig)) #intersect between type and treatment
x64 = intersect(rownames(type.sig), rownames(inter.sig)) #intersect between type and interaction
x54 = intersect(rownames(treat.sig), rownames(inter.sig)) #intersect between treatment and interaction
x4 = intersect(x24,x54) #intersect between type, treatment and interaction
x1 = rownames(type.sig)[which(is.na(match(rownames(type.sig),unique(c(x24,x64)))))] #genes unique to type 
x2 = x24[which(is.na(match(x24,x4)))] #intersect between type and treatment
x3 = rownames(treat.sig)[which(is.na(match(rownames(treat.sig),unique(c(x24,x54)))))] #genes unique to treatment
x5 = x54[which(is.na(match(x54,x4)))] #intersect between treatment and interaction
x6 = x64[which(is.na(match(x64,x4)))] #intersect between type and interaction
x7 = rownames(inter.sig)[which(is.na(match(rownames(inter.sig),unique(c(x54,x64)))))] #genes unique to the interaction

# =============================================================================
# Re-do the above intersection analysis, but now only with genes that show at 
# least a 1.5 fold change (for sweet and iron, not interaction)
# =============================================================================

#type.sig = type.sig[which(abs(type.sig$coef.TypeSweet) > 1.5),]
#treat.sig = treat.sig[which(abs(treat.sig$coef.TreatmentN) > 1.5),]

# =============================================================================
# code below was done to extract the data for descriptive statistics in 
# no_of_annotated_sig_genes_check.py file
# =============================================================================

# subset_treat.sig <- subset(treat.sig, select = "qvals.Treatment")
# write.table(subset_treat.sig, "result_qval_trt.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
# 
# subset_type.sig <- subset(type.sig, select = "qvals.Type")
# write.table(subset_type.sig, "result_qval_type.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
# 
# subset_int.sig <- subset(inter.sig, select = "qvals.Treatment.Type")
# write.table(subset_int.sig, "result_qval_int.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# 
# =============================================================================

