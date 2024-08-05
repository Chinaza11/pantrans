rm(list=ls())
ls()

library(topGO)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/Rio_reference_genome/")

# =============================================================================
# importing glmmseq result and preparing it
# =============================================================================

stats = read.table("glmmseq/glmm-allGenes-results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# get a list of gene names without v2 or v1
geneNames = gsub("\\.v2\\.1","", rownames(stats))

# =============================================================================
# read in the annotation file, and create a geneID2GO object
# =============================================================================

annot = read.table("GO_analysis/Rio_annotationInfo-reduced_deduplicate.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
x = strsplit(annot$GO, split=",")
names(x) = annot$locusName
head(x)

#------------------------------------------------------------------------------
# annotating all hybrid/merged genes
# hybrid/merged genes here means reads that mapped to more than one locus, these genes are joined by "--" e.g. SbRio.02G124800--SbRio.07G193200

# finding all rownames that has "--" and putting this in a list
hyb_index = grep("--", rownames(stats))

# create an empty dataframe to store the hybrid genes
hyb_genes = data.frame() 

# finding all hybrid genes and adding them to the new dataframe created
for (i in hyb_index) {
  selected_row = stats[i,] # select the row based on the index
  hyb_genes = rbind(hyb_genes, selected_row) # append the selected row to the new dataframe
}

# subseting dataframe to include only the q-vals for type, treatment and interactions
hyb_genes = hyb_genes[,c(22,23,24)]

# removing ".v2.1" at the end of the genes names
rownames(hyb_genes) = gsub("\\.v2\\.1","", rownames(hyb_genes)) 


# annotating all hybrid genes
for (i in 1:length(rownames(hyb_genes))) {
  temp = rownames(hyb_genes)[i] 
  genes = unlist(strsplit(temp, split="--"))
  m = match(genes, annot$locusName)
  hyb_genes$GO[i] = paste0(unique(annot$GO[m]), collapse=",")
}


#------------------------------------------------------------------------------
# adding the hybrid genes to the geneID2GO object

y = strsplit(hyb_genes$GO, split=",")
names(y) = rownames(hyb_genes)

x= c(x,y)

# =============================================================================
# FACTOR: Type
# =============================================================================

# predefining a list of interesting genes: assigning 0 and 1 based on p-value score of each gene from permutation result 
result_qval_type = read.table("glmmseq/permutation_analysis/result_qval_type.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

qval_type_gene_names = gsub("\\.v2\\.1","", result_qval_type$Gene)
qval_type = result_qval_type$perm_pval_type
names(qval_type) = qval_type_gene_names
sig_qval_type = qval_type[qval_type <= 0.05]

sig_qval_type_genes = names(sig_qval_type)
type_geneList <- factor(as.integer(geneNames %in% sig_qval_type_genes))
table(type_geneList)

names(type_geneList) <- geneNames
str(type_geneList)

#------------------------------------------------------------------------------
# create topGO objects
GOdata.type.MF = new("topGOdata", 
                     description="Sweet v Biomass", 
                     ontology="MF", 
                     allGenes=type_geneList, 
                     nodeSize=5, 
                     annot = annFUN.gene2GO, gene2GO = x)

#------------------------------------------------------------------------------
# enrichment test using Fisher's exact method
fisher.type.mf = runTest(GOdata.type.MF, algorithm="classic", statistic="fisher")

#------------------------------------------------------------------------------
# get table of test result
res.type.mf = GenTable(GOdata.type.MF, classicFisher=fisher.type.mf, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes=100)
res.type.mf$classicFisher = as.numeric(res.type.mf$classicFisher)
res.type.mf = res.type.mf[res.type.mf$classicFisher < 0.05,]

# write.table(res.type.mf, "GO_analysis/permutation/res.type.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)

# =============================================================================
# FACTOR: Interaction
# =============================================================================

# predefining a list of interesting genes: assigning 0 and 1 based on p-value score of each gene from permutation result 

result_qval_trt_type = read.table("glmmseq/permutation_analysis/result_qval_trt_type.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

qval_trt_type_gene_names = gsub("\\.v2\\.1","", result_qval_trt_type$Gene)
qval_trt_type = result_qval_trt_type$perm_pval_trt_type
names(qval_trt_type) = qval_trt_type_gene_names
sig_qval_trt_type = qval_trt_type[qval_trt_type <= 0.05]

sig_qval_trt_type_genes = names(sig_qval_trt_type)
trt_type_geneList <- factor(as.integer(geneNames %in% sig_qval_trt_type_genes))
table(trt_type_geneList)

names(trt_type_geneList) <- geneNames
str(trt_type_geneList)

#------------------------------------------------------------------------------
# create topGO objects
GOdata.int.MF = new("topGOdata", 
                     description="Interaction", 
                     ontology="MF", 
                     allGenes=trt_type_geneList, 
                     nodeSize=5,
                     annot = annFUN.gene2GO, gene2GO = x)

#------------------------------------------------------------------------------
# enrichment test using Fisher's exact method
fisher.int.mf = runTest(GOdata.int.MF, algorithm="classic", statistic="fisher")

#------------------------------------------------------------------------------
# get table of test result
res.int.mf = GenTable(GOdata.int.MF, classicFisher=fisher.int.mf, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes=100)
res.int.mf$classicFisher = as.numeric(res.int.mf$classicFisher)
res.int.mf = res.int.mf[res.int.mf$classicFisher < 0.05,]

# write.table(res.int.mf, "GO_analysis/permutation/res.int.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)
