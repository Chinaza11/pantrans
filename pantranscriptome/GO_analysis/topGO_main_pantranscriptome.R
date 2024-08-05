rm(list=ls())
ls()

library(topGO)

setwd("G:/My Drive/PhD/project/Iron_RNASeq_sorghum/data_analysis/data_and_result/pantranscriptome/")

# =============================================================================
# importing glmmseq result and preparing it
# =============================================================================

stats = read.table("glmmseq/glmm-allGenes-results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Get a list of gene names
geneNames = rownames(stats)
head(geneNames, n=10)

# Create 3 vectors of p-values: one for each factor
typeList = stats$qvals.Type
treatList = stats$qvals.Treatment
intList = stats$qvals.Treatment.Type

names(typeList) = geneNames
names(treatList) = geneNames
names(intList) = geneNames

head(typeList)

# =============================================================================
# Read in the annotation file, and create a geneID2GO object
# =============================================================================

annot = read.csv("GO_analysis/pantranscriptome_plus_GO_cleaned.csv")

x = strsplit(annot$GO_annotations_cleaned, split=",")
names(x) = annot$sequence_name
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

# annotating all hybrid genes
for (i in 1:length(rownames(hyb_genes))) {
  temp = rownames(hyb_genes)[i] 
  genes = unlist(strsplit(temp, split="--"))
  m = match(genes, annot$sequence_name)
  hyb_genes$GO[i] = paste0(unique(annot$GO_annotations_cleaned[m]), collapse=",")
}

# write.table(hyb_genes, "GO_analysis/hyb_genes.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)

#------------------------------------------------------------------------------
#adding the hybrid genes to the geneID2GO object

y = strsplit(hyb_genes$GO, split=",")
names(y) = rownames(hyb_genes)

x = c(x,y)

# =============================================================================
# Create topGO objects for the Treatment, Type, and Interaction effects
# =============================================================================

topDiffGenes = function (allScore){
  return(allScore < 0.05)
}

GOdata.type.MF = new("topGOdata", 
                     description="Sweet v Biomass", 
                     ontology="MF", 
                     allGenes=typeList, 
                     geneSel = topDiffGenes, 
                     nodeSize=5, 
                     annot = annFUN.gene2GO, gene2GO = x)

GOdata.treat.MF = new("topGOdata", 
                      description="Iron Treatment", 
                      ontology="MF", 
                      allGenes=treatList, 
                      geneSel = topDiffGenes, 
                      nodeSize=5, 
                      annot = annFUN.gene2GO, gene2GO = x)

GOdata.int.MF = new("topGOdata", 
                    description="Interaction", 
                    ontology="MF", 
                    allGenes=intList, 
                    geneSel = topDiffGenes, 
                    nodeSize=5, 
                    annot = annFUN.gene2GO, gene2GO = x)

# =============================================================================
# Perform enrichment tests with Fisher method for all 3 sets
# =============================================================================

fisher.type.mf = runTest(GOdata.type.MF, algorithm="classic", statistic="fisher")
fisher.treat.mf = runTest(GOdata.treat.MF, algorithm="classic", statistic="fisher")
fisher.int.mf = runTest(GOdata.int.MF, algorithm="classic", statistic="fisher")

# =============================================================================
# Get table of test result
# =============================================================================

res.type.mf = GenTable(GOdata.type.MF, classicFisher=fisher.type.mf, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes=100)
res.type.mf$classicFisher = gsub('< ', '', res.type.mf$classicFisher)
res.type.mf$classicFisher = as.numeric(res.type.mf$classicFisher)
res.type.mf = res.type.mf[res.type.mf$classicFisher < 0.05,]

res.treat.mf = GenTable(GOdata.treat.MF, classicFisher=fisher.treat.mf, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes=100)
res.treat.mf$classicFisher = as.numeric(res.treat.mf$classicFisher)
res.treat.mf <- res.treat.mf[res.treat.mf$classicFisher<0.05,]

res.int.mf = GenTable(GOdata.int.MF, classicFisher=fisher.int.mf, orderBy="classicFisher", ranksOf="classicFisher", topNodes=100)
res.int.mf$classicFisher = as.numeric(res.int.mf$classicFisher)
res.int.mf <- res.int.mf[res.int.mf$classicFisher<0.05,]

# write.table(res.type.mf, "GO_analysis/res.type.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)
# write.table(res.treat.mf, "GO_analysis/res.treat.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)
# write.table(res.int.mf, "GO_analysis/res.int.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)


# =============================================================================
# Checking the significant genes associated to an enriched GO term
# =============================================================================

df_long <- stack(x)
View(df_long)

stats$locusName = rownames(stats)

threshold <- 0.05
significant_int_genes <- names(intList)[intList < threshold]

no_of_annotated_genes = 0
for (i in 1:length(df_long$values)){
  if ("GO:0015145" %in% df_long$values[i]){
    no_of_annotated_genes = no_of_annotated_genes + 1
    for (j in significant_int_genes){
      if (j %in% df_long$ind[i]) {
        print(j)
      }
    }
  }
}

View(stats) # manually search here to confirm that the DEGs is significant
View(annot) # manually search here to get the locs tag that can be searched on NCBI

# =============================================================================
# Plot the results
# =============================================================================

# Interaction effect result
res.int.mf$classicFisher <- as.numeric(res.int.mf$classicFisher)
goEnrichment_int <- res.int.mf[res.int.mf$classicFisher<0.05,]
goEnrichment_int <- goEnrichment_int[,c("GO.ID","Term","classicFisher")]
goEnrichment_int$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment_int$Term)
goEnrichment_int$Term <- gsub("\\.\\.\\.$", "", goEnrichment_int$Term)
goEnrichment_int$Term <- paste(goEnrichment_int$GO.ID, goEnrichment_int$Term, sep=", ")
goEnrichment_int$Term <- factor(goEnrichment_int$Term, levels=rev(goEnrichment_int$Term))

require(ggplot2)

png(file="GO_analysis/int_fisher_all_terms.png", width=9, height=8, units="in", res=500)
ggplot(goEnrichment_int, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="steelblue", color="steelblue") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Interaction") +
  scale_y_continuous(limits = c(0, max(-log10(goEnrichment_int$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()

#------------------------------------------------------------------------------
# Treatment effect result
res.treat.mf$classicFisher <- as.numeric(res.treat.mf$classicFisher)
goEnrichment_treat <- res.treat.mf[res.treat.mf$classicFisher<0.05,]
goEnrichment_treat <- goEnrichment_treat[,c("GO.ID","Term","classicFisher")]
goEnrichment_treat$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment_treat$Term)
goEnrichment_treat$Term <- gsub("\\.\\.\\.$", "", goEnrichment_treat$Term)
goEnrichment_treat$Term <- paste(goEnrichment_treat$GO.ID, goEnrichment_treat$Term, sep=", ")
goEnrichment_treat$Term <- factor(goEnrichment_treat$Term, levels=rev(goEnrichment_treat$Term))

png(file="GO_analysis/treat_fisher_all_terms.png", width=9, height=8, units="in", res=500)
ggplot(goEnrichment_treat, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="orange", color="orange") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Treatment") +
  scale_y_continuous(limits = c(0, max(-log10(goEnrichment_treat$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()

#------------------------------------------------------------------------------
# Type effect result
res.type.mf$classicFisher[res.type.mf$classicFisher == '< 1e-30'] <- 1e-30
res.type.mf$classicFisher <- as.numeric(res.type.mf$classicFisher)
goEnrichment_type <- res.type.mf[res.type.mf$classicFisher<0.05,]
goEnrichment_type <- goEnrichment_type[,c("GO.ID","Term","classicFisher")]
goEnrichment_type$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment_type$Term)
goEnrichment_type$Term <- gsub("\\.\\.\\.$", "", goEnrichment_type$Term)
goEnrichment_type$Term <- paste(goEnrichment_type$GO.ID, goEnrichment_type$Term, sep=", ")
goEnrichment_type$Term <- factor(goEnrichment_type$Term, levels=rev(goEnrichment_type$Term))

png(file="GO_analysis/type_fisher_all_terms.png", width=9, height=8, units="in", res=500)
ggplot(goEnrichment_type, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="purple", color="purple") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Type") +
  scale_y_continuous(limits = c(0, max(-log10(goEnrichment_type$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()

# =============================================================================
# Trying different algorithms and statistical analysis
# =============================================================================

# Treatment
ks_classic.treat.mf = runTest(GOdata.treat.MF, algorithm="classic", statistic="ks")
fisher_weight01.treat.mf = runTest(GOdata.treat.MF, algorithm="weight01", statistic="fisher")
ks_weight01.treat.mf = runTest(GOdata.treat.MF, algorithm="weight01", statistic="ks")

all_res.treat.mf = GenTable(GOdata.treat.MF, classicFisher=fisher.treat.mf, orderBy="weight01KS", ranksOf="classicFisher", classicKS=ks_classic.treat.mf, weight01KS=ks_weight01.treat.mf, weight01Fisher=fisher_weight01.treat.mf, topNodes=100)

# write.table(all_res.treat.mf, "GO_analysis/all_res.treat.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)

#------------------------------------------------------------------------------
# Type
ks_classic.type.mf = runTest(GOdata.type.MF, algorithm="classic", statistic="ks")
fisher_weight01.type.mf = runTest(GOdata.type.MF, algorithm="weight01", statistic="fisher")
ks_weight01.type.mf = runTest(GOdata.type.MF, algorithm="weight01", statistic="ks")

all_res.type.mf = GenTable(GOdata.type.MF, classicFisher=fisher.type.mf, orderBy="weight01KS", ranksOf="classicFisher", classicKS=ks_classic.type.mf, weight01KS=ks_weight01.type.mf, weight01Fisher=fisher_weight01.type.mf, topNodes=100)

# write.table(all_res.type.mf, "GO_analysis/all_res.type.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)

#------------------------------------------------------------------------------
# Interaction
ks_classic.int.mf = runTest(GOdata.int.MF, algorithm="classic", statistic="ks")
fisher_weight01.int.mf = runTest(GOdata.int.MF, algorithm="weight01", statistic="fisher")
ks_weight01.int.mf = runTest(GOdata.int.MF, algorithm="weight01", statistic="ks")

all_res.int.mf = GenTable(GOdata.int.MF, classicFisher=fisher.int.mf, orderBy="weight01KS", ranksOf="classicFisher", classicKS=ks_classic.int.mf, weight01KS=ks_weight01.int.mf, weight01Fisher=fisher_weight01.int.mf, topNodes=100)

# write.table(all_res.int.mf, "GO_analysis/all_res.int.mf.txt", quote=FALSE, sep="\t",row.names=TRUE, col.names=TRUE)

# =============================================================================
# sub-setting result and writing it out so you can use REVIGO(http://revigo.irb.hr/) 
# to eliminate redundant terms. Zea mays was selected under species since it is the 
# closest relative to sorghum
# =============================================================================

sub_res.type.mf <- res.type.mf[,c("GO.ID","classicFisher")]
sub_res.treat.mf <- res.treat.mf[,c("GO.ID","classicFisher")]
sub_res.int.mf <- res.int.mf[,c("GO.ID","classicFisher")]

# write.table(sub_res.type.mf, "GO_analysis/sub_res.type.mf.txt", quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)
# write.table(sub_res.treat.mf, "GO_analysis/sub_res.treat.mf.txt", quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)
# write.table(sub_res.int.mf, "GO_analysis/sub_res.int.mf.txt", quote=FALSE, sep="\t",row.names=FALSE, col.names=TRUE)

#------------------------------------------------------------------------------
# Type

revigo_type = read.table("GO_analysis/Revigo_MF_Type_OnScreenTable.tsv", header=TRUE, sep="\t")

merge_type = merge(goEnrichment_type, revigo_type, by.x = "GO.ID", by.y = "TermID")

merge_type = subset(merge_type, Representative == "null") # remove non-representative GO terms

merge_type = merge_type[,c("GO.ID","Term", "classicFisher")]

png(file="GO_analysis/type_fisher_redundant_terms_removed.png", width=9, height=8, units="in", res=500)
ggplot(merge_type, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="purple", color="purple") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Type") +
  scale_y_continuous(limits = c(0, max(-log10(merge_type$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()

#------------------------------------------------------------------------------
# Treatment

revigo_treatment = read.table("GO_analysis/Revigo_MF_Treat_OnScreenTable.tsv", header=TRUE, sep="\t")

merge_treatment = merge(goEnrichment_treat, revigo_treatment, by.x = "GO.ID", by.y = "TermID")

merge_treatment = subset(merge_treatment, Representative == "null")

merge_treatment = merge_treatment[,c("GO.ID","Term", "classicFisher")]

png(file="GO_analysis/treatment_fisher_redundant_terms_removed.png", width=9, height=8, units="in", res=500)
ggplot(merge_treatment, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="orange", color="orange") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Treatment") +
  scale_y_continuous(limits = c(0, max(-log10(merge_treatment$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()

#------------------------------------------------------------------------------
# Interaction

revigo_int = read.table("GO_analysis/Revigo_MF_Int_OnScreenTable.tsv", header=TRUE, sep="\t")

merge_int = merge(goEnrichment_int, revigo_int, by.x = "GO.ID", by.y = "TermID")

merge_int = subset(merge_int, Representative == "null")

merge_int = merge_int[,c("GO.ID","Term", "classicFisher")]

png(file="GO_analysis/int_fisher_redundant_terms_removed.png", width=9, height=8, units="in", res=500)
ggplot(merge_int, aes(x=Term, y=-log10(classicFisher))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", fill="steelblue", color="steelblue") +
  xlab("") +
  ylab(bquote(Enrichment * " " ~(-Log[10]~ italic(P-value)))) +
  ggtitle("Interaction") +
  scale_y_continuous(limits = c(0, max(-log10(merge_int$classicFisher))+2)) + 
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=12, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=12),  #Text size
    title=element_text(size=12)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()
dev.off()
