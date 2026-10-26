########################################################################################
# HD_100K vs. HCL_100K
########################################################################################
setwd("~/Desktop/Emato")
load("scripts/RData/Filippo_DTU_100_HD_vs_HCL.RData")

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)

# select significant results (gene-level)
res_gene = res_gene[res_gene$qval < 0.01,]
dim(res_gene)

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[4:5], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,1]) != which.max(x[,2])
}))

round(100 * mean(res_gene_split_2), 2)
# 67.1

########################################################################################
# HD_20K vs. HCL_20K
########################################################################################
setwd("~/Desktop/Emato")
load("scripts/RData/Filippo_DTU_20_HD_vs_HCL.RData")

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)

# select significant results (gene-level)
res_gene = res_gene[res_gene$qval < 0.01,]
dim(res_gene)

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[4:5], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,1]) != which.max(x[,2])
}))
round(100 * mean(res_gene_split_2), 2)
# 56.82

########################################################################################
# HD_20K vs. HD_100K
########################################################################################
setwd("~/Desktop/Emato")
load("scripts/RData/Filippo_DTU_HD_100_vs_20.RData")

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)

# select significant results (gene-level)
res_gene = res_gene[res_gene$qval < 0.01,]
dim(res_gene)

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[4:5], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,1]) != which.max(x[,2])
}))

round(100 * mean(res_gene_split_2), 2)
# 56.27

########################################################################################
# HCL_20K vs. HCL_100K
########################################################################################
setwd("~/Desktop/Emato")
load("scripts/RData/Filippo_DTU_HCL_100_vs_20.RData")

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)

# select significant results (gene-level)
res_gene = res_gene[res_gene$qval < 0.01,]
dim(res_gene)

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[4:5], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,1]) != which.max(x[,2])
}))

round(100 * mean(res_gene_split_2), 2)
# 45.92
