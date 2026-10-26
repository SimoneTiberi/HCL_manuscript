setwd("~/Desktop/Emato")
########################################################################################
# HCL_100K vs. HCL_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HCL_100_vs_20.RData")
res_DGE = data.frame(gene = rownames(res$table), FDR_DGE = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DTU_HCL_100_vs_20.RData")
res_DTU = data.frame(gene = res_gene$gene, FDR_DTU = res_gene$qval ); rm(res); rm(res_gene)

RES = merge(res_DGE, res_DTU, by = "gene")

DF_upset = data.frame(DGE = ifelse(RES$FDR_DGE < 0.01, 1, 0),
                      DTU = ifelse(RES$FDR_DTU < 0.01, 1, 0) )

library(UpSetR)
upset(DF_upset)

########################################################################################
# HD_100K vs. HD_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HD_100_vs_20.RData")
res_DGE = data.frame(gene = rownames(res$table), FDR_DGE = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DTU_HD_100_vs_20.RData")
res_DTU = data.frame(gene = res_gene$gene, FDR_DTU = res_gene$qval ); rm(res); rm(res_gene)

RES = merge(res_DGE, res_DTU, by = "gene")

DF_upset = data.frame(DGE = ifelse(RES$FDR_DGE < 0.01, 1, 0),
                      DTU = ifelse(RES$FDR_DTU < 0.01, 1, 0) )

library(UpSetR)
upset(DF_upset)

########################################################################################
# HD_100K vs. HCL_100K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_100_HD_vs_HCL.RData")
res_DGE = data.frame(gene = rownames(res$table), FDR_DGE = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DTU_100_HD_vs_HCL.RData")
res_DTU = data.frame(gene = res_gene$gene, FDR_DTU = res_gene$qval ); rm(res); rm(res_gene)

RES = merge(res_DGE, res_DTU, by = "gene")

DF_upset = data.frame(DGE = ifelse(RES$FDR_DGE < 0.01, 1, 0),
                      DTU = ifelse(RES$FDR_DTU < 0.01, 1, 0) )

library(UpSetR)
upset(DF_upset)

########################################################################################
# HD_20K vs. HCL_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_20_HD_vs_HCL.RData")
res_DGE = data.frame(gene = rownames(res$table), FDR_DGE = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DTU_20_HD_vs_HCL.RData")
res_DTU = data.frame(gene = res_gene$gene, FDR_DTU = res_gene$qval ); rm(res); rm(res_gene)

RES = merge(res_DGE, res_DTU, by = "gene")

DF_upset = data.frame(DGE = ifelse(RES$FDR_DGE < 0.01, 1, 0),
                      DTU = ifelse(RES$FDR_DTU < 0.01, 1, 0) )

library(UpSetR)
upset(DF_upset)

