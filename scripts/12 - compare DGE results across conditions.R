setwd("~/Desktop/Emato")
########################################################################################
# Filippo
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HCL_100_vs_20.RData")
res_1 = data.frame(gene = rownames(res$table), HCL_100_vs_20 = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DGE_HD_100_vs_20.RData")
res_2 = data.frame(gene = rownames(res$table), HD_100_vs_20 = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DGE_100_HD_vs_HCL.RData")
res_3 = data.frame(gene = rownames(res$table), `100_HD_vs_HCL` = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DGE_20_HD_vs_HCL.RData")
res_4 = data.frame(gene = rownames(res$table), `20_HD_vs_HCL` = res$table$FDR ); rm(res); rm(lrt)

RES = merge(res_1, res_2, by = "gene")
RES = merge(RES, res_3, by = "gene")
RES = merge(RES, res_4, by = "gene")

DF_upset = data.frame(HCL20_vs_HCL100 = ifelse(RES$HCL_100_vs_20 < 0.01, 1, 0),
                      HD20_vs_HD100 = ifelse(RES$HD_100_vs_20 < 0.01, 1, 0),
                      HD100_vs_HCL100 = ifelse(RES$X100_HD_vs_HCL < 0.01, 1, 0),
                      HD20_vs_HCL20 = ifelse(RES$X20_HD_vs_HCL < 0.01, 1, 0))

library(UpSetR)
a = upset(DF_upset)
pdf("scripts/Images/Filippo/DGE/upsetR.pdf",
    width = 8,
    height = 8)
a
dev.off()

########################################################################################
# Filippo - 20 vs 100
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HCL_100_vs_20.RData")
res_1 = data.frame(gene = rownames(res$table), HCL_100_vs_20 = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DGE_HD_100_vs_20.RData")
res_2 = data.frame(gene = rownames(res$table), HD_100_vs_20 = res$table$FDR ); rm(res); rm(lrt)

RES = merge(res_1, res_2, by = "gene")

DF_upset = data.frame(HCL20_vs_HCL100 = ifelse(RES$HCL_100_vs_20 < 0.01, 1, 0),
                      HD20_vs_HD100 = ifelse(RES$HD_100_vs_20 < 0.01, 1, 0))

library(UpSetR)
a = upset(DF_upset)
pdf("scripts/Images/Filippo/DGE/upsetR_20_vs_100.pdf",
    width = 8,
    height = 8)
a
dev.off()

########################################################################################
# Filippo - HD vs HCL
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_100_HD_vs_HCL.RData")
res_1 = data.frame(gene = rownames(res$table), `100_HD_vs_HCL` = res$table$FDR ); rm(res); rm(lrt)
load("scripts/RData/Filippo_DGE_20_HD_vs_HCL.RData")
res_2 = data.frame(gene = rownames(res$table), `20_HD_vs_HCL` = res$table$FDR ); rm(res); rm(lrt)

RES = merge(res_1, res_2, by = "gene")

DF_upset = data.frame(HD100_vs_HCL100 = ifelse(RES$X100_HD_vs_HCL < 0.01, 1, 0),
                      HD20_vs_HCL20 = ifelse(RES$X20_HD_vs_HCL < 0.01, 1, 0))

library(UpSetR)
a = upset(DF_upset)
pdf("scripts/Images/Filippo/DGE/upsetR_HD_vs_HCL.pdf",
    width = 8,
    height = 8)
a
dev.off()
