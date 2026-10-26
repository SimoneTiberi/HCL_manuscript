setwd("~/Desktop/Emato")
library(ggplot2)

rm(list = ls())

# Replace with the actual path to your GFF file.
gff <- rtracklayer::import(
  "Kallisto/Homo_sapiens.GRCh38.112.gff3.gz",
  format = "gff3",
  colnames = c("ID", "Name", "description")
)
# Keep gene records, including noncoding genes and pseudogenes.
g <- gff[grepl("^gene:", gff$ID)]
# annotation:
annotation <- unique(data.frame(
  ensembl_gene_id = sub("^gene:", "", as.character(g$ID)),
  gene_name = as.character(g$Name),
  description = as.character(g$description),
  stringsAsFactors = FALSE
))

########################################################################################
# HCL_100K vs. HCL_20K
########################################################################################
# rm(list =ls())
load("scripts/RData/data_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HCL_100K", "HCL_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi_gene$counts[,sel]
normMat = txi_gene$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi_gene)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
set.seed(169612)
eff.lib = calcNormFactors(normCts) * colSums(normCts)
normMat = sweep(normMat, 2, eff.lib, "*")
normMat = log(normMat)
y = DGEList(cts, group = group)
y = scaleOffset(y, normMat)

# filtering
design = model.matrix(~ group)
design # design matrix of our model, the only covariate is the group (A or B)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep, ]
# norm.factors column - a gene and sample specific offset was created to account for the average transcript length.

y = estimateDisp(y, design)

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit = glmFit(y, design)
design
lrt = glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

# check the counts for the top genes:
res = topTags(lrt, sort.by="PValue", n = Inf)

table(res$table$FDR < 0.01)
# ~ 2k significant genes
mean(res$table$FDR < 0.01)
# 17% significant results.
hist(res$table$PValue)

plotBCV(y)
# plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DGE_HCL_100_vs_20.RData")

# match gene name and description:
gene = rownames(res$table)
res$table$gene = gene

# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids <- sub("\\.[0-9]+$", "", as.character(gene))

G_list <- annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes <- setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res$table); dim(G_list)
res$table = merge(res$table, G_list, by = "gene", all.x = TRUE)
dim(res$table)
res$table = res$table[order(res$table$PValue), ]

library(EnhancedVolcano)
EnhancedVolcano(res$table,
                lab = res$table$gene_name,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DGE/Volcano_HCL_100_vs_20.pdf",
       width = 8,
       height = 8)

# plot CPM:
CPM = cpm(y$counts, log = FALSE)
colSums(CPM)
par(mfrow = c(3,3))
ordering = order(group)
group[ordering]
for(i in 1:27){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ordering ]
  plot(y = y,
       x = factor(seq_along(samples)),
       ylab = "CPMs",
       xlab = "samples",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(v = sum(group[ordering] == group[ordering][1]) + 0.5,
         col = "gray")
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}

pdf("scripts/images/Filippo/DGE/boxplots_HCL_100_vs_20.pdf",         # File name
    width = 12, height = 12)          # Paper size
par(mfrow = c(5,5))
for(i in 1:25){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ]
  plot(y = y,
       x = group,
       ylab = "CPMs",
       xlab = "groups",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}
dev.off() 

### Save results:
DF = res$table
colnames(DF)[2] = "log2FC"

sel_A = group %in% levels(group)[1]
avg_cpm = data.frame( gene = rownames(CPM),
                      rowMeans(CPM[,sel_A]),
                      rowMeans(CPM[,!sel_A]),
                      CPM)
colnames(avg_cpm)[2:3] = paste0("Avg_CPM_", levels(group))
colnames(avg_cpm)[-c(1:3)] = paste0("CPM_", samples)
DF = merge(DF, avg_cpm, by = "gene")

DF = DF[order(DF$PValue), ]
head(DF)

DF = DF[, c(1, 7:8, 5:6,2,9:ncol(DF)) ]
head(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/DGE_HCL_100_vs_20.xlsx")


########################################################################################
# HD_100K vs. HD_20K
########################################################################################
# rm(list =ls())
load("scripts/RData/data_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HD_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi_gene$counts[,sel]
normMat = txi_gene$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi_gene)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
set.seed(169612)
eff.lib = calcNormFactors(normCts) * colSums(normCts)
normMat = sweep(normMat, 2, eff.lib, "*")
normMat = log(normMat)
y = DGEList(cts, group = group)
y = scaleOffset(y, normMat)

# filtering
design = model.matrix(~ group)
design # design matrix of our model, the only covariate is the group (A or B)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep, ]
# norm.factors column - a gene and sample specific offset was created to account for the average transcript length.

y = estimateDisp(y, design)

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit = glmFit(y, design)
design
lrt = glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

# check the counts for the top genes:
res = topTags(lrt, sort.by="PValue", n = Inf)

table(res$table$FDR < 0.01)
# ~ 5k significant genes
mean(res$table$FDR < 0.01)
# 42% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DGE_HD_100_vs_20.RData")

# match gene name and description:
gene = rownames(res$table)
res$table$gene = gene

# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids <- sub("\\.[0-9]+$", "", as.character(gene))

G_list <- annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes <- setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res$table); dim(G_list)
res$table = merge(res$table, G_list, by = "gene", all.x = TRUE)
dim(res$table)
res$table = res$table[order(res$table$PValue), ]

library(EnhancedVolcano)
EnhancedVolcano(res$table,
                lab = res$table$gene_name,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DGE/Volcano_HD_100_vs_20.pdf",
       width = 8,
       height = 8)

# plot CPM:
CPM = cpm(y$counts, log = FALSE)
colSums(CPM)
par(mfrow = c(3,3))
ordering = order(group)
group[ordering]
for(i in 1:27){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ordering ]
  plot(y = y,
       # x = group,
       x = factor(seq_along(samples)),
       ylab = "CPMs",
       xlab = "samples",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(v = sum(group[ordering] == group[ordering][1]) + 0.5,
         col = "gray")
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}

pdf("scripts/images/Filippo/DGE/boxplots_HD_100_vs_20.pdf",         # File name
    width = 12, height = 12)          # Paper size
par(mfrow = c(5,5))
for(i in 1:25){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ]
  plot(y = y,
       x = group,
       ylab = "CPMs",
       xlab = "groups",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}
dev.off() 

### Save results:
DF = res$table
colnames(DF)[2] = "log2FC"

sel_A = group %in% levels(group)[1]
avg_cpm = data.frame( gene = rownames(CPM),
                      rowMeans(CPM[,sel_A]),
                      rowMeans(CPM[,!sel_A]),
                      CPM)
colnames(avg_cpm)[2:3] = paste0("Avg_CPM_", levels(group))
colnames(avg_cpm)[-c(1:3)] = paste0("CPM_", samples)
DF = merge(DF, avg_cpm, by = "gene")

DF = DF[order(DF$PValue), ]
head(DF)

DF = DF[, c(1, 7:8, 5:6,2,9:ncol(DF)) ]
head(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/DGE_HD_100_vs_20.xlsx")

########################################################################################
# HD_100K vs. HCL_100K
########################################################################################
# rm(list =ls())
load("scripts/RData/data_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HCL_100K")
samples = samples[sel]; samples
group = factor(group_ALL[sel], levels = c("HD_100K", "HCL_100K")); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi_gene$counts[,sel]
normMat = txi_gene$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi_gene)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
set.seed(169612)
eff.lib = calcNormFactors(normCts) * colSums(normCts)
normMat = sweep(normMat, 2, eff.lib, "*")
normMat = log(normMat)
y = DGEList(cts, group = group)
y = scaleOffset(y, normMat)

# filtering
design = model.matrix(~ group)
design # design matrix of our model, the only covariate is the group (A or B)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep, ]
# norm.factors column - a gene and sample specific offset was created to account for the average transcript length.

y = estimateDisp(y, design)

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit = glmFit(y, design)
design
lrt = glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

# check the counts for the top genes:
res = topTags(lrt, sort.by="PValue", n = Inf)

table(res$table$FDR < 0.01)
# ~ 2k significant genes
mean(res$table$FDR < 0.01)
# 14-15% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DGE_100_HD_vs_HCL.RData")

# match gene name and description:
gene = rownames(res$table)
res$table$gene = gene

# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids <- sub("\\.[0-9]+$", "", as.character(gene))

G_list <- annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes <- setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res$table); dim(G_list)
res$table = merge(res$table, G_list, by = "gene", all.x = TRUE)
dim(res$table)
res$table = res$table[order(res$table$PValue), ]

library(EnhancedVolcano)
EnhancedVolcano(res$table,
                lab = res$table$gene_name,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DGE/Volcano_100_HD_vs_HC.pdf",
       width = 8,
       height = 8)

# plot CPM:
CPM = cpm(y$counts, log = FALSE)
colSums(CPM)
par(mfrow = c(3,3))
ordering = order(group)
group[ordering]
for(i in 1:27){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ordering ]
  plot(y = y,
       # x = group,
       x = factor(seq_along(samples)),
       ylab = "CPMs",
       xlab = "samples",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(v = sum(group[ordering] == group[ordering][1]) + 0.5,
         col = "gray")
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}

pdf("scripts/images/Filippo/DGE/boxplots_100_HD_vs_HC.pdf",         # File name
    width = 12, height = 12)          # Paper size
par(mfrow = c(5,5))
for(i in 1:25){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ]
  plot(y = y,
       x = group,
       ylab = "CPMs",
       xlab = "groups",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}
dev.off() 

### Save results:
DF = res$table
colnames(DF)[2] = "log2FC"

sel_A = group %in% levels(group)[1]
avg_cpm = data.frame( gene = rownames(CPM),
                      rowMeans(CPM[,sel_A]),
                      rowMeans(CPM[,!sel_A]),
                      CPM)
colnames(avg_cpm)[2:3] = paste0("Avg_CPM_", levels(group))
colnames(avg_cpm)[-c(1:3)] = paste0("CPM_", samples)
DF = merge(DF, avg_cpm, by = "gene")

DF = DF[order(DF$PValue), ]
head(DF)

DF = DF[, c(1, 7:8, 5:6,2,9:ncol(DF)) ]
head(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/DGE_100_HD_vs_HCL.xlsx")

########################################################################################
# HD_20K vs. HCL_20K
########################################################################################
# rm(list =ls())
load("scripts/RData/data_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_20K", "HCL_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel], levels = c("HD_20K", "HCL_20K")); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi_gene$counts[,sel]
normMat = txi_gene$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi_gene)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
set.seed(169612)
eff.lib = calcNormFactors(normCts) * colSums(normCts)
normMat = sweep(normMat, 2, eff.lib, "*")
normMat = log(normMat)
y = DGEList(cts, group = group)
y = scaleOffset(y, normMat)

# filtering
design = model.matrix(~ group)
design # design matrix of our model, the only covariate is the group (A or B)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep, ]
# norm.factors column - a gene and sample specific offset was created to account for the average transcript length.

y = estimateDisp(y, design)

# To perform likelihood ratio test, edgeR applies a glm model, keeping fixed the gene dispersions to the estimated values.
fit = glmFit(y, design)
design
lrt = glmLRT(fit, coef=2) # coef specifies the coefficient to test.
lrt

# check the counts for the top genes:
res = topTags(lrt, sort.by="PValue", n = Inf)

table(res$table$FDR < 0.01)
# ~ 1.4k significant genes
mean(res$table$FDR < 0.01)
# 13% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DGE_20_HD_vs_HCL.RData")

# match gene name and description:
gene = rownames(res$table)
res$table$gene = gene

# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids <- sub("\\.[0-9]+$", "", as.character(gene))

G_list <- annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes <- setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res$table); dim(G_list)
res$table = merge(res$table, G_list, by = "gene", all.x = TRUE)
dim(res$table)
res$table = res$table[order(res$table$PValue), ]

library(EnhancedVolcano)
EnhancedVolcano(res$table,
                lab = res$table$gene_name,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DGE/Volcano_20_HD_vs_HC.pdf",
       width = 8,
       height = 8)

# plot CPM:
CPM = cpm(y$counts, log = FALSE)
colSums(CPM)
par(mfrow = c(3,3))
ordering = order(group)
group[ordering]
for(i in 1:27){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ordering ]
  plot(y = y,
       # x = group,
       x = factor(seq_along(samples)),
       ylab = "CPMs",
       xlab = "samples",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(v = sum(group[ordering] == group[ordering][1]) + 0.5,
         col = "gray")
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}


pdf("scripts/images/Filippo/DGE/boxplots_20_HD_vs_HC.pdf",         # File name
    width = 12, height = 12)          # Paper size
par(mfrow = c(5,5))
for(i in 1:25){
  gene = res$table$gene[i]
  sel = which(rownames(CPM) == gene)
  y = CPM[sel, ]
  plot(y = y,
       x = group,
       ylab = "CPMs",
       xlab = "groups",
       ylim = c(0, max(y)),
       main = res$table$gene_name[i])
  abline(h = 0, col = "gray")
  # main represents the raking of the gene
}
dev.off() 

### Save results:
DF = res$table
colnames(DF)[2] = "log2FC"

sel_A = group %in% levels(group)[1]
avg_cpm = data.frame( gene = rownames(CPM),
                      rowMeans(CPM[,sel_A]),
                      rowMeans(CPM[,!sel_A]),
                      CPM)
colnames(avg_cpm)[2:3] = paste0("Avg_CPM_", levels(group))
colnames(avg_cpm)[-c(1:3)] = paste0("CPM_", samples)
DF = merge(DF, avg_cpm, by = "gene")

DF = DF[order(DF$PValue), ]
head(DF)

DF = DF[, c(1, 7:8, 5:6,2,9:ncol(DF)) ]
head(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/DGE_20_HD_vs_HCL.xlsx")
