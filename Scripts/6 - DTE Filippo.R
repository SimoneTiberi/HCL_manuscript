setwd("~/Desktop/Emato")
library(ggplot2)
########################################################################################
# HCL_100K vs. HCL_20K
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HCL_100K", "HCL_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi$counts[,sel]
rm(txi)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
y = DGEList(cts, group = group)
design = model.matrix(~ group)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep,,keep.lib.sizes=FALSE]
y = normLibSizes(y)
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
# 1.8 k transcripts
mean(res$table$FDR < 0.01)
# 17% significant results.
hist(res$table$PValue)

plotBCV(y)
# plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DTE_HCL_100_vs_20.RData")


########################################################################################
# HD_100K vs. HD_20K
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HD_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi$counts[,sel]
normMat = txi$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
y = DGEList(cts, group = group)
design = model.matrix(~ group)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep,,keep.lib.sizes=FALSE]
y = normLibSizes(y)
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
# 2.3 k transcripts
mean(res$table$FDR < 0.01)
# 42% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DTE_HD_100_vs_20.RData")




########################################################################################
# HD_100K vs. HCL_100K
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HCL_100K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi$counts[,sel]
normMat = txi$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
y = DGEList(cts, group = group)
design = model.matrix(~ group)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep,,keep.lib.sizes=FALSE]
y = normLibSizes(y)
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
# 0.8 k transcripts
mean(res$table$FDR < 0.01)
# 14-15% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DTE_100_HD_vs_HCL.RData")



########################################################################################
# HD_20K vs. HCL_20K
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_20K", "HCL_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

cts     = txi$counts[,sel]
normMat = txi$length[,sel]
normMat = normMat/exp(rowMeans(log(normMat)))
normCts = cts/normMat
rm(txi)
# We compute a normalization factor (to add as an offset in the NB) for the average transcript length.
# see: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Section: Use with downstream Bioconductor DGE packages

library(edgeR)
y = DGEList(cts, group = group)
design = model.matrix(~ group)
keep = filterByExpr(y, design)
mean(keep); sum(keep)
y = y[keep,,keep.lib.sizes=FALSE]
y = normLibSizes(y)
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
# 0.7 k transcripts
mean(res$table$FDR < 0.01)
# 13% significant results.
hist(res$table$PValue)

plotBCV(y)
#plotSmear(y)

save(lrt, res, file = "scripts/RData/Filippo_DTE_20_HD_vs_HCL.RData")

