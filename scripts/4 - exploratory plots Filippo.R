setwd("~/Desktop/Emato")
library(ggplot2)

rm(list =ls())
load("scripts/RData/data_Filippo.RData")
Gene_counts = txi_gene$counts
colnames(Gene_counts) = samples
length(samples)

sort(round(colSums(Gene_counts)/10^6, 2))
# HD 47 have almost 0 counts
# maybe remove the sample

# define colour scheme
cols_HD = ifelse(group_HD == "HCL", "red", "blue")
cols_20 = ifelse(group_20 == "100K", "red", "blue")
cols_ALL = ifelse(group_ALL == "HCL_100K", "red",
                  ifelse(group_ALL == "HCL_20K", "orange",
                         ifelse(group_ALL == "HD_100K", "blue",
                                "lightblue"
                         ) 
                    )
                  )

# 2D MDS plot based on Gene-level counts:
library(edgeR);
dd = DGEList(counts = Gene_counts, samples = samples)
dd = calcNormFactors(dd)
m3 = plotMDS(dd, labels = samples, main = "MDS")
m3 = plotMDS(dd, labels = group_HD, col = cols_HD, main = "MDS")
m3 = plotMDS(dd, labels = group_20, col = cols_20, main = "MDS")
m3 = plotMDS(dd, labels = group_ALL, col = cols_ALL, main = "MDS")

pdf("scripts/images/Filippo/Exploratory/MDS.pdf",         # File name
    width = 8, height = 8)          # Paper size
plotMDS(dd, labels = group_ALL, col = cols_ALL, main = "MDS")
# Closing the graphical device
dev.off() 

# PCA plot:
library(ggfortify)
logcounts = cpm(dd, log = TRUE)
colnames(logcounts) = samples
pca = prcomp(t(logcounts))
rownames(pca$x) = group_ALL
pdf("scripts/images/Filippo/Exploratory/PCA.pdf",         # File name
    width = 8, height = 8)          # Paper size
autoplot( pca, col = cols_ALL, 
          label = TRUE) + theme_bw()
dev.off() 

# k-means clustering:
km = kmeans(t(logcounts), 4, nstart = 100, iter.max = 1000)
table(km$cluster, group_ALL)
# 1 group: odd sample

km = kmeans(t(logcounts), 6, nstart = 100, iter.max = 1000)
table(km$cluster, group_ALL)
# 20 and 100 never mix in the same group (except for the 2 outliers)
# 20 vs. 100 K difference bigger than HCL vs. HD

# heatmap based on the most variable TRANSCRIPTS:
var_genes = apply(logcounts, 1, var)
# we select the 2*10^3 most variable log-cpm
select_var = names(sort(var_genes, decreasing=TRUE))[1:2000]
# Subset logcounts matrix (only keep the most variables log-cpm)
highly_variable_lcpm = logcounts[select_var,]
library(corrplot); library( gplots) 
# we apply a hierarchical clustering on the samples and on the rows.
colnames(highly_variable_lcpm) = group_ALL

pdf("scripts/images/Filippo/Exploratory/heatmap.pdf",         # File name
    width = 8, height = 8)          # Paper size
heatmap.2(highly_variable_lcpm, trace="none", 
          main="Top 2000 variable Genes")
dev.off() 
# 20 and 100 never mix in the same branch (except for the 2 outliers)
# 20 vs. 100 K difference bigger than HCL vs. HD

# scatterplot of pairs of covariates:
# pairs(highly_variable_lcpm)
# pairs(dd$counts)

# Correlation plots between samples's expression patterns:
colnames( logcounts ) = group_ALL
logcounts = logcounts[,order(group_ALL)]

pdf("scripts/images/Filippo/Exploratory/Corrrelation.pdf",         # File name
    width = 8, height = 8)          # Paper size
corrplot( cor(logcounts), method="color" )
dev.off() 


# AVERAGE logcounts ACROSS SAMPLESS
avg_logcounts = sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  rowMeans(logcounts[,sel])
})
colnames(avg_logcounts) = levels(group_ALL)

pdf("scripts/images/Filippo/Exploratory/Corrrelation_groups.pdf",         # File name
    width = 8, height = 8)          # Paper size
corrplot( cor(avg_logcounts),
          method = "pie")
dev.off() 


# 20 vs. 100 difference stronger in HD than in HCL
# Clear separation of the samples in the 4 group_ALL

