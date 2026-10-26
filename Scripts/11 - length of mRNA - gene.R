setwd("~/Desktop/Emato")

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
save(annotation, file = "gene_annotation.RData")

########################################################################################
# Dorian
########################################################################################
rm(list =ls())
load("scripts/RData/data_Dorian.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi_gene$abundance) > 0
txi_gene$length = txi_gene$length[sel,]
txi_gene$abundance = txi_gene$abundance[sel,]
txi_gene$counts = txi_gene$counts[sel,]

sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  mean(colSums(txi_gene$counts[,sel]))
})/10^6
# average million reads per sample in each group:


# compute length:
len = rowMeans(txi_gene$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi_gene$abundance; colSums(TPM)

pdf("scripts/Images/Dorian/Exploratory/gene_length.pdf",
    width = 8,
    height = 8)
par(mfrow = c(2,2))
name = "HD"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "AML_low"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "AML_int"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "AML_high"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
dev.off()

# NON-zero transcripts:
X = sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    mean(xx[ rep(len_groups, sum(sel)) == len] > 0)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#     Length Quartile HCL_100K HCL_20K HD_100K HD_20K
#[1,]               1     0.34    0.15    0.32   0.17
#[2,]               2     0.67    0.43    0.65   0.47
#[3,]               3     0.88    0.71    0.86   0.75
#[4,]               4     0.96    0.84    0.94   0.88

# 0.95 quantile:
X = sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    quantile(xx[ rep(len_groups, sum(sel)) == len], 0.95)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#    Length Quartile HCL_100K HCL_20K HD_100K HD_20K
#95%               1     5.65    2.15   11.82   2.99
#95%               2     4.29    2.92    7.63   3.22
#95%               3     4.03    4.00    5.50   3.55
#95%               4     2.97    2.73    4.15   2.53



# Save gene names and lengths:
DF = data.frame(gene = names(len_groups),
                TPM = rowMeans(txi_gene$abundance),
                length = len,
                length_quartile = len_groups)

DF = DF[ order(DF$length),]

library('biomaRt')
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list = getBM(filters= "ensembl_gene_id", 
               attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
               values=DF$gene,
               mart = mart)
colnames(G_list) = c("gene", "gene_name", "description")
dim(DF); dim(G_list)
DF = merge(DF, G_list, by = "gene", all.x = TRUE)
dim(DF)
DF = DF[ , c(1,5,2:4, 6)]

DF = DF[ order(DF$length),]
head(DF)
tail(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Dorian/Gene_length.xlsx")


########################################################################################
# Filippo
########################################################################################
rm(list =ls())
load("scripts/RData/data_Filippo.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi_gene$abundance) > 0
txi_gene$length = txi_gene$length[sel,]
txi_gene$abundance = txi_gene$abundance[sel,]
txi_gene$counts = txi_gene$counts[sel,]

sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  mean(colSums(txi_gene$counts[,sel]))
})/10^6
# average million reads per sample in each group:
#  HCL_100K   HCL_20K   HD_100K    HD_20K 
# 1.7201410 2.4812875 0.7346279 1.5755525 

# compute length:
len = rowMeans(txi_gene$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi_gene$abundance; colSums(TPM)

if(FALSE){
  library(ggplot2)
  name = "HCL_100K"
  sel = group_ALL == name
  DF = data.frame(TPM = c(TPM[,sel]),
                  length_quartile = factor(rep(len_groups, sum(sel))))
  ggplot(DF, aes(x=length_quartile, y=TPM)) +
    geom_violin(trim=TRUE, drop = TRUE,
                adjust = 1,
                bounds = c(0, 0.1))
  
  name = "HCL_20K"
  sel = group_ALL == name
  DF = data.frame(TPM = c(TPM[,sel]),
                  length_quartile = factor(rep(len_groups, sum(sel))))
  ggplot(DF, aes(x=length_quartile, y=TPM)) +
    geom_violin(trim=TRUE, drop = TRUE,
                adjust = 1,
                bounds = c(0, 0.1))
  
  name = "HD_100K"
  sel = group_ALL == name
  DF = data.frame(TPM = c(TPM[,sel]),
                  length_quartile = factor(rep(len_groups, sum(sel))))
  ggplot(DF, aes(x=length_quartile, y=TPM)) +
    geom_violin(trim=TRUE, drop = TRUE,
                adjust = 1,
                bounds = c(0, 0.1))
  
  name = "HD_20K"
  sel = group_ALL == name
  DF = data.frame(TPM = c(TPM[,sel]),
                  length_quartile = factor(rep(len_groups, sum(sel))))
  ggplot(DF, aes(x=length_quartile, y=TPM)) +
    geom_violin(trim=TRUE, drop = TRUE,
                adjust = 1,
                bounds = c(0, 0.1))
}

pdf("scripts/Images/Filippo/Exploratory/gene_length.pdf",
    width = 8,
    height = 8)
par(mfrow = c(2,2))
name = "HCL_100K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "HCL_20K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "HD_100K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
name = "HD_20K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name)
dev.off()

# NON-zero transcripts:
X = sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    mean(xx[ rep(len_groups, sum(sel)) == len] > 0)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#     Length Quartile HCL_100K HCL_20K HD_100K HD_20K
#[1,]               1     0.34    0.15    0.32   0.17
#[2,]               2     0.67    0.43    0.65   0.47
#[3,]               3     0.88    0.71    0.86   0.75
#[4,]               4     0.96    0.84    0.94   0.88

# 0.95 quantile:
X = sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    quantile(xx[ rep(len_groups, sum(sel)) == len], 0.95)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#    Length Quartile HCL_100K HCL_20K HD_100K HD_20K
#95%               1     5.65    2.15   11.82   2.99
#95%               2     4.29    2.92    7.63   3.22
#95%               3     4.03    4.00    5.50   3.55
#95%               4     2.97    2.73    4.15   2.53

# Save gene names and lengths:
DF = data.frame(gene = names(len_groups),
                TPM = rowMeans(txi_gene$abundance),
                length = len,
                length_quartile = len_groups)

DF = DF[ order(DF$length),]

load("gene_annotation.RData")

# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids <- sub("\\.[0-9]+$", "", as.character(DF$gene))

G_list <- annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes <- setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(DF); dim(G_list)
DF = merge(DF, G_list, by = "gene", all.x = TRUE)
dim(DF)
DF = DF[ , c(1,5,2:4, 6)]

DF = DF[ order(DF$length),]
head(DF)
tail(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/Gene_length.xlsx")


########################################################################################
# Ghazal
########################################################################################
rm(list =ls())
load("scripts/RData/data_Ghazal.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi_gene$abundance) > 0
txi_gene$length = txi_gene$length[sel,]
txi_gene$abundance = txi_gene$abundance[sel,]
txi_gene$counts = txi_gene$counts[sel,]

sapply( levels(group), function(gr){
  sel = group == gr
  mean(colSums(txi_gene$counts[,sel]))
})/10^6
# average million reads per sample in each group:
#    ET        HD        PV 
# 0.5024640 0.6408453 0.5759614 


# compute length:
len = rowMeans(txi_gene$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi_gene$abundance; colSums(TPM)

pdf("scripts/Images/Ghazal/Exploratory/gene_length.pdf",
    width = 8,
    height = 4)
par(mfrow = c(1,3))
name = "HD"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name,
         col = "blue")
name = "ET"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name,
         col = "darkgreen")
name = "PV"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Gene length quartiles", ylab = "TPMs", main = name,
         col = "pink")
dev.off()

# NON-zero transcripts:
X = sapply( levels(group), function(gr){
  sel = group == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    mean(xx[ rep(len_groups, sum(sel)) == len] > 0)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#     Length Quartile   ET   HD   PV
#[1,]               1 0.21 0.31 0.23
#[2,]               2 0.54 0.62 0.52
#[3,]               3 0.75 0.80 0.71
#[4,]               4 0.87 0.88 0.82

# 0.95 quantile:
X = sapply( levels(group), function(gr){
  sel = group == gr
  xx = c(TPM[,sel])
  sapply(1:4, function(len){
    quantile(xx[ rep(len_groups, sum(sel)) == len], 0.95)
  })
})
X = cbind(1:4, X)
colnames(X)[1] = "Length Quartile"
round(X, 2)
#   Length Quartile   ET   HD   PV
#95%               1 5.76 8.46 4.05
#95%               2 4.00 5.17 2.84
#95%               3 2.84 3.75 2.07
#95%               4 2.13 3.03 1.63



# Save gene names and lengths:
DF = data.frame(gene = names(len_groups),
                TPM = rowMeans(txi_gene$abundance),
                length = len,
                length_quartile = len_groups)

DF = DF[ order(DF$length),]

library('biomaRt')
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list = getBM(filters= "ensembl_gene_id", 
               attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
               values=DF$gene,
               mart = mart)
colnames(G_list) = c("gene", "gene_name", "description")
dim(DF); dim(G_list)
DF = merge(DF, G_list, by = "gene", all.x = TRUE)
dim(DF)
DF = DF[ , c(1,5,2:4, 6)]

DF = DF[ order(DF$length),]
head(DF)
tail(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Ghazal/Gene_length.xlsx")
