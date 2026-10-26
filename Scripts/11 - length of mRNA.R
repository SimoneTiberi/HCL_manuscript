setwd("~/Desktop/Emato")
########################################################################################
# Dorian
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Dorian.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi$abundance) > 0
txi$length = txi$length[sel,]
txi$abundance = txi$abundance[sel,]
txi$counts = txi$counts[sel,]

sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  mean(colSums(txi$counts[,sel]))
})/10^6
# average million reads per sample in each group:
#  HCL_100K   HCL_20K   HD_100K    HD_20K 
# 1.0283360 0.9744331 0.9243600 0.8133555 

# compute length:
len = rowMeans(txi$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi$abundance; colSums(TPM)

pdf("scripts/Images/Dorian/Exploratory/tr_length.pdf",
    width = 8,
    height = 8)
par(mfrow = c(2,2))
name = "HD"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "AML_low"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "AML_int"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "AML_high"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
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
#[1,]               1     0.30    0.18    0.29   0.20
#[2,]               2     0.36    0.25    0.33   0.26
#[3,]               3     0.43    0.33    0.38   0.34
#[4,]               4     0.52    0.42    0.47   0.43

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


########################################################################################
# Filippo
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Filippo.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi$abundance) > 0
txi$length = txi$length[sel,]
txi$abundance = txi$abundance[sel,]
txi$counts = txi$counts[sel,]

sapply( levels(group_ALL), function(gr){
  sel = group_ALL == gr
  mean(colSums(txi$counts[,sel]))
})/10^6
# average million reads per sample in each group:
#  HCL_100K   HCL_20K   HD_100K    HD_20K 
# 1.7201410 2.4812875 0.7346279 1.5755525 


# compute length:
len = rowMeans(txi$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi$abundance; colSums(TPM)

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

pdf("scripts/Images/Filippo/Exploratory/tr_length.pdf",
    width = 8,
    height = 8)
par(mfrow = c(2,2))
name = "HCL_100K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "HCL_20K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "HD_100K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
name = "HD_20K"
sel = group_ALL == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name)
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
#[1,]               1     0.30    0.18    0.29   0.20
#[2,]               2     0.36    0.25    0.33   0.26
#[3,]               3     0.43    0.33    0.38   0.34
#[4,]               4     0.52    0.42    0.47   0.43

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

########################################################################################
# Ghazal
########################################################################################
rm(list =ls())
load("scripts/RData/data_transcript_Ghazal.RData")
ls()

# remove transcripts which are ALWAYS zero:
sel = rowSums(txi$abundance) > 0
txi$length = txi$length[sel,]
txi$abundance = txi$abundance[sel,]
txi$counts = txi$counts[sel,]

sapply( levels(group), function(gr){
  sel = group == gr
  mean(colSums(txi$counts[,sel]))
})/10^6
# average million reads per sample in each group:
#    ET        HD        PV 
# 0.5024640 0.6408453 0.5759614 


# compute length:
len = rowMeans(txi$length)
quant = quantile(len, c(1/4, 1/2, 3/4)); quant
len_groups = ifelse(len < quant[1],
                    1, ifelse(len < quant[2], 2,
                              ifelse(len < quant[3], 3, 4)))

TPM = txi$abundance; colSums(TPM)

pdf("scripts/Images/Ghazal/Exploratory/tr_length.pdf",
    width = 8,
    height = 4)
par(mfrow = c(1,3))
name = "HD"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name,
         col = "blue")
name = "ET"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name,
         col = "darkgreen")
name = "PV"
sel = group == name
boxplot( c(TPM[,sel]) ~ rep(len_groups, sum(sel)), outline = FALSE,
         xlab = "Transcript length quartiles", ylab = "TPMs", main = name,
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
#Length Quartile   ET   HD   PV
#               1 0.22 0.28 0.22
#               2 0.26 0.32 0.25
#               3 0.31 0.36 0.29
#               4 0.40 0.44 0.38

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
