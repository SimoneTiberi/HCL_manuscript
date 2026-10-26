setwd("~/Desktop/Emato")

rm(list = ls())

# make gene-transript file:
library(txdbmaker)
gffFile = "Kallisto/Homo_sapiens.GRCh38.112.gff3.gz"
txdb = makeTxDbFromGFF(gffFile,
                        dataSource = "Ensembl",
                        organism="Homo sapiens")
txdb = as.list(txdb)
tx2gene = data.frame(transcript = txdb$transcripts$tx_name, gene = txdb$genes$gene_id)
save(tx2gene, file = "scripts/RData/tx2gene.RData")

########################################################################################
# Filippo:
########################################################################################
setwd("~/Desktop/Emato")
rm(list =ls())
load("scripts/RData/tx2gene.RData")

samples = list.files("Kallisto/Filippo/")
files = paste0("Kallisto/Filippo/", samples, "/abundance.h5")
files[ ! file.exists(files) ] 
length(samples)

# import quantification files with tximport:
library(tximport)
txi = tximport(files = files, txOut = TRUE, type = "kallisto")
str(txi)

# remove .1 at the end of transcript names:
tmp = strsplit(rownames(txi$counts), ".",  fixed = TRUE)
tmp = sapply(tmp, function(x){
  x[[1]]
})
rownames(txi$abundance) = tmp
rownames(txi$counts) = tmp
rownames(txi$length) = tmp

# summarize info at the Gene-level:
txi_gene = summarizeToGene(txi, tx2gene = tx2gene)
str(txi_gene)
Gene_counts = txi_gene$counts
dim(Gene_counts)
# 36 k genes and 33 samples
head(Gene_counts, 20) # prints the first 20 rows of counts
# The length matrix contains the average transcript length for each gene which can be used as an offset for gene-level analysis.
# million counts per sample:
hist(colSums(Gene_counts)/10^6)

sort(round(colSums(Gene_counts)/10^6, 1))
# 0.0 0.0 0.2 0.3 0.4 0.6 0.9 0.9 0.9 1.0 1.0 1.2 1.2 1.3
# 1.3 1.3 1.5 1.6 1.6 1.7 1.8 1.8 1.9 1.9 1.9 2.2 2.3 2.5 
# 2.5 2.5 2.7 3.1 3.1

n_samples = length(samples); n_samples
group_20 = rep("100K", n_samples)
group_20[grep("20K", samples)] = "20K"
group_20 = factor(group_20); table(group_20)

group_HD = rep("HCL", n_samples)
group_HD[grep("HD", samples)] = "HD"
group_HD = factor(group_HD); table(group_HD)

group_ALL = factor(paste(group_HD, group_20, sep = "_")); table(group_ALL)

save(txi_gene, samples, group_20, group_HD, group_ALL, 
     file = "scripts/RData/data_ALL_Filippo.RData")

save(txi, samples, group_20, group_HD, group_ALL, 
     file = "scripts/RData/data_ALL_transcript_Filippo.RData")

# remove HD_47_100K and HD_47_20K
# identified as outliers (see exploratory plots)
sel = !(samples %in% c("HD_47_100K", "HD_47_20K", "BO_LYMPHO_056_100K","BO_LYMPHO_056_20K"))
sum(sel); sum(!sel)
samples = samples[sel]
group_ALL = group_ALL[sel]
group_HD = group_HD[sel]
group_20 = group_20[sel]
txi$abundance = txi$abundance[,sel]
txi$counts = txi$counts[,sel]
txi$length = txi$length[,sel]
txi_gene$abundance = txi_gene$abundance[,sel]
txi_gene$counts = txi_gene$counts[,sel]
txi_gene$length = txi_gene$length[,sel]

save(txi_gene, samples, group_20, group_HD, group_ALL, 
     file = "scripts/RData/data_Filippo.RData")

save(txi, samples, group_20, group_HD, group_ALL, 
     file = "scripts/RData/data_transcript_Filippo.RData")
