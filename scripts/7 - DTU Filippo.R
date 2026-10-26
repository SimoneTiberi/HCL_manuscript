setwd("~/Desktop/Emato")
rm(list =ls())
########################################################################################
# count fraction of protein coding transcript reads:
########################################################################################
load("scripts/RData/data_transcript_Filippo.RData")

# % of protein_coding vs. non-protein_coding mRNA counts:
counts = txi$counts
counts_by_group = sapply(levels(group_ALL), function(gr){
  rowSums(counts[, group_ALL == gr])
})
counts_by_group = as.data.frame(counts_by_group)
counts_by_group$transcript = rownames(counts_by_group)

# match TRANSCRIPT name and description:
gff = rtracklayer::import(
  "Kallisto/Homo_sapiens.GRCh38.112.gff3.gz",
  format = "gff3",
  colnames = c("ID", "biotype")
)
# Keep gene records, including noncoding genes and pseudogenes.
tx = gff[grepl("^transcript:", gff$ID)]
# annotation:
tx_annotation = unique(data.frame(
  ensembl_transcript_id = sub("^transcript:", "", as.character(tx$ID)),
  transcript_biotype = as.character(tx$biotype),
  stringsAsFactors = FALSE
))
transcript_ids = sub(
  "\\.[0-9]+$", "",
  as.character(counts_by_group$transcript)
)

# Select annotations for the requested transcripts.
G_list = tx_annotation[
  tx_annotation$ensembl_transcript_id %in% transcript_ids,
  ,
  drop = FALSE
]

colnames(G_list) = c("transcript", "transcript_biotype")
counts_by_group = merge(counts_by_group, G_list, by = "transcript")

protein_coding = list()
for(gr in levels(group_ALL)){
  counts_split = split(counts_by_group[,gr], counts_by_group$transcript_biotype)
  counts_split = sapply(counts_split, sum)
  protein_coding[[gr]] = counts_split / sum(counts_split)
}
protein_coding = do.call(cbind, protein_coding)
colSums(protein_coding)

protein_coding["protein_coding",]

DF = data.frame(protein_coding = protein_coding["protein_coding",], 
                group = colnames(protein_coding) )
library(ggplot2)
p = ggplot(data=DF, aes(x=group, y=protein_coding)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal() + labs(y = "Percentage of protein coding transcripts")
p
ggsave("scripts/Images/Filippo/DTU/Protein_coding_transcripts.pdf",
       width = 8,
       height = 8)


protein_coding = protein_coding[ order(rowSums(protein_coding), decreasing= TRUE), ]
DF = round(100 * protein_coding, 2)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/Transcript_type.xlsx")
rm(DF)

#> protein_coding["protein_coding",]
# HCL_100K  HCL_20K  HD_100K   HD_20K 
# 66.2     74.1     59.8     71.9 

#> protein_coding["protein_coding",]
# ET        HD        PV 
# 0.5370792 0.5332337 0.4847119 
# 53.7      53.3      48.5 
########################################################################################
# load gene names
########################################################################################
rm(list =ls())
# Replace with the actual path to your GFF file.
gff = rtracklayer::import(
  "Kallisto/Homo_sapiens.GRCh38.112.gff3.gz",
  format = "gff3",
  colnames = c("ID", "Name", "description")
)
# Keep gene records, including noncoding genes and pseudogenes.
g = gff[grepl("^gene:", gff$ID)]
# annotation:
annotation = unique(data.frame(
  ensembl_gene_id = sub("^gene:", "", as.character(g$ID)),
  gene_name = as.character(g$Name),
  description = as.character(g$description),
  stringsAsFactors = FALSE
))

save(annotation, file = "annotation.RData")

########################################################################################
# load transcript biotype
########################################################################################
rm(list =ls())
gff_file = path.expand(
  "~/Desktop/Emato/Kallisto/Homo_sapiens.GRCh38.112.gff3.gz"
)

normalize_transcript = function(x) {
  sub("\\.[0-9]+$", "", sub("^transcript:", "", as.character(x)))
}

gff_attribute = function(attributes, key) {
  positions = regexpr(
    paste0("(?:^|;)", key, "=\\K[^;]*"), attributes, perl = TRUE
  )
  values = rep(NA_character_, length(attributes))
  values[positions > 0L] = regmatches(attributes, positions)
  values
}

read_transcript_biotypes = function(path) {
  con = gzfile(path, open = "rt")
  on.exit(close(con))
  chunks = list()
  repeat {
    lines = readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    if (any(lines == "##FASTA")) {
      lines = lines[seq_len(match("##FASTA", lines) - 1L)]
      last_chunk = TRUE
    } else {
      last_chunk = FALSE
    }
    # Include tutti i tipi di trascritto, anche lncRNA, pseudogeni, ecc.
    # Il biotipo viene dall'attributo biotype del TRASCRITTO, non del gene
    # e non dalla terza colonna del GFF (mRNA, lnc_RNA, ...).
    lines = lines[!startsWith(lines, "#")]
    lines = grep("(?:\\t|;)ID=transcript:", lines, value = TRUE, perl = TRUE)
    if (length(lines)) {
      fields = strsplit(lines, "\t", fixed = TRUE)
      if (any(lengths(fields) != 9L)) stop("Riga GFF non valida.")
      attributes = vapply(fields, `[`, character(1), 9L)
      ids = gff_attribute(attributes, "transcript_id")
      if (anyNA(ids)) stop("Attributo transcript_id mancante nel GFF.")
      chunks[[length(chunks) + 1L]] = data.frame(
        transcript = normalize_transcript(ids),
        transcript_biotype = gff_attribute(attributes, "biotype"),
        stringsAsFactors = FALSE
      )
    }
    if (last_chunk) break
  }
  if (!length(chunks)) stop("Nessun trascritto trovato nel GFF.")
  annotation = unique(do.call(rbind, chunks))
  if (anyNA(annotation$transcript_biotype) ||
      any(!nzchar(annotation$transcript_biotype))) {
    stop("Uno o piu' trascritti nel GFF non hanno biotype.")
  }
  if (anyDuplicated(annotation$transcript)) {
    stop("Uno stesso ID di trascritto ha biotipi discordanti nel GFF.")
  }
  annotation
}

transcript_annotation = read_transcript_biotypes(gff_file)
save(transcript_annotation, file = "transcript_annotation.RData")

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

counts = txi$counts[,sel]
# remove 0 count transcripts:
dim(counts)
counts = counts[rowSums(counts) > 0, ]
dim(counts)
colnames(counts) = samples
rm(txi)

load("scripts/RData/tx2gene.RData")
matches = match( rownames(counts), tx2gene$transcript)
gene_id = as.character( tx2gene$gene[matches] )

# remove NAs in matches:
counts = counts[!is.na(matches),]
gene_id = gene_id[!is.na(matches)]

# keep transcripts with >= 10 counts per condition (across all samples):
counts_split = split(counts, gene_id)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = group %in% levels(group)[1]
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
mean(sel)
genes_kept = unique(names(sel[sel]))
length(genes_kept)
# 23018
sel = gene_id %in% genes_kept
counts = counts[sel,]
gene_id = gene_id[sel]

dim(counts); length(gene_id)

# DEXSeq:
library(DEXSeq)
set.seed(169612)

# Load truth table:
design = data.frame( row.names = samples,
                     condition = group )

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData = design,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(counts),
                    groupID = gene_id, 
                    transcripts = rownames(counts))

library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd = testForDEU(dxd, reducedModel = ~sample + exon, BPPARAM = BPPARAM)

res = DEXSeqResults(dxd, independentFiltering = FALSE)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

sel_a = group %in% "HCL_20K"
transcripts = rownames(counts)
RES = lapply(res_gene$gene, function(gene){
  sel = which(gene_id  == gene)
  sums_A = rowSums(counts[sel, sel_a])
  sums_B = rowSums(counts[sel, !sel_a])
  
  data.frame(gene = gene, 
             transcript = transcripts[sel],
             HCL_20K = sums_A/sum(sums_A), 
             HCL_100K = sums_B/sum(sums_B))
})
RES = do.call(rbind, RES)

res_gene = merge(res_gene, RES, by = "gene")

hist(res_gene$qval)
sum(res_gene$qval < 0.01)
# 2 k

save(res, res_gene,
     file = "scripts/RData/Filippo_DTU_HCL_100_vs_20.RData")

# match gene name and description:
# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids = sub("\\.[0-9]+$", "", as.character(res_gene$gene))

load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes = setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res_gene); dim(G_list)
res_gene = merge(res_gene, G_list, by = "gene", all.x = TRUE)
dim(res_gene)

### Save results:
colnames(res_gene)[2] = "FDR"

res_gene = res_gene[order(res_gene$FDR), ]
rownames(res_gene) = NULL
head(res_gene)

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)

res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene = res_gene[,c(2,3,1,8,4:7)]

# match TRANSCRIPT name and description:
load("transcript_annotation.RData")
idx = match(res_gene$transcript,
             transcript_annotation$transcript)
res_gene$transcript_biotype = transcript_annotation$transcript_biotype[idx]

res_gene = res_gene[,c(1:4, 9, 5:8)]

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
openxlsx::write.xlsx(res_gene, file = "scripts/results/Filippo/Diff_Splic_HCL_100_vs_20.xlsx")

# Volcano plot:
RES_trans = data.frame(transcript = res$featureID, PValue = res$pvalue)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene$logFC = log2( (0.01+res_gene[,7])/(0.01+res_gene[,6]) )
res_gene$logFC

library(EnhancedVolcano)
EnhancedVolcano(res_gene,
                lab = res_gene$transcript,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DTU/Volcano_HCL_100_vs_20.pdf",
       width = 8,
       height = 8)


# % of times dominant transcript is protein_coding
res_gene$transcript_biotype[is.na(res_gene$transcript_biotype)] = "NA"
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = t(sapply(res_gene_split, function(x){
  c(x[ which.max(x[,2]), 1], x[ which.max(x[,3]), 1])
}))
res_gene_split_2 = do.call(rbind, res_gene_split_2)
colnames(res_gene_split_2) = colnames(res_gene)[6:7]
colMeans(res_gene_split_2 == "protein_coding")
round(100 *colMeans(res_gene_split_2 == "protein_coding"),1)
# HCL_20K  HCL_100K 
# 0.8385878 0.7796841 
# HCL_20K HCL_100K 
# 83.9     78.0 

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,2]) != which.max(x[,3])
}))
mean(res_gene_split_2)
# 0.3291752

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
source("scripts/7 - DTU function plot.R")
ggp = PLOT_props(res_gene$gene[1], res_gene)
ggsave("scripts/Images/Filippo/DTU/Proportions_HCL_100_vs_20.pdf",
       width = 8,
       height = 8)


########################################################################################
# HD_100K vs. HD_20K
########################################################################################

load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HD_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

counts = txi$counts[,sel]
# remove 0 count transcripts:
dim(counts)
counts = counts[rowSums(counts) > 0, ]
dim(counts)
colnames(counts) = samples
rm(txi)

load("scripts/RData/tx2gene.RData")
matches = match( rownames(counts), tx2gene$transcript)
gene_id = as.character( tx2gene$gene[matches] )

# remove NAs in matches:
counts = counts[!is.na(matches),]
gene_id = gene_id[!is.na(matches)]

# keep transcripts with >= 10 counts per condition (across all samples):
counts_split = split(counts, gene_id)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = group %in% levels(group)[1]
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
mean(sel)
genes_kept = unique(names(sel[sel]))
length(genes_kept)
# 23018
sel = gene_id %in% genes_kept
counts = counts[sel,]
gene_id = gene_id[sel]

dim(counts); length(gene_id)

# DEXSeq:
library(DEXSeq)
set.seed(169612)

# Load truth table:
design = data.frame( row.names = samples,
                     condition = group )

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData = design,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(counts),
                    groupID = gene_id, 
                    transcripts = rownames(counts))

library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd = testForDEU(dxd, reducedModel = ~sample + exon, BPPARAM = BPPARAM)

res = DEXSeqResults(dxd, independentFiltering = FALSE)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

sel_a = group %in% "HD_20K"
transcripts = rownames(counts)
RES = lapply(res_gene$gene, function(gene){
  sel = which(gene_id  == gene)
  sums_A = rowSums(counts[sel, sel_a])
  sums_B = rowSums(counts[sel, !sel_a])
  
  data.frame(gene = gene, 
             transcript = transcripts[sel],
             HD_20K = sums_A/sum(sums_A), 
             HD_100K = sums_B/sum(sums_B))
})
RES = do.call(rbind, RES)

res_gene = merge(res_gene, RES, by = "gene")

hist(res_gene$qval)
sum(res_gene$qval < 0.01)
# 1.1 k

save(res, res_gene,
     file = "scripts/RData/Filippo_DTU_HD_100_vs_20.RData")

# match gene name and description:
# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids = sub("\\.[0-9]+$", "", as.character(res_gene$gene))

load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes = setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res_gene); dim(G_list)
res_gene = merge(res_gene, G_list, by = "gene", all.x = TRUE)
dim(res_gene)

### Save results:
colnames(res_gene)[2] = "FDR"

res_gene = res_gene[order(res_gene$FDR), ]
rownames(res_gene) = NULL
head(res_gene)

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)

res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene = res_gene[,c(2,3,1,8,4:7)]

# match TRANSCRIPT name and description:
load("transcript_annotation.RData")
idx = match(res_gene$transcript,
            transcript_annotation$transcript)
res_gene$transcript_biotype = transcript_annotation$transcript_biotype[idx]

res_gene = res_gene[,c(1:4, 9, 5:8)]

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
openxlsx::write.xlsx(res_gene, file = "scripts/results/Filippo/Diff_Splic_HD_100_vs_20.xlsx")

# Volcano plot:
RES_trans = data.frame(transcript = res$featureID, PValue = res$pvalue)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene$logFC = log2( (0.01+res_gene[,7])/(0.01+res_gene[,6]) )
res_gene$logFC

library(EnhancedVolcano)
EnhancedVolcano(res_gene,
                lab = res_gene$transcript,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DTU/Volcano_HD_100_vs_20.pdf",
       width = 8,
       height = 8)


# % of times dominant transcript is protein_coding
res_gene$transcript_biotype[is.na(res_gene$transcript_biotype)] = "NA"
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = t(sapply(res_gene_split, function(x){
  c(x[ which.max(x[,2]), 1], x[ which.max(x[,3]), 1])
}))
res_gene_split_2 = do.call(rbind, res_gene_split_2)
colnames(res_gene_split_2) = colnames(res_gene)[6:7]
colMeans(res_gene_split_2 == "protein_coding")
round(100 *colMeans(res_gene_split_2 == "protein_coding"),1)
# HD_20K   HD_100K 
# 0.8160572 0.7134488 
# 81.6    71.3 

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,2]) != which.max(x[,3])
}))
mean(res_gene_split_2)
# 0.4148208

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
source("scripts/7 - DTU function plot.R")
ggp = PLOT_props(res_gene$gene[1], res_gene)
ggsave("scripts/Images/Filippo/DTU/Proportions_HD_100_vs_20.pdf",
       width = 8,
       height = 8)


########################################################################################
# HD_100K vs. HCL_100K
########################################################################################

load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_100K", "HCL_100K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

counts = txi$counts[,sel]
# remove 0 count transcripts:
dim(counts)
counts = counts[rowSums(counts) > 0, ]
dim(counts)
colnames(counts) = samples
rm(txi)

load("scripts/RData/tx2gene.RData")
matches = match( rownames(counts), tx2gene$transcript)
gene_id = as.character( tx2gene$gene[matches] )

# remove NAs in matches:
counts = counts[!is.na(matches),]
gene_id = gene_id[!is.na(matches)]

# keep transcripts with >= 10 counts per condition (across all samples):
counts_split = split(counts, gene_id)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = group %in% levels(group)[1]
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
mean(sel)
genes_kept = unique(names(sel[sel]))
length(genes_kept)
# 23018
sel = gene_id %in% genes_kept
counts = counts[sel,]
gene_id = gene_id[sel]

dim(counts); length(gene_id)

# DEXSeq:
library(DEXSeq)
set.seed(169612)

# Load truth table:
design = data.frame( row.names = samples,
                     condition = group )

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData = design,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(counts),
                    groupID = gene_id, 
                    transcripts = rownames(counts))

library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd = testForDEU(dxd, reducedModel = ~sample + exon, BPPARAM = BPPARAM)

res = DEXSeqResults(dxd, independentFiltering = FALSE)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

sel_a = group %in% "HD_100K"
transcripts = rownames(counts)
RES = lapply(res_gene$gene, function(gene){
  sel = which(gene_id  == gene)
  sums_A = rowSums(counts[sel, sel_a])
  sums_B = rowSums(counts[sel, !sel_a])
  
  data.frame(gene = gene, 
             transcript = transcripts[sel],
             HD_100K = sums_A/sum(sums_A), 
             HCL_100K = sums_B/sum(sums_B))
})
RES = do.call(rbind, RES)

res_gene = merge(res_gene, RES, by = "gene")

hist(res_gene$qval)
sum(res_gene$qval < 0.01)
# 0.2 k

save(res, res_gene,
     file = "scripts/RData/Filippo_DTU_100_HD_vs_HCL.RData")

# match gene name and description:
# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids = sub("\\.[0-9]+$", "", as.character(res_gene$gene))

load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes = setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res_gene); dim(G_list)
res_gene = merge(res_gene, G_list, by = "gene", all.x = TRUE)
dim(res_gene)

### Save results:
colnames(res_gene)[2] = "FDR"

res_gene = res_gene[order(res_gene$FDR), ]
rownames(res_gene) = NULL
head(res_gene)

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)

res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene = res_gene[,c(2,3,1,8,4:7)]

# match TRANSCRIPT name and description:
load("transcript_annotation.RData")
idx = match(res_gene$transcript,
            transcript_annotation$transcript)
res_gene$transcript_biotype = transcript_annotation$transcript_biotype[idx]

res_gene = res_gene[,c(1:4, 9, 5:8)]

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
openxlsx::write.xlsx(res_gene, file = "scripts/results/Filippo/Diff_Splic_100_HD_vs_HCL.xlsx")

# Volcano plot:
RES_trans = data.frame(transcript = res$featureID, PValue = res$pvalue)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene$logFC = log2( (0.01+res_gene[,7])/(0.01+res_gene[,6]) )
res_gene$logFC

library(EnhancedVolcano)
EnhancedVolcano(res_gene,
                lab = res_gene$transcript,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DTU/Volcano_100_HD_vs_HCL.pdf",
       width = 8,
       height = 8)


# % of times dominant transcript is protein_coding
res_gene$transcript_biotype[is.na(res_gene$transcript_biotype)] = "NA"
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = t(sapply(res_gene_split, function(x){
  c(x[ which.max(x[,2]), 1], x[ which.max(x[,3]), 1])
}))
res_gene_split_2 = do.call(rbind, res_gene_split_2)
colnames(res_gene_split_2) = colnames(res_gene)[6:7]
colMeans(res_gene_split_2 == "protein_coding")
round(100 *colMeans(res_gene_split_2 == "protein_coding"),1)
# HD_100K   HCL_100K 
# 0.7146889 0.7786452 
# 71.5      77.9 


# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,2]) != which.max(x[,3])
}))
mean(res_gene_split_2)
# 0.3272061

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
source("scripts/7 - DTU function plot.R")
ggp = PLOT_props(res_gene$gene[1], res_gene)
ggsave("scripts/Images/Filippo/DTU/Proportions_100_HD_vs_HCL.pdf",
       width = 8,
       height = 8)

########################################################################################
# HD_20K vs. HCL_20K
########################################################################################

load("scripts/RData/data_transcript_Filippo.RData")
ls()

# select group to test:
group_ALL
sel = group_ALL %in% c("HD_20K", "HCL_20K")
samples = samples[sel]; samples
group = factor(group_ALL[sel]); group
rm(group_ALL); rm(group_20); rm(group_HD)

counts = txi$counts[,sel]
# remove 0 count transcripts:
dim(counts)
counts = counts[rowSums(counts) > 0, ]
dim(counts)
colnames(counts) = samples
rm(txi)

load("scripts/RData/tx2gene.RData")
matches = match( rownames(counts), tx2gene$transcript)
gene_id = as.character( tx2gene$gene[matches] )

# remove NAs in matches:
counts = counts[!is.na(matches),]
gene_id = gene_id[!is.na(matches)]

# keep transcripts with >= 10 counts per condition (across all samples):
counts_split = split(counts, gene_id)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = group %in% levels(group)[1]
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
mean(sel)
genes_kept = unique(names(sel[sel]))
length(genes_kept)
# 23018
sel = gene_id %in% genes_kept
counts = counts[sel,]
gene_id = gene_id[sel]

dim(counts); length(gene_id)

# DEXSeq:
library(DEXSeq)
set.seed(169612)

# Load truth table:
design = data.frame( row.names = samples,
                     condition = group )

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData = design,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(counts),
                    groupID = gene_id, 
                    transcripts = rownames(counts))

library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)

dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd = testForDEU(dxd, reducedModel = ~sample + exon, BPPARAM = BPPARAM)

res = DEXSeqResults(dxd, independentFiltering = FALSE)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

sel_a = group %in% "HD_20K"
transcripts = rownames(counts)
RES = lapply(res_gene$gene, function(gene){
  sel = which(gene_id  == gene)
  sums_A = rowSums(counts[sel, sel_a])
  sums_B = rowSums(counts[sel, !sel_a])
  
  data.frame(gene = gene, 
             transcript = transcripts[sel],
             HD_20K = sums_A/sum(sums_A), 
             HCL_20K = sums_B/sum(sums_B))
})
RES = do.call(rbind, RES)

res_gene = merge(res_gene, RES, by = "gene")

hist(res_gene$qval)
sum(res_gene$qval < 0.01)
# 88

save(res, res_gene,
     file = "scripts/RData/Filippo_DTU_20_HD_vs_HCL.RData")

# match gene name and description:
# Remove version suffixes, if present: ENSG00000141510.18 -> ENSG00000141510
gene_ids = sub("\\.[0-9]+$", "", as.character(res_gene$gene))

load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% gene_ids,
  ,
  drop = FALSE
]

# Check whether any requested IDs were absent.
missing_genes = setdiff(gene_ids, annotation$ensembl_gene_id)
length(missing_genes)

colnames(G_list) = c("gene", "gene_name", "description")
dim(res_gene); dim(G_list)
res_gene = merge(res_gene, G_list, by = "gene", all.x = TRUE)
dim(res_gene)

### Save results:
colnames(res_gene)[2] = "FDR"

res_gene = res_gene[order(res_gene$FDR), ]
rownames(res_gene) = NULL
head(res_gene)

# ADD transcript FDR 
RES_trans = data.frame(transcript = res$featureID, FDR_transcript = res$padj)

res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene = res_gene[,c(2,3,1,8,4:7)]

# match TRANSCRIPT name and description:
load("transcript_annotation.RData")
idx = match(res_gene$transcript,
            transcript_annotation$transcript)
res_gene$transcript_biotype = transcript_annotation$transcript_biotype[idx]

res_gene = res_gene[,c(1:4, 9, 5:8)]

openxlsx::write.xlsx(res_gene, file = "scripts/results/Filippo/Diff_Splic_20_HD_vs_HCL.xlsx")

# Volcano plot:
RES_trans = data.frame(transcript = res$featureID, PValue = res$pvalue)
res_gene = merge(res_gene, RES_trans, by = "transcript", all.x = TRUE)
res_gene$logFC = log2( (0.01+res_gene[,7])/(0.01+res_gene[,6]) )
res_gene$logFC

library(EnhancedVolcano)
EnhancedVolcano(res_gene,
                lab = res_gene$transcript,
                x = 'logFC',
                y = 'PValue')
ggsave("scripts/Images/Filippo/DTU/Volcano_20_HD_vs_HCL.pdf",
       width = 8,
       height = 8)


# % of times dominant transcript is protein_coding
res_gene$transcript_biotype[is.na(res_gene$transcript_biotype)] = "NA"
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = t(sapply(res_gene_split, function(x){
  c(x[ which.max(x[,2]), 1], x[ which.max(x[,3]), 1])
}))
res_gene_split_2 = do.call(rbind, res_gene_split_2)
colnames(res_gene_split_2) = colnames(res_gene)[6:7]
colMeans(res_gene_split_2 == "protein_coding")
round(100 *colMeans(res_gene_split_2 == "protein_coding"),1)
# HCL_20K  HCL_100K 
# 0.8169216 0.8415880 
# 81.7    84.2 

# % of times dominant transcript changes between conditions:
res_gene_split = split(res_gene[5:7], res_gene$gene)
res_gene_split_2 = unlist(sapply(res_gene_split, function(x){
  which.max(x[,2]) != which.max(x[,3])
}))
mean(res_gene_split_2)
# 0.2791303

res_gene = res_gene[order(res_gene$FDR, res_gene$gene, res_gene$FDR_transcript), ]
source("scripts/7 - DTU function plot.R")
ggp = PLOT_props(res_gene$gene[1], res_gene)
ggsave("scripts/Images/Filippo/DTU/Proportions_20_HD_vs_HCL.pdf",
       width = 8,
       height = 8)
