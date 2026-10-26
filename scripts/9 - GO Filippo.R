setwd("~/Desktop/Emato")
library(ggplot2)
########################################################################################
# HCL_100K vs. HCL_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HCL_100_vs_20.RData")

genes = rownames(res$table)[res$table$FDR<0.01]
length(genes)

# map Ensembl ids to Entrez ids
library(clusterProfiler)
set.seed(169612)
if(FALSE){
  GGO = groupGO(gene     = genes,
                OrgDb    = 'org.Hs.eg.db',
                keyType  = 'ENSEMBL',
                ont      = "CC",
                level    = 3,
                readable = TRUE)
  head(GGO@result[,-5], 100)
}

EGO = enrichGO(gene         = genes,
               OrgDb        = 'org.Hs.eg.db',
               keyType      = 'ENSEMBL',
               ont          = "ALL",
               pvalueCutoff = 1)

head(EGO@result[,-which(colnames(EGO@result) == "geneID")])

sum(EGO@result$p.adjust < 0.01)
# 130

EGO@result = EGO@result[order(EGO@result$p.adjust), ]

library(DOSE)
dotplot(EGO, showCategory=20, x = "p.adjust")
ggsave("scripts/Images/Filippo/ORA/ORA_HCL_100_vs_20.pdf",
       width = 8,
       height = 12)

# emapplot(EGO, showCategory = 10)

DF = EGO@result[,-which(colnames(EGO@result) == "geneID")]
class(DF)
head(DF)

DF = data.frame(ONTOLOGY = DF$ONTOLOGY,
                ID = DF$ID,
                Description = DF$Description, 
                pvalue = DF$pvalue, 
                FDR = DF$p.adjust)
dim(DF)

# add gene ids for each othology
gene = sapply(EGO@result$geneID, strsplit, "/")

gene_unique = data.frame(gene = unique(unlist(gene)))

# save gene names for each onthology.
load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% unique(unlist(gene)),
  ,
  drop = FALSE
]

G_list = G_list[,1:2]
colnames(G_list) = c("gene", "gene_name")
dim(gene_unique); dim(G_list)
gene_unique = merge(gene_unique, G_list, by = "gene", all.x = TRUE)
dim(gene_unique)
head(gene_unique)

gene_name = sapply(gene, function(gg){
  gene_unique$gene_name[match(gg, gene_unique$gene)]
})

DF$gene = gene
DF$gene_name = gene_name

DF = DF[order(DF$pvalue),]
openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/GO_ORA_HCL_100_vs_20.xlsx")

########################################################################################
# HD_100K vs. HD_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_HD_100_vs_20.RData")

genes = rownames(res$table)[res$table$FDR<0.01]
length(genes)

# map Ensembl ids to Entrez ids
library(clusterProfiler)
set.seed(169612)
if(FALSE){
  GGO = groupGO(gene     = genes,
                OrgDb    = 'org.Hs.eg.db',
                keyType  = 'ENSEMBL',
                ont      = "CC",
                level    = 3,
                readable = TRUE)
  head(GGO@result[,-5], 100)
}

EGO = enrichGO(gene         = genes,
               OrgDb         = 'org.Hs.eg.db',
               keyType       = 'ENSEMBL',
               ont          = "ALL",
               pvalueCutoff = 1)
head(EGO@result[,-which(colnames(EGO@result) == "geneID")])

sum(EGO@result$p.adjust < 0.01)
# 845

EGO@result = EGO@result[order(EGO@result$p.adjust), ]

library(DOSE)
dotplot(EGO, showCategory=20, x = "p.adjust")
ggsave("scripts/Images/Filippo/ORA/ORA_HD_100_vs_20.pdf",
       width = 8,
       height = 12)

DF = EGO@result[,-which(colnames(EGO@result) == "geneID")]
class(DF)
head(DF)

DF = data.frame(ONTOLOGY = DF$ONTOLOGY,
                ID = DF$ID,
                Description = DF$Description, 
                pvalue = DF$pvalue, 
                FDR = DF$p.adjust)
dim(DF)

# add gene ids for each othology
gene = sapply(EGO@result$geneID, strsplit, "/")

gene_unique = data.frame(gene = unique(unlist(gene)))

# save gene names for each onthology.
load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% unique(unlist(gene)),
  ,
  drop = FALSE
]

G_list = G_list[,1:2]
colnames(G_list) = c("gene", "gene_name")
dim(gene_unique); dim(G_list)
gene_unique = merge(gene_unique, G_list, by = "gene", all.x = TRUE)
dim(gene_unique)
head(gene_unique)

gene_name = sapply(gene, function(gg){
  gene_unique$gene_name[match(gg, gene_unique$gene)]
})

DF$gene = gene
DF$gene_name = gene_name

DF = DF[order(DF$pvalue),]

flatten_list_cols <- function(df, sep = "/") {
  is_list <- vapply(df, is.list, logical(1))
  df[is_list] <- lapply(df[is_list], function(col)
    vapply(col, function(v) { v <- unlist(v); paste(replace(v, is.na(v), "NA"), collapse = sep) },
           character(1), USE.NAMES = FALSE))
  df
}
DF = flatten_list_cols(DF)

openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/GO_ORA_HD_100_vs_20.xlsx")

########################################################################################
# HD_100K vs. HCL_100K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_100_HD_vs_HCL.RData")

genes = rownames(res$table)[res$table$FDR<0.01]
length(genes)

# map Ensembl ids to Entrez ids
library(clusterProfiler)
set.seed(169612)
if(FALSE){
  GGO = groupGO(gene     = genes,
                OrgDb    = 'org.Hs.eg.db',
                keyType  = 'ENSEMBL',
                ont      = "CC",
                level    = 3,
                readable = TRUE)
  head(GGO@result[,-5], 100)
}

EGO = enrichGO(gene         = genes,
               OrgDb         = 'org.Hs.eg.db',
               keyType       = 'ENSEMBL',
               ont          = "ALL",
               pvalueCutoff = 1)
head(EGO@result[,-which(colnames(EGO@result) == "geneID")])

sum(EGO@result$p.adjust < 0.01)
# 459

EGO@result = EGO@result[order(EGO@result$p.adjust), ]

library(DOSE)
dotplot(EGO, showCategory=20, x = "p.adjust")
ggsave("scripts/Images/Filippo/ORA/ORA_100_HD_vs_HCL.pdf",
       width = 8,
       height = 12)

DF = EGO@result[,-which(colnames(EGO@result) == "geneID")]
class(DF)
head(DF)

DF = data.frame(ONTOLOGY = DF$ONTOLOGY,
                ID = DF$ID,
                Description = DF$Description, 
                pvalue = DF$pvalue, 
                FDR = DF$p.adjust)
dim(DF)


# add gene ids for each othology
gene = sapply(EGO@result$geneID, strsplit, "/")

gene_unique = data.frame(gene = unique(unlist(gene)))

# save gene names for each onthology.
load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% unique(unlist(gene)),
  ,
  drop = FALSE
]

G_list = G_list[,1:2]
colnames(G_list) = c("gene", "gene_name")
dim(gene_unique); dim(G_list)
gene_unique = merge(gene_unique, G_list, by = "gene", all.x = TRUE)
dim(gene_unique)
head(gene_unique)

gene_name = sapply(gene, function(gg){
  gene_unique$gene_name[match(gg, gene_unique$gene)]
})

DF$gene = gene
DF$gene_name = gene_name

DF = DF[order(DF$pvalue),]
openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/GO_ORA_100_HD_vs_HCL.xlsx")


########################################################################################
# HD_20K vs. HCL_20K
########################################################################################
rm(list = ls())
load("scripts/RData/Filippo_DGE_20_HD_vs_HCL.RData")

genes = rownames(res$table)[res$table$FDR<0.01]
length(genes)

library(clusterProfiler)
set.seed(169612)
if(FALSE){
  GGO = groupGO(gene     = genes,
                OrgDb    = 'org.Hs.eg.db',
                keyType  = 'ENSEMBL',
                ont      = "CC",
                level    = 3,
                readable = TRUE)
  head(GGO@result[,-5], 100)
}

EGO = enrichGO(gene         = genes,
               OrgDb         = 'org.Hs.eg.db',
               keyType       = 'ENSEMBL',
               ont          = "ALL",
               pvalueCutoff = 1)
head(EGO@result[,-which(colnames(EGO@result) == "geneID")])

sum(EGO@result$p.adjust < 0.01)
# 271

EGO@result = EGO@result[order(EGO@result$p.adjust), ]

library(DOSE)
dotplot(EGO, showCategory=20, x = "p.adjust")
ggsave("scripts/Images/Filippo/ORA/ORA_20_HD_vs_HCL.pdf",
       width = 8,
       height = 12)

DF = EGO@result[,-which(colnames(EGO@result) == "geneID")]
class(DF)
head(DF)

DF = data.frame(ONTOLOGY = DF$ONTOLOGY,
                ID = DF$ID,
                Description = DF$Description, 
                pvalue = DF$pvalue, 
                FDR = DF$p.adjust)
dim(DF)


# add gene ids for each othology
gene = sapply(EGO@result$geneID, strsplit, "/")

gene_unique = data.frame(gene = unique(unlist(gene)))

# save gene names for each onthology.
load("annotation.RData")
G_list = annotation[
  annotation$ensembl_gene_id %in% unique(unlist(gene)),
  ,
  drop = FALSE
]

G_list = G_list[,1:2]
colnames(G_list) = c("gene", "gene_name")
dim(gene_unique); dim(G_list)
gene_unique = merge(gene_unique, G_list, by = "gene", all.x = TRUE)
dim(gene_unique)
head(gene_unique)

gene_name = sapply(gene, function(gg){
  gene_unique$gene_name[match(gg, gene_unique$gene)]
})

DF$gene = gene
DF$gene_name = gene_name

DF = DF[order(DF$pvalue),]
openxlsx::write.xlsx(DF, file = "scripts/results/Filippo/GO_ORA_20_HD_vs_HCL.xlsx")


if(FALSE){
  # map Ensembl ids to Entrez ids
  DF = data.frame(logFC = res$table$logFC,
                  ensembl_gene_id = rownames(res$table))
  
  library('biomaRt')
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list = getBM(filters= "ensembl_gene_id", 
                 attributes= c("ensembl_gene_id",
                               "entrezgene_id"),
                 values=DF$ensembl_gene_id,
                 mart= mart)
  
  G_list = merge(DF, G_list, by = "ensembl_gene_id")
  G_list = G_list[!is.na(G_list$entrezgene_id),]
  
  ranks <- G_list$logFC
  names(ranks) <- G_list$entrezgene_id
  head(ranks)
  barplot(sort(ranks, decreasing = T))
  ?fgsea
}