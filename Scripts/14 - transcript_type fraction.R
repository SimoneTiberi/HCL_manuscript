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

# group labels for the figures (the order of this vector = left-to-right order):
group_labels = c(HD_100K  = "HD S-EVs",
                 HCL_100K = "HCL S-EVs",
                 HD_20K   = "HD L-EVs",
                 HCL_20K  = "HCL L-EVs")

DF = data.frame(protein_coding = protein_coding["protein_coding",], 
                group = colnames(protein_coding) )
DF$group = factor(group_labels[DF$group], levels = group_labels)
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

########################################################################################
# Stacked bar plot of transcript biotypes (% of counts) per group:
########################################################################################
# DF: biotypes (rows, sorted by overall abundance) x groups (columns), in %.
# Keep the n_top most abundant biotypes; pool the remaining ones into "Others".
n_top = min(6, nrow(DF))
DF_plot = DF[seq_len(n_top), , drop = FALSE]
if(nrow(DF) > n_top){
  DF_plot = rbind(DF_plot,
                  Others = colSums(DF[-seq_len(n_top), , drop = FALSE]))
}

# long format:
biotypes = gsub("_", " ", rownames(DF_plot))
DF_long = data.frame(biotype = factor(rep(biotypes, ncol(DF_plot)), levels = biotypes),
                     group = factor(group_labels[rep(colnames(DF_plot), each = nrow(DF_plot))],
                                    levels = group_labels),
                     percentage = c(DF_plot))

# colours (colour-blind safe; extend the vector if n_top > 7):
cols = c("#2A78D6", "#EB6834", "#1BAF7A", "#EDA100",
         "#E87BA4", "#008300", "#6250D6")
cols = setNames(cols[seq_len(n_top)], biotypes[seq_len(n_top)])
if("Others" %in% biotypes) cols["Others"] = "#B4B2A9" # neutral grey for pooled biotypes

p_types = ggplot(DF_long, aes(x = group, y = percentage, fill = biotype)) +
  geom_col(width = 0.95, position = position_stack(reverse = TRUE),
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Percentage (%)", fill = NULL,
       title = "Transcript biotypes") +
  theme_classic(base_size = 13) +
  theme(axis.line = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.9, "cm"),
        legend.key.width = unit(0.9, "cm"),
        legend.key = element_rect(colour = "white", linewidth = 1.5))
p_types
ggsave("scripts/Images/Filippo/DTU/Transcript_type_fractions.pdf",
       plot = p_types,
       width = 7,
       height = 5)

rm(DF)