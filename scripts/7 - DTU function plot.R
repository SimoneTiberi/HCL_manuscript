PLOT_props = function(gene, counts){
  res = res_gene[res_gene$gene == gene, ]
  means = res[,6:7]
  means = means/colSums(means)
  
  prop_samp = data.frame(transcript_it = factor( rep(res$transcript, 2)), 
                         proportion = unlist(c(means)),
                         group = rep(colnames(means), each = nrow(means)),
                         stringsAsFactors = FALSE)
  
  # Plot the estimated average proportions of each groups:
  ggp = ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "transcript_it", y = "proportion", 
                                          fill = "group"),
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(paste("Gene:", gene, sep=" ")) +
    xlab("Transcripts") +
    ylab("Proportions") +
    geom_point(data = prop_samp, 
               aes_string(x = "transcript_it", y = "proportion", 
                          group = "group", fill = "group"), 
               position = position_dodge(width = 0.9), size = 3, shape = 1, 
               alpha = 0.75)
  ggp
}