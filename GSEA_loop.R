library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr) 

timepoints <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")

#Lists to store results
GSEA_results <- list()      
GSEA_top10   <- list()      
ridge_plots  <- list()      

# Loop over each timepoint
for (tp in timepoints) {
  

  file_path <- paste0("~/Documents/Master Thesis/", tp, "AV_diff_expression_symbols.csv")
  
  
  df <- read.csv(file_path, header = TRUE) %>%
    drop_na(log2FoldChange) %>%
    drop_na(stat) %>%
    dplyr::rename(gene = X)
  
 
  geneList <- df$stat
  names(geneList) <- df$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  #GO Terms
  gsea_result <- gseGO(
    geneList      = geneList,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    minGSSize     = 10,
    maxGSSize     = 500,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    verbose       = FALSE
  )
  
  
  GSEA_results[[tp]] <- gsea_result
  
  # Top 10 terms based on adjusted p-value
  top10_df <- as.data.frame(gsea_result@result) %>%
    arrange(p.adjust) %>%
    head(10)
  GSEA_top10[[tp]] <- top10_df
  
  # Create and store a ridge plot for this timepoint
  plot <- ridgeplot(gsea_result,
                  showCategory = 10,
                  fill         = "p.adjust") +
    labs(title = paste("Ridge Plot -", tp),
         x     = "Gene Ranking",
         y     = "GO Term")
  
  ridge_plots[[tp]] <- plot
  
  
  print(plot)
}
