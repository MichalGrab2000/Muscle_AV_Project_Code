#########Gene Set Enrchment Analysis code###############################################


library(clusterProfiler)
library(org.Hs.eg.db)

################################################################################################################################################
########################################################      LOAD LISTS   #####################################################################
################################################################################################################################################

bulk_post_vs_pre <- read.csv ("~/Documents/Master Thesis/Diff_expression_Post_vs_Pre.csv", header = T) %>%
    drop_na(log2FoldChange) %>%
    drop_na(xiao_score) %>%
    dplyr::rename(gene = X)

bulk_rec_vs_pre <- read.csv ("~/Documents/Master Thesis/Diff_expression_Recovery_vs_Pre.csv", header = T) %>%
    drop_na(log2FoldChange) %>%
    drop_na(xiao_score) %>%
    dplyr::rename(gene = X)


################################################################################################################################################
##################################################      CREATE NAMED VECTOR   ##################################################################
################################################################################################################################################
#Filter the Diff Expression results for xiao_score and LFC values

DE_post_up <- bulk_post_vs_pre %>% dplyr::filter(xiao_score < 0.05) %>% dplyr::filter(log2FoldChange > 0) 
DE_post_down <- bulk_post_vs_pre %>% dplyr::filter(xiao_score < 0.05) %>% dplyr::filter(log2FoldChange < 0)

DE_rec_up <- bulk_rec_vs_pre %>% dplyr::filter(xiao_score < 0.05) %>% dplyr::filter(log2FoldChange > 0)
DE_rec_down <- bulk_rec_vs_pre %>% dplyr::filter(xiao_score < 0.05) %>% dplyr::filter(log2FoldChange < 0)

################################################################################################################################################
######################################################       GSEA FOR POST UP      #############################################################
################################################################################################################################################

# Perform GSEA -------------------------------------------------------------
GSEA_post_up <- enrichGO(
    gene          = DE_post_up$gene,
    universe      = bulk_post_vs_pre$gene,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)

# Simplify GO terms redundancy --------------------------------------------
GSEA_post_up_simplify <- simplify(GSEA_post_up, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
GSEA_post_up_simplify_df <- as.data.frame(GSEA_post_up_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
GSEA_post_up_simplify_df <- mutate(GSEA_post_up_simplify_df, foldEnrich =
                                            (as.numeric(sub("/\\d+", "", GSEA_post_up_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", GSEA_post_up_simplify_df$GeneRatio))) /
                                            (as.numeric(sub("/\\d+", "", GSEA_post_up_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", GSEA_post_up_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(GSEA_post_up_simplify_df,
          file="~/Documents/Master Thesis/GSEA_results_post_up.csv")

################################################################################################################################################
######################################################       GSEA FOR POST down      #############################################################
################################################################################################################################################

# Perform GSEA -------------------------------------------------------------
GSEA_post_down <- enrichGO(
    gene          = DE_post_down$gene,
    universe      = bulk_post_vs_pre$gene,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)

# NO SIGNIFICANT GO TERMS


################################################################################################################################################
######################################################       GSEA FOR REC UP      #############################################################
################################################################################################################################################

# Perform GSEA -------------------------------------------------------------
GSEA_rec_up <- enrichGO(
    gene          = DE_rec_up$gene,
    universe      = bulk_rec_vs_pre$gene,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)

# Simplify GO terms redundancy --------------------------------------------
GSEA_rec_up_simplify <- simplify(GSEA_rec_up, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
GSEA_rec_up_simplify_df <- as.data.frame(GSEA_rec_up_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
GSEA_rec_up_simplify_df <- mutate(GSEA_rec_up_simplify_df, foldEnrich =
                                       (as.numeric(sub("/\\d+", "", GSEA_rec_up_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", GSEA_rec_up_simplify_df$GeneRatio))) /
                                       (as.numeric(sub("/\\d+", "", GSEA_rec_up_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", GSEA_rec_up_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(GSEA_rec_up_simplify_df,
          file="~/Documents/Master Thesis/GSEA_results_rec_up.csv")


################################################################################################################################################
######################################################       GSEA FOR REC DOWN      ############################################################
################################################################################################################################################

# Perform GSEA -------------------------------------------------------------
GSEA_rec_down <- enrichGO(
    gene          = DE_rec_down$gene,
    universe      = bulk_rec_vs_pre$gene,
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)

# Simplify GO terms redundancy --------------------------------------------
GSEA_rec_down_simplify <- simplify(GSEA_rec_down, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
GSEA_rec_down_simplify_df <- as.data.frame(GSEA_rec_down_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
GSEA_rec_down_simplify_df <- mutate(GSEA_rec_down_simplify_df, foldEnrich =
                                      (as.numeric(sub("/\\d+", "", GSEA_rec_down_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", GSEA_rec_down_simplify_df$GeneRatio))) /
                                      (as.numeric(sub("/\\d+", "", GSEA_rec_down_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", GSEA_rec_down_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(GSEA_rec_down_simplify_df,
          file="~/Documents/Master Thesis/GSEA_results_rec_down.csv")


#PLOT results
barplot(GSEA_rec_up)

dotplot(GSEA_rec_down, showCategory=30) + ggtitle("dotplot for GSEA")
GSEA_rec_down_DF <- GSEA_rec_down@result


heatplot(GSEA_post_up)
