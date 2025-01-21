

library("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
library(readxl)
library (ggrepel)
library("dplyr")


counts <- read.table("/Users/MichalGrabowski/Documents/Master Thesis/bulk_counts.txt", header = TRUE, sep = "\t")

counts[1:5,1:5]

# renaming rownames, formatting table for RNA number as header 
rownames(counts) <- counts$GeneSymbol

counts <- counts[,-1]

counts[1:5,1:5]

meta <- read_excel("/Users/MichalGrabowski/Documents/Master Thesis/bulk_Metadata_3.xlsx")



# I am only interested in the placebo condition, therefore filter condition for 'P'
meta_flt <- meta1 %>% 
  filter(Condition == 'P')

# factorise timepoint
meta_flt$Time <- factor(meta_flt$Time, levels = c("Pre", "Post", "Recovery"))

meta_flt$Subject <- factor(meta_flt$Subject)

#Match columns of count data to rows of metadata 
rownames(meta_flt) <- meta_flt$RNA_number

#Extract RNA numbers for the placebo condition 
RNA_Placebo <- select(meta_flt, RNA_number, Condition)
rownames(RNA_Placebo) <- RNA_Placebo$RNA_number

#Filter for Placebo (P) corresponding RNAnumbers in counts data
RNA_placebo_rows <- rownames(RNA_Placebo)
filtered_counts <- counts[, RNA_placebo_rows]

#Check 
colnames(filtered_counts) == rownames(meta_flt)


#Data quality assessment
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
PCAplot <- plotPCA(vsd, intgroup=c("Condition", "Subject"))
PCAplot+
  theme_bw(base_size = 15)


############################################################################################
###################     Manual PCA plot       ##############################################
############################################################################################
#Explore plot deeper, greater manipulation, identify what factors drive each principal component 

vsd_counts <- assay(vsd) %>% as.data.frame()

dim(vsd_counts)

vsd_counts$var <- apply(vsd_counts, 1, var)

#filter for top 500 most variable genes
vsd_counts_flt <- vsd_counts %>% 
  arrange(-var) %>% 
  head(500) %>% 
  dplyr::select(-var)

# PCA analysis

pca <- prcomp(t(vsd_counts_flt), scale. = TRUE)

pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta_flt, by = "RNA_number") %>% 
  ggplot(aes(x = PC1, y = PC2))+
    geom_point(size = 3, alpha = 1, aes(colour = Subject)) + theme_classic(base_size = 15) + 
        geom_label_repel(label = colnames(vsd_counts_flt), max.overlaps = Inf, size = 3)
    #Identified RNA_number 021436 as the outlier 

pca$x %>%
    as.data.frame() %>% 
    rownames_to_column(var = "RNA_number") %>% 
    merge(.,meta_flt, by = "RNA_number") %>% 
    mutate(Label = paste0(Subject, "_", Time)) %>% 
    ggplot(aes(x = PC1, y = PC2))+
    geom_point(size = 3, alpha = 1, aes(colour = Subject)) + theme_classic(base_size = 15) + 
    #geom_label(label = colnames(vsd_counts_flt),nudge_x = -5, nudge_y = 1)
    geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)

# PCA time
pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta_flt, by = "RNA_number") %>% 
  mutate(Label = paste0(Subject, "_", Time)) %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Time)) + theme_classic(base_size = 15) 
  

# what drives PC1
pca$rotation %>% 
  as.data.frame() %>% 
  arrange(-PC1) %>% 
  dplyr::select(PC1) %>% 
  rownames_to_column(var = "gene") %>% pull(gene)

############################################################################################
###################     Manual PCA plot       ##############################################
############################################################################################

#Remove identified RNA_number 021436 outlier (a column in filtered counts and a row (42) in meta_flt)
filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA021436))
filtered_counts_noOutliers
meta_flt_noOutliers <- meta_flt[-42,]
meta_flt_noOutliers


#DESeq data set construction
dds <- DESeqDataSetFromMatrix(countData = filtered_counts_noOutliers,
                              colData = meta_flt_noOutliers,
                              design = ~ Subject + Time)

dim(dds)

#Pre Filtering, only keep counts above 10
keep <- rowSums(counts(dds) >= 10) >= 41 #number of rows (we just removed one) 
dds <- dds[keep,]

dim(dds)

dds$condition <- factor(dds$Time, levels = c("Pre","Post", "Recovery"))

summary(meta_flt_noOutliers)



#Differential expression analysis 
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
Diff_expression_Recovery_vs_Pre <- results(dds, contrast=c("Time","Recovery","Pre"))
Diff_expression_Recovery_vs_Pre
Diff_expression_Recovery_vs_PreOrdered <- Diff_expression_Recovery_vs_Pre[order(res$pvalue),]  #Ordered by the smallest p-value 
#Add XIAO score column (Xiao et al., 2014) for enhanced gene ranking in preparation of GSEA
Diff_expression_Recovery_vs_PreOrdered$xiao_score = Diff_expression_Recovery_vs_PreOrdered$log2FoldChange * -log(Diff_expression_Recovery_vs_PreOrdered$pvalue)
Diff_expression_Recovery_vs_PreOrdered
summary(res)

Diff_expression_Recovery_vs_Post <- results(dds, contrast=c("Time","Recovery","Post"))
Diff_expression_Recovery_vs_Post
Diff_expression_Recovery_vs_PostOrdered <- Diff_expression_Recovery_vs_Post[order(res$pvalue),] 
#Add XIAO score column (Xiao et al., 2014) for enhanced gene ranking in preparation of GSEA
Diff_expression_Recovery_vs_PostOrdered$xiao_score = Diff_expression_Recovery_vs_PostOrdered$log2FoldChange * -log(Diff_expression_Recovery_vs_PostOrdered$pvalue)
Diff_expression_Recovery_vs_PostOrdered

Diff_expression_Post_vs_Pre <- results(dds, contrast=c("Time","Post","Pre"))
Diff_expression_Post_vs_Pre  #Purely the effect of exercise
Diff_expression_Post_vs_PreOrdered <- Diff_expression_Post_vs_Pre[order(res$pvalue),]
#Add XIAO score column (Xiao et al., 2014) for enhanced gene ranking in preparation of GSEA
Diff_expression_Post_vs_PreOrdered$xiao_score = Diff_expression_Post_vs_PreOrdered$log2FoldChange * -log(Diff_expression_Post_vs_PreOrdered$pvalue)
Diff_expression_Post_vs_PreOrdered


#LogFoldChange (i.e. effect size) --> useful for visualization and ranking of genes
resLFC <- lfcShrink(dds, coef="Time_Recovery_vs_Pre", type="apeglm")
resLFC 
# Shrinkage is the reduction in effects of sampling variation (useful for lowly expressed genes?)


#MA plots 
plotMA(Diff_expression_Recovery_vs_Pre, ylim=c(-2,2))

#MA plot for shrinkage: removes noise associated with fold-changes among low count genes 
plotMA(resLFC, ylim=c(-2,2))

#Idenitfy genes by clicking directly on graph 
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]



#Export Diff_expression results 
write.csv(as.data.frame(Diff_expression_Recovery_vs_PreOrdered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expression_Recovery_vs_Pre.csv")
write.csv(as.data.frame(Diff_expression_Recovery_vs_PostOrdered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expression_Recovery_vs_Post.csv")
write.csv(as.data.frame(Diff_expression_Post_vs_PreOrdered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expression_Post_vs_Pre.csv")


#################################### VOLCANO PLOTS ######################################################

path <- "/Users/MichalGrabowski/Documents/Master Thesis"

Recovery_vs_Pre_Results <- read.csv(paste0(path, '/Diff_expression_Recovery_vs_Pre.csv'), row.names = 1)

#Add explicit gene_symbol column for easier processing
Recovery_vs_Pre_Results$Gene_symbol <- rownames(Recovery_vs_Pre_Results)

#ggplot(data = Recovery_vs_Pre_Results, aes(x = log2FoldChange, y = -log10(pvalue))) +
  #geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = 'dashed') +  #Addition of threshold lines
  #geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +  #Horizontal line for p-value threshold 
  #geom_point() + theme_minimal()


# add a column of NAs
Recovery_vs_Pre_Results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Recovery_vs_Pre_Results$diffexpressed[Recovery_vs_Pre_Results$log2FoldChange > 0.6 & Recovery_vs_Pre_Results$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Recovery_vs_Pre_Results$diffexpressed[Recovery_vs_Pre_Results$log2FoldChange < -0.6 & Recovery_vs_Pre_Results$pvalue < 0.05] <- "DOWN"


p <- ggplot(data = Recovery_vs_Pre_Results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  geom_point() + theme_minimal()

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p2 <- p + scale_colour_manual(values = mycolors)
p2

Recovery_vs_Pre_Results$label <- NA
Recovery_vs_Pre_Results$label[Recovery_vs_Pre_Results$diffexpressed != "NO"] <- Recovery_vs_Pre_Results$Gene_symbol[Recovery_vs_Pre_Results$diffexpressed != "NO"]
#Recovery_vs_Pre_Results$delabel <- ifelse(Recovery_vs_Pre_Results$Gene_symbol %in% head(Recovery_vs_Pre_Results[order(Recovery_vs_Pre_Results$padj), "Gene_symbol"], 30), Recovery_vs_Pre_Results$Gene_symbol, NA)


ggplot(data = Recovery_vs_Pre_Results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("red", "black", "green")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="blue", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="blue", linetype = 'dashed') +
  geom_text_repel(max.overlaps = 20)
  
################# HEATMAP ######################################

# get top 20 DEGs post vs. pre

top20_degs_post_vs_pre <- Diff_expression_Post_vs_Pre %>% 
  as.data.frame() %>% 
  arrange(-abs(log2FoldChange)) %>% 
  head(20) %>% 
  rownames_to_column(var = "Gene_id") %>% 
  pull(Gene_id)

norm_counts <- counts(dds,normalized=TRUE)
norm_counts[top20_degs_post_vs_pre,]

#Filter RNA numbers for just PRE condition 
samples_pre <- meta_flt %>% 
  filter(Time == "Pre") %>% 
  pull(RNA_number)

# Norm counts PRE
pre_counts <- norm_counts[top20_degs_post_vs_pre, samples_pre]

# Rowmeans
pre_means <- rowMeans(pre_counts) %>% 
  as.data.frame() %>% 
  dplyr::select(pre = ".")

#Fitler RNA numbers for POST condition
samples_post <- meta_flt %>%
  filter(Time == "Post") %>%
  pull(RNA_number)

#Norm counts POST
post_counts <- norm_counts[top20_degs_post_vs_pre, samples_post]

#Rowmeans
post_means <- rowMeans(post_counts) %>%
  as.data.frame() %>% 
  dplyr::select(post = ".")

#Fitler RNA numbers for REC condition
samples_recovery <- meta_flt_noOutliers %>% 
  filter(Time == "Recovery") %>% 
  pull(RNA_number)

#Norm counts RECOVERY
recovery_counts <- norm_counts[top20_degs_post_vs_pre, samples_recovery]

#Rowmeans
recovery_means <- rowMeans(recovery_counts) %>%
  as.data.frame() %>% 
  dplyr::select(recovery = ".")

#Combine all timepoint means into same data frame
PrePostRecoveryMeans <- cbind(pre_means, post_means, recovery_means)

#Convert dataframe to Long Form for heatmap creation
PrePostRecoveryMeansLongForm <- PrePostRecoveryMeans %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)

Heatmap <- ggplot(data = PrePostRecoveryMeansLongForm, aes(x = rowname, y = colname, fill = log(value))) +
  xlab("Gene") +
  ylab("Time") + 
  geom_tile() +
  theme_minimal()
Heatmap

################# INDIVIDUAL GENE PLOTS ################################################

#Format norm counts
norm_counts_formatted <- norm_counts %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Gene_id")

#Only interested in NR4A3 gene
NR4A3only <- norm_counts_formatted %>% 
  filter(Gene_id == "NR4A3") 
 
#Transpose to long form
NR4A3onlyTransposed <- t(NR4A3only) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  rename(c("V1" = "NR4A3")) %>% 
  filter(RNA_number != "Gene_id")

#Merge with meta data for timepoints, convert to numeric form 
NR4A3_meta <- merge(NR4A3onlyTransposed, meta_flt_noOutliers) %>% 
  transform(NR4A3num = as.numeric(NR4A3))

NR4A3_meta %>% str()

#Boxplot 
ggplot(NR4A3_meta, aes(x = Time, y = NR4A3num)) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group=Subject)) 














ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

library("colorspace")
library("pheatmap")

norm_counts <- counts(dds,normalized=TRUE)





df <- as.data.frame(colData(dds)[,c("Time")])

pheatmap(assay(vsd)[PrePostRecoveryMeans,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

heatmap(df, scale = "none")
