library("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
library(readxl)
library (ggrepel)
library("dplyr")
library(clusterProfiler)
library(org.Hs.eg.db)

#Combine all individual sample txt. count files into 1
file_path <- "/Users/MichalGrabowski/Documents/Master Thesis/AV data"
file_list <- list.files(path = file_path, pattern = ".txt", full.names = TRUE)

i = 1
samples <- list()
for (i in 1:length(file_list)) {
  
  name <- str_split_i(file_list[i], "_", 1)
  name <- str_split_i(name, "/", 7)
  
  count <- read.delim(file_list[i], header = F)
  
  colnames(count) <- c("Gene_ID", name)
  
  samples[[name]] <- count
  
  print(i)
  
}

raw_counts <- samples[[1]]

# merge samples into table
for (i in 2:length(samples)) {
  raw_counts <- merge(raw_counts, samples[[i]], by = "Gene_ID")
  
  print(i)
  
}

#############################################################

#Convert gene_id to symbol 
library(AnnotationDbi)
library(org.Hs.eg.db) # Hs = human

gene_ids <- counts %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
print(gene_symbols)

#Convert the symbols to data frame to merge with counts
symbols_data_frame <- read.table(text = gene_symbols)
counts$gene_symbol <- symbols_data_frame$V1
counts <- counts %>% 
  relocate(gene_symbol, .before = everything())

counts <- counts %>%     ##### Spills error as some gene symbols are NA --> results in duplicate rownames 
  column_to_rownames(var = "gene_symbol")

#Use gene_id as rownames for now 
rownames(counts) <- counts$gene_id

##############################################################

meta <- read_excel("/Users/MichalGrabowski/Documents/Master Thesis/Meta_AV_data.xlsx")

# factorise timepoint
meta$Time <- factor(meta$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))

# Factors
meta$Subject <- factor(meta$Subject)
meta$Subject <- factor(meta$Subject)
meta$Condition <- factor(meta$Condition, levels = c("A", "V"))


rownames(meta) <- meta$RNA_number

raw_counts <- raw_counts %>%     
  column_to_rownames(var = "Gene_ID")


#Cut rows that are not ENSG.... Gene ID
filtered_counts <- raw_counts[-(1:5), ]
filtered_counts <- raw_counts[-(62783:62883), ]

n_row <- nrow(raw_counts)  # Get the total number of rows
filtered_counts <- raw_counts[-c(1:5, (n_row-173):n_row), ]



#DESeq2..............................................................................................................

#DESeq data set construction
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = meta,
                              design = ~ Subject + Condition + Time + Condition:Time)

#DESeq data set construction
dds1 <- DESeqDataSetFromMatrix(countData = baseline_counts,
                              colData = meta_baseline,
                              design = ~ Subject + combined)

#DESeq
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                               colData = meta,
                               design = ~ Subject + Time*Condition)

# variance stbilized counts
vsd <- vst(dds, blind = T)

vsd_counts <- assay(vsd) %>% as.data.frame()

vsd_counts$var <- apply(vsd_counts, 1, var)

# get the top 500 most variable point
vsd_counts_filtered <- vsd_counts %>% 
  arrange(-var) %>% 
  head(500) %>% 
  dplyr::select(-var)

#MANUAL PCA

pca <- prcomp(t(vsd_counts_filtered), scale. = TRUE)


pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta, by = "RNA_number") %>% 
  mutate(Label = paste0(Subject, "_", Time)) %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) + 
  #geom_label(label = colnames(vsd_counts_flt),nudge_x = -5, nudge_y = 1)
  geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)


pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta, by = "RNA_number") %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) +
  geom_label_repel(label = colnames(vsd_counts_filtered), max.overlaps = Inf, size = 3)

#The two outliers are RNA031403 and RNA031398
#Remove identified RNA_number outliers 
filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA031403, RNA031398))
remove_RNA <- c("RNA031403", "RNA031398")
meta_noOutliers <- meta[!rownames(meta) %in% remove_RNA, ]
rownames(meta_noOutliers) <- meta_noOutliers$RNA_number

#PCA again wthout Outliers
dds_clean <- DESeqDataSetFromMatrix(countData = filtered_counts_noOutliers,
                              colData = meta_noOutliers,
                              design = ~ Subject + Condition + Time + Condition:Time)

dim(dds_clean)

# variance stbilized counts
vsd_2 <- vst(dds_clean, blind = T)

vsd_counts_2 <- assay(vsd_2) %>% as.data.frame()

vsd_counts_2$var <- apply(vsd_counts_2, 1, var)

# get the top 500 most variable point
vsd_counts_filtered_2 <- vsd_counts_2 %>% 
  arrange(-var) %>% 
  head(500) %>% 
  dplyr::select(-var)

#MANUAL PCA AGAIN 

pca <- prcomp(t(vsd_counts_filtered_2), scale. = TRUE)

pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta_noOutliers, by = "RNA_number") %>% 
  mutate(Label = paste0(Subject, "_", Time)) %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Time)) + theme_classic(base_size = 15)+
  geom_label(label = colnames(vsd_counts_filtered_2),nudge_x = -5, nudge_y = 1)
  #geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)

#Four more Outliers (RNA031468, RNA031438, RNA031454, RNA031400) --> mainly different along PCA2
#RNA031400 has a high duplication count and low reads overall --> can be removed on that matter
# what drives PC2
PCA2_genes <- pca$rotation %>% 
  as.data.frame() %>% 
  arrange(-PC2) %>% 
  head(100) %>% 
  dplyr::select(PC2) %>% 
  rownames_to_column(var = "gene") %>% 
  pull(gene)

#Convert IDs to symbols 
gene_symbols <- mapIds(org.Hs.eg.db, keys = PCA2_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_ids_PCA2 <- mapIds(org.Hs.eg.db, keys = PCA2_genes, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first") %>% as.data.frame()
gene_ids_PCA2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(gene_ids_PCA2), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
print(gene_symbols)

colnames(gene_ids_PCA2) <- c("Entrezid", "symbol")


#Overrepresentation analysis 
go_enrich <- enrichGO(gene = gene_ids_PCA2$Entrezid, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "ALL",  # "BP" for Biological Process, "MF" for Molecular Function, or "CC" for Cellular Component
                      pvalueCutoff = 0.05)

head(go_enrich)
barplot(go_enrich, showCategory = 10)

#Pre Filtering, only keep counts above 10
keep <- rowSums(counts(dds_clean) >= 10) >= 90 #number of rows 
dds_keep <- dds_clean[keep,]

dim(dds_keep)

#Later: >10 counts for 80% of samples at each timepoint? 

filtered_counts_10 <- counts(dds_keep, normalized = FALSE)



plot(density(filtered_counts_noOutliers[,1]), 
     main = "Density Plot of Counts per Sample", 
     xlab = "Counts", 
     ylab = "Density", 
     col = 1,  # Color for the first sample
     lwd = 2, 
     xlim = c(0,50000),
     ylim = c(0, 1)) 
    

# Add lines for other samples
for (i in 2:ncol(filtered_counts_noOutliers)) {
  lines(density(filtered_counts_noOutliers[,i]), col = i, lwd = 2)  # Different color for each sample
}



# Add legend
legend("topright", legend = colnames(data), col = 1:ncol(data), lty = 1, lwd = 2)




























dds <- DESeq(dds)
results(dds)
Diff_expression_AV_baseline <- results(dds, contrast=c("Condition","A","V"))
Diff_expression_AV_baseline
Diff_expression_AV_baselineOrdered <- Diff_expression_AV_baseline[order(res$pvalue),]

#Diff expression for A vs V across all time points --> not necessarily useful
Diff_expression <- results(dds1, contrast=list("Condition_V_vs_A"))
Diff_expression

resultsNames(dds)


########################## LIMMMMMMMMAAAAAA ###################################
library(limma)
library(edgeR)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

design <- model.matrix(~ Subject + Time * Condition, data = meta)
colnames(design)
colnames(design) <- make.names(colnames(design))

dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge)
voom_output <- voom(dge, design, plot = TRUE)

# Extract normalized expression matrix
norm_counts <- voom_output$E

fit <- lmFit(voom_output, design)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(
  Ex2_vs_Ex1 = ((TimeEx2 - TimeEx2.ConditionV) - 
                  (TimeEx1 - TimeEx1.ConditionV)),
  #Baseline compared to Exercise
  Baseline_vs_Ex2 = ((TimeEx2 - TimeEx2.ConditionV) -
                       (X.Intercept. - ConditionV)),
  #Baseline compared to Recovery
  Baseline_vs_Rec2 = ((TimeRec2 - TimeRec2.ConditionV) - 
                        (X.Intercept. - ConditionV)),
  #Exercise compared to Recovery
  Ex2_vs_Rec2 = ((TimeRec2 - TimeRec2.ConditionV) - 
                   (TimeEx2 - TimeEx2.ConditionV)),

  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top DE genes for the Ex2_vs_Ex1 contrast

top_genes <- topTable(fit2, coef = "Baseline_vs_Ex2", adjust = "fdr", number = 20)

#How many significant genes? 
sig_genes <- topTable(fit2, coef = "Baseline_vs_Ex2", number = Inf, p.value = 0.05)

#WHAT DO WE DO WITH THE FLOW RATE NOW?
#Multiplt raw counts by the flow rate before running LIMMA

#keep <- rowSums(counts >= 10) 


pca_result <- prcomp(t(norm_counts), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(norm_counts),
  Condition = meta$Condition 
)

# Create PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")





dds <- DESeqDataSetFromMatrix(countData = norm_expr, 
                              colData = meta, 
                              design = ~ Condition)

#Volcano Plot
ggplot(top_genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot: Ex2_vs_Ex1", x = "Log Fold Change", y = "-Log10 P-Value")








#Data quality assessment
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
PCAplot <- plotPCA(vsd, intgroup=c("Condition", "Subject"))
PCAplot+
  theme_bw(base_size = 15)




















#Create new combined column in meta
meta <- meta %>% 
  mutate(combined = paste0(Condition,"_",Time))

#Combined as a factor
meta_baseline$combined <- factor(meta_baseline$combined, levels = c("A_Baseline", "V_Baseline"))
meta_baseline$combined <- relevel(meta_baseline$combined, ref = "A_Baseline")


meta_baseline <- meta[c(Baseline_RNAs),]
rownames(meta_baseline) <- meta_baseline$RNA_number

rownames(meta) <- meta$RNA_number

#Fish out RNA_numbers in meta data, filter for those in counts data (some are not assigned)
RNA_Numbers <- meta %>% 
  pull(RNA_number)
filtered_counts <- counts[, c(RNA_Numbers)]

#Extract RNAs for baseline condition
Baseline_RNAs <- meta %>% 
  filter(Time == "Baseline") %>% 
  pull(RNA_number)

baseline_counts <- counts[, c(Baseline_RNAs)]














#Make dataframes according to Time point, split into A and V with the same subject in each numbered column 
#Baseline..........................................................................................

#Extract RNA numbers for the Baseline_A Condition 
RNA_Numbers_Baseline_A <- meta %>% 
  filter(Time == "Baseline") %>% 
  filter(Condition == 'A') %>% 
  pull(RNA_number)

#Extract RNA numbers for the Baseline_V Condition 
RNA_Numbers_Baseline_V <- meta %>% 
  filter(Time == "Baseline") %>% 
  filter(Condition == 'V') %>% 
  pull(RNA_number)

#Baseline A counts, one column is one subject
Baseline_A_counts <- counts[, c(RNA_Numbers_Baseline_A)]

#Baseline V counts, one column is one subject
Baseline_V_counts <- counts[, c(RNA_Numbers_Baseline_V)]

#Subtract both count data frames to get the AV difference
Baseline_AV_diff <- Baseline_A_counts - Baseline_V_counts


#Passive............................................................................................

#Extract RNA numbers for the Passive_A Condition 
RNA_Numbers_Passive_A <- meta %>% 
  filter(Time == "Passive") %>% 
  filter(Condition == 'A') %>% 
  pull(RNA_number)

#Extract RNA numbers for the Passive_V Condition 
RNA_Numbers_Passive_V <- meta %>% 
  filter(Time == "Passive") %>% 
  filter(Condition == 'V') %>% 
  pull(RNA_number)

#Baseline A counts, one column is one subject
Passive_A_counts <- counts[, c(RNA_Numbers_Passive_A)]

#Baseline V counts, one column is one subject
Passive_V_counts <- counts[, c(RNA_Numbers_Passive_V)]

#Subtract both count data frames to get the AV difference
Passive_AV_diff <- Passive_V_counts - Passive_A_counts
