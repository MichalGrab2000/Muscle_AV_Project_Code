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
library(ggpubr)

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

#############################################################
#Code to convert IDs to symbols --> best to do later, after DGE output 
library(AnnotationDbi)
library(org.Hs.eg.db) # Hs = human

filtered_counts_geneid <- rownames_to_column(filtered_counts, var = "gene_id")

gene_ids <- filtered_counts_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
print(gene_symbols)

#Convert the symbols to data frame to merge with counts
symbols_data_frame <- read.table(text = gene_symbols)
filtered_counts$gene_symbol <- symbols_data_frame$V1
filtered_counts_gene_symbol <- filtered_counts %>% 
  relocate(gene_symbol, .before = everything())

rownames(filtered_counts_gene_symbol) <- NULL

counts_gene_symbol <- filtered_counts_gene_symbol %>%     ##### Spills error as some gene symbols are NA and there are some duplciate genes --> results in duplicate rownames 
  column_to_rownames(var = "gene_symbol")

#Use gene_id as rownames for now 
rownames(counts) <- counts$gene_id

##############################################################



#DESeq2..............................................................................................................

#DESeq data set construction
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = meta,
                              design = ~ Subject + Condition + Time + Condition:Time)

#DESeq data set construction
dds1 <- DESeqDataSetFromMatrix(countData = baseline_counts,
                              colData = meta_baseline,
                              design = ~ Subject + combined)

# variance stbilized counts
vsd <- vst(dds, blind = T)

vsd_counts <- assay(vsd) %>% as.data.frame()

vsd_counts$var <- apply(vsd_counts, 1, var)

# get the top 500 most variable points
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
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) 
  geom_label(label = colnames(vsd_counts_flt),nudge_x = -5, nudge_y = 1) +
  geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)


pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta, by = "RNA_number") %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) +
  geom_label_repel(label = colnames(vsd_counts_filtered), max.overlaps = Inf, size = 3)

#The two outliers are RNA031403 and RNA031398 ---> same subject, different Time Points 
#Remove identified RNA_number outliers 
filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA031403, RNA031398))
remove_RNA <- c("RNA031403", "RNA031398")
meta_noOutliers <- meta[!rownames(meta) %in% remove_RNA, ]
rownames(meta_noOutliers) <- meta_noOutliers$RNA_number

filtered_counts_noOutliers_gene_symbol <- filtered_counts_noOutliers
filtered_counts_noOutliers_gene_symbol$gene_symbol <- symbols_data_frame$V1



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

pca2 <- prcomp(t(vsd_counts_filtered_2), scale. = TRUE)

pca2$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta_noOutliers, by = "RNA_number") %>% 
  mutate(Label = paste0(Subject, "_", Time)) %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) +
  geom_label(label = colnames(vsd_counts_filtered_2),nudge_x = -5, nudge_y = 1) 
  #geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)

#Four more Outliers (RNA031468, RNA031438, RNA031454, RNA031400) --> mainly different along PCA2
#RNA031400 has a high duplication count and low reads overall --> can be removed on that matter?
remove_RNA <- c("RNA031403", "RNA031398", "RNA031400")
meta_noOutliers <- meta[!rownames(meta) %in% remove_RNA, ]

filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA031403, RNA031398,RNA031400))

# what drives PC2
PCA2_genes <- pca2$rotation %>% 
  as.data.frame() %>% 
  arrange(PC2) %>% 
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
keep <- rowSums(counts(dds_clean) >= 10) >= 90 #number of columns 
dds_keep <- dds_clean[keep,]

dim(dds_keep)

#Later: >10 counts for 80% of samples at each timepoint? 

filtered_counts_10 <- counts(dds_keep, normalized = FALSE)



plot(density(filtered_counts_10[,1]), 
     main = "Density Plot of Counts per Sample", 
     xlab = "Counts", 
     ylab = "Density", 
     col = 1,  # Color for the first sample
     lwd = 2, 
     xlim = c(0,50000),
     ylim = c(0, 0.1)) 
    

# Add lines for other samples
for (i in 2:ncol(filtered_counts_10)) {
  lines(density(filtered_counts_10[,i]), col = i, lwd = 2)  # Different color for each sample
}

# Add legend
legend("topright", legend = colnames(filtered_counts_10), col = 1:ncol(data), lty = 1, lwd = 2)


#Filter raw count data for the spikes
counts_spikes <- raw_counts[-(1:5), ]


row_groups <- substr(rownames(counts_spikes), 1, 2)
unique_groups <- unique(row_groups)

summed_data <- data.frame(
  Group = unique_groups, 
  t(sapply(unique_groups, function(group) {
    rows_in_group <- rownames(counts_spikes)[row_groups == group]
    colSums(counts_spikes[rows_in_group, ], na.rm = TRUE)  # Sum all columns for this group
  }))
)

summed_long <- pivot_longer(
  summed_data,
  cols = -Group,
  names_to = "RNA_number",
  values_to = "Spike_sum"
) %>% filter(!RNA_number %in% c("RNA031403", "RNA031398"))


ggplot(summed_long, aes(x = RNA_number, y = Spike_sum, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  labs(title = "Grouped Bar Chart by RNA ID",
       x = "RNA ID",
       y = "Sum Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#Now group by A/V and time point --> to ensure spikes can normalize for the flow rate 


#Extract RNA numbers Based on time point 
Extract_RNA_numbers_time <- function(meta, Timepoint) {
  filtered_metadata_time <- meta[meta$Time %in% Timepoint, ]
  return(filtered_metadata_time$RNA_number)
}

RNA_numbers_baseline <- Extract_RNA_numbers_time(meta, c("Baseline"))
RNA_numbers_passive <- Extract_RNA_numbers_time(meta, c("Passive"))
RNA_numbers_Ex1 <- Extract_RNA_numbers_time(meta, c("Ex1"))
RNA_numbers_Ex2 <- Extract_RNA_numbers_time(meta, c("Ex2"))
RNA_numbers_Ex3 <- Extract_RNA_numbers_time(meta, c("Ex3"))
RNA_numbers_Rec1 <- Extract_RNA_numbers_time(meta, c("Rec1"))
RNA_numbers_Rec2<- Extract_RNA_numbers_time(meta, c("Rec2"))

#Extract RNA numbers Based on Time & Condition 
Extract_rna_numbers_time_condition <- function(meta, Timepoint, Conditions) {
  filtered_metadata <- subset(meta, 
                              Time %in% Timepoint & 
                                Condition %in% Conditions)
  return(filtered_metadata$RNA_number)
}

RNA_numbers_Baseline_A <- Extract_rna_numbers_time_condition(meta, c("Baseline"), c("A"))
RNA_numbers_Baseline_V <- Extract_rna_numbers_time_condition(meta, c("Baseline"), c("V"))
#..........

#Merge summarized data with meta data
merged_data <- summed_long %>%
  left_join(meta, by = "RNA_number")

summarized_data <- merged_data %>%
  group_by(Group, Condition, Time) %>%
  summarize(Sum_Value = sum(Spike_sum, na.rm = TRUE)) %>% 
  mutate(combined = paste0(Time, "_", Condition))

ggplot(summarized_data, aes(x = combined, y = Sum_Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  #facet_wrap(~ Time) +  # Facet by group (EN, ER, etc.)
  labs(
    title = "Sum of counts and spikes across Timepoint and Condition",
    x = "Timepoint",
    y = "log10(Sum of Values)",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(summarized_data %>% filter(Time %in% c("Baseline", "Ex2")), aes(x = combined, y = Sum_Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  #facet_wrap(~ Time) +  # Facet by group (EN, ER, etc.)
  labs(
    title = "Grouped Bar Plot by Condition and Time",
    x = "Timepoint",
    y = "Sum of Values",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#We would expect the Spike counts to go down at the Exercising timepoints as they should be more diluted due to increased flow rate

#Remove counts to display spike data more accurately without exponential distortion 
Only_spikes <- summarized_data[summarized_data$Group != "EN", ]
Only_spikes$Time <- factor(Only_spikes$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))


ggplot(Only_spikes %>% filter(Time %in% c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")), aes(x = combined, y = Sum_Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  #facet_wrap(~ Time) +  # Facet by group (EN, ER, etc.)
  labs(
    title = "Sum of spikes by Time and Cndition",
    x = "Timepoint",
    y = "Sum",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Sum and Flow rate correlation graphs

flow_data <- read_excel("/Users/MichalGrabowski/Documents/Master Thesis/AV_Flow_data.xlsx")

merged_with_flowrate <- summed_long %>%
  left_join(flow_data, by = "RNA_number")

merged_with_flowrate_spikes <- merged_with_flowrate[!merged_with_flowrate$Group %in% c("EN", "RN"), ]



ggplot(merged_with_flowrate_spikes, aes(x = Flow, y = Value, color = Group)) +
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE) +  #trendline 
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  geom_text(aes(label = RNA_number), vjust = -0.5, size = 3) +
  labs(title = "Flow vs Sum of spikes by spike type ",
       x = "Flow",
       y = "Sum of spikes") +
  theme_minimal()
#One notable outlier in ERCC spikes --> RNA031400 


#Spike count bar chart by timepoint 
merged_with_flowrate_spikes$Time <- factor(
  merged_with_flowrate_spikes$Time,
  levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
)

ggplot(na.omit(merged_with_flowrate_spikes), aes(x = Time, y = Value, fill = Group)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(), alpha = 0.7) +
  geom_jitter(aes(color = Group), width = 0.2, size = 3) +
  #geom_line(aes(group = Subject, color = Group), size = 1) +
  labs(
    title = "Spike counts across Timepoints ",
    x = "Timepoint",
    y = "Spike sum",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal()




#SPIKE PCA plot

just_spike_counts <- counts_spikes[!substr(rownames(counts_spikes), 1, 3) %in% c("ENS"), ]

ERCC_spikes <- just_spike_counts[substr(rownames(just_spike_counts), 1, 3) %in% c("ERC"), ]
R1_spikes <- just_spike_counts[substr(rownames(just_spike_counts), 1, 3) %in% c("R1"), ]
R2_spikes <- just_spike_counts[substr(rownames(just_spike_counts), 1, 3) %in% c("R2"), ]


pca_ERCC <- prcomp(t(ERCC_spikes), scale. = FALSE)

pca_ERCC$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta, by = "RNA_number") %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) +
  geom_label_repel(label = colnames(ERCC_spikes), max.overlaps = Inf, size = 3)


#Filter for genes that appear in at least 80% of the samples at given time point  
#MAKE IT A FUNCTION
# Define a function to process each timepoint
process_timepoint <- function(timepoint, meta, filtered_counts) {
  # Extract RNA numbers for the specific timepoint
  timepoint_RNAs <- meta %>%
    filter(Time == timepoint) %>%
    pull(RNA_number)
  
  # Subset counts data
  timepoint_counts <- filtered_counts[, c(timepoint_RNAs)]
  
  # Subset metadata
  meta_timepoint <- meta[c(timepoint_RNAs),]
  rownames(meta_timepoint) <- meta_timepoint$RNA_number
  
  # Convert columns to factors
  meta_timepoint$Subject <- factor(meta_timepoint$Subject)
  meta_timepoint$Condition <- factor(meta_timepoint$Condition, levels = c("A", "V"))
  
  # Create DESeqDataSet
  dds_timepoint <- DESeqDataSetFromMatrix(countData = timepoint_counts,
                                          colData = meta_timepoint,
                                          design = ~ Subject + Condition)
  
  # Filter out lowly expressed genes
  timepoint_keep <- rowSums(counts(dds_timepoint) >= 10) >= (0.8 * ncol(timepoint_counts)) 
  dds_timepoint_keep <- dds_timepoint[timepoint_keep,]
  
  # Return the processed DESeqDataSet object
  return(dds_timepoint_keep)
}

# List of timepoints
timepoints <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")

# Initialize a list to store the results for each timepoint
dds_list <- list()

# Loop over the timepoints and apply the function
for (timepoint in timepoints) {
  dds_list[[timepoint]] <- process_timepoint(timepoint, meta, filtered_counts)
  # print the dimensions of the filtered dds to see the results
  print(paste(timepoint, "dimensions:", dim(dds_list[[timepoint]])))
}

#Check structure of individual timepoints
dds_Ex1 <- dds_list[["Ex1"]]
dim(dds_Ex1)
colData(dds_Ex1)
counts(dds_Ex1)

#Compare the DESeq2 object across timepoint to check for duplicate &/or unique genes 

# Extract gene names for each timepoint
gene_lists <- lapply(dds_list, rownames)

distinct_genes <- unlist(gene_lists) %>% 
  as.data.frame() %>% 
  rename("Ensembl" = ".") %>% 
  distinct(Ensembl) %>%
  pull(Ensembl)   #These are all the distinct genes that appear in at least 80% of individual timepoints 

filtered_counts_noOutliers80pct <- filtered_counts_noOutliers[rownames(filtered_counts_noOutliers) %in% distinct_genes,]


#Common genes across timepoints
common_genes <- Reduce(intersect, gene_lists)
print(length(common_genes))  # Number of common genes
 
#Unique genes across all timepoints 
all_genes <- Reduce(union, gene_lists)
print(length(all_genes))  # Total unique genes across all timepoints

#Unique genes across specific timepoint 
ex2_unique_genes <- setdiff(gene_lists[["Ex2"]], unlist(gene_lists[names(gene_lists) != "Ex2"]))
print(length(ex2_unique_genes))  

#Overlaps between two timepoints 
baseline_ex2_overlap <- intersect(gene_lists[["Baseline"]], gene_lists[["Ex2"]])
print(length(baseline_ex2_overlap))  # Number of overlapping genes


#SUMMARY TABLE TO SHOW UNIQUE AND OVERLAPPING GENES

#empty data frame to store results
summary_table <- data.frame(Timepoint = character(), 
                            Unique_Genes = numeric(), 
                            Common_Genes = numeric(), 
                            stringsAsFactors = FALSE)

# Calculate unique and common genes for each timepoint
for (timepoint in names(gene_lists)) {
  unique_genes <- setdiff(gene_lists[[timepoint]], unlist(gene_lists[names(gene_lists) != timepoint]))
  common_genes <- intersect(gene_lists[[timepoint]], unlist(gene_lists[names(gene_lists) != timepoint]))
  
  summary_table <- rbind(summary_table, 
                         data.frame(Timepoint = timepoint, 
                                    Unique_Genes = length(unique_genes), 
                                    Common_Genes = length(common_genes)))
}

print(summary_table)

library(ggVennDiagram)

#Gene sets for Venn Diagram comparison 
gene_sets <- list(
  Baseline = gene_lists[["Baseline"]],
  Exercise = unique(c(gene_lists[["Ex1"]], gene_lists[["Ex2"]], gene_lists[["Ex3"]])),
  Recovery = unique(c(gene_lists[["Rec1"]], gene_lists[["Rec2"]]))
)

ggVennDiagram(gene_sets) +
  scale_fill_gradient(low = "white", high = "blue", name = "Number of Genes") +
  theme_void()



#Normalization with spike ins. Each raw count has to be multiplied by sum of ERCC spikes? 

spike_sums <- setNames(summed_long$Spike_sum, summed_long$RNA_number)

filtered_counts_noOutliers80pct_NORM <- sweep(filtered_counts_noOutliers80pct, 2, spike_sums[colnames(filtered_counts_noOutliers80pct)], "*")


#Dynamic range density plot of genes 

gene_medians <- apply(filtered_counts_noOutliers, 1, median)
gene_medians
filtered_counts_medians <- filtered_counts_noOutliers_gene_symbol
filtered_counts_medians$RowMedian <- gene_medians
str(filtered_counts_medians)
filtered_counts_medians$RowMedian <- as.integer(filtered_counts_medians$RowMedian)
any(is.na(filtered_counts_medians$RowMedian))
unique(filtered_counts_medians$RowMedian)
filtered_counts_medians$RowMedian[is.na(filtered_counts_medians$RowMedian)] <- 0

total_counts_sum <- sum(filtered_counts_medians$RowMedian)
total_counts_sum

filtered_counts_medians$MedianRatio <- (filtered_counts_medians$RowMedian / total_counts_sum) * 100

gene_medians <- filtered_counts_medians[, c("gene_symbol", "RowMedian", "MedianRatio")]
gene_medians$gene_rank <- rank(-gene_medians$MedianRatio, ties.method = "average")
gene_medians$adjusted_rank <- ave(gene_medians$gene_rank, gene_medians$MedianRatio, 
                                  FUN = function(x) seq(min(x), max(x), length.out = length(x)))

gene_medians <- gene_medians %>% filter(MedianRatio != 0)
gene_medians$Log10MedianRatio <- log10(gene_medians$MedianRatio)
#62,000 genes in raw counts, only 14,000 without a Median >0 

gene_medians$gene_symbol[rownames(gene_medians) == "ENSG00000210082"] <- "MT-RNR2"
gene_medians$gene_symbol[rownames(gene_medians) == "ENSG00000211459"] <- "MT-RNR1"


ggplot(gene_medians, aes(x = gene_rank, y = Log10MedianRatio)) +
  geom_point(color = "blue") +  # Plot the points
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = Inf, y = 0, label = "1%", vjust = -0.5, hjust = 1.1, color = "black") +
  geom_hline(yintercept = -1, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = Inf, y = -1, label = "0.1%", vjust = -0.5, hjust = 1.1, color = "black") +
  geom_hline(yintercept = 0.7, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = Inf, y = 0.7, label = "5%", vjust = -0.5, hjust = 1.1, color = "black") +
  geom_hline(yintercept = -2, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = Inf, y = -2, label = "0.01%", vjust = -0.5, hjust = 1.1, color = "black") +
  geom_label_repel(data = gene_medians[order(gene_medians$gene_rank), ][1:5, ], aes(label = gene_symbol), box.padding = 1.3) +  
  geom_ribbon(data = gene_medians[gene_medians$gene_rank <= 9122, ], aes(x = gene_rank, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.1) +
  labs(x = "Gene Rank", y = "Log10 (% Total Counts)") +
  xlim (-100, 15000) +
  ylim (-5, 2) + 
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", size = 12),      
    axis.title = element_text(face = "bold", size = 14)     
  )



###Count number of genes per sample

genes_with_count10 <- apply(filtered_counts, 2, function(x) sum(x > 10))
genes_with_count10

df_genesabove10 <- data.frame(Sample = names(genes_with_count10), Count = genes_with_count10)

ggplot(df_genesabove10, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +  # Reduce width for spacing
  theme_classic() +
  labs(title = "Number of Genes with Count > 10", x = "Samples", y = "Gene Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = df_genesabove10$Sample[seq(1, nrow(df_genesabove10), by = 2)])



#Box plot for each sequin row --> look for inidivual outliers
just_spike_counts_noOutliers <- subset(just_spike_counts, select = -c(RNA031403, RNA031398, RNA031400))
sequin_counts <- just_spike_counts_noOutliers[grepl("^R[12]", rownames(just_spike_counts_noOutliers)), ]

for (i in 1:nrow(sequin_counts)) {
  # Extract and reshape the row data
  sequin_long <- data.frame(Value = as.numeric(sequin_counts[i, ]))  # Convert row to a long format
  
  # Calculate mean and 3 standard deviations
  row_mean <- mean(sequin_long$Value)
  row_sd <- sd(sequin_long$Value)
  upper_limit <- row_mean + (3 * row_sd)  # Mean + 3SD
  
  # Create the box plot with individual dots
  p <- ggplot(sequin_long, aes(x = "", y = Value, label = colnames(sequin_counts))) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue") +  # Box plot without default outliers
    geom_point(width = 0.1, size = 3, color = "red", alpha = 0.7) +  # Individual dots
    geom_text(aes(y = Value), vjust = -1, size = 4, color = "black") +  # Add RNA labels
    geom_hline(yintercept = upper_limit, linetype = "dashed", color = "blue", size = 1) +  # +3SD line
    labs(y = "Values", title = paste("Box Plot for", rownames(sequin_counts)[i])) +
    theme_minimal()
  
  # Print the plot
  print(p)
  
}


#Normalization factor = (sum of sequin for sample 1 / mean of sequin sum of all samples) x flow rate
sequin_sample_sums <- colSums(sequin_counts) %>%  as.data.frame() %>% t() 
rownames(sequin_sample_sums) <- c("Sequin sums")

sequin_mean_sum <- mean(as.numeric(sequin_sample_sums[1, ]))
sequin_mean_sum_row <- rep(sequin_mean_sum, ncol(sequin_sample_sums))

flow_data_repeated <- read_excel("/Users/MichalGrabowski/Documents/Master Thesis/AV_Flow_data_repeated.xlsx") %>% t()
flow <- flow_data_repeated[2, ]

rm(norm_factor_df)
norm_factor_df <- rbind(sequin_sample_sums, sequin_mean_sum_row)
norm_factor_df <- rbind(norm_factor_df, flow) %>% as.data.frame()
norm_factor_df[] <- lapply(norm_factor_df, function(x) as.numeric(as.character(x)))

norm_factor <- (norm_factor_df[1, ] / norm_factor_df[2, ])

# Append the new row to the dataframe
norm_factor_df <- rbind(norm_factor_df, norm_factor)
rownames(norm_factor_df)[nrow(norm_factor_df)] <- "Norm_factor"


########################### DESeq2 for within timepoint AV comparison ##################################

###########EX2

Ex2_RNAs <- meta_noOutliers %>% 
  filter(Time == "Ex2") %>% 
  pull(RNA_number)

Ex2_counts <- filtered_counts_noOutliers80pct[, c(Ex2_RNAs)]

meta_Ex2 <- meta[c(Ex2_RNAs),]
rownames(meta_Ex2) <- meta_Ex2$RNA_number

norm_factors_Ex2 <-norm_factor_df["Norm_factor", Ex2_RNAs] 
norm_factors_Ex2


# Factors 
meta_Ex2$Subject <- factor(meta_Ex2$Subject)
meta_Ex2$Condition <- factor(meta_Ex2$Condition, levels = c("A", "V"))

dds_Ex2 <- DESeqDataSetFromMatrix(countData = Ex2_counts,
                                       colData = meta_Ex2,
                                       design = ~ Subject + Condition)


sizeFactors(dds_Ex2) <- norm_factors_Ex2
sizeFactors(dds_Ex2)


Ex2AV <- DESeq(dds_Ex2)
res <- results(Ex2AV)
head(res)
Diff_expressionEx2_A_vs_V <- results(Ex2AV, contrast=c("Condition","A","V"))
Diff_expressionEx2_A_vs_V
Diff_expressionEx2_A_vs_Vordered <- Diff_expressionEx2_A_vs_V[order(res$pvalue),]
Diff_expressionEx2_A_vs_Vordered 

write.csv(as.data.frame(Diff_expressionEx2_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv")



######VOLCANO PLOTS EX2


Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)

#Add explicit gene_symbol column for easier processing
Ex2AV_results$Gene_symbol <- rownames(Ex2AV_results)


# add a column of NAs
Ex2AV_results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Ex2AV_results$diffexpressed[Ex2AV_results$log2FoldChange > 0.6 & Ex2AV_results$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Ex2AV_results$diffexpressed[Ex2AV_results$log2FoldChange < -0.6 & Ex2AV_results$pvalue < 0.05] <- "DOWN"


Ex2AV_results$label <- NA
Ex2AV_results_top_genes <- Ex2AV_results[order(-abs(Ex2AV_results$log2FoldChange), Ex2AV_results$pvalue), ][1:10, ]
Ex2AV_results_bottom_genes <- Ex2AV_results[order(Ex2AV_results$log2FoldChange, Ex2AV_results$pvalue), ][1:5, ]

# Add labels only for the top 10 genes
Ex2AV_results$label <- ifelse(
  rownames(Ex2AV_results) %in% rownames(Ex2AV_results_top_genes) | 
    rownames(Ex2AV_results) %in% rownames(Ex2AV_results_bottom_genes), 
  Ex2AV_results$Gene_symbol, 
  NA
)



ggplot(data = Ex2AV_results, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "black", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    title = "Volcano Plot: Ex2 A_V", 
    x = "Log2 Fold Change", 
    y = "-Log10 P-Value"
  ) +
  geom_text_repel(
    aes(fill = diffexpressed),
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
    color = "black",
    fontface = "bold",
    size = 3.5,
    min.segment.length = 0
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )
































################################################### LIMMMMMMMMAAAAAA ###########################################################################




design <- model.matrix(~ Subject + Time * Condition, data = meta_noOutliers)
colnames(design)
colnames(design) <- make.names(colnames(design))

dge <- DGEList(counts = filtered_counts_noOutliers80pct)
dge <- calcNormFactors(dge)
voom_output <- voom(dge, design, plot = TRUE)

# Extract normalized expression matrix
norm_counts <- voom_output$E

fit <- lmFit(voom_output, design)
fit <- eBayes(fit)


contrast.matrix <- makeContrasts(
  Ex2_vs_Ex1 = ((TimeEx2.ConditionV - TimeEx2) - 
                  (TimeEx1.ConditionV - TimeEx1)),
  #Baseline compared to Exercise
  Baseline_vs_Ex2 = ((TimeEx2.ConditionV - TimeEx2) -
                       (ConditionV - X.Intercept.)),
  #Baseline compared to Recovery
  Baseline_vs_Rec2 = ((TimeRec2.ConditionV - TimeRec2) - 
                        (ConditionV - X.Intercept.)),
  #Exercise compared to Recovery
  Ex2_vs_Rec2 = ((TimeRec2.ConditionV - TimeRec2) - 
                   (TimeEx2.ConditionV - TimeEx2)),
  
  levels = design
)


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get top DE genes for the contrast

Baseline_vs_Ex2_LimmaResults <- topTable(fit2, coef = "Baseline_vs_Ex2", adjust = "BH", number = Inf)
Baseline_vs_Ex2_LimmaResults

# Order by p-value?
Baseline_vs_Ex2_LimmaResults_ordered <- Baseline_vs_Ex2_LimmaResults[order(Baseline_vs_Ex2_LimmaResults$P.Value), ]
Baseline_vs_Ex2_LimmaResults_ordered

#How many significant genes? 
sig_genes <- topTable(fit2, coef = "Baseline_vs_Ex2", number = Inf, p.value = 0.05)
sig_genes



#Convert gene IDS to symbols now for easier processing

#Rownames to column for easier data manipulation 
Baseline_vs_Ex2_LimmaResults_geneid<- rownames_to_column(Baseline_vs_Ex2_LimmaResults_ordered, var = "gene_id")

gene_ids_Baseline_Ex2 <- Baseline_vs_Ex2_LimmaResults_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Baseline_Ex2 <- mapIds(org.Hs.eg.db, keys = gene_ids_Baseline_Ex2, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Baseline_Ex2 <- data.frame(
  gene_id = names(gene_symbols_Baseline_Ex2),
  gene_symbol = ifelse(is.na(gene_symbols_Baseline_Ex2), names(gene_symbols_Baseline_Ex2), gene_symbols_Baseline_Ex2),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Baseline_vs_Ex2_LimmaResults_merged <- Baseline_vs_Ex2_LimmaResults_geneid %>%
  left_join(symbols_data_frame_Baseline_Ex2, by = "gene_id")

# Reorder columns to place gene_symbol first
Baseline_vs_Ex2_LimmaResults_ordered_gene_symbol <- Baseline_vs_Ex2_LimmaResults_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Baseline_vs_Ex2_LimmaResults_ordered_gene_symbol <- Baseline_vs_Ex2_LimmaResults_ordered_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
Baseline_vs_Ex2_gene_symbol <- Baseline_vs_Ex2_LimmaResults_ordered_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
Baseline_vs_Ex2_gene_symbol <- Baseline_vs_Ex2_gene_symbol[, setdiff(names(Baseline_vs_Ex2_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(Baseline_vs_Ex2_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/Baseline_vs_Ex2_diff_expression.csv", row.names = TRUE)


#################################### VOLCANO PLOTS ######################################################


Baseline_vs_Ex2_Results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Baseline_vs_Ex2_diff_expression.csv", row.names = 1)

#Add explicit gene_symbol column for easier processing
Baseline_vs_Ex2_Results$Gene_symbol <- rownames(Baseline_vs_Ex2_Results)


# add a column of NAs
Baseline_vs_Ex2_Results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Baseline_vs_Ex2_Results$diffexpressed[Baseline_vs_Ex2_Results$logFC > 0.6 & Baseline_vs_Ex2_Results$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Baseline_vs_Ex2_Results$diffexpressed[Baseline_vs_Ex2_Results$logFC < -0.6 & Baseline_vs_Ex2_Results$P.Value < 0.05] <- "DOWN"


p <- ggplot(data = Baseline_vs_Ex2_Results, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  geom_point() + theme_minimal()

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p2 <- p + scale_colour_manual(values = mycolors)
p2

Baseline_vs_Ex2_Results$label <- NA
Baseline_vs_Ex2_Results_top_genes <- Baseline_vs_Ex2_Results[order(-abs(Baseline_vs_Ex2_Results$logFC), Baseline_vs_Ex2_Results$P.Value), ][1:10, ]
Baseline_vs_Ex2_Results_bottom_genes <- Baseline_vs_Ex2_Results[order(Baseline_vs_Ex2_Results$logFC, Baseline_vs_Ex2_Results$P.Value), ][1:5, ]

# Add labels only for the top 10 genes
Baseline_vs_Ex2_Results$label <- ifelse(
  rownames(Baseline_vs_Ex2_Results) %in% rownames(Baseline_vs_Ex2_Results_top_genes) | 
    rownames(Baseline_vs_Ex2_Results) %in% rownames(Baseline_vs_Ex2_Results_bottom_genes), 
  Baseline_vs_Ex2_Results$Gene_symbol, 
  NA
)


  
ggplot(data = Baseline_vs_Ex2_Results, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "black", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    title = "Volcano Plot: Baseline_vs_Ex2", 
    x = "Log2 Fold Change", 
    y = "-Log10 P-Value"
  ) +
  geom_text_repel(
    aes(fill = diffexpressed),
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.5,
    color = "black",
    fontface = "bold",
    size = 3.5,
    min.segment.length = 0
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )
  




















######Verify model AV behaviour by plotting some of the top genes AV difference per timepoint 

str(meta_noOutliers)
design <- model.matrix(~ Subject + Time * Condition, data = meta_noOutliers)
colnames(design)
colnames(design) <- make.names(colnames(design))

dge <- DGEList(counts = filtered_counts_noOutliers80pct)
dge <- calcNormFactors(dge)
voom_output <- voom(dge, design, plot = TRUE)

# Extract normalized expression matrix
norm_counts <- voom_output$E

fit <- lmFit(voom_output, design)
fit <- eBayes(fit)



contrast.matrix_timepoint <- makeContrasts(
  Baseline = ((X.Intercept.+ ConditionV) - X.Intercept.),
  Ex2 = ((X.Intercept.+ TimeEx2.ConditionV) - (X.Intercept.+ TimeEx2)),
  Rec2 = ((X.Intercept.+ TimeRec2.ConditionV) - (X.Intercept.+ TimeRec2)),
  levels = design
)


fit2 <- contrasts.fit(fit, contrast.matrix_timepoint)
fit2 <- eBayes(fit2)


#Extract results for each timepoint, fitler for only logFC and the top 5 genes 
top5_genes_BaselineEx2 <- c("ENSG00000210082", "ENSG00000276168", 
  "ENSG00000211459", "ENSG00000198888", 
  "ENSG00000205542")


BaselineAV_LimmaResults <- topTable(fit2, coef = "Baseline", adjust = "BH", number = Inf)
BaselineAV_LimmaLFC <- subset(BaselineAV_LimmaResults, select = "logFC")
BaselineAV_LimmaLFC_top5 <- BaselineAV_LimmaLFC[top5_genes_BaselineEx2, , drop = FALSE]


Ex2AV_LimmaResults <- topTable(fit2, coef = "Ex2", adjust = "BH", number = Inf)
Ex2AV_LimmaLFC <- subset(Ex2AV_LimmaResults, select = "logFC")
Ex2AV_LimmaLFC_top5 <- Ex2AV_LimmaLFC[top5_genes_BaselineEx2, , drop = FALSE]

# explore baseline samples
BaselineAV_LimmaResults %>% 
  filter(adj.P.Val < 0.05) %>%
  arrange(logFC)

voom_output$E %>% as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  filter(gene_id == "ENSG00000210082") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::select(voom_counts = 1) %>% 
  rownames_to_column(var = "RNA_number") %>% 
  filter(RNA_number != "gene_id") %>% 
  merge(., meta, by = "RNA_number") %>% 
  mutate(voom_counts = as.numeric(voom_counts)) %>% 
  group_by(Time, Condition) %>% 
  summarize(mean  = mean(voom_counts),
            sd = sd(voom_counts)) %>% 
  ggplot(aes(x = Time, y = mean, fill = Condition))+
  geom_bar(stat = "identity", position = position_dodge(1))+
  geom_errorbar(aes(ymax = mean+sd, ymin = mean-sd), position = position_dodge(1), width = 0.2)


# Combine the data into one data frame
BaselineAV_LimmaLFC_top5$TimePoint <- "Baseline"
Ex2AV_LimmaLFC_top5$TimePoint <- "Ex2"

# Add gene names as a column
BaselineAV_LimmaLFC_top5$Gene <- rownames(BaselineAV_LimmaLFC_top5)
Ex2AV_LimmaLFC_top5$Gene <- rownames(Ex2AV_LimmaLFC_top5)

# Merge the two datasets
BaselineAV_Ex2AV <- rbind(BaselineAV_LimmaLFC_top5, Ex2AV_LimmaLFC_top5)

# Bar chart
ggplot(BaselineAV_Ex2AV, aes(x = Gene, y = logFC, fill = TimePoint)) +
  geom_bar(stat = "identity", position = "dodge") +  # Grouped bars
  theme_minimal() +
  labs(title = "Log-Fold Change Across Time Points",
       x = "Gene",
       y = "Log-Fold Change (LFC)",
       fill = "Time Point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









#Go even further back... Plot just log normalized counts for A and V seperately at each timepoint 
#Work on just norm counts


norm_counts_top5genes <- norm_counts[top5_genes_BaselineEx2, , drop = FALSE]

norm_counts_top5genes_BaselineA <- norm_counts_top5genes[, colnames(norm_counts_top5genes) %in% RNA_numbers_Baseline_A]
norm_counts_top5genes_BaselineV <- norm_counts_top5genes[, colnames(norm_counts_top5genes) %in% RNA_numbers_Baseline_V]



top5_genes_BaselineEx2 <- c("ENSG00000210082", "ENSG00000276168", 
                            "ENSG00000211459", "ENSG00000198888", 
                            "ENSG00000205542")

ENSG00000210082 <- filtered_counts_noOutliers80pct["ENSG00000210082", ]
mean(as.numeric(ENSG00000210082))
sd(as.numeric(ENSG00000210082))

ENSG00000198888 <- filtered_counts_noOutliers80pct["ENSG00000198888", ]
mean(as.numeric(ENSG00000198888))
sd(as.numeric(ENSG00000198888))

ENSG00000205542 <- filtered_counts_noOutliers80pct["ENSG00000205542", ]
mean(as.numeric(ENSG00000205542))
sd(as.numeric(ENSG00000205542))





























#WHAT DO WE DO WITH THE FLOW RATE NOW?
#Multiplt raw counts by the flow rate before running LIMMA


#Volcano Plot
ggplot(Baseline_vs_Ex2_LimmaResults, aes(x = logFC, y = -log10(P.Value))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot: Baseline_vs_Ex2", x = "Log Fold Change", y = "-Log10 P-Value")

















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







#BASELINE

#Extract RNAs for baseline condition
Baseline_RNAs <- meta %>% 
  filter(Time == "Baseline") %>% 
  pull(RNA_number)

baseline_counts <- filtered_counts[, c(Baseline_RNAs)]

meta_baseline <- meta[c(Baseline_RNAs),]
rownames(meta_baseline) <- meta_baseline$RNA_number

# Factors baseline
meta_baseline$Subject <- factor(meta_baseline$Subject)
meta_baseline$Condition <- factor(meta_baseline$Condition, levels = c("A", "V"))

dds_baseline <- DESeqDataSetFromMatrix(countData = baseline_counts,
                                       colData = meta_baseline,
                                       design = ~ Subject + Condition)

baseline_keep <- rowSums(counts(dds_baseline) >= 10) >= (0.8*ncol(baseline_counts)) 
dds_baseline_keep <- dds_baseline[baseline_keep,]

dim(dds_baseline_keep)


