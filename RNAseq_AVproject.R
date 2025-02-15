library("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", Force = TRUE)

BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
library(readxl)
library (ggrepel)
library("dplyr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

#Combine all individual sample txt count files into 1
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

# Factors
meta$Time <- factor(meta$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))
meta$Subject <- factor(meta$Subject)
meta$Subject <- factor(meta$Subject)
meta$Condition <- factor(meta$Condition, levels = c("A", "V"))


rownames(meta) <- meta$RNA_number

raw_counts <- raw_counts %>%     
  column_to_rownames(var = "Gene_ID")


#Filtering: Cut rows that are not ENSG.... Gene ID
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

#Initial QC: PCA plot, outlier detection 

#DESeq data set construction
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = meta,
                              design = ~ Subject + Condition + Time + Condition:Time)



# Variance stbilized counts
vsd <- vst(dds, blind = T)

vsd_counts <- assay(vsd) %>% as.data.frame()

vsd_counts$var <- apply(vsd_counts, 1, var)

# Get the top 500 most variable points
vsd_counts_filtered <- vsd_counts %>% 
  arrange(-var) %>% 
  head(500) %>% 
  dplyr::select(-var)

#MANUAL PCA plot 
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

###TWO OUTLIER IDENTIFIED
#The two outliers are RNA031403 and RNA031398 ---> same subject, different Time Points 
#Remove identified RNA_number outliers 
filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA031403, RNA031398))
remove_RNA <- c("RNA031403", "RNA031398")
meta_noOutliers <- meta[!rownames(meta) %in% remove_RNA, ]
rownames(meta_noOutliers) <- meta_noOutliers$RNA_number



#REDO PCA AGAIN WITHOUT OUTLIERS --> DOUBLE CHECK WITH QC DATA (Deduplication Counts, Sequence length Distribution etc... )
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
#RNA031400 has a high duplication count and low reads overall --> can be removed on that matter? YES
remove_RNA <- c("RNA031403", "RNA031398", "RNA031400")
meta_noOutliers <- meta[!rownames(meta) %in% remove_RNA, ]

filtered_counts_noOutliers <- subset(filtered_counts, select = -c(RNA031403, RNA031398,RNA031400))

# what drives PC2 variability? 
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


#Overrepresentation analysis of PCA2 drivers 
go_enrich <- enrichGO(gene = gene_ids_PCA2$Entrezid, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "ALL",  # "BP" for Biological Process, "MF" for Molecular Function, or "CC" for Cellular Component
                      pvalueCutoff = 0.05)

head(go_enrich)
barplot(go_enrich, showCategory = 10)

#Pre Filtering, remove genes with low number of counts only keep counts above 10
keep <- rowSums(counts(dds_clean) >= 10) >= 90 #number of columns 
dds_keep <- dds_clean[keep,]
dim(dds_keep)


###Count number of genes per sample - to confirm poor sequencing depth issues on outlines 

genes_with_count10 <- apply(filtered_counts, 2, function(x) sum(x > 10))
genes_with_count10

df_genesabove10 <- data.frame(Sample = names(genes_with_count10), Count = genes_with_count10)

ggplot(df_genesabove10, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +  # Reduce width for spacing
  theme_classic() +
  labs(title = "Number of Genes with Count > 10", x = "Samples", y = "Gene Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = df_genesabove10$Sample[seq(1, nrow(df_genesabove10), by = 2)])


######QC graphs on spike-ins

#Filter raw count data for the Spikes 
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
#Spikes counts appear to be appropriately consistent across all samples 


#Now group by A/V and time point --> can spikes (somewhwat) normalize for the flow rate ? 
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


#Now group by A/V and time point --> can spikes (somewhwat) normalize for the flow rate ? 

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

#Closer inspection between two timepoints 
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

#We would expect the Spike counts to go down at the Exercising time points as they should be more diluted due to increased flow rate

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
#One notable outlier in ERCC spikes --> RNA031400, was removed previously. 


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


######### GENE COUNT MATRIX FILTERING #################################################

#Filter for genes that appear in at least 80% of the samples at given time point  
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
  
  return(dds_timepoint_keep)
}

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

#Venn Diagram showing unique/overlapping genes by Timepoint

library(ggVennDiagram)

#Gene sets for Venn Diagram comparison 
gene_sets <- list(
  Baseline = gene_lists[["Baseline"]],
  Exercise = unique(c(gene_lists[["Ex1"]], gene_lists[["Ex2"]], gene_lists[["Ex3"]])),
  Recovery = unique(c(gene_lists[["Rec1"]], gene_lists[["Rec2"]]))
)

ggVennDiagram(gene_sets) +
  scale_fill_gradient(low = "white", high = "red", name = "Number of Genes") +
  theme_void()


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



#Box plot for each sequin "row" --> are there any major individual outliers? 
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

#####NORMALIZATION FACTOR CALCULATION################################

#Normalization factor = (sum of sequin for sample # / mean of sequin sum of all samples)
sequin_sample_sums <- colSums(sequin_counts) %>%  as.data.frame() %>% t() 
rownames(sequin_sample_sums) <- c("Sequin sums")

sequin_mean_sum <- mean(as.numeric(sequin_sample_sums[1, ]))
sequin_mean_sum_row <- rep(sequin_mean_sum, ncol(sequin_sample_sums))

#Analyse sequin_sum distribution
hist(sequin_sample_sums)

sequin_median_sum <- median(as.numeric(sequin_sample_sums[1, ]))
sequin_median_sum_row <- rep(sequin_median_sum, ncol(sequin_sample_sums))

flow_data_repeated <- read_excel("/Users/MichalGrabowski/Documents/Master Thesis/AV_Flow_data_repeated.xlsx") %>% t()
flow <- flow_data_repeated[2, ]

rm(norm_factor_df)
norm_factor_df <- rbind(sequin_sample_sums, sequin_median_sum_row)
norm_factor_df <- rbind(norm_factor_df, flow) %>% as.data.frame()
norm_factor_df[] <- lapply(norm_factor_df, function(x) as.numeric(as.character(x)))

norm_factor <- (norm_factor_df[1, ] / norm_factor_df[2, ])

# Append the new row to the dataframe
norm_factor_df <- rbind(norm_factor_df, norm_factor)
rownames(norm_factor_df)[nrow(norm_factor_df)] <- "Norm_factor"

norm_factor <- as.matrix(norm_factor)
hist(norm_factor)
#One sample with an abnormally high norm_factor identified RNA031418
max(norm_factor)


###########################################################################################################
########################### DESeq2 for within timepoint AV comparison ##################################

###########BASELINE################################################

Baseline_RNAs <- meta_noOutliers %>% 
  filter(Time == "Baseline") %>% 
  pull(RNA_number)


Baseline_counts <- filtered_counts_noOutliers[, c(Baseline_RNAs)]
Baseline80_keep <- rowSums(Baseline_counts >= 10) >= (0.8 * ncol(Baseline_counts)) 
Baseline_counts80 <- Baseline_counts[Baseline80_keep,]


meta_Baseline <- meta[c(Baseline_RNAs),]
rownames(meta_Baseline) <- meta_Baseline$RNA_number

norm_factors_Baseline <-norm_factor_df["Norm_factor", Baseline_RNAs] 
norm_factors_Baseline
norm_factors_Baseline <- as.vector(as.numeric(norm_factors_Baseline[1, ]))


# Factors 
meta_Baseline$Subject <- factor(meta_Baseline$Subject)
meta_Baseline$Condition <- factor(meta_Baseline$Condition, levels = c("A", "V"))

dds_Baseline <- DESeqDataSetFromMatrix(countData = Baseline_counts80,
                                       colData = meta_Baseline,
                                       design = ~ Subject +  Condition)


sizeFactors(dds_Baseline) <- norm_factors_Baseline
sizeFactors(dds_Baseline)

BaselineAV <- DESeq(dds_Baseline)

#Spike normalized
Baseline_sequin_norm_counts <- counts(BaselineAV, normalized=TRUE)

res_baseline <- results(BaselineAV)
head(res_baseline)
Diff_expressionBaseline_A_vs_V <- results(BaselineAV, contrast=c("Condition","V","A"))
Diff_expressionBaseline_A_vs_V
Diff_expressionBaseline_A_vs_Vordered <- Diff_expressionBaseline_A_vs_V[order(res_baseline$pvalue),]
Diff_expressionBaseline_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionBaseline_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionBaseline_A_vs_Vordered$log2FoldChange > 0, 
                                                            10^-(sqrt(log10(1 / (Diff_expressionBaseline_A_vs_Vordered$pvalue^Diff_expressionBaseline_A_vs_Vordered$log2FoldChange))^2)), 
                                                            NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionBaseline_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionBaseline_A_vs_Vordered$log2FoldChange < 0, 
                                                            10^-(sqrt(log10(1 / (Diff_expressionBaseline_A_vs_Vordered$pvalue^Diff_expressionBaseline_A_vs_Vordered$log2FoldChange))^2)), 
                                                            NA)


write.csv(as.data.frame(Diff_expressionBaseline_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionBaseline_A_vs_V.csv")


BaselineAV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionBaseline_A_vs_V.csv", row.names = 1)
BaselineAV_results <- as.data.frame (Baseline_sequin_norm_counts)




########ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
BaselineAV_results_geneid<- rownames_to_column(BaselineAV_results, var = "gene_id")

gene_ids_Baseline <- BaselineAV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Baseline <- mapIds(org.Hs.eg.db, keys = gene_ids_Baseline, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Baseline <- data.frame(
  gene_id = names(gene_symbols_Baseline),
  gene_symbol = ifelse(is.na(gene_symbols_Baseline), names(gene_symbols_Baseline), gene_symbols_Baseline),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
BaselineAV_results_merged <- BaselineAV_results_geneid %>%
  left_join(symbols_data_frame_Baseline, by = "gene_id")

# Reorder columns to place gene_symbol first
BaselineAV_results_gene_symbol <- BaselineAV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
BaselineAV_results_gene_symbol <- BaselineAV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
BaselineAV_gene_symbol <- BaselineAV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
BaselineAV_gene_symbol <- BaselineAV_gene_symbol[, setdiff(names(BaselineAV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(BaselineAV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/BaselineAV_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS BASELINE


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
BaselineAV_gene_symbol 


#Add explicit gene_symbol column for easier processing
BaselineAV_gene_symbol$Gene_symbol <- rownames(BaselineAV_gene_symbol)


# add a column of NAs
BaselineAV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
BaselineAV_gene_symbol$diffexpressed[BaselineAV_gene_symbol$log2FoldChange > 0.6 & BaselineAV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
BaselineAV_gene_symbol$diffexpressed[BaselineAV_gene_symbol$log2FoldChange < -0.6 & BaselineAV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN "

#Get tup gene labels based on xiao score
BaselineAV_gene_symbol$label <- NA
BaselineAV_results_top_genes <- BaselineAV_gene_symbol[order(BaselineAV_gene_symbol$xiao_score_V), ][1:5, ]
BaselineAV_results_bottom_genes <- BaselineAV_gene_symbol[order(BaselineAV_gene_symbol$xiao_score_A), ][1:5, ]

# Add labels only for the top 10 genes
BaselineAV_gene_symbol$label <- ifelse(
  rownames(BaselineAV_gene_symbol) %in% rownames(BaselineAV_results_top_genes) | 
    rownames(BaselineAV_gene_symbol) %in% rownames(BaselineAV_results_bottom_genes), 
  BaselineAV_gene_symbol$Gene_symbol, 
  NA
)


BaselineAV_volcano <- ggplot(data = BaselineAV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    x = "Log2 Fold Change", 
    y = "-Log10 P-Value"
  ) +
  geom_text_repel(
    aes(color = diffexpressed),
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )

ggsave(filename = "BaselineAV_VOLCANO.png", plot = BaselineAV_volcano, 
       width = 6, height = 5, dpi = 600)



###FILTER FOR DIFF EXPRESSED GENES

BaselineAV_padj <- subset(BaselineAV_gene_symbol, padj < 0.05)
nrow(BaselineAV_padj)

BaselineAV_pval_V <- BaselineAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(BaselineAV_pval_V)

BaselineAV_pval_A <- BaselineAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(BaselineAV_pval_A)

BaselineAV_pval_LFC_V <- BaselineAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(BaselineAV_pval_LFC_V)

BaselineAV_pval_LFC_A <- BaselineAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(BaselineAV_pval_LFC_A)

BaselineAV_xiao_V <- subset(BaselineAV_gene_symbol, xiao_score_V < 0.05)
nrow(BaselineAV_xiao_V)

BaselineAV_xiao_A <- subset(BaselineAV_gene_symbol, xiao_score_A < 0.05)
nrow(BaselineAV_xiao_A)

##############################################################################
##############################################################
####################################################
###########PASSIVE##################

Passive_RNAs <- meta_noOutliers %>% 
  filter(Time == "Passive") %>% 
  pull(RNA_number)

Passive_counts <- filtered_counts_noOutliers[, c(Passive_RNAs)]
Passive80_keep <- rowSums(Passive_counts >= 10) >= (0.8 * ncol(Passive_counts)) 
Passive_counts80 <- Passive_counts[Passive80_keep,]


meta_Passive <- meta[c(Passive_RNAs),]
rownames(meta_Passive) <- meta_Passive$RNA_number

norm_factors_Passive <-norm_factor_df["Norm_factor", Passive_RNAs] 
norm_factors_Passive
norm_factors_Passive <- as.vector(as.numeric(norm_factors_Passive[1, ]))


# Factors 
meta_Passive$Subject <- factor(meta_Passive$Subject)
meta_Passive$Condition <- factor(meta_Passive$Condition, levels = c("A", "V"))

dds_Passive <- DESeqDataSetFromMatrix(countData = Passive_counts80,
                                       colData = meta_Passive,
                                       design = ~ Subject + Condition)


sizeFactors(dds_Passive) <- norm_factors_Passive
sizeFactors(dds_Passive)


PassiveAV <- DESeq(dds_Passive)

#Spike normalized
Passive_sequin_norm_counts <- counts(PassiveAV, normalized=TRUE)
PassiveAV@colData$sizeFactor


res_Passive <- results(PassiveAV)
head(res_Passive)
Diff_expressionPassive_A_vs_V <- results(PassiveAV, contrast=c("Condition","V","A"))
Diff_expressionPassive_A_vs_V
Diff_expressionPassive_A_vs_Vordered <- Diff_expressionPassive_A_vs_V[order(res_Passive$pvalue),]
Diff_expressionPassive_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionPassive_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionPassive_A_vs_Vordered$log2FoldChange > 0, 
                                                            10^-(sqrt(log10(1 / (Diff_expressionPassive_A_vs_Vordered$pvalue^Diff_expressionPassive_A_vs_Vordered$log2FoldChange))^2)), 
                                                            NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionPassive_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionPassive_A_vs_Vordered$log2FoldChange < 0, 
                                                            10^-(sqrt(log10(1 / (Diff_expressionPassive_A_vs_Vordered$pvalue^Diff_expressionPassive_A_vs_Vordered$log2FoldChange))^2)), 
                                                            NA)

write.csv(as.data.frame(Diff_expressionPassive_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionPassive_A_vs_V.csv")


PassiveAV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionPassive_A_vs_V.csv", row.names = 1)
PassiveAV_results <- as.data.frame (Passive_sequin_norm_counts)



PassiveAV_results_sequin <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionPassive_A_vs_V.csv", row.names = 1)

#Rownames to column for easier data manipulation 
PassiveAV_results_geneid<- rownames_to_column(PassiveAV_results, var = "gene_id")

PassiveAV_results_sequin<- rownames_to_column(PassiveAV_results_sequin, var = "gene_id")

merged_Passive <- merge(PassiveAV_results_geneid, PassiveAV_results_sequin, by = "gene_id", suffixes = c("_1","_2"))

cor_value <- cor(merged_Passive$log2FoldChange_1, merged_Passive$log2FoldChange_2)

ggplot(merged_Passive, aes(x = log2FoldChange_1, y = log2FoldChange_2)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "log2FoldChange - Dataset 1",
    y = "log2FoldChange - Dataset 2",
    title = paste("Correlation of log2FoldChange (r =", round(cor_value, 2), ")")
  ) +
  theme_minimal() +
  geom_abline()


###Take it one step back: pick a few genes and see how the AV difference behaves in spike and DESeq normalization 
dds_norm_Passive
dds_norm_Passive_sequin

gene_of_interest <- "ENSG00000149932"

dds_norm_Passive <- dds_norm_Passive[gene_of_interest,]
dds_norm_Passive_sequin <- dds_norm_Passive_sequin[gene_of_interest,]

rownames(meta_Passive) <- meta_Passive$RNA_number

df_plot1 <- data.frame(
  sample = names(dds_norm_Passive),           # the sample names (colnames of df1)
  counts = as.numeric(dds_norm_Passive)       # the count values for geneOfInterest
)


# Merge with metadata to get Condition (and other info if needed)
df_plot1$Condition <- meta_Passive[df_plot1$sample,"Condition",  drop = TRUE]


df_plot2 <- data.frame(
  sample = names(dds_norm_Passive_sequin),
  counts = as.numeric(dds_norm_Passive_sequin)
)

df_plot2$Condition <- meta_Passive[df_plot2$sample, "Condition", drop = TRUE]

df_plot1$Normalization <- "DESeq"
df_plot2$Normalization <- "Sequin"

df_combined <- rbind(df_plot1, df_plot2)

df_summarized <- df_combined %>%
  group_by(Normalization, Condition) %>%
  summarize(
    MeanCount = mean(counts),
    .groups = "drop"
  )

ggplot() +
  # 1) Bars: one per (Normalization, Condition) from df_summarized
  geom_bar(
    data = df_summarized,
    aes(x = Condition, y = MeanCount, fill = NA),
    stat = "identity",
    width = 0.6,
    color = "blue"
  ) +
  # 2) Individual points: from the original df_combined
  geom_jitter(
    data = df_combined,
    aes(x = Condition, y = counts, color = Condition),
    position = position_jitter(width = 0.15, height = 0),
    size = 2,
    alpha = 0.7
  ) +
  # Create one panel for Norm1, another for Norm2
  facet_wrap(~ Normalization, nrow = 1) +
  theme_bw() +
  labs(
    title = "Mean Counts ",
    x = "Condition",
    y = "Count"
  ) +
  theme(legend.position = "none")




########ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
PassiveAV_results_geneid<- rownames_to_column(PassiveAV_results, var = "gene_id")

gene_ids_Passive <- PassiveAV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Passive <- mapIds(org.Hs.eg.db, keys = gene_ids_Passive, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Passive <- data.frame(
  gene_id = names(gene_symbols_Passive),
  gene_symbol = ifelse(is.na(gene_symbols_Passive), names(gene_symbols_Passive), gene_symbols_Passive),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
PassiveAV_results_merged <- PassiveAV_results_geneid %>%
  left_join(symbols_data_frame_Passive, by = "gene_id")

# Reorder columns to place gene_symbol first
PassiveAV_results_gene_symbol <- PassiveAV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
PassiveAV_results_gene_symbol <- PassiveAV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
PassiveAV_gene_symbol <- PassiveAV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
PassiveAV_gene_symbol <- PassiveAV_gene_symbol[, setdiff(names(PassiveAV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(PassiveAV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/PassiveAV_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS PASSIVE


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
PassiveAV_gene_symbol 


#Add explicit gene_symbol column for easier processing
PassiveAV_gene_symbol$Gene_symbol <- rownames(PassiveAV_gene_symbol)


# add a column of NAs
PassiveAV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
PassiveAV_gene_symbol$diffexpressed[PassiveAV_gene_symbol$log2FoldChange > 0.6 & PassiveAV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
PassiveAV_gene_symbol$diffexpressed[PassiveAV_gene_symbol$log2FoldChange < -0.6 & PassiveAV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN"


PassiveAV_gene_symbol$label <- NA
PassiveAV_results_top_genes <- PassiveAV_gene_symbol[order(PassiveAV_gene_symbol$xiao_score_V), ][1:5, ]
PassiveAV_results_bottom_genes <- PassiveAV_gene_symbol[order(PassiveAV_gene_symbol$xiao_score_A), ][1:3, ]

# Add labels only for the top 10 genes
PassiveAV_gene_symbol$label <- ifelse(
  rownames(PassiveAV_gene_symbol) %in% rownames(PassiveAV_results_top_genes) | 
    rownames(PassiveAV_gene_symbol) %in% rownames(PassiveAV_results_bottom_genes), 
  PassiveAV_gene_symbol$Gene_symbol, 
  NA
)




PassiveAV_Volcano <- ggplot(data = PassiveAV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    x = "Log2 Fold Change", 
    y = "-Log10 P-Value"
  ) +
  geom_text_repel(
    aes(color = diffexpressed),
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )


ggsave(filename = "PassiveAV_VOLCANO.png", plot = PassiveAV_Volcano, 
       width = 6, height = 5, dpi = 600)



###FILTER FOR DIFF EXPRESSED GENES

PassiveAV_padj <- subset(PassiveAV_gene_symbol, padj < 0.05)
nrow(PassiveAV_padj)

PassiveAV_pval_V <- PassiveAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(PassiveAV_pval_V)

PassiveAV_pval_A <- PassiveAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(PassiveAV_pval_A)

PassiveAV_pval_LFC_V <- PassiveAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(PassiveAV_pval_LFC_V)

PassiveAV_pval_LFC_A <- PassiveAV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(PassiveAV_pval_LFC_A)

PassiveAV_xiao_V <- subset(PassiveAV_gene_symbol, xiao_score_V < 0.05)
nrow(PassiveAV_xiao_V)

PassiveAV_xiao_A <- subset(PassiveAV_gene_symbol, xiao_score_A < 0.05)
nrow(PassiveAV_xiao_A)



##################################################
###########################################
####################################
###########EX1##################

EX1_RNAs <- meta_noOutliers %>% 
  filter(Time == "Ex1") %>% 
  pull(RNA_number)

Ex1_counts <- filtered_counts_noOutliers[, c(EX1_RNAs)]
Ex1_80_keep <- rowSums(Ex1_counts >= 10) >= (0.8 * ncol(Ex1_counts)) 
Ex1_counts80 <- Ex1_counts[Ex1_80_keep,]

meta_EX1 <- meta[c(EX1_RNAs),]
rownames(meta_EX1) <- meta_EX1$RNA_number

norm_factors_EX1 <-norm_factor_df["Norm_factor", EX1_RNAs] 
norm_factors_EX1
norm_factors_EX1 <- as.vector(as.numeric(norm_factors_EX1[1, ]))


# Factors 
meta_EX1$Subject <- factor(meta_EX1$Subject)
meta_EX1$Condition <- factor(meta_EX1$Condition, levels = c("A", "V"))

dds_EX1 <- DESeqDataSetFromMatrix(countData = Ex1_counts80,
                                      colData = meta_EX1,
                                      design = ~ Subject + Condition)


sizeFactors(dds_EX1) <- norm_factors_EX1
sizeFactors(dds_EX1)


EX1AV <- DESeq(dds_EX1)

#Spike normalized
Ex1_sequin_norm_counts <- counts(EX1AV, normalized=TRUE)


res_EX1 <- results(EX1AV)
head(res_EX1)
Diff_expressionEX1_A_vs_V <- results(EX1AV, contrast=c("Condition","V","A"))
Diff_expressionEX1_A_vs_V
Diff_expressionEX1_A_vs_Vordered <- Diff_expressionEX1_A_vs_V[order(res_EX1$pvalue),]
Diff_expressionEX1_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionEX1_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionEX1_A_vs_Vordered$log2FoldChange > 0, 
                                                           10^-(sqrt(log10(1 / (Diff_expressionEX1_A_vs_Vordered$pvalue^Diff_expressionEX1_A_vs_Vordered$log2FoldChange))^2)), 
                                                           NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionEX1_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionEX1_A_vs_Vordered$log2FoldChange < 0, 
                                                           10^-(sqrt(log10(1 / (Diff_expressionEX1_A_vs_Vordered$pvalue^Diff_expressionEX1_A_vs_Vordered$log2FoldChange))^2)), 
                                                           NA)

write.csv(as.data.frame(Diff_expressionEX1_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEX1_A_vs_V.csv")


EX1AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEX1_A_vs_V.csv", row.names = 1)
Ex1AV_results <- as.data.frame (Ex1_sequin_norm_counts)


########ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
EX1AV_results_geneid<- rownames_to_column(EX1AV_results, var = "gene_id")

gene_ids_EX1 <- EX1AV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_EX1 <- mapIds(org.Hs.eg.db, keys = gene_ids_EX1, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_EX1 <- data.frame(
  gene_id = names(gene_symbols_EX1),
  gene_symbol = ifelse(is.na(gene_symbols_EX1), names(gene_symbols_EX1), gene_symbols_EX1),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
EX1AV_results_merged <- EX1AV_results_geneid %>%
  left_join(symbols_data_frame_EX1, by = "gene_id")

# Reorder columns to place gene_symbol first
EX1AV_results_gene_symbol <- EX1AV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
EX1AV_results_gene_symbol <- EX1AV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
EX1AV_gene_symbol <- EX1AV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
EX1AV_gene_symbol <- EX1AV_gene_symbol[, setdiff(names(EX1AV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(EX1AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/EX1AV_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS EX1


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
EX1AV_gene_symbol 


#Add explicit gene_symbol column for easier processing
EX1AV_gene_symbol$Gene_symbol <- rownames(EX1AV_gene_symbol)


# add a column of NAs
EX1AV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
EX1AV_gene_symbol$diffexpressed[EX1AV_gene_symbol$log2FoldChange > 0.6 & EX1AV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
EX1AV_gene_symbol$diffexpressed[EX1AV_gene_symbol$log2FoldChange < -0.6 & EX1AV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN "


EX1AV_gene_symbol$label <- NA
EX1AV_results_top_genes <- EX1AV_gene_symbol[order(EX1AV_gene_symbol$xiao_score_V), ][1:5, ]
EX1AV_results_bottom_genes <- EX1AV_gene_symbol[order(EX1AV_gene_symbol$xiao_score_A), ][1:5, ]
bottom_5 <- rownames(EX1AV_results_bottom_genes)[1:5]


# Add labels only for the top 10 genes
EX1AV_gene_symbol$label <- ifelse(
  test = rownames(EX1AV_gene_symbol) %in% bottom_5,
  yes = EX1AV_gene_symbol$Gene_symbol,   # or whichever column holds gene names
  no  = NA
)

EX1AV_gene_symbol$label <- ifelse(
  rownames(EX1AV_gene_symbol) %in% rownames(EX1AV_results_top_genes) | 
    rownames(EX1AV_gene_symbol) %in% rownames(EX1AV_results_bottom_genes), 
  PassiveAV_gene_symbol$Gene_symbol, 
  NA
)




Ex1AV_Volcano <- ggplot(data = EX1AV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    x = "Log2 Fold Change", 
    y = "-Log10 P-Value"
  ) +
  geom_text_repel(
    aes(color = diffexpressed),
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )

ggsave(filename = "Ex1AV_VOLCANO.png", plot = Ex1AV_Volcano, 
       width = 6, height = 5, dpi = 600)



###FILTER FOR DIFF EXPRESSED GENES BY CUTOFF

Ex1AV_padj <- subset(EX1AV_gene_symbol, padj < 0.05)
nrow(Ex1AV_padj)

Ex1AV_pval_V <- EX1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(Ex1AV_pval_V)

Ex1AV_pval_A <- EX1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(Ex1AV_pval_A)

Ex1AV_pval_LFC_V <- EX1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(Ex1AV_pval_LFC_V)

Ex1AV_pval_LFC_A <- EX1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(Ex1AV_pval_LFC_A)

Ex1AV_xiao_V <- subset(EX1AV_gene_symbol, xiao_score_V < 0.05)
nrow(Ex1AV_xiao_V)

Ex1AV_xiao_A <- subset(EX1AV_gene_symbol, xiao_score_A < 0.05)
nrow(Ex1AV_xiao_A)




##################################################
###########################################
####################################
###########EX2##################

Ex2_RNAs <- meta_noOutliers %>% 
  filter(Time == "Ex2") %>% 
  pull(RNA_number)

Ex2_counts <- filtered_counts_noOutliers[, c(Ex2_RNAs)]
Ex2_80_keep <- rowSums(Ex2_counts >= 10) >= (0.8 * ncol(Ex2_counts)) 
Ex2_counts80 <- Ex2_counts[Ex2_80_keep,]


meta_Ex2 <- meta[c(Ex2_RNAs),]
rownames(meta_Ex2) <- meta_Ex2$RNA_number

norm_factors_Ex2 <-norm_factor_df["Norm_factor", Ex2_RNAs] 
norm_factors_Ex2
norm_factors_Ex2 <- as.vector(as.numeric(norm_factors_Ex2[1, ]))


# Factors 
meta_Ex2$Subject <- factor(meta_Ex2$Subject)
meta_Ex2$Condition <- factor(meta_Ex2$Condition, levels = c("A", "V"))

dds_Ex2 <- DESeqDataSetFromMatrix(countData = Ex2_counts80,
                                       colData = meta_Ex2,
                                       design = ~ Subject + Condition)

sizeFactors(dds_Ex2) <- norm_factors_Ex2

Ex2AV <- DESeq(dds_Ex2)

#Spike normalized
Ex2_sequin_norm_counts <- counts(Ex2AV, normalized=TRUE)


res <- results(Ex2AV)
head(res)
Diff_expressionEx2_A_vs_V <- results(Ex2AV, contrast=c("Condition","V","A"))
Diff_expressionEx2_A_vs_V
Diff_expressionEx2_A_vs_Vordered <- Diff_expressionEx2_A_vs_V[order(res$pvalue),]
Diff_expressionEx2_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionEx2_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionEx2_A_vs_Vordered$log2FoldChange > 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionEx2_A_vs_Vordered$pvalue^Diff_expressionEx2_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionEx2_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionEx2_A_vs_Vordered$log2FoldChange < 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionEx2_A_vs_Vordered$pvalue^Diff_expressionEx2_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)

write.csv(as.data.frame(Diff_expressionEx2_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv")

Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
Ex2AV_results <- as.data.frame (Ex2_sequin_norm_counts)


########ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
Ex2AV_results_geneid<- rownames_to_column(Ex2AV_results, var = "gene_id")

gene_ids_Ex2 <- Ex2AV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Ex2 <- mapIds(org.Hs.eg.db, keys = gene_ids_Ex2, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Ex2 <- data.frame(
  gene_id = names(gene_symbols_Ex2),
  gene_symbol = ifelse(is.na(gene_symbols_Ex2), names(gene_symbols_Ex2), gene_symbols_Ex2),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Ex2AV_results_merged <- Ex2AV_results_geneid %>%
  left_join(symbols_data_frame_Ex2, by = "gene_id")

# Reorder columns to place gene_symbol first
Ex2AV_results_gene_symbol <- Ex2AV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Ex2AV_results_gene_symbol <- Ex2AV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
Ex2AV_gene_symbol <- Ex2AV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
Ex2AV_gene_symbol <- Ex2AV_gene_symbol[, setdiff(names(Ex2AV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
rownames(Ex2AV_gene_symbol)[rownames(Ex2AV_gene_symbol) == "ENSG00000173209"] <- "AHSA2"



#Save to CSV file
write.csv(Ex2AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/EX2AV_diff_expression_symbols.txt", row.names = TRUE)



######VOLCANO PLOTS EX2


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
Ex2AV_gene_symbol 


#Add explicit gene_symbol column for easier processing
Ex2AV_gene_symbol$Gene_symbol <- rownames(Ex2AV_gene_symbol)


# add a column of NAs
Ex2AV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Ex2AV_gene_symbol$diffexpressed[Ex2AV_gene_symbol$log2FoldChange > 0.6 & Ex2AV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Ex2AV_gene_symbol$diffexpressed[Ex2AV_gene_symbol$log2FoldChange < -0.6 & Ex2AV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN"


Ex2AV_gene_symbol$label <- NA
Ex2AV_results_top_genes <- Ex2AV_gene_symbol[order(Ex2AV_gene_symbol$xiao_score_V), ][1:5, ]
Ex2AV_results_bottom_genes <- Ex2AV_gene_symbol[order(Ex2AV_gene_symbol$xiao_score_A), ][1:5, ]

# Add labels only for the top 10 genes
Ex2AV_gene_symbol$label <- ifelse(
  rownames(Ex2AV_gene_symbol) %in% rownames(Ex2AV_results_top_genes) | 
    rownames(Ex2AV_gene_symbol) %in% rownames(Ex2AV_results_bottom_genes), 
  Ex2AV_gene_symbol$Gene_symbol, 
  NA
)



Ex2AV_Volcano <- ggplot(data = Ex2AV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
   # title = "Volcano Plot: Ex2 A_V", 
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )

ggsave(filename = "Ex2AV_VOLCANO.png", plot = Ex2AV_Volcano, 
       width = 6, height = 5, dpi = 600)

###FILTER FOR DIFF EXPRESSED GENES BY CUTOFF

Ex2AV_padj <- subset(Ex2AV_gene_symbol, padj < 0.05)
nrow(Ex2AV_padj)

Ex2AV_pval_V <- Ex2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(Ex2AV_pval_V)

Ex2AV_pval_A <- Ex2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(Ex2AV_pval_A)

Ex2AV_pval_LFC_V <- Ex2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(Ex2AV_pval_LFC_V)

Ex2AV_pval_LFC_A <- Ex2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(Ex2AV_pval_LFC_A)

Ex2AV_xiao_V <- subset(Ex2AV_gene_symbol, xiao_score_V < 0.05)
nrow(Ex2AV_xiao_V)

Ex2AV_xiao_A <- subset(Ex2AV_gene_symbol, xiao_score_A < 0.05)
nrow(Ex2AV_xiao_A)



##################################################
###########################################
####################################
###########EX3##################

EX3_RNAs <- meta_noOutliers %>% 
  filter(Time == "Ex3") %>% 
  pull(RNA_number)

Ex3_counts <- filtered_counts_noOutliers[, c(EX3_RNAs)]
Ex3_80_keep <- rowSums(Ex3_counts >= 10) >= (0.8 * ncol(Ex3_counts)) 
Ex3_counts80 <- Ex3_counts[Ex3_80_keep,]

meta_EX3 <- meta[c(EX3_RNAs),]
rownames(meta_EX3) <- meta_EX3$RNA_number

norm_factors_EX3 <-norm_factor_df["Norm_factor", EX3_RNAs] 
norm_factors_EX3
norm_factors_EX3 <- as.vector(as.numeric(norm_factors_EX3[1, ]))


# Factors 
meta_EX3$Subject <- factor(meta_EX3$Subject)
meta_EX3$Condition <- factor(meta_EX3$Condition, levels = c("A", "V"))

dds_EX3 <- DESeqDataSetFromMatrix(countData = Ex3_counts80,
                                  colData = meta_EX3,
                                  design = ~ Subject + Condition)


sizeFactors(dds_EX3) <- norm_factors_EX3
sizeFactors(dds_EX3)


EX3AV <- DESeq(dds_EX3)

#Spike normalized
Ex3_sequin_norm_counts <- counts(EX3AV, normalized=TRUE)

res_EX3 <- results(EX3AV)
head(res_EX3)
Diff_expressionEX3_A_vs_V <- results(EX3AV, contrast=c("Condition","V","A"))
Diff_expressionEX3_A_vs_V
Diff_expressionEX3_A_vs_Vordered <- Diff_expressionEX3_A_vs_V[order(res_EX3$pvalue),]
Diff_expressionEX3_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionEX3_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionEX3_A_vs_Vordered$log2FoldChange > 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionEX3_A_vs_Vordered$pvalue^Diff_expressionEX3_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionEX3_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionEX3_A_vs_Vordered$log2FoldChange < 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionEX3_A_vs_Vordered$pvalue^Diff_expressionEX3_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)

write.csv(as.data.frame(Diff_expressionEX3_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEX3_A_vs_V.csv")


EX3AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEX3_A_vs_V.csv", row.names = 1)
Ex3AV_results <- as.data.frame (Ex3_sequin_norm_counts)

########ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
Ex3AV_results_geneid<- rownames_to_column(EX3AV_results, var = "gene_id")

gene_ids_Ex3 <- Ex3AV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Ex3 <- mapIds(org.Hs.eg.db, keys = gene_ids_Ex3, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Ex3 <- data.frame(
  gene_id = names(gene_symbols_Ex3),
  gene_symbol = ifelse(is.na(gene_symbols_Ex3), names(gene_symbols_Ex3), gene_symbols_Ex3),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Ex3AV_results_merged <- Ex3AV_results_geneid %>%
  left_join(symbols_data_frame_Ex3, by = "gene_id")

# Reorder columns to place gene_symbol first
Ex3AV_results_gene_symbol <- Ex3AV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Ex3AV_results_gene_symbol <- Ex3AV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
Ex3AV_gene_symbol <- Ex3AV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
Ex3AV_gene_symbol <- Ex3AV_gene_symbol[, setdiff(names(Ex3AV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000260742"] <- "ITPRID2-DT"
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000237984"] <- "PTENP1"
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000290791"] <- "ENSG...290791"
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000267469"] <- "ENSG...267469"


#Save to CSV file
write.csv(Ex3AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/EX3AV_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS EX3


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
Ex3AV_gene_symbol 


#Add explicit gene_symbol column for easier processing
Ex3AV_gene_symbol$Gene_symbol <- rownames(Ex3AV_gene_symbol)


# add a column of NAs
Ex3AV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Ex3AV_gene_symbol$diffexpressed[Ex3AV_gene_symbol$log2FoldChange > 0.6 & Ex3AV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Ex3AV_gene_symbol$diffexpressed[Ex3AV_gene_symbol$log2FoldChange < -0.6 & Ex3AV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN"


Ex3AV_gene_symbol$label <- NA
Ex3AV_results_top_genes <- Ex3AV_gene_symbol[order(Ex3AV_gene_symbol$xiao_score_V), ][1:7, ]
Ex3AV_results_bottom_genes <- Ex3AV_gene_symbol[order(Ex3AV_gene_symbol$xiao_score_A), ][1:7, ]

# Add labels only for the top 10 genes
Ex3AV_gene_symbol$label <- ifelse(
  rownames(Ex3AV_gene_symbol) %in% rownames(Ex3AV_results_top_genes) | 
    rownames(Ex3AV_gene_symbol) %in% rownames(Ex3AV_results_bottom_genes), 
  Ex3AV_gene_symbol$Gene_symbol, 
  NA
)



Ex3AV_Volcano <- ggplot(data = Ex3AV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    # title = "Volcano Plot: Ex2 A_V", 
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )


ggsave(filename = "Ex3AV_VOLCANO.png", plot = Ex3AV_Volcano, 
       width = 6, height = 5, dpi = 600)


###FILTER FOR DIFF EXPRESSED GENES BY CUTOFF

EX3AV_padj <- subset(Ex3AV_gene_symbol, padj < 0.05)
nrow(EX3AV_padj)

Ex3AV_pval_V <- Ex3AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(Ex3AV_pval_V)

Ex3AV_pval_A <- Ex3AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(Ex3AV_pval_A)

Ex3AV_pval_LFC_V <- Ex3AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(Ex3AV_pval_LFC_V)

Ex3AV_pval_LFC_A <- Ex3AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(Ex3AV_pval_LFC_A)

Ex3AV_xiao_V <- subset(Ex3AV_gene_symbol, xiao_score_V < 0.05)
nrow(Ex3AV_xiao_V)

Ex3AV_xiao_A <- subset(Ex3AV_gene_symbol, xiao_score_A < 0.05)
nrow(Ex3AV_xiao_A)

##################################################
###########################################
####################################
###########REC1##################

Rec1_RNAs <- meta_noOutliers %>% 
  filter(Time == "Rec1") %>% 
  pull(RNA_number)

Rec1_counts <- filtered_counts_noOutliers[, c(Rec1_RNAs)]
Rec1_80_keep <- rowSums(Rec1_counts >= 10) >= (0.8 * ncol(Rec1_counts)) 
Rec1_counts80 <- Rec1_counts[Rec1_80_keep,]

meta_Rec1 <- meta[c(Rec1_RNAs),]
rownames(meta_Rec1) <- meta_Rec1$RNA_number

norm_factors_Rec1 <-norm_factor_df["Norm_factor", Rec1_RNAs] 
norm_factors_Rec1
norm_factors_Rec1 <- as.vector(as.numeric(norm_factors_Rec1[1, ]))


# Factors 
meta_Rec1$Subject <- factor(meta_Rec1$Subject)
meta_Rec1$Condition <- factor(meta_Rec1$Condition, levels = c("A", "V"))

dds_Rec1 <- DESeqDataSetFromMatrix(countData = Rec1_counts80,
                                  colData = meta_Rec1,
                                  design = ~ Subject + Condition)


sizeFactors(dds_Rec1) <- norm_factors_Rec1
sizeFactors(dds_Rec1)


Rec1AV <- DESeq(dds_Rec1)

#Spike normalized
Rec1_sequin_norm_counts <- counts(Rec1AV, normalized=TRUE)

res_Rec1 <- results(Rec1AV)
head(res_Rec1)
Diff_expressionRec1_A_vs_V <- results(Rec1AV, contrast=c("Condition","V","A"))
Diff_expressionRec1_A_vs_V
Diff_expressionRec1_A_vs_Vordered <- Diff_expressionRec1_A_vs_V[order(res_Rec1$pvalue),]
Diff_expressionRec1_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionRec1_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionRec1_A_vs_Vordered$log2FoldChange > 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionRec1_A_vs_Vordered$pvalue^Diff_expressionRec1_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionRec1_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionRec1_A_vs_Vordered$log2FoldChange < 0, 
                                                       10^-(sqrt(log10(1 / (Diff_expressionRec1_A_vs_Vordered$pvalue^Diff_expressionRec1_A_vs_Vordered$log2FoldChange))^2)), 
                                                       NA)


write.csv(as.data.frame(Diff_expressionRec1_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionRec1_A_vs_V.csv")


Rec1AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionRec1_A_vs_V.csv", row.names = 1)
Rec1AV_results <- as.data.frame (Rec1_sequin_norm_counts)


#######ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
Rec1AV_results_geneid<- rownames_to_column(Rec1AV_results, var = "gene_id")

gene_ids_Rec1 <- Rec1AV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Rec1 <- mapIds(org.Hs.eg.db, keys = gene_ids_Rec1, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Rec1 <- data.frame(
  gene_id = names(gene_symbols_Rec1),
  gene_symbol = ifelse(is.na(gene_symbols_Rec1), names(gene_symbols_Rec1), gene_symbols_Rec1),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Rec1AV_results_merged <- Rec1AV_results_geneid %>%
  left_join(symbols_data_frame_Rec1, by = "gene_id")

# Reorder columns to place gene_symbol first
Rec1AV_results_gene_symbol <- Rec1AV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Rec1AV_results_gene_symbol <- Rec1AV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
Rec1AV_gene_symbol <- Rec1AV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
Rec1AV_gene_symbol <- Rec1AV_gene_symbol[, setdiff(names(Rec1AV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(Rec1AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/Rec1AV_diff_expression_symbols.csv", row.names = TRUE)



######VOLCANO PLOTS REC1


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
Rec1AV_gene_symbol 


#Add explicit gene_symbol column for easier processing
Rec1AV_gene_symbol$Gene_symbol <- rownames(Rec1AV_gene_symbol)


# add a column of NAs
Rec1AV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Rec1AV_gene_symbol$diffexpressed[Rec1AV_gene_symbol$log2FoldChange > 0.6 & Rec1AV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Rec1AV_gene_symbol$diffexpressed[Rec1AV_gene_symbol$log2FoldChange < -0.6 & Rec1AV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN "


Rec1AV_gene_symbol$label <- NA
Rec1AV_results_top_genes <- Rec1AV_gene_symbol[order(Rec1AV_gene_symbol$xiao_score_V), ][1:5, ]
Rec1AV_results_bottom_genes <- Rec1AV_gene_symbol[order(Rec1AV_gene_symbol$xiao_score_A), ][1:5, ]

# Add labels only for the top 10 genes
Rec1AV_gene_symbol$label <- ifelse(
  rownames(Rec1AV_gene_symbol) %in% rownames(Rec1AV_results_top_genes) | 
    rownames(Rec1AV_gene_symbol) %in% rownames(Rec1AV_results_bottom_genes), 
  Rec1AV_gene_symbol$Gene_symbol, 
  NA
)



Rec1AV_volcano <- ggplot(data = Rec1AV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    # title = "Volcano Plot: Ex2 A_V", 
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )

ggsave(filename = "Rec1AV_VOLCANO.png", plot = Rec1AV_volcano, 
       width = 6, height = 5, dpi = 600)

###FILTER FOR DIFF EXPRESSED GENES BY CUTOFF

Rec1AV_padj <- subset(Rec1AV_gene_symbol, padj < 0.05)
nrow(Rec1AV_padj)

Rec1AV_pval_V <- Rec1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(Rec1AV_pval_V)

Rec1AV_pval_A <- Rec1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(Rec1AV_pval_A)

Rec1AV_pval_LFC_V <- Rec1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(Rec1AV_pval_LFC_V)

Rec1AV_pval_LFC_A <- Rec1AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(Rec1AV_pval_LFC_A)

Rec1AV_xiao_V <- subset(Rec1AV_gene_symbol, xiao_score_V < 0.05)
nrow(Rec1AV_xiao_V)

Rec1AV_xiao_A <- subset(Rec1AV_gene_symbol, xiao_score_A < 0.05)
nrow(Rec1AV_xiao_A)





##################################################
###########################################
####################################
###########REC2##################

Rec2_RNAs <- meta_noOutliers %>% 
  filter(Time == "Rec2") %>% 
  pull(RNA_number)

Rec2_counts <- filtered_counts_noOutliers[, c(Rec2_RNAs)]
Rec2_80_keep <- rowSums(Rec2_counts >= 10) >= (0.8 * ncol(Rec2_counts)) 
Rec2_counts80 <- Rec2_counts[Rec2_80_keep,]

meta_Rec2 <- meta[c(Rec2_RNAs),]
rownames(meta_Rec2) <- meta_Rec2$RNA_number

norm_factors_Rec2 <-norm_factor_df["Norm_factor", Rec2_RNAs] 
norm_factors_Rec2
norm_factors_Rec2 <- as.vector(as.numeric(norm_factors_Rec2[1, ]))


# Factors 
meta_Rec2$Subject <- factor(meta_Rec2$Subject)
meta_Rec2$Condition <- factor(meta_Rec2$Condition, levels = c("A", "V"))

dds_Rec2 <- DESeqDataSetFromMatrix(countData = Rec2_counts80,
                                   colData = meta_Rec2,
                                   design = ~ Subject + Condition)


sizeFactors(dds_Rec2) <- norm_factors_Rec2
sizeFactors(dds_Rec2)


Rec2AV <- DESeq(dds_Rec2)

#Spike normalized
Rec2_sequin_norm_counts <- counts(Rec2AV, normalized=TRUE)


res_Rec2 <- results(Rec2AV)
head(res_Rec2)
Diff_expressionRec2_A_vs_V <- results(Rec2AV, contrast=c("Condition","V","A"))
Diff_expressionRec2_A_vs_V
Diff_expressionRec2_A_vs_Vordered <- Diff_expressionRec2_A_vs_V[order(res_Rec2$pvalue),]
Diff_expressionRec2_A_vs_Vordered 
#Add XIAO score column just for positive LFC (Venous or UP) 
Diff_expressionRec2_A_vs_Vordered$xiao_score_V = ifelse(Diff_expressionRec2_A_vs_Vordered$log2FoldChange > 0, 
                                                        10^-(sqrt(log10(1 / (Diff_expressionRec2_A_vs_Vordered$pvalue^Diff_expressionRec2_A_vs_Vordered$log2FoldChange))^2)), 
                                                        NA)
#Add XIAO score column just for negative LFC (Arterial or DOWN ) 
Diff_expressionRec2_A_vs_Vordered$xiao_score_A = ifelse(Diff_expressionRec2_A_vs_Vordered$log2FoldChange < 0, 
                                                        10^-(sqrt(log10(1 / (Diff_expressionRec2_A_vs_Vordered$pvalue^Diff_expressionRec2_A_vs_Vordered$log2FoldChange))^2)), 
                                                        NA)

write.csv(as.data.frame(Diff_expressionRec2_A_vs_Vordered), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionRec2_A_vs_V.csv")


Rec2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionRec2_A_vs_V.csv", row.names = 1)
Rec2AV_results <- as.data.frame (Rec2_sequin_norm_counts)


#######ADD GENE SYMBOLS ################################


#Rownames to column for easier data manipulation 
Rec2AV_results_geneid<- rownames_to_column(Rec2AV_results, var = "gene_id")

gene_ids_Rec2 <- Rec2AV_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Rec2 <- mapIds(org.Hs.eg.db, keys = gene_ids_Rec2, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Rec2 <- data.frame(
  gene_id = names(gene_symbols_Rec2),
  gene_symbol = ifelse(is.na(gene_symbols_Rec2), names(gene_symbols_Rec2), gene_symbols_Rec2),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Rec2AV_results_merged <- Rec2AV_results_geneid %>%
  left_join(symbols_data_frame_Rec2, by = "gene_id")

# Reorder columns to place gene_symbol first
Rec2AV_results_gene_symbol <- Rec2AV_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Rec2AV_results_gene_symbol <- Rec2AV_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
Rec2AV_gene_symbol <- Rec2AV_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
Rec2AV_gene_symbol <- Rec2AV_gene_symbol[, setdiff(names(Rec2AV_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000210082"] <- "MT-RNR2"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"


#Save to CSV file
write.csv(Rec2AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/Rec2AV_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS REC2


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
Rec2AV_gene_symbol 


#Add explicit gene_symbol column for easier processing
Rec2AV_gene_symbol$Gene_symbol <- rownames(Rec2AV_gene_symbol)


# add a column of NAs
Rec2AV_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Rec2AV_gene_symbol$diffexpressed[Rec2AV_gene_symbol$log2FoldChange > 0.6 & Rec2AV_gene_symbol$pvalue < 0.05] <- "VENOUS HIGHER / UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Rec2AV_gene_symbol$diffexpressed[Rec2AV_gene_symbol$log2FoldChange < -0.6 & Rec2AV_gene_symbol$pvalue < 0.05] <- "ARTERIAL HIGHER / DOWN"


Rec2AV_gene_symbol$label <- NA
Rec2AV_results_top_genes <- Rec2AV_gene_symbol[order(Rec2AV_gene_symbol$xiao_score_V), ][1:5, ]
Rec2AV_results_bottom_genes <- Rec2AV_gene_symbol[order(Rec2AV_gene_symbol$xiao_score_A), ][1:5, ]

# Add labels only for the top 10 genes
Rec2AV_gene_symbol$label <- ifelse(
  rownames(Rec2AV_gene_symbol) %in% rownames(Rec2AV_results_top_genes) | 
    rownames(Rec2AV_gene_symbol) %in% rownames(Rec2AV_results_bottom_genes), 
  Rec2AV_gene_symbol$Gene_symbol, 
  NA
)



Rec2AV_volcano <- ggplot(data = Rec2AV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("red", "grey", "green"),
    name = "Differential expression"  # Add legend title
  ) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed', size = 0.2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed', size = 0.2) +
  labs(
    # title = "Volcano Plot: Ex2 A_V", 
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
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black", size = 0.8),  # Darker axis lines
    axis.ticks = element_line(color = "black", size = 0.8)  # Darker axis ticks
  )



ggsave(filename = "Rec2AV_VOLCANO.png", plot = Rec2AV_volcano, 
       width = 6, height = 5, dpi = 600)

###FILTER FOR DIFF EXPRESSED GENES BY CUTOFF

Rec2AV_padj <- subset(Rec2AV_gene_symbol, padj < 0.05)
nrow(Rec2AV_padj)

Rec2AV_pval_V <- Rec2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0)
nrow(Rec2AV_pval_V)

Rec2AV_pval_A <- Rec2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < 0)
nrow(Rec2AV_pval_A)

Rec2AV_pval_LFC_V <- Rec2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange > 0.6)
nrow(Rec2AV_pval_LFC_V)

Rec2AV_pval_LFC_A <- Rec2AV_gene_symbol %>% filter(pvalue < 0.05 & log2FoldChange < -0.6)
nrow(Rec2AV_pval_LFC_A)

Rec2AV_xiao_V <- subset(Rec2AV_gene_symbol, xiao_score_V < 0.05)
nrow(Rec2AV_xiao_V)

Rec2AV_xiao_A <- subset(Rec2AV_gene_symbol, xiao_score_A < 0.05)
nrow(Rec2AV_xiao_A)





#EXTRACT COMMON GENES ACROSS TIMEPOINTS 

DEG_list <- list(BaselineAV_gene_symbol, PassiveAV_gene_symbol, EX1AV_gene_symbol, Ex2AV_gene_symbol, Ex3AV_gene_symbol, Rec1AV_gene_symbol, Rec2AV_gene_symbol)


gene_lists <- lapply(DEG_list, function(df) {
  rownames(subset(df, log2FoldChange > 0.6 & pvalue < 0.05))
})

# Combine all gene identifiers into one vector.
all_genes <- unlist(gene_lists)

# Count how many timepoints each gene appears in.
gene_freq <- table(all_genes)

# Select genes that appear in at least 3 timepoints.
selected_genes <- names(gene_freq)[gene_freq >= 2]

# Display the selected genes.
print(selected_genes)


timepoints <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")



#Single gene graphs across timepoints 

R3HCC1_list <- lapply(seq_along(DEG_list), function(i) {
  df <- DEG_list[[i]]
  # Check if TSR3 is in the rownames; if not, assign NA
  if ("CHST8" %in% rownames(df)) {
    lfc <- df["CHST8", "log2FoldChange"]
  } else {
    lfc <- NA
  }
  # Create a small data frame
  data.frame(
    Timepoint = timepoints[i],
    log2FoldChange = lfc
  )
})

# Combine into one data frame
R3HCC1_data <- do.call(rbind, R3HCC1_list)

# Make sure Timepoint is a factor in the correct order
R3HCC1_data$Timepoint <- factor(R3HCC1_data$Timepoint, levels = timepoints)

ggplot(R3HCC1_data, aes(x = Timepoint, y = log2FoldChange)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "TSR3 log2FoldChange across Timepoints (Boxplot)",
    x = "Timepoint",
    y = "log2FoldChange"
  )



##### EXTRACT GENES WHICH ONLY INCREASE DURING EXERCISE

# Genes passing threshold at Baseline
Baseline_sig_genes <- rownames(subset(BaselineAV_gene_symbol, 
                                      log2FoldChange > 0.6 & pvalue < 0.05))


Passive_sig_genes <- rownames(subset(PassiveAV_gene_symbol, 
                                     log2FoldChange > 0.6 & pvalue < 0.05))

# Genes passing threshold at Ex1
Ex1_sig_genes <- rownames(subset(EX1AV_gene_symbol, 
                                 log2FoldChange > 0.6 & pvalue < 0.05))

# Genes passing threshold at Ex2
Ex2_sig_genes <- rownames(subset(Ex2AV_gene_symbol, 
                                 log2FoldChange > 0.6 & pvalue < 0.05))

# Genes passing threshold at Ex3
Ex3_sig_genes <- rownames(subset(Ex3AV_gene_symbol, 
                                 log2FoldChange > 0.6 & pvalue < 0.05))

Excluded_genes <- union(Baseline_sig_genes, Passive_sig_genes)


Ex_all <- union(Ex1_sig_genes, union(Ex2_sig_genes, Ex3_sig_genes))

Only_Ex_genes <- setdiff(Ex_all, Excluded_genes)

length(Only_Ex_genes)
head(Only_Ex_genes)

Ex_Rec <- intersect(Rec_all, Ex_all)

setdiff(Ex_Rec, Excluded_genes)


##### EXTRACT GENES WHICH ONLY INCREASE DURING RECOVERY


Rec1_sig_genes <- rownames(subset(Rec1AV_gene_symbol, 
                                 log2FoldChange > 0.6 & pvalue < 0.05))

Rec2_sig_genes <- rownames(subset(Rec2AV_gene_symbol, 
                                 log2FoldChange > 0.6 & pvalue < 0.05))

Rec_all <- union(Rec1_sig_genes, Rec2_sig_genes)

Common_rec_genes <- intersect(Rec1_sig_genes, Rec2_sig_genes)

Excluded_genes <- union(Ex_all, union(Baseline_sig_genes, Passive_sig_genes))


Only_Rec_genes <- setdiff(Rec_all, Excluded_genes)

length(Only_Rec_genes)
head(Only_Rec_genes)

rownames

for (gene in Only_Ex_genes) {
  
  # Extract log2FoldChange for the current gene at each timepoint
  gene_data_list <- lapply(seq_along(DEG_list), function(i) {
    df <- DEG_list[[i]]
    if (gene %in% rownames(df)) {
      lfc <- df[gene, "log2FoldChange"]
    } else {
      lfc <- NA  # or 0, or some default if gene not found
    }
    data.frame(
      Timepoint = timepoints[i],
      log2FoldChange = lfc
    )
  })
  
  # Combine into a single data frame
  gene_data <- do.call(rbind, gene_data_list)
  
  # Set the factor level ordering for timepoints
  gene_data$Timepoint <- factor(gene_data$Timepoint, levels = timepoints)
  
  # Create the plot
  p <- ggplot(gene_data, aes(x = Timepoint, y = log2FoldChange)) +
    geom_boxplot(fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(
      title = paste(gene, "log2FoldChange across Timepoints (Boxplot)"),
      x = "Timepoint",
      y = "log2FoldChange"
    )
  
  # Display the plot
  print(p)
  

}


#UTP20 - LARGE INCREASE IN RECOVERY
#PHF5A
#IL1B
#VPS18
#MRPL40
#IPO10, IL1R1, MRPL54 - ONLY INCREASED (OR EVEN PRESENT) DURING RECOVERY 
#FLYWCH1



#Heatmap of the Only Ex genes 

ExGenes_mat <- cbind(
  Baseline = BaselineAV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Passive  = PassiveAV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Ex1      = EX1AV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Ex2      = Ex2AV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Ex3      = Ex3AV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Rec1     = Rec1AV_gene_symbol[Only_Ex_genes, "log2FoldChange"],
  Rec2     = Rec2AV_gene_symbol[Only_Ex_genes, "log2FoldChange"]
)
rownames(ExGenes_mat) <- Only_Ex_genes


Time_order <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
ExGenes_mat <- ExGenes_mat[, Time_order]
head(rownames(ExGenes_mat))


library(pheatmap)


# Define a color scale from red (negative) to white (zero) to darkgreen (positive)
my_palette <- colorRampPalette(c("red", "white", "darkgreen"))(100)

pheatmap(
  ExGenes_mat,
  color = my_palette, 
  scale = "none",            # Optionally scale each gene's row to highlight relative changes
  cluster_rows = FALSE,      # or FALSE, depending on whether you want clustering
  cluster_cols = FALSE, 
  na_col        = "grey",
  show_rownames = TRUE,
  angle_col     = 45 # Turn off if you have many genes
)


#Heatmap of the Only Rec genes 

RecGenes_mat <- cbind(
  Baseline = BaselineAV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Passive  = PassiveAV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Ex1      = EX1AV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Ex2      = Ex2AV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Ex3      = Ex3AV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Rec1     = Rec1AV_gene_symbol[Only_Rec_genes, "log2FoldChange"],
  Rec2     = Rec2AV_gene_symbol[Only_Rec_genes, "log2FoldChange"]
)
rownames(RecGenes_mat) <- Only_Rec_genes


#ExGenes_mat[is.na(ExGenes_mat)] <- 0
Time_order <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
RecGenes_mat <- RecGenes_mat[, Time_order]
head(rownames(RecGenes_mat))

pheatmap(
  RecGenes_mat,
  color = my_palette, 
  scale = "none",            # Optionally scale each gene's row to highlight relative changes
  cluster_rows = FALSE,      # or FALSE, depending on whether you want clustering
  cluster_cols = FALSE, 
  na_col        = "grey",
  show_rownames = TRUE,
  angle_col     = 45 # Turn off if you have many genes
)



###VENN DIAGRAM OF SIGNIFICANCE GENES ACROSS BASELINE EXERCISE RECOVERY 


sig_DEG_list <- list(Baseline_sig_genes, Ex1_sig_genes, Ex2_sig_genes, Ex3_sig_genes, Rec1_sig_genes, Rec2_sig_genes)



sig_venn_list <- list(
  "Baseline" = Baseline_sig_genes,
  "Exercise" = Ex_all,
  "Recovery" = Rec_all
)


ggVennDiagram(sig_venn_list) +
  scale_fill_gradient(low="white", high="blue",  name = "Number of differentially 
 expressed RNAs") +
  theme_void()



#Graph summarizing significant gene numbers 

Sig_results <- read_xlsx("/Users/MichalGrabowski/Documents/Master Thesis/DiffExGenes.xlsx")

Sig_long <- pivot_longer(
  Sig_results,
  cols = c("VENOUS HIGHER / UP", "ARTERIAL HIGHER / DOWN"),
  names_to = "Category",
  values_to = "Value"
)

Sig_results  <- pivot_longer(Sig_results, cols = c("Total_genes", "VENOUS HIGHER / UP", "ARTERIAL HIGHER / DOWN"), names_to = "Category",
                            values_to = "Value")
Sig_long$Time <- factor(Sig_long$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))



ggplot(Sig_long, aes(x = Time, y = Value, group = Category, color = Category)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_text(data = Sig_results, aes(x = Time, y = max(Sig_long$Value) + 200, label = Total_genes), inherit.aes = FALSE, size = 4) + 
  theme_minimal() +
  labs(
    x = "Timepoint",
    y = "Number of Genes",
    color = NULL
  ) +
  scale_color_manual(values = c("VENOUS HIGHER / UP" = "green", "ARTERIAL HIGHER / DOWN" = "red")) +
  scale_y_continuous(breaks = seq(0, max(Sig_long$Value), by = 200)) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold",angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),                      # Remove gridlines
    axis.line = element_line(size = 1.2, color = "black"),  # Thicker axis lines
    axis.ticks = element_line(size = 1.2, color = "black")  # Thicker tick marks
  ) +
  geom_hline(yintercept = max(Sig_long$Value) + 150, linetype = "dashed", color = "gray") +
  annotate("text", x = 3.5, y = max(Sig_long$Value) + 350, label = "Total Genes", size = 5, fontface = "bold")


##GSEA BASELINE

Baseline_GSEA <- read.csv ("~/Documents/Master Thesis/BaselineAV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Baseline_GSEA$stat  
names(geneList) <- Baseline_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_BASELINE <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  verbose      = FALSE
)


GSEA_Baseline_df <- as.data.frame(GSEA_BASELINE@result) %>%
  arrange(p.adjust) %>%  # Sort by adjusted p-value
  head(10)

ridgeplot(GSEA_BASELINE, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")



##GENE SET ENRICHMENT ANALYSIS ON PASSIVE 

library(clusterProfiler)
library(org.Hs.eg.db)

Passive_GSEA <- read.csv ("~/Documents/Master Thesis/PassiveAV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Passive_GSEA$stat  
names(geneList) <- Passive_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_PASSIVE <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)


ridgeplot(GSEA_PASSIVE, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")






Passive_GSEA <- read.csv ("~/Documents/Master Thesis/PassiveAV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(xiao_score_V) %>%
  dplyr::rename(gene = X)

Passive_UP <- Passive_GSEA %>% dplyr::filter(xiao_score_V < 0.05) 

GSEA_Passive_UP <- enrichGO(
  gene          = Passive_UP$gene,
  universe      = Passive_GSEA$gene,
  keyType       = "SYMBOL",
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE)


barplot(GSEA_Passive_UP)
heatplot (GSEA_Passive_UP)


Passive_GSEA_plotdata <- GSEA_Passive_UP %>%
  as.data.frame() %>%
  mutate(
    GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),  # Convert GeneRatio
    LogQ = -log10(p.adjust)  # Negative log10 of adj pvalue
  )

ggplot(Passive_GSEA_plotdata, aes(x = Count, y = reorder(Description, Count), size = GeneRatio, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "purple", high = "red", name = "p.adjust") +
  labs(
    x = "Gene Count",
    y = "Description",
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),  # Slightly darker grey plot background
    axis.text.y = element_text(size = 12, face = "bold"),# Adjust y-axis text size
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )+ 
coord_fixed(ratio = 8)

GSEA_passive_UP_df <- as.data.frame(GSEA_Passive_UP)


write_csv(GSEA_passive_UP_df,
          file="~/Documents/Master Thesis/GSEA_PassiveUP.csv")



#Very few significant enriched terms, simplification not necessary 


# Simplify GO terms redundancy --------------------------------------------
GSEA_passive_UP_simplify <- simplify(GSEA_Passive_UP, cutoff=0.7, by="p.adjust", select_fun=min)

# Get ORA result in dataframe --------------------------------------------
GSEA_passive_UP_simplify_df <- as.data.frame(GSEA_passive_UP_simplify)

# Fold enrichment:  ratio of the frequency of input genes annotated in a term to the frequency of all genes annotated to that term --------------------------------------------
GSEA_passive_UP_simplify_df <- mutate(GSEA_passive_UP_simplify_df, foldEnrich =
                                    (as.numeric(sub("/\\d+", "", GSEA_passive_UP_simplify_df$GeneRatio)) / as.numeric(sub(".*/", "", GSEA_passive_UP_simplify_df$GeneRatio))) /
                                    (as.numeric(sub("/\\d+", "", GSEA_passive_UP_simplify_df$BgRatio)) / as.numeric(sub(".*/", "", GSEA_passive_UP_simplify_df$BgRatio)))
)

# Save results ------------------------------------------------------------
write_csv(GSEA_Passive_UP,
          file="~/Documents/Master Thesis/GSEA_PassiveUP.csv")

barplot(GSEA_passive_UP_simplify)


##GSEA EX1

Ex1_GSEA <- read.csv ("~/Documents/Master Thesis/EX1AV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Ex1_GSEA$stat  
names(geneList) <- Ex1_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_EX1 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)


ridgeplot(GSEA_EX1, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")



#### GSEA ON EX2

Ex2_GSEA <- read.csv ("~/Documents/Master Thesis/EX2AV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)

Ex2_DOWN <- Ex2_GSEA %>% dplyr::filter(xiao_score_A < 0.05) 

GSEA_Ex2_DOWN <- enrichGO(
  gene          = Ex2_DOWN$gene,
  universe      = Ex2_GSEA$gene,
  keyType       = "SYMBOL",
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)



Ex2_GSEA_plotdata <- GSEA_Ex2_DOWN %>%
  as.data.frame() %>%
  mutate(
    GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),  # Convert GeneRatio
    LogQ = -log10(p.adjust)  # Negative log10 of adj pvalue
  )

ggplot(Ex2_GSEA_plotdata, aes(x = Count, y = reorder(Description, Count), size = GeneRatio, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "purple", high = "red", name = "p.adjust") +
  labs(
    x = "Gene Count",
    y = "Description",
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),  # Slightly darker grey plot background
    axis.text.y = element_text(size = 12, face = "bold"),# Adjust y-axis text size
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )+ 
  coord_fixed(ratio = 8)






geneList <- Ex2_GSEA$stat  
names(geneList) <- Ex2_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)

anyDuplicated(names(geneList))  



GSEA_Ex2_DOWN <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.3,
  verbose      = FALSE
)


GSEA_EX2_df <- as.data.frame(GSEA_Ex2_DOWN@result) %>%
  arrange(p.adjust) %>%  # Sort by adjusted p-value
  head(10)

GSEA_EX2_df$Count <- sapply(strsplit(GSEA_EX2_df$core_enrichment, "/"), length)
GSEA_EX2_df$GeneRatio <- GSEA_EX2_df$Count / GSEA_EX2_df$setSize

GSEA_EX2_df$mean_LFC <- sapply(strsplit(GSEA_EX2_df$core_enrichment, "/"), function(genes) {
  mean(Ex2_GSEA$log2FoldChange[Ex2_GSEA$gene %in% genes], na.rm = TRUE)  # Match genes & compute mean LFC
})

# Check if values are correctly assigned
head(GSEA_EX2_df[, c("Description", "mean_LFC")])





ggplot(GSEA_EX2_df, aes(x = Count, y = reorder(Description, Count), size = GeneRatio, color = p.adjust)) +
  geom_point() +  # Create the dot plot
  scale_color_gradient(low = "purple", high = "red") +  # Adjust color scale
  labs(title = "Top Enriched GO Terms in Exercised Muscle",
       x = "Gene Count",
       y = "GO Term",
       color = "p.adjust",
       size = "GeneRatio") +
  theme_minimal() +  # Clean theme
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))




ridgeplot(GSEA_Ex2_DOWN, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")



##GSEA EX3

Ex3_GSEA <- read.csv ("~/Documents/Master Thesis/EX3AV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Ex3_GSEA$stat  
names(geneList) <- Ex3_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_EX3 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)


ridgeplot(GSEA_EX3, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")

#GSEA REC1 

Rec1_GSEA <- read.csv ("~/Documents/Master Thesis/Rec1AV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Rec1_GSEA$stat  
names(geneList) <- Rec1_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


set.seed(1)

GSEA_Rec1 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  verbose      = FALSE
)


ridgeplot(GSEA_Rec1, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")


#GSEA REC2

Rec2_GSEA <- read.csv ("~/Documents/Master Thesis/Rec2AV_diff_expression_symbols.csv", header = T) %>%
  drop_na(log2FoldChange) %>%
  drop_na(stat) %>%
  dplyr::rename(gene = X)


geneList <- Rec2_GSEA$stat  
names(geneList) <- Rec2_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_Rec2 <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 10,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 1.0,
  verbose      = FALSE
)


ridgeplot(GSEA_Rec2, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")



######## IMMUNE DECONVOLUTION ################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GSVA")
BiocManager::install("xCell")
BiocManager::install("sva")
BiocManager::install("preprocessCore")

devtools::install_github("IOBR/IOBR")
BiocManager::install("ComplexHeatmap")
library(IOBR)
library(reshape2)
library(car)

Ex2AV_gene_symbol 

# Compute CPM
dge <- DGEList(counts = Ex2AV_gene_symbol)
cpm_matrix <- cpm(dge)

write.table(cpm_matrix, file = "/Users/MichalGrabowski/Documents/Master Thesis/EX2_norm_counts_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


data("lm22")
lm22

EX2AV_tpm <- count2tpm(
  Ex2AV_results,  # Now using Ensembl IDs
  idType = "Ensembl",
  org = "hsa",
  source = "local",
  effLength = NULL,  # Uses built-in IOBR gene lengths
  id = "row.names",  # Ensembl IDs are now row names
  check_data = FALSE
)


EX2AV_cibersort <- CIBERSORT(
  sig_matrix = lm22,
  EX2AV_tpm,
  perm = 100,
  QN = TRUE,
  absolute = FALSE,
  abs_method = "sig.score"
)

EX2AV_cibersort <- EX2AV_cibersort %>% as.data.frame()

group_lymphocytes <- grep("^(B cells|T cells|Plasma cells|NK cells)", colnames(EX2AV_cibersort), value = TRUE)
group_monocytes <- grep("^(Mono|Macro)", colnames(EX2AV_cibersort), value = TRUE)
group_DendriticMast <- grep("^(Dend|Mast)", colnames(EX2AV_cibersort), value = TRUE)

EX2AV_cibersort[, group_lymphocytes] <- lapply(EX2AV_cibersort[, group_lymphocytes], as.numeric)
EX2AV_cibersort[, group_monocytes] <- lapply(EX2AV_cibersort[, group_monocytes], as.numeric)
EX2AV_cibersort[, group_DendriticMast] <- lapply(EX2AV_cibersort[, group_DendriticMast], as.numeric)


EX2AV_cibersort$Lymphocytes <- rowSums(EX2AV_cibersort[, group_lymphocytes, drop = FALSE], na.rm = TRUE)
EX2AV_cibersort$Monocytes <- rowSums(EX2AV_cibersort[, group_monocytes, drop = FALSE], na.rm = TRUE)
EX2AV_cibersort$Dendritic_Mast <- rowSums(EX2AV_cibersort[, group_DendriticMast, drop = FALSE], na.rm = TRUE)

EX2AV_cibersort <- EX2AV_cibersort[, !(colnames(EX2AV_cibersort) %in% c(group_lymphocytes, group_monocytes, group_DendriticMast))]



timepoints <- c("BaselineAV", "PassiveAV", "Ex1AV", "Ex2AV", "Ex3AV", "Rec1AV", "Rec2AV")

# Create an empty list to store results
cibersort_results <- list()

# Loop through each timepoint
for (tp in timepoints) {
  
  # Dynamically get the count matrix variable using `get()`
  count_matrix <- get(paste0(tp, "_results"))
  
  # Convert count matrix to TPM
  tpm_matrix <- count2tpm(
    count_matrix, 
    idType = "Ensembl",
    org = "hsa",
    source = "local",
    effLength = NULL,  # Uses built-in IOBR gene lengths
    id = "row.names",
    check_data = FALSE
  )
  
  # Run CIBERSORT analysis
  cibersort_output <- CIBERSORT(
    sig_matrix = lm22,
    mixture_file = tpm_matrix,
    perm = 500,
    QN = TRUE,
    absolute = FALSE,
    abs_method = "sig.score"
  )
  
  # Convert to data frame
  cibersort_output <- as.data.frame(cibersort_output)
  
  # Identify immune cell groups
  group_lymphocytes <- grep("^(B cells|T cells|Plasma cells|NK cells)", colnames(cibersort_output), value = TRUE)
  group_monocytes <- grep("^(Mono|Macro)", colnames(cibersort_output), value = TRUE)
  group_DendriticMast <- grep("^(Dend|Mast)", colnames(cibersort_output), value = TRUE)
  
  # Convert relevant columns to numeric
  cibersort_output[, group_lymphocytes] <- lapply(cibersort_output[, group_lymphocytes], as.numeric)
  cibersort_output[, group_monocytes] <- lapply(cibersort_output[, group_monocytes], as.numeric)
  cibersort_output[, group_DendriticMast] <- lapply(cibersort_output[, group_DendriticMast], as.numeric)
  
  # Compute row sums for immune groups
  cibersort_output$Lymphocytes <- rowSums(cibersort_output[, group_lymphocytes, drop = FALSE], na.rm = TRUE)
  cibersort_output$Total_Monocytes <- rowSums(cibersort_output[, group_monocytes, drop = FALSE], na.rm = TRUE)
  cibersort_output$Dendritic_Mast <- rowSums(cibersort_output[, group_DendriticMast, drop = FALSE], na.rm = TRUE)
  
  # Remove individual immune cell columns after summing
  cibersort_output <- cibersort_output[, !(colnames(cibersort_output) %in% c(group_lymphocytes, group_monocytes, group_DendriticMast))]
  
  # Store results in the list
  cibersort_results[[tp]] <- cibersort_output
  
  print(paste("Completed CIBERSORT for:", tp))  # Print progress
}


cibersort_results$Ex2AV

combined_results <- do.call(rbind, cibersort_results)



# Convert the results list into a single data frame
combined_results <- do.call(rbind, lapply(names(cibersort_results), function(tp) {
  data <- cibersort_results[[tp]]
  data$Timepoint <- tp  # Add timepoint column
  return(data)
}))

# Keep only the relevant columns
plot_data <- combined_results[, c("Timepoint", "Lymphocytes", "Total_Monocytes", "Dendritic_Mast", "Eosinophils", "Neutrophils")]

# Convert to long format for ggplot
plot_data_long <- melt(plot_data, id.vars = "Timepoint", variable.name = "CellType", value.name = "Proportion")


#Box plot for immune cell groups across Timepoints AandV combined

ggplot(plot_data_long, aes(x = Timepoint, y = Proportion, fill = CellType)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Immune Cell Distributions Across Timepoints",
       x = "Timepoints",
       y = "Proportion",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





timepoint_dfs <- list()

# Iterate over each timepoint in cibersort_results
for (tp in names(cibersort_results)) {
  
  # Get the CIBERSORT results for the current timepoint
  data <- cibersort_results[[tp]]
  
  # Add RNA numbers as a column (they are currently row names)
  data$RNA_number <- rownames(data)
  
  # Merge with metadata to get Condition (A/V)
  data <- merge(data, meta_noOutliers, by = "RNA_number", all.x = TRUE)
  
  # Store the updated dataframe
  timepoint_dfs[[tp]] <- data
}

# Combine all timepoints into one dataframe
combined_results <- do.call(rbind, timepoint_dfs)

plot_data <- combined_results[, c("Time", "Condition", "Lymphocytes", "Total_Monocytes", "Dendritic_Mast", "Eosinophils")]

# Convert to long format
plot_data_long <- melt(plot_data, id.vars = c("Time", "Condition"), 
                       variable.name = "CellType", value.name = "Proportion")

ggplot(plot_data_long, aes(x = Time, y = Proportion, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.7)) + 
  #geom_jitter(aes(color = Condition), size = 2, alpha = 0.6, position = position_dodge(width = 0.7)) +  # Individual points
  theme_bw(base_size = 15) +
  scale_y_continuous(limits = c(0,1))+
  labs(title = "Immune Cell Distributions Across Timepoints (A vs V)",
       x = "Timepoints",
       y = "Proportion",
       fill = "Condition (A/V)") +
  facet_wrap(~CellType, scales = "free") +  # Creates separate plots for each cell type
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




  
plot_data_long %>%
    group_by(Time, Condition, CellType) %>%  # Ensures unique rows
    summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%  # Take mean to avoid list columns
    pivot_wider(names_from = CellType, values_from = Proportion) %>% 
  ggplot(aes(x = Eosinophils, y = Lymphocytes))+
  geom_point()+
  geom_smooth(method = "lm")
  

plot_data_wide <- plot_data_long %>%
  group_by(Time, Condition, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = CellType, values_from = Proportion)

# View structure to ensure proper format
str(plot_data_wide)
head(plot_data_wide)


celltype_data <- plot_data_wide %>%
  dplyr::select(-Time, -Condition)  # Remove non-numeric columns

# Create an empty matrix to store R² values
r_squared_matrix <- matrix(NA, ncol = ncol(celltype_data), nrow = ncol(celltype_data),
                           dimnames = list(colnames(celltype_data), colnames(celltype_data)))

# Compute pairwise R² values using linear regression
for (i in 1:ncol(celltype_data)) {
  for (j in 1:ncol(celltype_data)) {
    if (i != j) {
      # Perform linear regression
      model <- lm(celltype_data[[i]] ~ celltype_data[[j]])
      # Extract R-squared value
      r_squared_matrix[i, j] <- summary(model)$r.squared
    }
  }
}

# Convert matrix to data frame for easier visualization
r_squared_df <- as.data.frame(as.table(r_squared_matrix))


# View R² values
print(r_squared_df)


##Append Cell Type proportions to Metadata for each Time point 

















#ANALYSE THE VOLCANO PLOTS FURTHER: IS THE MUSCLE UPTAKING ANY IMMUNE CELLS?
#MONOCYTES

monocyte_genes <- c("ITGAM", "LYZ", "S100A9", "S100A8", 
                    "ITGAL", "VCAN", "FCN1","MNDA", "CD163", "CX3CR1", "HLA-DRA")

#c("ITGAM", "LYZ", "S100A9", "S100A8", "CSF3R", "LYN", "FCGR3A", 
  #"CSF1R", "ITGAL", "VCAN", "FCN1", "SELL","MNDA", "CD163", "CX3CR1", "SPN", "SIGLEC10", "HLA-DRA")



baseline_values_monocytes <- numeric(length(monocyte_genes))
passive_values_monocytes <- numeric(length(monocyte_genes))
ex1_values_monocytes <- numeric(length(monocyte_genes))
ex2_values_monocytes <- numeric(length(monocyte_genes))
ex3_values_monocytes <- numeric(length(monocyte_genes))
rec1_values_monocytes <- numeric(length(monocyte_genes))
rec2_values_monocytes <- numeric(length(monocyte_genes))



for (i in seq_along(monocyte_genes)) {
  gene <- monocyte_genes[i]
  
  # Extract log2FoldChange values, handling missing data
  if (!is.na(BaselineAV_gene_symbol[gene, "log2FoldChange"])) {
    baseline_values_monocytes[i] <- BaselineAV_gene_symbol[gene, "log2FoldChange"]
  } else {
    baseline_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(PassiveAV_gene_symbol[gene, "log2FoldChange"])) {
    passive_values_monocytes[i] <- PassiveAV_gene_symbol[gene, "log2FoldChange"]
  } else {
    passive_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(EX1AV_gene_symbol[gene, "log2FoldChange"])) {
    ex1_values_monocytes[i] <- EX1AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex1_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Ex2AV_gene_symbol[gene, "log2FoldChange"])) {
    ex2_values_monocytes[i] <- Ex2AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex2_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Ex3AV_gene_symbol[gene, "log2FoldChange"])) {
    ex3_values_monocytes[i] <- Ex3AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex3_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Rec1AV_gene_symbol[gene, "log2FoldChange"])) {
    rec1_values_monocytes[i] <- Rec1AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    rec1_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Rec2AV_gene_symbol[gene, "log2FoldChange"])) {
    rec2_values_monocytes[i] <- Rec2AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    rec2_values_monocytes[i] <- NA  # Assign NA if value is missing
  }
  
}

log2FoldChange_monocytes <- data.frame(
  Gene = rep(monocyte_genes, 7),  # Repeat genes for both conditions
  Time = c(rep("Baseline", length(monocyte_genes)), rep("Passive", length(monocyte_genes)), rep("Ex1", length(monocyte_genes)), 
           rep("Ex2", length(monocyte_genes)), rep("Ex3", length(monocyte_genes)), rep("Rec1", length(monocyte_genes)), rep("Rec2", length(monocyte_genes))),
  log2FoldChange = c(baseline_values_monocytes, passive_values_monocytes, ex1_values_monocytes, ex2_values_monocytes, ex3_values_monocytes, rec1_values_monocytes, rec2_values_monocytes)
)

log2FoldChange_monocytes$Time <- factor(
  log2FoldChange_monocytes$Time,
  levels = c("Baseline", "Passive",  "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
)

ggplot(log2FoldChange_monocytes, aes(x = Gene, y = log2FoldChange, fill = Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
    labs(title = "Comparison of monocyte gene log2FoldChange Across Time",
       x = "Gene",
       y = "log2FoldChange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("BaselineAV" = "blue", "PassiveAV" = "red", "Ex1AV" = "green", "Ex2AV" = "yellow","Ex3AV" = "orange", "Rec1AV" = "black" ,  "Rec2AV" = "purple"))


ggplot(log2FoldChange_monocytes, aes(x = Time, y = log2FoldChange)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + # Box plot with transparency
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_jitter(aes(color = Gene), width = 0.2, size = 2, alpha = 0.8) +  # Jitter for individual genes
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(
       x = "Timepoint",
       y = "log2FoldChange",
       fill = "Timepoint") +
  scale_color_manual(values = rainbow(length(unique(log2FoldChange_monocytes$Gene)))) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1)) 



#LYMPHOCYTES

lymphocyte_genes <- c("CD4", "IGHG1", "GNLY", "GZMB", 
                    "PRF1", "IL4", "IL2","IGK", "IL10")


baseline_values_lymphocyte <- numeric(length(lymphocyte_genes))
passive_values_lymphocyte <- numeric(length(lymphocyte_genes))
ex1_values_lymphocyte <- numeric(length(lymphocyte_genes))
ex2_values_lymphocyte <- numeric(length(lymphocyte_genes))
ex3_values_lymphocyte <- numeric(length(lymphocyte_genes))
rec1_values_lymphocyte <- numeric(length(lymphocyte_genes))
rec2_values_lymphocyte <- numeric(length(lymphocyte_genes))



for (i in seq_along(lymphocyte_genes)) {
  gene <- lymphocyte_genes[i]
  
  # Extract log2FoldChange values, handling missing data
  if (!is.na(BaselineAV_gene_symbol[gene, "log2FoldChange"])) {
    baseline_values_lymphocyte[i] <- BaselineAV_gene_symbol[gene, "log2FoldChange"]
  } else {
    baseline_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(PassiveAV_gene_symbol[gene, "log2FoldChange"])) {
    passive_values_lymphocyte[i] <- PassiveAV_gene_symbol[gene, "log2FoldChange"]
  } else {
    passive_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(EX1AV_gene_symbol[gene, "log2FoldChange"])) {
    ex1_values_lymphocyte[i] <- EX1AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex1_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Ex2AV_gene_symbol[gene, "log2FoldChange"])) {
    ex2_values_lymphocyte[i] <- Ex2AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex2_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Ex3AV_gene_symbol[gene, "log2FoldChange"])) {
    ex3_values_lymphocyte[i] <- Ex3AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    ex3_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Rec1AV_gene_symbol[gene, "log2FoldChange"])) {
    rec1_values_lymphocyte[i] <- Rec1AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    rec1_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
  if (!is.na(Rec2AV_gene_symbol[gene, "log2FoldChange"])) {
    rec2_values_lymphocyte[i] <- Rec2AV_gene_symbol[gene, "log2FoldChange"]
  } else {
    rec2_values_lymphocyte[i] <- NA  # Assign NA if value is missing
  }
  
}

log2FoldChange_lymphocytes <- data.frame(
  Gene = rep(lymphocyte_genes, 7),  # Repeat genes for both conditions
  Time = c(rep("Baseline", length(lymphocyte_genes)), rep("Passive", length(lymphocyte_genes)), rep("Ex1", length(lymphocyte_genes)), 
           rep("Ex2", length(lymphocyte_genes)), rep("Ex3", length(lymphocyte_genes)), rep("Rec1", length(lymphocyte_genes)), rep("Rec2", length(lymphocyte_genes))),
  log2FoldChange = c(baseline_values_lymphocyte, passive_values_lymphocyte, ex1_values_lymphocyte, ex2_values_lymphocyte, ex3_values_lymphocyte, rec1_values_lymphocyte, rec2_values_lymphocyte)
)

log2FoldChange_lymphocytes$Time <- factor(
  log2FoldChange_lymphocytes$Time,
  levels = c("Baseline", "Passive",  "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
)


ggplot(log2FoldChange_lymphocytes, aes(x = Time, y = log2FoldChange)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + # Box plot with transparency
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_jitter(aes(color = Gene), width = 0.2, size = 2, alpha = 0.8) +  # Jitter for individual genes
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(
    x = "Timepoint",
    y = "log2FoldChange",
    fill = "Timepoint") +
  scale_color_manual(values = rainbow(length(unique(log2FoldChange_monocytes$Gene)))) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1)) 



PassiveAV_gene_symbol["S100A8", ]
Ex2AV_gene_symbol["IL1B", ]
Rec1AV_gene_symbol["CXCL1", ]















#Monocyte analysis yielded quite convicing results. Try cell deconvolution 

library(immunedeconv)

deconv_results <- deconvolute(norm_counts, method = "cibersort")









#BETWEEN TIMEPOINTS COMPARISON 

counts_long <- as.data.frame(filtered_counts_noOutliers80pct) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "RNA_number",
    values_to = "count"
  )

counts_long_merge <- counts_long %>%
  left_join(meta_noOutliers, by = "RNA_number")

counts_avdiff <- counts_long_merge %>%
  group_by(Gene, Subject, Time) %>%
  summarize(
    # venous count minus arterial count:
    av_diff = count[Condition == "V"] - count[Condition == "A"],
    flow = Flow[Condition == "A"],  # or FlowRate[Site=="Arterial"] if it's the same
    .groups = "drop"
  )

counts_avdiff$flow <- as.numeric(counts_avdiff$flow)
counts_avdiff$flux <- as.numeric(counts_avdiff$flux)


counts_avdiff <- counts_avdiff %>%
  mutate(
    flux = av_diff * flow
  )


flux_wide <- counts_avdiff %>%
  tidyr::unite("SampleLabel", Subject, Time, remove = FALSE) %>%
  select(Gene, SampleLabel, flux) %>%
  pivot_wider(names_from = "SampleLabel", values_from = "flux")


flux_matrix <- as.matrix(flux_wide[,-1])  # Drop Gene column
rownames(flux_matrix) <- flux_wide$Gene

for (j in seq_along(named_norm_factor)) {
  flux_matrix[, j] <- flux_matrix[, j] / spike_in_norm[j]
}


######LIMMA: REPEATED MEASURES ACROSS TIME 

metadata_flux <- counts_avdiff %>%
  distinct(Subject, Time) %>%
  tidyr::unite("SampleLabel", Subject, Time, remove=FALSE)

design <- model.matrix(~ Time, data = metadata_flux)

colnames(design)
colnames(design) <- make.names(colnames(design))

corfit <- duplicateCorrelation(flux_matrix, design, 
                               block = metadata_flux$Subject)
cat("Estimated within-subject correlation: ", corfit$consensus, "\n")

fit <- lmFit(flux_matrix, design, 
             block = metadata_flux$Subject,
             correlation = corfit$consensus)
fit <- eBayes(fit)

contr.matrix <- makeContrasts(
  Ex3_vs_Baseline = TimeEx3 - X.Intercept.,
  levels = design
)

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

# Extract top results for each contrast
res_ex3 <- topTable(fit2, coef = "Ex3_vs_Baseline", number = Inf)



###Take it one step back: pick a few genes and see how the AV difference behaves, and whether the Limma LFC output makes numerical sense 

#Take the normlization factor dataframe
norm_factor_df 
norm_factor

#Apply normalization factor directly on fitlered counts 
rownames(norm_factor) <- gsub("Sequin sums", "Factor", rownames(norm_factor))

named_norm_factor <- setNames(as.numeric(as.matrix(norm_factor)), rep(colnames(norm_factor), each = nrow(norm_factor)))


filtered_counts_noOutliers80pct_NORM <- sweep(filtered_counts_noOutliers80pct, 2, named_norm_factor, "*")

#Apply log transform the norm_counts matrix 
log_norm_counts <- log1p(filtered_counts_noOutliers80pct_NORM)


#Filter for top few genes

select_genes <- log_norm_counts[rownames(log_norm_counts) %in% c("ENSG00000276168", "ENSG00000101162", "ENSG00000100345"), ]

Baseline_Ex3RNAs <- meta_noOutliers$RNA_number[meta_noOutliers$Time %in% c("Baseline", "Ex3")]

filtered_log_counts <- select_genes[, colnames(select_genes) %in% Baseline_Ex3RNAs, drop = FALSE]

Baseline_Ex3_meta <- meta_noOutliers[meta_noOutliers$RNA_number %in% Baseline_Ex3RNAs, ]

#One subject is missing a Baseline pair (due to being an Outlier) --> Remove it completely 

filtered_log_counts <- filtered_log_counts %>% select(-"RNA031399")

Baseline_Ex3_meta <- Baseline_Ex3_meta %>% filter(RNA_number != "RNA031399")


Baseline_Ex3_meta_wide <- Baseline_Ex3_meta %>%
  select(RNA_number, Subject, Time, Condition) %>%
  tidyr::pivot_wider(
    id_cols  = c("Subject", "Time"),
    names_from = Condition,  # this will create "A" and "V" columns
    values_from = RNA_number
  )


subject_time_labels <- paste(Baseline_Ex3_meta_wide$Subject, Baseline_Ex3_meta_wide$Time, sep="_")
difference_mat <- matrix(
  nrow = nrow(filtered_log_counts), 
  ncol = nrow(Baseline_Ex3_meta_wide),
  dimnames = list(rownames(filtered_log_counts), subject_time_labels)
)

for(i in seq_len(nrow(Baseline_Ex3_meta_wide))) {
  a_col <- Baseline_Ex3_meta_wide$A[i]
  v_col <- Baseline_Ex3_meta_wide$V[i]
  
  difference_mat[, i] <- filtered_log_counts[, v_col] - filtered_log_counts[, a_col]
}

diff_long <- difference_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(
    cols      = -gene_id,
    names_to  = "Subject_Time",
    values_to = "difference"
  )


diff_long <- diff_long %>%
  separate(Subject_Time, into = c("Subject", "Time"), sep = "_")


ggplot(diff_long, aes(x = Time, y = difference)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightblue", alpha = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(aes(color = Subject), width = 0.1) +
  facet_wrap(~ gene_id, scales = "free_y") +
  theme_minimal()



##PCA WITH SPIKE NORMALIZED COUNTS


norm_factor_all <- unname(named_norm_factor)



dds <- DESeqDataSetFromMatrix(countData = filtered_counts_noOutliers,
                              colData = meta_noOutliers,
                              design = ~ Subject + Condition + Time + Condition:Time)

sizeFactors(dds) <- named_norm_factor



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
  #geom_label(label = colnames(vsd_counts_flt),nudge_x = -5, nudge_y = 1) +
  #geom_label_repel(aes(label = Label), max.overlaps = Inf, size = 3)


pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var = "RNA_number") %>% 
  merge(.,meta, by = "RNA_number") %>% 
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(size = 3, alpha = 1, aes(colour = Condition)) + theme_classic(base_size = 15) +
  geom_label_repel(label = colnames(vsd_counts_filtered), max.overlaps = Inf, size = 3)


















































################################################### LIMMMMMMMMAAAAAA ###########################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", force = TRUE)

#Append flow rate to Meta to use as a covariate in model 
flow <- as.vector(t(flow))
meta_noOutliers$flow <- flow
meta_noOutliers$flow <- as.numeric(as.character(meta_noOutliers$flow))


meta_noOutliers$Subject <- factor(meta_noOutliers$Subject)
meta_noOutliers$Condition <- factor(meta_noOutliers$Condition, levels = c("A", "V"))
meta_noOutliers$Time <- factor(meta_noOutliers$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))
meta_noOutliers$logFlowRate <- log(meta_noOutliers$flow)

rownames(meta_noOutliers) <- meta_noOutliers$RNA_number

#Spike-in normalization factors 
spike_in_norm <- norm_factor

#Offset approach 
offset_mat <- matrix(
 rep(log(spike_in_norm), each = nrow(filtered_counts_noOutliers80pct)),
 nrow = nrow(filtered_counts_noOutliers80pct),
ncol = ncol(filtered_counts_noOutliers80pct)
)


dge <- DGEList(counts = filtered_counts_noOutliers80pct)
#dge <- calcNormFactors(dge) Do we want this standard TMM normalization? 
#dge <- calcNormFactors(dge)

#Library size adjustment approach
#dge$samples$lib.size <- dge$samples$lib.size * spike_in_norm
#dge <- calcNormFactors(dge, method = "none")


design <- model.matrix(~ logFlowRate + Time * Condition, data = meta_noOutliers)
colnames(design)
colnames(design) <- make.names(colnames(design))

v <- voom(dge, design, plot = TRUE)
v$E <- v$E + offset_mat




#Within subject correlation via duplicateCorrelation
corfit <- duplicateCorrelation(v, design, block = meta_noOutliers$Subject)
cat("Estimated within-subject correlation:", corfit$consensus, "\n")

fit <- lmFit(v, design, block = meta_noOutliers$Subject, correlation = corfit$consensus)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(
 
  #Baseline compared to Exercise
  Baseline_vs_Ex2 = ((TimeEx2.ConditionV - TimeEx2)) - 
                       ((X.Intercept.+ ConditionV) - X.Intercept.),
  
  levels = design
)


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


Baseline_vs_Ex2_LimmaResults <- topTable(fit2, coef = "Baseline_vs_Ex2", adjust.method = "BH", number = Inf)
head(Baseline_vs_Ex2_LimmaResults)

hist(Baseline_vs_Ex2_LimmaResults$P.Value, breaks=50, main="Raw P-values for ConditionB", xlab="P-value")

hist(Baseline_vs_Ex2_LimmaResults$adj.P.Val, breaks=50, main="Adjusted P-values for ConditionB", xlab="Adjusted P-value")


ggplot(Baseline_vs_Ex2_LimmaResults, aes(x = logFC, y = -log(P.Value))) +
  geom_point(alpha = 0.6) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(P-value)") +
  theme_minimal() +
  theme(legend.title = element_blank())

















contrast.matrix <- makeContrasts(
  Baseline = ((X.Intercept.+ ConditionV) - X.Intercept.),
  Ex2 = ((X.Intercept.+ TimeEx2.ConditionV) - (X.Intercept.+ TimeEx2)),
  Rec2 = ((X.Intercept.+ TimeRec2.ConditionV) - (X.Intercept.+ TimeRec2)),
  levels = design
)





colnames(design)
colnames(design) <- make.names(colnames(design))


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


