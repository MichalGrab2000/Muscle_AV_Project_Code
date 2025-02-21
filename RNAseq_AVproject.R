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

###TWO OUTLIERS IDENTIFIED
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

#We would expect the Spike counts to go down at the Exercising time points as they should be more diluted due to increased flow rate?

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


################################## GENE COUNT MATRIX FILTERING #################################################

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



meta_noOutliers <- meta_noOutliers %>%
  column_to_rownames(var = "RNA_number")



###########################################################################################################
########################### DESeq2 for within timepoint AV comparison ##################################

###########BASELINE################################################

Baseline_RNAs <- meta_noOutliers %>% 
  filter(Time == "Baseline") %>% 
  pull(RNA_number)


Baseline_counts <- filtered_counts_noOutliers[, c(Baseline_RNAs)]
Baseline80_keep <- rowSums(Baseline_counts >= 10) >= (0.8 * ncol(Baseline_counts)) 
Baseline_counts80 <- Baseline_counts[Baseline80_keep,]


meta_Baseline <- meta_noOutliers[c(Baseline_RNAs),]
rownames(meta_Baseline) <- meta_Baseline$RNA_number

norm_factors_Baseline <-norm_factor_df["Norm_factor", Baseline_RNAs] 
norm_factors_Baseline
norm_factors_Baseline <- as.vector(as.numeric(norm_factors_Baseline[1, ]))


# Factors 
meta_Baseline$Subject <- factor(meta_Baseline$Subject)
meta_Baseline$Condition <- factor(meta_Baseline$Condition, levels = c("A", "V"))

dds_Baseline <- DESeqDataSetFromMatrix(countData = Baseline_counts80,
                                       colData = meta_Baseline,
                                       design = ~ Subject + Neutrophils + Lymphocytes + Condition)


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
rownames(BaselineAV_gene_symbol)[rownames(BaselineAV_gene_symbol) == "ENSG00000255478"] <- "ENSG...255478 (lnc)"
rownames(BaselineAV_gene_symbol)[rownames(BaselineAV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
rownames(BaselineAV_gene_symbol)[rownames(BaselineAV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
       width = 8, height = 6, dpi = 800)



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


meta_Passive <- meta_noOutliers[c(Passive_RNAs),]
rownames(meta_Passive) <- meta_Passive$RNA_number

norm_factors_Passive <-norm_factor_df["Norm_factor", Passive_RNAs] 
norm_factors_Passive
norm_factors_Passive <- as.vector(as.numeric(norm_factors_Passive[1, ]))


# Factors 
meta_Passive$Subject <- factor(meta_Passive$Subject)
meta_Passive$Condition <- factor(meta_Passive$Condition, levels = c("A", "V"))

dds_Passive <- DESeqDataSetFromMatrix(countData = Passive_counts80,
                                       colData = meta_Passive,
                                       design = ~ Subject + Eosinophils + Neutrophils + Condition)


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
  sample = names(dds_norm_Passive),           
  counts = as.numeric(dds_norm_Passive)      
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
rownames(PassiveAV_gene_symbol)[rownames(PassiveAV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"
rownames(PassiveAV_gene_symbol)[rownames(PassiveAV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
PassiveAV_results_top_genes <- PassiveAV_gene_symbol[order(PassiveAV_gene_symbol$xiao_score_V), ][1:7, ]
PassiveAV_results_bottom_genes <- PassiveAV_gene_symbol[order(PassiveAV_gene_symbol$xiao_score_A), ][1:3, ]

# Add labels only for the top 10 genes
PassiveAV_gene_symbol$label <- ifelse(
  rownames(PassiveAV_gene_symbol) %in% rownames(PassiveAV_results_top_genes), #| 
   # rownames(PassiveAV_gene_symbol) %in% rownames(PassiveAV_results_bottom_genes), 
  PassiveAV_gene_symbol$Gene_symbol, 
  NA
)




PassiveAV_Volcano <- ggplot(data = PassiveAV_gene_symbol, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=label)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(
    values = c("grey", "green"),
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
       width = 8, height = 6, dpi = 800)



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

meta_EX1 <- meta_noOutliers[c(EX1_RNAs),]
rownames(meta_EX1) <- meta_EX1$RNA_number

norm_factors_EX1 <-norm_factor_df["Norm_factor", EX1_RNAs] 
norm_factors_EX1
norm_factors_EX1 <- as.vector(as.numeric(norm_factors_EX1[1, ]))


# Factors 
meta_EX1$Subject <- factor(meta_EX1$Subject)
meta_EX1$Condition <- factor(meta_EX1$Condition, levels = c("A", "V"))

dds_EX1 <- DESeqDataSetFromMatrix(countData = Ex1_counts80,
                                      colData = meta_EX1,
                                      design = ~ Subject + Total_Monocytes + Eosinophils + Condition)

#Colinearity / Rank troubleshooting 

colData_df <- as.data.frame(colData(dds_EX1))

modelMatrix <- model.matrix(
  design(dds_EX1),  
  data = colData_df
)

qr_mm <- qr(modelMatrix)
rank_mm <- qr_mm$rank
ncol_mm <- ncol(modelMatrix)

cat("Rank:", rank_mm, "\n")
cat("Number of columns:", ncol_mm, "\n")
if (rank_mm < ncol_mm) {
  cat("Warning: The model matrix is not full rank!\n")
}



sizeFactors(dds_EX1) <- norm_factors_EX1
sizeFactors(dds_EX1)


dds_EX1 <- estimateDispersions(dds_EX1, fitType="parametric")
dds_EX1 <- nbinomWaldTest(dds_EX1, maxit = 1000)
sum(mcols(dds_EX1)$betaConv == FALSE)



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
rownames(EX1AV_gene_symbol)[rownames(EX1AV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"
rownames(EX1AV_gene_symbol)[rownames(EX1AV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
  EX1AV_gene_symbol$Gene_symbol, 
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
       width = 8, height = 6, dpi = 800)



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


meta_Ex2 <- meta_noOutliers[c(Ex2_RNAs),]
rownames(meta_Ex2) <- meta_Ex2$RNA_number

norm_factors_Ex2 <-norm_factor_df["Norm_factor", Ex2_RNAs] 
norm_factors_Ex2
norm_factors_Ex2 <- as.vector(as.numeric(norm_factors_Ex2[1, ]))


# Factors 
meta_Ex2$Subject <- factor(meta_Ex2$Subject)
meta_Ex2$Condition <- factor(meta_Ex2$Condition, levels = c("A", "V"))

dds_Ex2 <- DESeqDataSetFromMatrix(countData = Ex2_counts80,
                                       colData = meta_Ex2,
                                       design = ~ Subject + Eosinophils  +  Condition)

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
rownames(Ex2AV_gene_symbol)[rownames(Ex2AV_gene_symbol) == "ENSG00000272369"] <- "ENSG...272369"
rownames(Ex2AV_gene_symbol)[rownames(Ex2AV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
rownames(Ex2AV_gene_symbol)[rownames(Ex2AV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"



#Save to CSV file
write.csv(Ex2AV_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/EX2AV_diff_expression_symbols.csv", row.names = TRUE)



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
Ex2AV_results_top_genes <- Ex2AV_gene_symbol[order(Ex2AV_gene_symbol$xiao_score_V), ][1:4, ]
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
       width = 8, height = 8, dpi = 800)

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

meta_EX3 <- meta_noOutliers[c(EX3_RNAs),]
rownames(meta_EX3) <- meta_EX3$RNA_number

norm_factors_EX3 <-norm_factor_df["Norm_factor", EX3_RNAs] 
norm_factors_EX3
norm_factors_EX3 <- as.vector(as.numeric(norm_factors_EX3[1, ]))


# Factors 
meta_EX3$Subject <- factor(meta_EX3$Subject)
meta_EX3$Condition <- factor(meta_EX3$Condition, levels = c("A", "V"))

dds_EX3 <- DESeqDataSetFromMatrix(countData = Ex3_counts80,
                                  colData = meta_EX3,
                                  design = ~ Subject + Dendritic_Mast + Total_Monocytes + Condition)


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
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG...267469"] <- "ENSG...267469 (lncRNA)"
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
rownames(Ex3AV_gene_symbol)[rownames(Ex3AV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
       width = 8, height = 6, dpi = 800)


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

meta_Rec1 <- meta_noOutliers[c(Rec1_RNAs),]
rownames(meta_Rec1) <- meta_Rec1$RNA_number

norm_factors_Rec1 <-norm_factor_df["Norm_factor", Rec1_RNAs] 
norm_factors_Rec1
norm_factors_Rec1 <- as.vector(as.numeric(norm_factors_Rec1[1, ]))


# Factors 
meta_Rec1$Subject <- factor(meta_Rec1$Subject)
meta_Rec1$Condition <- factor(meta_Rec1$Condition, levels = c("A", "V"))

dds_Rec1 <- DESeqDataSetFromMatrix(countData = Rec1_counts80,
                                  colData = meta_Rec1,
                                  design = ~ Subject + Dendritic_Mast + Lymphocytes + Condition)


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
rownames(Rec1AV_gene_symbol)[rownames(Rec1AV_gene_symbol) == "ENSG00000259790"] <- "ANP32BP1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000211459"] <- "MT-RNR1"
#rownames(Baseline_vs_Ex2_gene_symbol)[rownames(Baseline_vs_Ex2_gene_symbol) == "ENSG00000262413"] <- "lncRNA AC145207.2"
rownames(Rec1AV_gene_symbol)[rownames(Rec1AV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
rownames(Rec1AV_gene_symbol)[rownames(Rec1AV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
Rec1AV_results_top_genes <- Rec1AV_gene_symbol[order(Rec1AV_gene_symbol$xiao_score_V), ][1:7, ]
Rec1AV_results_bottom_genes <- Rec1AV_gene_symbol[order(Rec1AV_gene_symbol$xiao_score_A), ][1:7, ]

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
       width = 8, height = 6, dpi = 800)

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

meta_Rec2 <- meta_noOutliers[c(Rec2_RNAs),]
rownames(meta_Rec2) <- meta_Rec2$RNA_number

norm_factors_Rec2 <-norm_factor_df["Norm_factor", Rec2_RNAs] 
norm_factors_Rec2
norm_factors_Rec2 <- as.vector(as.numeric(norm_factors_Rec2[1, ]))


# Factors 
meta_Rec2$Subject <- factor(meta_Rec2$Subject)
meta_Rec2$Condition <- factor(meta_Rec2$Condition, levels = c("A", "V"))

dds_Rec2 <- DESeqDataSetFromMatrix(countData = Rec2_counts80,
                                   colData = meta_Rec2,
                                   design = ~ Subject + Eosinophils + Lymphocytes + Condition)


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
rownames(Rec2AV_gene_symbol)[rownames(Rec2AV_gene_symbol) == "ENSG00000290918"] <- "ENSG...290918 (lncRNA)"
rownames(Rec2AV_gene_symbol)[rownames(Rec2AV_gene_symbol) == "ENSG00000203395"] <- "ENSG...203395 (lncRNA)"


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
Rec2AV_results_top_genes <- Rec2AV_gene_symbol[order(Rec2AV_gene_symbol$xiao_score_V), ][1:7, ]
Rec2AV_results_bottom_genes <- Rec2AV_gene_symbol[order(Rec2AV_gene_symbol$xiao_score_A), ][1:7, ]

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
       width = 8, height = 6, dpi = 800)

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
  if ("SVIL-AS1" %in% rownames(df)) {
    lfc <- df["SVIL-AS1", "log2FoldChange"]
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

Excluded_genes <- unique(c(Baseline_sig_genes, Passive_sig_genes, Rec_all))

Ex_all <- unique(c(Ex1_sig_genes, Ex2_sig_genes, Ex3_sig_genes))

Only_Ex_genes <- setdiff(Ex_all, Excluded_genes)






Ex_all <- union(Ex1_sig_genes, union(Ex2_sig_genes, Ex3_sig_genes))

Excluded_genes <- union(Rec_all, union(Baseline_sig_genes, Passive_sig_genes))

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

Z <- Ex3AV_gene_symbol[rownames(Ex3AV_gene_symbol) %in% Only_Ex_genes, ]



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


sig_DEG_list <- list(Baseline_sig_genes, Passive_sig_genes, Ex1_sig_genes, Ex2_sig_genes, Ex3_sig_genes, Rec1_sig_genes, Rec2_sig_genes)



sig_venn_list <- list(
  "Baseline" = Baseline_sig_genes,
  "Passive" = Passive_sig_genes,
  "Exercise" = Ex_all,
  "Recovery" = Rec_all
)

library(ggVennDiagram)

ggVennDiagram(sig_venn_list) +
  scale_fill_gradient(low="white", high="blue",  name = "Number of diff.
expressed genes") +
  theme_void()



#Graph summarizing significant gene numbers 

Sig_results <- read_xlsx("/Users/MichalGrabowski/Documents/Master Thesis/DiffExGenes.xlsx")

Sig_long <- pivot_longer(
  Sig_results,
  cols = c("VENOUS HIGHER / UP (ImmuneCell Corrected)", "VENOUS HIGHER / UP", "ARTERIAL HIGHER / DOWN (ImmuneCell Corrected)", "ARTERIAL HIGHER / DOWN"),
  names_to = "Category",
  values_to = "Value"
)

Sig_results  <- pivot_longer(Sig_results, cols = c("VENOUS HIGHER / UP (ImmuneCell Corrected)", "VENOUS HIGHER / UP", "ARTERIAL HIGHER / DOWN (ImmuneCell Corrected)", "ARTERIAL HIGHER / DOWN"), names_to = "Category",
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
  scale_color_manual(values = c("VENOUS HIGHER / UP" = "darkgreen", "VENOUS HIGHER / UP (ImmuneCell Corrected)" = "green" , "ARTERIAL HIGHER / DOWN" = "darkred", "ARTERIAL HIGHER / DOWN (ImmuneCell Corrected)" = "red" )) +
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
  pvalueCutoff = 0.8,
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



GSEA_Ex2 <- gseGO(
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

Ex2_GSEA_df <- GSEA_Ex2 %>%
  as.data.frame() 


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
  geom_point() +  
  scale_color_gradient(low = "purple", high = "red") +  # Adjust color scale
  labs(title = "Top Enriched GO Terms in Exercised Muscle",
       x = "Gene Count",
       y = "GO Term",
       color = "p.adjust",
       size = "GeneRatio") +
  theme_minimal() +  
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))




ridgeplot(GSEA_Ex2, showCategory = 10, fill = "p.adjust") +
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
  pvalueCutoff = 0.4,
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
  pvalueCutoff = 0.6,
  verbose      = FALSE
)


ridgeplot(GSEA_Rec2, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")


##########GSEA COMBINED HEATMAP##############################



# Combine GSEA results into a single data frame
list_of_gsea <- list(
  "Baseline" = GSEA_BASELINE,
  "Passive" = GSEA_PASSIVE,
  "Ex1" = GSEA_EX1,
  "Ex2" = GSEA_Ex2,
  "Ex3" = GSEA_EX3,
  "Rec1" = GSEA_Rec1,
  "Rec2" = GSEA_Rec2
)


is_gsea_result <- sapply(list_of_gsea, function(x) inherits(x, "gseaResult"))

pval_threshold <- 0.5
sig_threshold <- 0.05

extract_gsea_data <- function(gsea_obj, timepoint) {
  df <- as.data.frame(gsea_obj) %>%
    dplyr::filter(p.adjust < sig_threshold) %>%  # Keep only significant pathways
    dplyr::select(Description, NES) %>%  # Use pathway descriptions instead of GO ID
    dplyr::rename(pathway = Description, !!timepoint := NES)  # Rename columns
  return(df)
}


gsea_dfs <- lapply(names(list_of_gsea), function(tp) extract_gsea_data(list_of_gsea[[tp]], tp))

merged_gsea <- Reduce(function(x, y) full_join(x, y, by = "pathway"), gsea_dfs)


merged_gsea <- merged_gsea[rowSums(!is.na(merged_gsea[,-1])) > 1, ]

# Convert to matrix format for heatmap
nes_matrix <- as.matrix(merged_gsea[,-1])  
rownames(nes_matrix) <- merged_gsea$pathway

# Handle missing values 
nes_matrix[is.na(nes_matrix)] <- 0  # Change if needed

color_breaks <- seq(-1.5, 1.5, length.out = 50)  # Adjust scale
color_palette <- colorRampPalette(c("red", "white", "green"))(50)

nes_matrix <- log1p(abs(nes_matrix)) * sign(nes_matrix)

timepoint_order <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
nes_matrix <- nes_matrix[, timepoint_order] 


rownames(nes_matrix) <- gsub("adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", 
                             "Adaptive immune response (recomb. of immune receptors)", 
                             rownames(nes_matrix))


GSEA_heatmap <- pheatmap(nes_matrix, 
         color = color_palette,
         breaks = color_breaks,  # Use custom color scale
         cluster_rows = TRUE,  
         cluster_cols = FALSE,  
         #treeheight_row = 0,  # Remove dendrogram lines
         #treeheight_col = 0,
         scale = "none",  
         display_numbers = FALSE,  
         fontsize_row = 12,  
         fontsize_col = 12, 
         angle_col = 45
  )

ggsave(filename = "GSEA_heatmap.png", plot = GSEA_heatmap, 
       width = 12, height = 10, dpi = 800)




# Define thresholds
pval_inclusion_threshold <- 0.2  # include terms with p.adjust < 0.5
sig_threshold <- 0.05            # highlight cells with p.adjust < 0.05

# Extraction functions that filter for inclusion threshold
extract_nes_data <- function(gsea_obj, timepoint) {
  as.data.frame(gsea_obj) %>%
    dplyr::filter(p.adjust < pval_inclusion_threshold) %>% 
    dplyr::select(Description, NES) %>%
    dplyr::rename(pathway = Description, !!timepoint := NES)
}

extract_padj_data <- function(gsea_obj, timepoint) {
  as.data.frame(gsea_obj) %>%
    dplyr::filter(p.adjust < pval_inclusion_threshold) %>% 
    dplyr::select(Description, p.adjust) %>%
    dplyr::rename(pathway = Description, !!timepoint := p.adjust)
}

# Apply functions for each timepoint in the list_of_gsea
nes_list <- lapply(names(list_of_gsea), function(tp) {
  extract_nes_data(list_of_gsea[[tp]], tp)
})
padj_list <- lapply(names(list_of_gsea), function(tp) {
  extract_padj_data(list_of_gsea[[tp]], tp)
})

# Merge the data frames across timepoints by pathway
merged_nes <- Reduce(function(x, y) full_join(x, y, by = "pathway"), nes_list)
merged_padj <- Reduce(function(x, y) full_join(x, y, by = "pathway"), padj_list)

# Convert merged data frames to matrices (rows: pathways, columns: timepoints)
nes_matrix <- as.matrix(merged_nes[,-1])
padj_matrix <- as.matrix(merged_padj[,-1])
rownames(nes_matrix) <- merged_nes$pathway
rownames(padj_matrix) <- merged_nes$pathway

# Replace missing values as needed (here, we set missing NES to 0 and missing p.adjust to 1)
nes_matrix[is.na(nes_matrix)] <- 0
padj_matrix[is.na(padj_matrix)] <- 1

# Create an annotation matrix: add "*" for cells where p.adjust < 0.05, blank otherwise.
annotation_matrix <- ifelse(padj_matrix < sig_threshold, "*", "")

# (Optional) Transform NES values (for example, a log transformation)
nes_matrix <- log1p(abs(nes_matrix)) * sign(nes_matrix)

# Order columns by your specified timepoints if desired
timepoint_order <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
nes_matrix <- nes_matrix[, timepoint_order, drop = FALSE]
annotation_matrix <- annotation_matrix[, timepoint_order, drop = FALSE]

# Define the heatmap color scale
color_breaks <- seq(-1.5, 1.5, length.out = 50)
color_palette <- colorRampPalette(c("red", "white", "green"))(50)

# Plot the heatmap with cell annotations (asterisks for p.adjust < 0.05)
GSEA_heatmap <- pheatmap(nes_matrix, 
                         color = color_palette,
                         breaks = color_breaks,
                         cluster_rows = FALSE,  
                         cluster_cols = FALSE,  
                         scale = "none",  
                         display_numbers = annotation_matrix,  
                         fontsize_row = 12,  
                         fontsize_col = 12, 
                         angle_col = 45)


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
  cibersort_output$Dendritic_MastCells <- rowSums(cibersort_output[, group_DendriticMast, drop = FALSE], na.rm = TRUE)
  
  # Remove individual immune cell columns after summing
  cibersort_output <- cibersort_output[, !(colnames(cibersort_output) %in% c(group_lymphocytes, group_monocytes, group_DendriticMast))]
  
  # Store results in the list
  cibersort_results[[tp]] <- cibersort_output
  
  print(paste("Completed CIBERSORT for:", tp))  # Print progress
}



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

plot_data <- combined_results[, c("Time", "Condition", "Lymphocytes", "Total_Monocytes", "Dendritic_Mast", "Eosinophils", "Neutrophils")]

# Convert to long format
plot_data_long <- melt(plot_data, id.vars = c("Time", "Condition"), 
                       variable.name = "CellType", value.name = "Proportion")

immune_cellsplot <- ggplot(plot_data_long, aes(x = Time, y = Proportion, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.7)) + 
  #geom_jitter(aes(color = Condition), size = 2, alpha = 0.6, position = position_dodge(width = 0.7)) +  # Individual points
  theme_bw(base_size = 15) +
  scale_y_continuous(limits = c(0,1))+
  labs(
       x = "Timepoints",
       y = "Proportion",
       fill = "Condition (A/V)") +
  facet_wrap(~CellType, scales = "free") +  # Creates separate plots for each cell type
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        panel.grid = element_blank())

ggsave(filename = "immune_cellsplot.png", plot = immune_cellsplot, 
       width = 12, height = 10, dpi = 800)

  
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

immunecell_data <- combined_results %>%
  dplyr::select(RNA_number, Lymphocytes, Total_Monocytes, Dendritic_Mast, Eosinophils, Neutrophils)

meta_noOutliers <- meta_noOutliers %>%
  left_join(immunecell_data, by = "RNA_number")

meta_noOutliers <- meta_noOutliers %>%
  column_to_rownames(var = "RNA_number")

meta_noOutliers$RNA_number <- rownames(meta_noOutliers)



#Statistical test: are the immune cell proprotions actually different by timepoint in A and V

timepoints <- c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2")
immune_vars <- c("Lymphocytes", "Total_Monocytes", "Dendritic_Mast", "Eosinophils", "Neutrophils")

# Initialize an empty list to store test results
test_results <- list()

# Loop through each timepoint
for(time in timepoints){
  # Subset your data for the current timepoint
  data_subset <- subset(combined_results, Time == time)
  
  # Initialize a sub-list for the current timepoint
  test_results[[time]] <- list()
  
  cat("Timepoint:", time, "\n")
  
  # Loop through each immune cell variable
  for(var in immune_vars){
    cat("  Testing variable:", var, "\n")
    
    # Construct the formula dynamically
    formula <- as.formula(paste(var, "~ Condition"))
    
    # Run the Wilcoxon test
    test <- wilcox.test(formula, data = data_subset)
    
    # Save the test result in the list
    test_results[[time]][[var]] <- test
    
    # Optionally print the result to the console
    print(test)
  }
  
  cat("\n")
}


Baseline_results <- test_results[["Baseline"]]
print(Baseline_results)

Passive_results <- test_results[["Passive"]]
print(Passive_results)

Ex2_results <- test_results[["Ex2"]]
print(Ex2_results)

Ex3_results <- test_results[["Ex3"]]
print(Ex3_results)

Rec1_results <- test_results[["Rec1"]]
print(Rec1_results)

Rec2_results <- test_results[["Rec2"]]
print(Rec2_results)



######################################################COMBINE WITH BULK MUSCLE DATA##################################

##RECOVERY / REC1 2

Muscle_rec <- read.csv ("~/Documents/Master Thesis/Diff_expression_Recovery_vs_Pre.csv", header = T)
rownames(Muscle_rec) <- Muscle_rec$X

Muscle_rec_up <- subset(Muscle_rec, log2FoldChange > 0.6 & pvalue < 0.05)


Rec1UP <- subset(Rec1AV_gene_symbol, log2FoldChange > 0.6 & pvalue < 0.05)
Rec2UP <- subset(Rec2AV_gene_symbol, log2FoldChange > 0.6 & pvalue < 0.05)

RecUPall <- union(rownames(Rec1UP), rownames(Rec2UP))

Rec1DOWN <- subset(Rec1AV_gene_symbol, log2FoldChange < -0.6 & pvalue < 0.05)
Rec2DOWN <- subset(Rec2AV_gene_symbol, log2FoldChange < -0.6 & pvalue < 0.05)

RecDOWNall <- union(rownames(Rec1DOWN), rownames(Rec2DOWN))


common_genes_Rec <- intersect(rownames(Muscle_rec_up), RecUPall)

common_genes_RecDown <- intersect(rownames(Muscle_rec_up), RecDOWNall)

length(common_genes_Rec)




common_in_Muscle_rec_up <- Muscle_rec_up[common_genes_Rec, ]
common_in_Rec2UP <- Rec2UP[common_genes_Rec, ]


#POST / EXERCISE 1 2 3 

Muscle_post <- read.csv ("~/Documents/Master Thesis/Diff_expression_Post_vs_Pre.csv", header = T)
rownames(Muscle_post) <- Muscle_post$X

Muscle_post_up <- subset(Muscle_post, log2FoldChange > 0.6 & pvalue < 0.05)

Ex1UP <- subset(EX1AV_gene_symbol, log2FoldChange > 0.6 & pvalue < 0.05)
Ex2UP <- subset(Ex2AV_gene_symbol, log2FoldChange > 0.6 & pvalue < 0.05)
Ex3UP <- subset(Ex3AV_gene_symbol, log2FoldChange > 0.6 & pvalue < 0.05)

ExUPall <- union(rownames(Ex1UP), union(rownames(Ex2UP), rownames(Ex3UP)))


Ex1DOWN <- subset(EX1AV_gene_symbol, log2FoldChange < -0.6 & pvalue < 0.05)
Ex2DOWN <- subset(Ex2AV_gene_symbol, log2FoldChange < -0.6 & pvalue < 0.05)
Ex3DOWN <- subset(Ex3AV_gene_symbol, log2FoldChange < -0.6 & pvalue < 0.05)

ExDOWNall <- union(rownames(Ex1DOWN), union(rownames(Ex2DOWN), rownames(Ex3DOWN)))

common_genes_Exup <- intersect(rownames(Muscle_post_up), ExUPall)

common_genes_Exdown <- intersect(rownames(Muscle_post_up), ExDOWNall)

common_genes_Exdown <- intersect(RecUPall, ExUPall)

library(VennDiagram)


AVmuscle_gene_list <- list(
  Muscle_Recovery   = rownames(Muscle_rec_up),
  V_Recovery        = RecUPall,
  Muscle_Exercise  = rownames(Muscle_post_up),
  V_Exercise        = ExUPall
)

ggVennDiagram(AVmuscle_gene_list) + 
  scale_fill_gradient(low = "white", high = "red") + 
  labs(title = "4-circle Venn Diagram of Gene Overlaps") +
  theme_minimal(base_size = 14)

library(ggvenn)

ggvenn(
  data = AVmuscle_gene_list, 
  columns = names(AVmuscle_gene_list),
  fill_color = c("cornflowerblue", "yellow", "green", "orchid"),
  stroke_size = 0.5,  # outline thickness
  set_name_size = 5
) + 
  ggtitle("4-circle Venn Diagram of Gene Overlaps") +
  theme_minimal(base_size = 14)



keep <- rowSums(filtered_counts_noOutliers >= 10) >= 89
fitlered_flt_counts_noOutliers <- filtered_counts_noOutliers[keep,]






####Between Timepoint Comparison#####################
##########################################################
#############LIMMAAAAAAA####################
library(edgeR)

meta_noOutliers$Subject <- factor(meta_noOutliers$Subject)
meta_noOutliers$Condition <- factor(meta_noOutliers$Condition, levels = c("A", "V"))
meta_noOutliers$Time <- factor(meta_noOutliers$Time, levels = c("Baseline", "Passive", "Ex1", "Ex2", "Ex3", "Rec1", "Rec2"))
meta_noOutliers$logFlowRate <- log(meta_noOutliers$flow)

#Merge all Ex timepoints into a single Exercise timepoint, likewise for recovery
meta_noOutliers$Time_combined <- meta_noOutliers$Time
levels(meta_noOutliers$Time_combined)[levels(meta_noOutliers$Time_combined) %in% c("Ex1","Ex2","Ex3")] <- "Exercise"
levels(meta_noOutliers$Time_combined)[levels(meta_noOutliers$Time_combined) %in% c("Rec1","Rec2")] <- "Recovery"
meta_noOutliers$Time_combined <- relevel(meta_noOutliers$Time_combined, ref = "Baseline")


y <- DGEList(counts = filtered_counts_noOutliers)

lib_sizes <- colSums(y$counts)  
print(lib_sizes)  

design <- model.matrix(~Lymphocytes 
                       + Total_Monocytes 
                       + Dendritic_Mast 
                       + Eosinophils  
                       + Time_combined * Condition, data = meta_noOutliers)

keep <- filterByExpr(y,
                     design = design)
                    

y <- y[keep, , keep.lib.sizes=FALSE]

spike_in_norm <- as.numeric(spike_in_norm)
y$samples$norm.factors <- spike_in_norm


colnames(design)
colnames(design) <- make.names(colnames(design))


v <- voom(y, design = design, plot = FALSE)

min_ave_expr <- 1
keep_genes <- rowMeans(v$E) > min_ave_expr
v <- v[keep_genes, ]



#Within subject correlation via duplicateCorrelation , Subject as a random effect 
corfit <- duplicateCorrelation(v, design, block = meta_noOutliers$Subject)
cat("Estimated within-subject correlation:", corfit$consensus, "\n")

fit <- lmFit(v, design, block = meta_noOutliers$Subject, correlation = corfit$consensus)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(
  
  #Recovery vs Exercise
  Recovery_vs_Exercise = Time_combinedRecovery.ConditionV - Time_combinedExercise.ConditionV,
  
  ExerciseAV = ConditionV + Time_combinedExercise.ConditionV,
  
  RecoveryAV = ConditionV + Time_combinedRecovery.ConditionV,
  
  levels = design
)


fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


BaselineAV <- topTable(fit, coef = "ConditionV", adjust.method = "BH", number = Inf)

ExerciseAV <- topTable(fit2, coef = "ExerciseAV", adjust.method = "BH", number = Inf)

RecoveryAV <- topTable(fit2, coef = "RecoveryAV", adjust.method = "BH", number = Inf)


Exercise_effect <- topTable(fit, coef = "Time_combinedExercise.ConditionV", adjust.method = "BH", number = Inf)
#Add XIAO score column 
Exercise_effect$xiao_score =  10^-(sqrt(log10(1 / (Exercise_effect$P.Value^Exercise_effect$logFC))^2))
                                                        
Recovery_effect <- topTable(fit, coef = "Time_combinedRecovery.ConditionV", adjust.method = "BH", number = Inf)
Recovery_effect$xiao_score =  10^-(sqrt(log10(1 / (Recovery_effect$P.Value^Recovery_effect$logFC))^2))

Recovery_Exercise <- topTable(fit2, coef = "Recovery_vs_Exercise", adjust.method = "BH", number = Inf)
Recovery_Exercise$xiao_score =  10^-(sqrt(log10(1 / (Recovery_Exercise$P.Value^Recovery_Exercise$logFC))^2))


write.csv(as.data.frame(Exercise_effect), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_Exercise_vs_Baseline.csv")

write.csv(as.data.frame(Recovery_effect), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_Recovery_vs_Baseline.csv")

write.csv(as.data.frame(Recovery_Exercise), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_Recovery_vs_Exercise.csv")

write.csv(as.data.frame(BaselineAV), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_BaselineAV.csv")

write.csv(as.data.frame(ExerciseAV), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_ExerciseAV.csv")

write.csv(as.data.frame(RecoveryAV), 
          file="/Users/MichalGrabowski/Documents/Master Thesis/Limma_RecoveryAV.csv")

##############

Exercise_Baseline_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Limma_Recovery_vs_Exercise.csv", row.names = 1)

#Rownames to column for easier data manipulation 
Exercise_Baseline_results_geneid<- rownames_to_column(Exercise_Baseline_results, var = "gene_id")

gene_ids_Exercise_Baseline <- Exercise_Baseline_results_geneid %>% 
  dplyr::select(gene_id) %>% 
  pull(gene_id)

gene_symbols_Exercise_Baseline <- mapIds(org.Hs.eg.db, keys = gene_ids_Exercise_Baseline, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

#Convert the symbols to data frame to merge with counts
symbols_data_frame_Exercise <- data.frame(
  gene_id = names(gene_symbols_Exercise_Baseline),
  gene_symbol = ifelse(is.na(gene_symbols_Exercise_Baseline), names(gene_symbols_Exercise_Baseline), gene_symbols_Exercise_Baseline),
  stringsAsFactors = FALSE
)

#Merge gene symbols with the original data
Exercise_results_merged <- Exercise_Baseline_results_geneid %>%
  left_join(symbols_data_frame_Exercise, by = "gene_id")

# Reorder columns to place gene_symbol first
Exercise_results_gene_symbol <- Exercise_results_merged %>%
  relocate(gene_symbol, .before = everything())

# Handle duplicate rownames by appending gene_id to duplicates
Exercise_results_gene_symbol <- Exercise_results_gene_symbol %>%
  mutate(
    gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(gene_symbol, "_", gene_id), gene_symbol)
  )

# Use the unique gene_symbol column as rownames
ExerciseLimma_gene_symbol <- Exercise_results_gene_symbol %>%
  column_to_rownames(var = "gene_symbol_unique")

#Remove symbol column
ExerciseLimma_gene_symbol <- ExerciseLimma_gene_symbol[, setdiff(names(ExerciseLimma_gene_symbol), c("gene_symbol", "gene_id"))]

#Rename those at the top still with IDs manually 
rownames(ExerciseLimma_gene_symbol)[rownames(ExerciseLimma_gene_symbol) == "ENSG00000187536"] <- "TPM3P7"
rownames(ExerciseLimma_gene_symbol)[rownames(ExerciseLimma_gene_symbol) == "ENSG00000267278"] <- "MAP3K14-AS1"
rownames(ExerciseLimma_gene_symbol)[rownames(ExerciseLimma_gene_symbol) == "ENSG00000290450"] <- "ENSG...290450 (lncRNA)"
rownames(ExerciseLimma_gene_symbol)[rownames(ExerciseLimma_gene_symbol) == "ENSG00000288398"] <- "ENSG...288398 (lncRNA)"

#Save to CSV file
write.csv(ExerciseLimma_gene_symbol, file = "/Users/MichalGrabowski/Documents/Master Thesis/RecoveryExerciseLimma_diff_expression_symbols.csv", row.names = TRUE)


######VOLCANO PLOTS BASELINE


#Ex2AV_results <- read.csv("/Users/MichalGrabowski/Documents/Master Thesis/Diff_expressionEx2_A_vs_V.csv", row.names = 1)
#Already labelled with symbols 
ExerciseLimma_gene_symbol 


#Add explicit gene_symbol column for easier processing
ExerciseLimma_gene_symbol$Gene_symbol <- rownames(ExerciseLimma_gene_symbol)


# add a column of NAs
ExerciseLimma_gene_symbol$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
ExerciseLimma_gene_symbol$diffexpressed[ExerciseLimma_gene_symbol$logFC > 0.6 & ExerciseLimma_gene_symbol$P.Value < 0.05] <- "UPREGULATED"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
ExerciseLimma_gene_symbol$diffexpressed[ExerciseLimma_gene_symbol$logFC < -0.6 & ExerciseLimma_gene_symbol$P.Value < 0.05] <- "DOWNREGULATED"

#Get up gene labels based on LFC
ExerciseLimma_gene_symbol$label <- NA
positive_genes <- ExerciseLimma_gene_symbol[ExerciseLimma_gene_symbol$logFC > 0, ]
ExerciseLimma_top_genes <- positive_genes[order(positive_genes$xiao_score), ][1:5, ]
negative_genes <- ExerciseLimma_gene_symbol[ExerciseLimma_gene_symbol$logFC < 0, ]
ExerciseLimma_bottom_genes <- negative_genes[order(negative_genes$xiao_score), ][1:5, ]

# Add labels only for the top 10 genes
ExerciseLimma_gene_symbol$label <- ifelse(
  rownames(ExerciseLimma_gene_symbol) %in% rownames(ExerciseLimma_top_genes) | 
    rownames(ExerciseLimma_gene_symbol) %in% rownames(ExerciseLimma_bottom_genes), 
  ExerciseLimma_gene_symbol$Gene_symbol, 
  NA
)


ExerciseLimma_Volcano <- ggplot(data = ExerciseLimma_gene_symbol, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=label)) + 
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

ggsave(filename = "Rec_Exercise.png", plot = ExerciseLimma_Volcano, 
       width = 8, height = 6, dpi = 800)


ExerciseLimma_gene_symbol["TCEAL3",]
ExerciseLimma_gene_symbol["CDK7",]
ExerciseLimma_gene_symbol["MBTPS1",]
ExerciseLimma_gene_symbol["FPGS",]
ExerciseLimma_gene_symbol["CC2D1A",]

ExerciseLimma_gene_symbol["UQCRQ",]
ExerciseLimma_gene_symbol["CYTH3",]
ExerciseLimma_gene_symbol["DNAJC17",]
ExerciseLimma_gene_symbol["LAMP5",]
ExerciseLimma_gene_symbol["GPKOW",]









###GSEA ON LIMMA RESULTS
#########################


BiocManager::install("pathview")
library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation database
library(DOSE)           # For GSEA analysis
library(enrichplot)     # Visualization
library(pathview)       # KEGG pathway visualization


LimmaExercise_GSEA <- read.csv ("~/Documents/Master Thesis/RecoveryLimma_diff_expression_symbols.csv", header = T) %>% 
drop_na(logFC) %>%
  dplyr::rename(gene = X)


geneList <- LimmaExercise_GSEA$logFC  
names(geneList) <- LimmaExercise_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)


GSEA_EXERCISE_AV <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,  # Human database
  keyType      = "SYMBOL",
  ont          = "ALL",
  minGSSize    = 5,  # Minimum genes per category
  maxGSSize    = 500, # Maximum genes per category
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  verbose      = FALSE
)



ridgeplot(GSEA_EXERCISE_AV, showCategory = 10, fill = "p.adjust") +
  labs(x = "Gene Ranking",
       y = "GO Term")



###KEGG
################

LimmaExercise_GSEA <- read.csv ("~/Documents/Master Thesis/RecoveryExerciseLimma_diff_expression_symbols.csv", header = T) %>% 
  drop_na(logFC) %>%
  dplyr::rename(gene = X)

geneList <- LimmaExercise_GSEA$logFC  
names(geneList) <- LimmaExercise_GSEA$gene
geneList <- sort(geneList, decreasing = TRUE)

gene_df <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gene_df <- na.omit(gene_df)

# Remove duplicate mappings: keep the first occurrence
gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]

# Subset geneList to keep only mapped symbols
geneList <- geneList[names(geneList) %in% gene_df$SYMBOL]

# Match geneList index to gene_df
matched_indices <- match(names(geneList), gene_df$SYMBOL)

# Rename geneList with Entrez IDs
names(geneList) <- gene_df$ENTREZID[matched_indices]

# Remove any non-finite values (NA, NaN, Inf)
geneList <- geneList[is.finite(geneList)]

# Re-sort in descending order
geneList <- sort(geneList, decreasing = TRUE)

# Quick check
cat("Number of genes in geneList after filtering: ", length(geneList), "\n")

head(names(geneList))



kegg_gsea <- gseKEGG(
  geneList     = geneList,
  organism     = "hsa",  # Human
  minGSSize    = 10,
  maxGSSize    = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.3
)

ridgeplot(kegg_gsea, showCategory = 15, fill = "p.adjust") + 
  labs(x = "Gene Ranking", y = "KEGG Pathway") + 
  ggtitle("KEGG Pathway GSEA")



significant_genes <- LimmaExercise_GSEA %>%
  # 1) Filter by adjusted p-value
  filter(P.Value < 0.05) %>%
  # 2) Filter by log fold change threshold
  filter(abs(logFC) > 0.6) %>%
  # 3) Extract only the gene column
  pull(gene)

gene_df <- bitr(significant_genes, 
                fromType = "SYMBOL", 
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# Remove any rows that failed to map
gene_df <- na.omit(gene_df)
# Remove duplicates
gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]

entrez_genes <- gene_df$ENTREZID
length(entrez_genes)
head(entrez_genes)


kegg_ora <- enrichKEGG(
  gene         = entrez_genes,  # the Entrez IDs
  organism     = "hsa",         # human
  keyType      = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)


ridgeplot(kegg_gsea, showCategory = 15, fill = "p.adjust") + 
  labs(x = "Gene Ranking", y = "KEGG Pathway") + 
  ggtitle("KEGG Pathway GSEA")


dotplot(kegg_ora, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment (ORA)") +
  theme_minimal()


####ORA

LimmaExercise_ORA <- read.csv(
  "~/Documents/Master Thesis/RecoveryExerciseLimma_diff_expression_symbols.csv",
  header = TRUE
) %>%
  drop_na(logFC) %>%
  dplyr::rename(gene = X)

up_genes <- LimmaExercise_ORA %>%
  filter(logFC > 0 & P.Value <0.05 ) %>%
  pull(gene)

background_genes <- unique(LimmaExercise_ORA$gene)

ORA_EXERCISE_AV <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",      # can be "BP", "MF", "CC", or "ALL"
  universe      = background_genes,  
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500
)


dotplot(ORA_EXERCISE_AV, showCategory = 10)










































































































