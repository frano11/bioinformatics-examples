# Load required libraries
library(tidyverse)
library(readr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("count", "dplyr")

# Define paths to FeatureCounts result folders
paths <- "/absolute/path/to/results/FeatureCounts"

# Find all *_counts.txt files recursively
count_files <- list.files(path = paths, pattern = "*_counts.txt$", full.names = TRUE)

# Initialize list for storing each sample's counts
count_list <- list()

# Loop through each count file
for (file in count_files) {
  # Extract clean sample name
  sample_name <- sub(".fq_counts.txt$", "", basename(file))
  
  # Read featureCounts file (skip first comment line)
  df <- read_tsv(file, comment = "#", col_types = cols())
  
  # Extract only Geneid and the last column (counts)
  df_sub <- df %>%
    dplyr::select(Geneid, count = last_col()) %>%  # Dynamically selects the last column
    rename(!!sample_name := count)  # Rename count column to sample name
  
  # Add to list
  count_list[[sample_name]] <- df_sub
}

# Merge all data frames by Geneid
count_matrix <- reduce(count_list, left_join, by = "Geneid")

# Optional: Set rownames and remove Geneid column if preferred by DESeq2
# count_matrix_rownames <- column_to_rownames(count_matrix, var = "Geneid")

# Renaming "Geneid"
count_matrix <- count_matrix %>% 
  rename(Gene_Id = Geneid)

# Save output
write_csv(count_matrix, "/absolute/path/to/R_count_matrix/count_matrix.csv")


#################################################
## ADDING THE COLUMN "Gene_Name"
#################################################

# Load biomaRt for annotation
library(biomaRt)

# Use Ensembl's mouse dataset
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract just the Ensembl IDs (remove version if needed)
ensembl_ids <- gsub("\\..*", "", count_matrix$Gene_Id)  # remove version like .5

# Map Ensembl IDs to gene symbols
annotations <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Rename for merging
colnames(annotations) <- c("Ensembl_ID", "Gene_name")

# Add stripped ID column for merging
count_matrix$Ensembl_ID <- ensembl_ids

# Merge by Ensembl ID
count_matrix_annotated <- left_join(count_matrix, annotations, by = "Ensembl_ID")

# Reorder columns: Gene_name, Geneid, then samples
count_matrix_annotated <- count_matrix_annotated %>%
  dplyr::select(Gene_name, Gene_Id, everything(), -Ensembl_ID)

# (Optional) Rename Geneid to Gene_Id
# count_matrix_annotated <- count_matrix_annotated %>%
#   rename(Gene_Id = Geneid)

# Save final annotated matrix
write_csv(count_matrix_annotated, "/absolute/path/to/R_count_matrix/count_matrix_annotated.csv")

##############################################
### VALIDATION OF GENE NAMES WITH ENSEMBL NAMES
##############################################

getBM(
  attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "mgi_symbol", "gene_biotype"),
  filters = "mgi_symbol",
  values = "Gapdh",
  mart = mart
)

getBM(
  attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "mgi_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = "ENSMUSG00000057666",
  mart = mart
)

#####################################################
### FILTERING GENES: NA, (starting with) GM, and Rik
#####################################################

count_matrix_filtered <- count_matrix_annotated %>%
  filter(
    !is.na(Gene_name),
    !grepl("^Gm", Gene_name),
    !grepl("Rik$", Gene_name)
  )


write_csv(
  count_matrix_filtered,
  "/absolute/path/to/R_count_matrix/count_matrix_annotated_clean.csv"
)

# preview what you're removing
removed_genes <- count_matrix_annotated %>%
  filter(
    is.na(Gene_name) | grepl("^Gm", Gene_name) | grepl("Rik$", Gene_name)
  )

View(removed_genes)


#####################################################
### DESeq2 ANALYSIS
#####################################################
# Install packages
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "pheatmap"))

# Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
# library(EnhancedVolcano)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

######## PART 1: PCA First

# Load count matrix
count_matrix <- read_csv("/absolute/path/to/R_count_matrix/count_matrix_annotated_clean.csv")

# Convert to data frame and set Gene_Id as rownames
count_df <- count_matrix %>%
  column_to_rownames(var = "Gene_Id") %>%
  dplyr::select(-Gene_name)

# Ensure it's integer matrix
count_df <- count_df %>%
  mutate(across(everything(), ~ as.integer(round(.))))  # Safely round if needed

# After coercion, check for NAs:
if (any(is.na(count_df))) {
  warning("Some counts became NA after coercion! Please check the input data carefully.")
}

summary(as.vector(as.matrix(count_df)))

####################################################
## INTERROGATING / MANIPULATING THE MATRIX
## FINDING THE GENE AND SAMPLE SHOWING MAX COUNTING (OPTIONAL)
####################################################
# Option 1

# Find where the maximum count is
max_value <- max(count_df, na.rm = TRUE)

# Find which gene and sample had that value
which_max <- which(as.matrix(count_df) == max_value, arr.ind = TRUE)

# Extract gene and sample
gene_with_max <- rownames(count_df)[which_max[,"row"]]
sample_with_max <- colnames(count_df)[which_max[,"col"]]

# Display
cat("Gene with max count:", gene_with_max, "\n")
cat("Sample with max count:", sample_with_max, "\n")
cat("Max count value:", max_value, "\n")

# Load your annotated matrix
count_matrix_annotated <- read_csv("/absolute/path/to/R_count_matrix/count_matrix_annotated_clean.csv")

# Search for your Ensembl ID (removing version if necessary)
gene_of_interest <- count_matrix_annotated %>%
  filter(Gene_Id == "ENSMUSG00000029368.10") %>%
  pull(Gene_name)

gene_of_interest

# Option 2

# Load biomaRt
library(biomaRt)

# Connect to Ensembl (mouse dataset)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Remove version (".10") part from Ensembl ID
ensembl_id_clean <- gsub("\\..*", "", "ENSMUSG00000029368.10")

# Query
gene_symbol <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_id_clean,
  mart = mart
)

gene_symbol # Alb: Albumin, massively expressed in the Tissue1

####################################################
####################################################

# Create metadata
sample_info <- data.frame(
  sample = colnames(count_df),
  tissue = c("Tissue1", "Tissue2", "Tissue1", "Tissue2", "Tissue1", "Tissue2", "Tissue1", "Tissue2"),  # Adjust to match your samples
  row.names = colnames(count_df)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_df, colData = sample_info, design = ~ tissue)

# Filter low-expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Variance Stabilizing Transformation
vsd <- vst(dds, blind = TRUE) # "blind = TRUE" is when you want an unbiased look at your sample relationships (for PCA)

head(assay(vsd))

# PCA Plot
pca_data <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# PCA plot rendered and saved
ggplot(pca_data, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_minimal() +
  ggtitle("PCA - Variance Stabilizing Transformation")

# Create the plot and save it
p <- ggplot(pca_data, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  ggtitle("PCA - Variance Stabilizing Transformation")

# Save the figure with white background
# ggsave(
#   filename = "/absolute/path/to/R_count_matrix/pca_plot_pre_deseq2.png",
#   plot = p,
#   width = 7, height = 7, dpi = 300,
#   bg = "white"
# )

#########################################
## Plot Raw vs VST Distributions (OPTIONAL)
#########################################

# Extract raw counts
raw_counts <- counts(dds, normalized = FALSE)

# Extract VST transformed counts
vst_counts <- assay(vsd)

# Turn both matrices into long format for ggplot
raw_long <- as.data.frame(raw_counts) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Count") %>%
  mutate(Data = "Raw counts (log10 + 1)")

vst_long <- as.data.frame(vst_counts) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Count") %>%
  mutate(Data = "VST counts")

# Combine
combined <- bind_rows(raw_long, vst_long)

# Apply log10 transformation to raw counts only
combined <- combined %>%
  mutate(Count = ifelse(Data == "Raw counts (log10 + 1)", log10(Count + 1), Count))

# Plot
ggplot(combined, aes(x = Count, fill = Data)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Data, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Raw vs VST Counts",
       x = "Expression value",
       y = "Density") +
  scale_fill_manual(values = c("Raw counts (log10 + 1)" = "#66c2a5", "VST counts" = "#fc8d62")) +
  theme(legend.position = "none")


##################################################
## GENERATING Mean-Variance Plot 
# for Raw Counts vs VST Counts (OPTIONAL)
##################################################

# Extract counts
raw_counts <- counts(dds, normalized = FALSE)
vst_counts <- assay(vsd)

# Compute mean and variance for each gene (row)
raw_mean_var <- data.frame(
  mean = rowMeans(raw_counts),
  variance = apply(raw_counts, 1, var),
  type = "Raw counts"
)

vst_mean_var <- data.frame(
  mean = rowMeans(vst_counts),
  variance = apply(vst_counts, 1, var),
  type = "VST counts"
)

# Combine into one dataframe
mean_var_df <- bind_rows(raw_mean_var, vst_mean_var)

# Plot mean vs variance
ggplot(mean_var_df, aes(x = mean, y = variance, color = type)) +
  geom_point(alpha = 0.5, size = 1) + # it can be used also "geom_smooth"
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ type, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean-Variance Relationship Before and After VST",
    x = "Mean expression (log10 scale)",
    y = "Variance (log10 scale)"
  ) +
  scale_color_manual(values = c("Raw counts" = "#1f78b4", "VST counts" = "#33a02c")) +
  theme(legend.position = "none")

##################################################

######## PART 2: DESeq

# Now run DESeq2
dds <- DESeq(dds)

# (Optional: get results to inspect next)
res <- results(dds)
summary(res)

view(res)

# === PREPARE RESULTS ===
# Already done: 
# dds <- DESeq(dds)
# res <- results(dds)

# === MA Plot ===
png("ma_plot.png", width = 800, height = 600)
plotMA(res, ylim = c(-5, 5), main = "MA Plot - Tissue1 vs Tissue2")
dev.off()
# ✅ X-axis = "mean normalized counts" 
# the x-axis represents the:
# log10-scaled (mean normalized Counts + 1) or the baseMean, which means the normalized counts across all samples for each gene.
# ✅ Y-axis = "log2 fold change" (LFC)

# Creating a prettier MA plot
# Create data frame from results
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "Gene_Id") %>%
  mutate(
    mean_expression = baseMean,
    log2FC = log2FoldChange,
    significant = case_when(
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

view(res_df)

library(ggrepel)

# Join Gene_name into res_df
res_df <- res_df %>%
  left_join(count_matrix %>% select(Gene_name, Gene_Id), by = "Gene_Id") %>%
  relocate(Gene_name, .before = Gene_Id)  # Move Gene_name to the left

View(res_df)

# Select top 10 most significant genes (lowest adjusted p-values)
top10_genes <- res_df %>%
  filter(!is.na(padj)) %>%          # Only keep genes with padj calculated
  arrange(padj) %>%                 # Sort by increasing padj
  slice(1:10)                       # Pick top 10

# Make the MA plot
MA <- ggplot(res_df, aes(x = log10(mean_expression + 1), y = log2FC)) +
  geom_point(aes(color = significant), alpha = 0.5, size = 1.5) +
  geom_text_repel(
    data = top10_genes,
    aes(label = Gene_name),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not significant" = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(ylim = c(-15,15)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Custom MA Plot - Tissue1 vs Tissue2",
    x = "log10(Mean Normalized Counts + 1)",
    y = "log2 Fold Change",
    color = "DE status"
  )

MA

# Save plot
# ggsave(
#   filename = "/absolute/path/to/R_count_matrix/ma_plot_prettier.png",
#   plot = MA,
#   width = 9,
#   height = 7,
#   dpi = 300,
#   bg = "white"
# )


# === Volcano Plot: Reference Tissue: Tissue1 ===

# Load ggrepel
library(ggrepel)

# Merge DESeq2 results with Gene_name
res_annot <- res_df %>%
  mutate(
    volcano_score = -log10(pvalue + 1e-300) * abs(log2FoldChange),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "Neutral"
      ),
    # Cap the -log10(pvalue) for plotting purposes
    plot_y = pmin(-log10(pvalue + 1e-300), 400) # We are putting y-axis a limit of 400, this will drag those few genes extremely upregulated, otherwise these genes are not possible to see in the volcano plot. In ggplot, 'plot_y' has to be in argument 'y ='
    )

view(res_annot)

# COMMENT: when calculating volcano_score, you can be a bit safer against tiny p-values by adding a +1e-300 to avoid possible -Inf from log10(0).
# This is especially useful if you have very very small p-values that could numerically reach 0 (floating point limitations).

# Select top 10 upregulated and top 10 downregulated genes
top_up <- res_annot %>%
  filter(direction == "Up") %>%
  arrange(desc(volcano_score)) %>%
  slice_head(n = 10)

top_down <- res_annot %>%
  filter(direction == "Down") %>%
  arrange(desc(volcano_score)) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

# Add label column: only top genes labeled
res_annot <- res_annot %>%
  mutate(label = ifelse(Gene_Id %in% top_genes$Gene_Id, Gene_name, NA))

# Plot volcano
volcano <- ggplot(res_annot, aes(x = log2FoldChange, y = plot_y)) + # y = -log10(pvalue)
  geom_point(aes(color = direction), alpha = 0.5, size = 2.2, stroke = 1, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # p-value threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +       # LFC thresholds
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Neutral" = "grey70")) +
  geom_point(aes(color = direction), alpha = 0.5, size = 1.5, stroke = 1, shape = 21) + # to highlight the significant genes
  geom_text_repel(
    data = subset(res_annot, !is.na(label)),
    aes(label = label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Tissue1 vs Tissue2 (Tissue1 as Reference)",
    subtitle = "Dashed lines: p-value < 0.05, |log2FC| > 1",
    x = "Log2 Fold Change",
    y = "-Log10 p-value",
    color = "Direction"
  ) +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, 5)) + # Adjust x-axis scale
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 100))  # Adjust y-axis scale

# Print
print(volcano)


# Save plot
# ggsave(
#   filename = "/absolute/path/to/R_count_matrix/volcano_plot_top20_balanced_Tissue2.png",
#   plot = volcano,
#   width = 9,
#   height = 7,
#   dpi = 600,
#   bg = "white"
# )

# === Volcano Plot: (changing) Reference Tissue: Tissue2 ===

# To check which tissue upregulates and which downregulates genes
levels(dds$tissue) # output: "Tissue1"  "Tissue2"

# COMMENT: The output indicates that the last tissue, "Tissue2" (the last level alphabetically), is compared to "Tissue1" (the first level alphabetically)
# COMMENT: Therefore:
# Positive log2FoldChange (right side of the volcano plot) indicates genes that are upregulated in Tissue2 compared to Tissue1, based on 'volcano_score'. This means these genes have higher expression levels in your Tissue2 samples.
# Negative log2FoldChange (left side of the volcano plot) indicates genes that are downregulated in Tissue2 compared to Tissue1, also based on 'volcano_score'. This also means these genes are upregulated in Tissue1 compared to Tissue2 (they have higher expression levels in your Tissue1 samples).
# Therefore: the downregulated genes (negative LFC, on the left side of your volcano plot) correspond to genes that have higher expression in Tissue1 compared to Tissue2. 

# If you want Tissue1 as the upregulated direction (that shows LFC > 0 upregulated genes), you need to relevel your tissue factor to make Tissue2 the reference instead:
dds$tissue <- relevel(dds$tissue, ref = "Tissue2")

# After releveling, run again DESeq 
dds <- DESeq(dds)

# (Optional: extract the new results an inspect)
res <- results(dds)
summary(res)

# Rebuild your res_df table (now res_df has the correct flipped LFC directions)
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "Gene_Id") %>%
  mutate(
    mean_expression = baseMean,
    log2FC = log2FoldChange,
    significant = case_when(
      padj < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(count_matrix %>% select(Gene_name, Gene_Id), by = "Gene_Id") %>%
  relocate(Gene_name, .before = Gene_Id)

# Volcano plot adapted for Tissue1 vs Tissue2
# Create res_annot
res_annot <- res_df %>%
  mutate(
    volcano_score = -log10(pvalue + 1e-300) * abs(log2FoldChange),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",    # Up = upregulated in Tissue1
      padj < 0.05 & log2FoldChange < 0 ~ "Down",  # Down = upregulated in Tissue2
      TRUE ~ "Neutral"
    ),
    plot_y = pmin(-log10(pvalue + 1e-300), 400)  # Cap y at 400
  )

# Select top 10 up and down
top_up <- res_annot %>%
  filter(direction == "Up") %>%
  arrange(desc(volcano_score)) %>%
  slice_head(n = 10)

top_down <- res_annot %>%
  filter(direction == "Down") %>%
  arrange(desc(volcano_score)) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

# Label top genes
res_annot <- res_annot %>%
  mutate(label = ifelse(Gene_Id %in% top_genes$Gene_Id, Gene_name, NA))

# Volcano plot
volcano_Tissue1_vs_Tissue2 <- ggplot(res_annot, aes(x = log2FoldChange, y = plot_y)) +
  geom_point(aes(color = direction), alpha = 0.5, size = 1.5, stroke = 1, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Neutral" = "grey70")) +
  geom_text_repel(
    data = subset(res_annot, !is.na(label)),
    aes(label = label),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Tissue1 vs Tissue2 (Tissue2 as Reference)",
    subtitle = "Dashed lines: p-value < 0.05, |log2FC| > 1",
    x = "Log2 Fold Change (positive = Tissue1 up)",
    y = "-Log10 p-value",
    color = "Direction"
  ) +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-20, 20, 5)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 100))

# Print
print(volcano_Tissue1_vs_Tissue2)

# Save plot
# ggsave(
#   filename = "/absolute/path/to/R_count_matrix/volcano_plot_top20_balanced_Tissue1.png",
#   plot = volcano_Tissue1_vs_Tissue2,
#   width = 9,
#   height = 7,
#   dpi = 600,
#   bg = "white"
# )

# COMMENT: Why is 'Alb' has the highest raw count of all Tissue1 genes it is not in the top ten most upregulated genes? To answer this question, run this code:
res_annot %>%
  filter(Gene_name %in% c("Alb","Akr1c14", "Iqgap2", 
                        "Akr1c12", "Rida", "Cp", 
                        "Slc25a13","Rad51b", "Id2",
                        "Gstt3", "Nrn1")) %>%
  select(Gene_name, baseMean, log2FoldChange, pvalue, padj, volcano_score) %>%
  arrange(desc(volcano_score))
# here you can see that the ranking of top upregulated genes and the labeling of these in the volcano plot is based on "volcano_score", which a calculation based on LFC and pvalue. The gene Alb, having actually the highest LFC, has its p-value relatively higher (though still very small), thus it has a lower 'volcano plot value.

library(patchwork)

combined_volcano <- volcano + volcano_Tissue1_vs_Tissue2 +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Volcano Plots: Tissue2 vs Tissue1")

# ggsave(
#   filename = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/R_count_matrix/combined_volcano_plots.pdf",
#   plot = combined_volcano,
#   width = 16,    # wide enough for two plots
#   height = 8,    # decent height
#   dpi = 300,     # good resolution
#   device = "pdf" # important: save as PDF
# )



# === Export Significant DEGs ===
# 1. Identify significant DEGs
res_sig <- res %>% 
  as.data.frame() %>%
  rownames_to_column("Gene_Id") %>%
  left_join(count_matrix %>% dplyr::select(Gene_name, Gene_Id), by = "Gene_Id") %>% # adding column "Gene_name"
  dplyr::select(Gene_name, Gene_Id, everything()) %>%
  filter(padj < 0.05) %>% # Keep only genes whose  padj is below 0.05. ✅ These are called Differentially Expressed Genes (DEGs).
  arrange(padj) # Genes with the smallest padj (most significant) appear first.

# write_csv(res_sig, "/absolute/path/to/R_count_matrix/DE_genes_padj_0.05.csv")

# 2. Separate Upregulated (Tissue1) and Downregulated (Tissue2)

# PLEASE REMEMBER THIS: 
# > levels(dds$tissue) 
# [1] "Tissue2" "Tissue1" # The first level ("Tissue2") is the reference.
# Interpretation:
# Positive LFC (> 0) means that Gene is upregulated in Tissue1 compared to Tissue2.
# Negative LFC (< 0) means that Gene is downregulated in Tissue1 (i.e., higher in Tissue2).

# Upregulated in Tissue1
res_up_Tissue1 <- res_sig %>%
  filter(log2FoldChange > 0) %>%
  arrange(padj)

# Upregulated in Tissue2
res_up_Tissue2 <- res_sig %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj)

# 3. Export to CSV
# write_csv(res_up_Tissue1, "/absolute/path/to/R_count_matrix/DEGs_up_in_Tissue1_padj0.05.csv")
# write_csv(res_up_Tissue2, "/absolute/path/to/R_count_matrix/DEGs_up_in_Tissue2_padj0.05.csv")

# COMMENT:
# Filter padj < 0.05	Keep only statistically significant genes
# Filter log2FoldChange > 0	Genes expressed more in Tissue1 than Tissue2
# Filter log2FoldChange < 0	Genes expressed more in Tissue2 than Tissue1
# Gene Behavior: Alb has log2FC > 0 -> Alb is upregulated in Tissue1
# Gene Behavior: Myh1 has log2FC < 0 -> Myh1 is upregulated in Tissue2

# 3. Count how many DEGs are upregulated in Tissue1 and in Tissue2
res_sig %>%
  mutate(direction = case_when(
    log2FoldChange > 0 ~ "Up_in_Tissue1",
    log2FoldChange < 0 ~ "Up_in_Tissue2"
  )) %>%
  count(direction) %>%
  mutate(percent = round(100 * n / sum(n), 1))


# === PCA again (now post-DESeq2) ===
vsd <- vst(dds, blind = FALSE)  # same as before, but not "blind" to design

pca_data <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = tissue)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_rect(fill = "white")) +
  ggtitle("PCA After DESeq2 Normalization")

ggsave("pca_post_deseq2.png", plot = p, width = 7, height = 7)

# === Heatmap of top 30 DE genes ===
# Select top DE genes
top_genes <- res_sig %>%
  slice_min(padj, n = 30) %>%
  pull(Gene_Id)

# Extract the variance-stabilized counts (VST matrix)
vsd_mat <- assay(vsd)

# Subset to top genes
heatmap_mat <- vsd_mat[top_genes, ]

# Map Gene_Id to Gene_name
gene_name_map <- res_sig %>% select(Gene_Id, Gene_name)
rownames(heatmap_mat) <- gene_name_map$Gene_name[match(rownames(heatmap_mat), gene_name_map$Gene_Id)]

# === Optional: Reorder sample columns manually ===

# Check actual column names
actual_samples <- colnames(heatmap_mat)

# Match manually — ensure sample names exist
ordered_samples <- c(
  "L1",  "M29",
  "L2",  "M30",
  "L3",  "M31",
  "L4",  "M32"
)

# Filter only those present
ordered_samples <- ordered_samples[ordered_samples %in% actual_samples]

# Reorder columns
heatmap_mat <- heatmap_mat[, ordered_samples]

# Reorder sample annotation
sample_annot <- data.frame(Tissue = colData(vsd)$tissue)
rownames(sample_annot) <- colnames(vsd)
sample_annot <- sample_annot[ordered_samples, , drop = FALSE]

# === Plot and save at high resolution ===
# output_path <- "/absolute/path/to/R_count_matrix/heatmap_top30_named_reordered_600dpi.png"

# Modify colors of "Tissue" annotation
annotation_colors <- list(
  Tissue = c(Tissue1 = "#90EE90",  # light green
             Tissue2 = "#AB82FF") # lavender
)

png(filename = output_path, width = 8, height = 10, units = "in", res = 600)

pheatmap(
  heatmap_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = sample_annot,
  annotation_colors = annotation_colors, 
  scale = "row",
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = viridis(100)
)

dev.off()

# palette options:
# color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100) # default
# color <- colorRampPalette(rev(brewer.pal(n = 9, name = "PRGn")))(100)
# color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
# custom palette
# color <- colorRampPalette(c("black", "yellow", "red"))(100)
# color <- colorRampPalette(c("#6A00B0", "white", "#FF6F00"))(100)
# library(viridis)
# color <- viridis(100)  # Options: "viridis", "magma", "plasma", "inferno"
# library(RColorBrewer)
# display.brewer.all()
















#####################################################
### GSEA (Gene Set Enrichment Analysis)
#####################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db") # Mouse gene annotations

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Making the ranking of genes based on LFC
# Option 1: Prepare ranked gene list (e.g., for GSEA)
# gene_list <- res_annot %>%
#   drop_na(log2FoldChange, padj) %>%
#   arrange(desc(log2FoldChange)) %>%
#   distinct(Gene_name, .keep_all = TRUE) %>%  # Avoid duplicate symbols
#   pull(log2FoldChange, name = Gene_name)

# Option 2: Use results with gene symbols and log2FoldChange (recommended)
# ranked_genes <- res_annot %>%
#   drop_na(Gene_name, log2FoldChange) %>%
#   distinct(Gene_name, .keep_all = TRUE) %>%
#   arrange(desc(log2FoldChange)) %>%
#   pull(log2FoldChange, name = Gene_name)

# Option 3: Use the Wald 'stat'-ranked gene list
# Assuming res_annot is already created with Gene_name and stat
gene_stat_ranks <- res_annot %>%
  drop_na(Gene_name, stat) %>%
  distinct(Gene_name, .keep_all = TRUE) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = Gene_name)

head(gene_stat_ranks) # Assure gene_stat_ranks is a named numeric vector
# Pygm     Pfkm     Eno3      Srl  Tmem38a    Aldoa 
# 47.74839 38.35730 37.05232 36.75979 36.70443 36.49116 

# Run GSEA using GO Biological Processes (BP)
gsea_go <- gseGO(
  geneList = gene_stat_ranks,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",     # Because we're using Gene_name / gene symbols
  ont = "BP",             # Choose GO ontology: "BP", "MF", or "CC"
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Dotplot of top enriched terms
dotplot(gsea_go, showCategory = 15, title = "GSEA - GO BP Enrichment")

# Enrichment map (network-style)
gsea_go <- pairwise_termsim(gsea_go)
# Select most significant ones (e.g., padj < 0.01 and NES > 1.5)
top_terms <- gsea_go@result %>%
  dplyr::filter(p.adjust < 0.01, abs(NES) > 1.5) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 15)
emapplot(gsea_go, showCategory = top_terms$ID)
# Alternative
gsea_go_simple <- simplify(gsea_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
emapplot(pairwise_termsim(gsea_go_simple))

# Ridge plot of enrichment scores
ridgeplot(gsea_go, showCategory = 20)

# GSEA Plot for a Specific Term
gseaplot2(gsea_res, geneSetID = gsea_res$ID[1], title = gsea_res$Description[1])

# Saving GSEA analysis as csv
# Convert to data frame
gsea_df <- as.data.frame(gsea_go)

# Write to CSV
# write_csv(gsea_df, "/absolute/path/to/R_count_matrix/GSEA_GO_BP_results.csv")

# Use filtered significant genes (Gene_name)
gene_sig <- res_sig$Gene_name[!is.na(res_sig$Gene_name)]


# GO Overrepresentation Analysis (Optional but Recommended)
go_enrich <- enrichGO(
  gene         = gene_sig,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",       # Can use "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

dotplot(go_enrich, showCategory = 15, title = "GO Enrichment - Overrepresentation")

# GO Enrichment
go_df <- as.data.frame(go_enrich)

# Write to CSV
# write_csv(go_df, "/absolute/path/to/R_count_matrix/GO_Enrichment_BP_results.csv")


#======== Performing GSEA with KEGG and Reactome pathways ========
# library(clusterProfiler)
# library(org.Mm.eg.db)

# Run GSEA using KEGG pathways
# Convert gene symbols to Entrez IDs
gene_df <- bitr(names(gene_stat_ranks), fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Create a named vector of stats with Entrez IDs
gene_list <- gene_stat_ranks[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID

# Run GSEA
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = 'mmu',
                     nPerm = 1000,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = TRUE)


# Dotplot of top enriched terms
dotplot(gsea_kegg, showCategory = 15, title = "GSEA - KEGG Enrichment")

# run GSEA using Reactome pathways
BiocManager::install("ReactomePA")
library(ReactomePA)

gsea_reactome <- gsePathway(geneList = gene_list,
                            organism = "mouse",
                            nPerm = 1000,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 0.05,
                            verbose = TRUE)

# Dotplot of top enriched terms
dotplot(gsea_reactome, showCategory = 15, title = "GSEA - Reactome pathways analysis")

# ======== Visualizing Enrichment results by tissue ===========

# 1. Genes upregulated in Tissue1
gene_list_Tissue1 <- res_annot %>%
  filter(log2FoldChange > 0) %>%               # Positive LFC = Tissue1 up
  drop_na(Gene_name, stat) %>%
  distinct(Gene_name, .keep_all = TRUE) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = Gene_name)


# 2. Genes upregulated in Tissue2
gene_list_Tissue2 <- res_annot %>%
  filter(log2FoldChange < 0) %>%               # Negative LFC = Tissue2 up
  drop_na(Gene_name, stat) %>%
  distinct(Gene_name, .keep_all = TRUE) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = Gene_name)



# List of your gene symbols
symbols_Tissue1 <- names(gene_list_Tissue1)

# Check how many are valid mouse symbols
valid_symbols <- symbols_Tissue1[symbols_Tissue1 %in% keys(org.Mm.eg.db, keytype = "SYMBOL")]

length(valid_symbols)         # How many are valid
length(symbols_Tissue1)         # Total

# Also confirm gene names are recognized
valid_symbols_Tissue1 <- names(gene_list_Tissue1)[names(gene_list_Tissue1) %in% keys(org.Mm.eg.db, keytype = "SYMBOL")]
length(valid_symbols_Tissue1)

# Optional: See a few missing ones
setdiff(symbols_Tissue1, valid_symbols)[1:10]


table(res_annot$direction)

dds$tissue <- relevel(dds$tissue, ref = "Tissue2")







# GSEA for Tissue1
gsea_Tissue1 <- gseGO(
  geneList = gene_list_Tissue1,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# GSEA for Tissue2
gsea_Tissue2 <- gseGO(
  geneList = gene_list_Tissue2,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

dotplot(gsea_Tissue1, showCategory = 15, title = "Tissue1 - GO BP Enrichment")

dotplot(gsea_Tissue2, showCategory = 15, title = "Tissue2 - GO BP Enrichment")


# library(clusterProfiler)

# Define gene lists for each tissue
gene_list_Tissue1 <- res_Tissue1 %>%
  filter(padj < 0.05) %>%
  pull(Gene_name)

gene_list_Tissue2 <- res_Tissue2 %>%
  filter(padj < 0.05) %>%
  pull(Gene_name)

# Create a list of gene sets
gene_clusters <- list(Tissue1 = gene_list_Tissue1, Tissue2 = gene_list_Tissue2)

# Perform GO enrichment analysis
go_results <- compareCluster(geneClusters = gene_clusters,
                             fun = "enrichGO",
                             OrgDb = org.Mm.eg.db,
                             keyType = "SYMBOL",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)

# Visualize the results
dotplot(go_results, showCategory = 20) + ggtitle("GO Enrichment by Tissue")
ridgeplot(go_results, showCategory = 20) + ggtitle("GO Enrichment by Tissue")


