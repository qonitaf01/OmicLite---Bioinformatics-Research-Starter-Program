# ==========================================================
# DIFFERENTIAL EXPRESSION ANALYSIS : TUBERCULOSIS VS HEALTHY
# ==========================================================

# Workflow:
# 1. Data acquisition & preprocessing
#    - Download dataset from GEO
#    - Extract expression matrix
#    - Perform log2 transformation if needed
#
# 2. Sample grouping & metadata cleaning
#    - Convert group names into valid R format using make.names()
#    - Convert categorical variables into factors
#
# 3. Differential expression analysis
#    - Create design matrix
#    - Fit linear model using limma (lmFit)
#    - Apply contrasts and compute DEGs using eBayes and topTable
#
# 4. Annotation & visualization
#    - Map probes to gene annotations
#    - Visualize results using boxplot, density plot,
#      UMAP, volcano plot, and heatmap
# =====================================================


# Load the required libraries
library(GEOquery)     # to download dataset and reference directly into R
library(limma)        # to perform DEG analysis
library(pheatmap)     # to generate clustered heatmap of the top significant genes
library(ggplot2)      # to generate volcano plot
library(ggrepel)      # to label genes on volcano plot without overlapping
library(dplyr)        # for data cleaning
library(AnnotationDbi)# to map IDs to gene symbols
library(umap)         # for sample clustering

# ====== 1. DATA ACQUISITION & PRE-PROCESSING ==========

# Download GEO dataset and load it as an ExpressionSet object
# GSEMatrix = TRUE -> load the processed expression matrix
# AnnotGPL = TRUE -> includes platform annotation information
gset <- getGEO("GSE313408", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Extract the gene expression matrix from the dataset
ex <- exprs(gset)

# Examine expression value distribution using quantiles
# qx contains: min, Q1, median, Q3, 99th percentile, max
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

# Determine if log2 transformation is needed
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
# Apply the log2 transformation if needed
if (LogTransform) {
  ex <- log2(ex + 1)
  message("Log2 transformation applied.")
} else {
  message("Log2 transformation not needed.")
}

# ====== 2. SAMPLE GROUPING & METADATA CLEANING ==========

# Get the information about the experiment groups from metadata
group_info <- pData(gset)[["title"]]

# Metadata cleaning: remove the extra text/number to get clean categories
clean_names <- gsub("_sample [0-9]+", "", group_info)
clean_names <- gsub("PBMCs_", "", clean_names)
# Verify the result
table(clean_names)

# Convert the categories name into a valid format for R
groups_rformat <- make.names(clean_names)
# Assign it to the data frame as a factor
gset$group <- factor(groups_rformat)
# Verify the result
table(gset$group)

# Get the levels of the factor
group_name <- levels(gset$group)
print(group_name)

# ======= 3. DIFFERENTIAL EXPRESSION ANALYSIS ========

# Create Design Matrix: Defines the groups in the study
design <- model.matrix(~0 + gset$group)
# colnames(): create a new column name
colnames(design) <- levels(gset$group)

# Define comparison groups : TB Patient vs. Healthy Control
group_patient <- "TB.patient"
group_healthy <- "healthy.control"
contrast_formula <- paste(group_patient, "-", group_healthy, sep = "")
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# Print contrast used in analysis
print(paste("Contrast Analyzed:", contrast_formula))

# ----------------------------
# Step-by-step limma pipeline
# ----------------------------

fit  <- lmFit(ex, design)              # 1. Fit linear model for each gene
fit2 <- contrasts.fit(fit, contrast_matrix) # 2. Apply the contrast (comparison)
fit2 <- eBayes(fit2)                   # 3. Empirical Bayes for variance stability

# -------------------------
# Extract significant DEGs
# -------------------------
# adjust = "fdr" : multiple testing correction using False Discovery Rate
# sort.by = "B"  : rank genes by log-odds of differential expression
# number = Inf   : return all genes
# p.value = 0.05 : significance threshold

topTableResults <- topTable(
  fit2,
  adjust  = "fdr",
  sort.by = "B",
  number  = Inf,
  p.value = 0.05
)

head(topTableResults)

# Create a separate full result table (no p-value filter) for volcano plot
allGenes <- topTable(
  fit2,
  adjust  = "fdr",
  sort.by = "B",
  number  = Inf       # No p.value filter — returns ALL genes
)

# ========= ANNOTATION & VISUALIZATION ============

# -----------------
# Gene Annotation
# -----------------

# Get the platform annotation
# for the specific GEO platform (GPL31250)
gpl_info       <- getGEO("GPL31250")
annotation_map <- Table(gpl_info)

# Match the Differential Expression results (significant only)
# with the annotation table to obtain gene names
final_results <- merge(
  topTableResults,
  annotation_map,
  by.x = "row.names",  # probe IDs from DEG results
  by.y = "ID"          # probe IDs from the platform annotation
)

# Merge ALL genes with annotation for volcano plot
final_results_all <- merge(
  allGenes,
  annotation_map,
  by.x = "row.names",
  by.y = "ID"
)

# Check the resulting annotated table
head(final_results)

# ---------------
# Visualization
# ---------------

# 1. BOXPLOT OF EXPRESSION VALUE DISTRIBUTION
# Boxplot is used to:
# - Check the distribution of expression values across samples
# - Detect potential batch effects
# - Evaluate whether normalization/log transformation looks reasonable

# Set colors based on group
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col     = group_colors,
  las     = 2,
  outline = FALSE,
  main    = "Expression Value Distribution per Sample",
  ylab    = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill   = unique(group_colors),
  cex    = 0.8
)


# 2. EXPRESSION VALUE DISTRIBUTION (DENSITY PLOT)
# Density plot shows the global distribution of gene expression values
# Used to:
# - Check the effect of log transformation
# - Compare distributions between groups

# Combine expression values and group labels into a data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group      = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Gene Expression Value Distribution",
    x     = "Expression Value (log2)",
    y     = "Density"
  )

# 3. UMAP (LOW-DIMENSIONAL VISUALIZATION)

# UMAP is used to:
# - Reduce thousands of genes into 2 dimensions
# - Observe global sample separation

# Transpose the expression matrix
# UMAP expects observations = samples
umap_input <- t(ex)

# Run UMAP
umap_result <- umap(umap_input)

# Save results to a data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: TB Patient vs Healthy Control",
    x     = "UMAP 1",
    y     = "UMAP 2"
  )


# 4. VOLCANO PLOT VISUALIZATION

# Volcano plot combines:
# - Log fold change (biological effect)
# - Statistical significance

# FIX #2 (continued): Use final_results_all (ALL genes) for volcano
# so that non-significant genes (grey dots) are properly shown
volcano_data <- data.frame(
  logFC     = final_results_all$logFC,
  adj.P.Val = final_results_all$adj.P.Val,
  Gene      = final_results_all$GENE_SYMBOL
)

# Classify gene status
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC >  1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

# Visualization
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot: TB Patient vs Healthy Control")


# 5. HEATMAP VISUALIZATION

# Heatmap is used to visualize gene expression patterns
# across samples based on the most significant genes

# Sort significant DEGs by adjusted p-value
sorted_DEG <- final_results[order(final_results$adj.P.Val), ]

# Removes rows where the GENE_SYMBOL column is empty
sorted_DEG_clean <- sorted_DEG %>%
  filter(!is.na(GENE_SYMBOL) & GENE_SYMBOL != "") %>%
  mutate(label = GENE_SYMBOL) %>%
  group_by(label) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(adj.P.Val)

top50_clean <- head(sorted_DEG_clean, 50)

# Extract expression matrix for selected genes
mat_heatmap <- ex[top50_clean$Row.names, ]

gene_label<- make.unique(top50_clean$label)
rownames(mat_heatmap) <- gene_label

# Data cleaning: removes rows where expression value is missing (NA)
mat_heatmap  <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
# Remove genes with 0 variance
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap  <- mat_heatmap[gene_variance > 0, ]

# Column annotation (sample groups)
annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale                    = "row",
  annotation_col           = annotation_col,
  show_colnames            = FALSE,
  show_rownames            = TRUE,
  fontsize_row             = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  main                     = "Top 50 Differentially Expressed Genes"
)


# ========= SAVE RESULTS ===================

#Significant DEGs with full annotation
write.csv(final_results, "GSE313408_DEG_significant_annotated.csv", row.names = FALSE)

message("Analysis completed. Results saved successfully.")