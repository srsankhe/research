#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis
# Generic script for pseudo-bulk analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript deseq2_analysis.R <config.json>\n")
  cat("\nConfig JSON format:\n")
  cat("{\n")
  cat('  "dataset_id": "GSE123456",\n')
  cat('  "csv_dir": "csv_output",\n')
  cat('  "output_dir": "deseq2_output",\n')
  cat('  "samples": ["sample1", "sample2", "sample3"],\n')
  cat('  "conditions": ["control", "control", "treatment"],\n')
  cat('  "comparison": {"numerator": "treatment", "denominator": "control"},\n')
  cat('  "additional_metadata": {"cell_type": ["CD45+", "CD45-", ...]} (optional)\n')
  cat("}\n")
  quit(status = 1)
}

config_file <- args[1]
cat(paste0("Loading configuration from: ", config_file, "\n"))
config <- fromJSON(config_file)

# Extract configuration
dataset_id <- config$dataset_id
csv_dir <- config$csv_dir
output_dir <- config$output_dir
samples <- config$samples
conditions <- config$conditions
comparison <- config$comparison
additional_metadata <- config$additional_metadata

# Validate inputs
if (length(samples) != length(conditions)) {
  stop("Number of samples must match number of conditions")
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set working directory to csv output
setwd(csv_dir)

cat("\n=== Loading Data ===\n")
count_matrices <- list()

for (sample in samples) {
  cat(paste0("  Loading ", sample, "...\n"))
  matrix_file <- paste0(sample, "_matrix.csv")
  
  if (!file.exists(matrix_file)) {
    stop(paste0("Matrix file not found: ", matrix_file))
  }
  
  # Read matrix
  mat <- read.csv(matrix_file, check.names = FALSE)
  gene_names <- mat[, 1]
  mat <- mat[, -1]
  
  # Sum across all cells to get pseudo-bulk counts
  pseudo_bulk <- rowSums(mat)
  names(pseudo_bulk) <- gene_names
  count_matrices[[sample]] <- pseudo_bulk
}

# Combine into a single count matrix
cat("\nCreating pseudo-bulk count matrix...\n")
count_matrix <- do.call(cbind, count_matrices)
colnames(count_matrix) <- samples

cat(paste0("  Matrix dimensions: ", nrow(count_matrix), " genes x ", 
           ncol(count_matrix), " samples\n"))

# Remove genes with zero counts
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
cat(paste0("  After filtering zero counts: ", nrow(count_matrix), " genes\n"))

# Create sample metadata
cat("\nCreating sample metadata...\n")
coldata <- data.frame(
  sample = samples,
  condition = conditions,
  row.names = samples
)

# Add additional metadata if provided
if (!is.null(additional_metadata)) {
  for (col_name in names(additional_metadata)) {
    coldata[[col_name]] <- additional_metadata[[col_name]]
  }
}

# Convert to factors
coldata$condition <- factor(coldata$condition, 
                           levels = c(comparison$denominator, comparison$numerator))

print(coldata)

# Save pseudo-bulk count matrix and metadata
cat("\nSaving data...\n")
setwd("..")  # Go back to parent directory
write.csv(count_matrix, file.path(output_dir, "pseudobulk_counts.csv"))
write.csv(coldata, file.path(output_dir, "sample_metadata.csv"))

cat("\n=== Running DESeq2 Analysis ===\n")

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

# Run DESeq2
cat("Running differential expression analysis...\n")
dds <- DESeq(dds)

# Get results
cat("\nExtracting results...\n")
comparison_name <- paste0(comparison$numerator, "_vs_", comparison$denominator)
res <- results(dds, contrast = c("condition", comparison$numerator, comparison$denominator))
summary(res)

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res_ordered),
          file.path(output_dir, paste0(comparison_name, "_results.csv")))

# Get significant genes
sig_genes <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
cat(paste0("\nSignificant genes (padj < 0.05, |log2FC| > 1): ", nrow(sig_genes), "\n"))
cat(paste0("  Upregulated in ", comparison$numerator, ": ", 
           sum(sig_genes$log2FoldChange > 0), "\n"))
cat(paste0("  Downregulated in ", comparison$numerator, ": ", 
           sum(sig_genes$log2FoldChange < 0), "\n"))

write.csv(as.data.frame(sig_genes),
          file.path(output_dir, paste0(comparison_name, "_significant_genes.csv")))

# Generate plots
cat("\n=== Generating Plots ===\n")

# MA plot
cat("Creating MA plot...\n")
png(file.path(output_dir, "MA_plot.png"), width = 800, height = 600)
plotMA(res, ylim = c(-5, 5),
       main = paste0("MA Plot: ", comparison_name))
dev.off()

# Volcano plot
cat("Creating volcano plot...\n")
volcano_data <- as.data.frame(res)
volcano_data$significant <- ifelse(volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1,
                                   "Significant", "Not significant")

png(file.path(output_dir, "volcano_plot.png"), width = 1000, height = 800)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = paste0("Volcano Plot: ", comparison_name),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "top")
dev.off()

# PCA plot
cat("Creating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)

png(file.path(output_dir, "PCA_plot.png"), width = 800, height = 600)
p <- plotPCA(vsd, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA Plot: Samples")
print(p)
dev.off()

# Heatmap of top 50 DE genes
cat("Creating heatmap...\n")
top_genes <- head(rownames(res_ordered), 50)
mat_top <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(
  Condition = coldata$condition
)
rownames(annotation_col) <- colnames(mat_scaled)

png(file.path(output_dir, "heatmap_top50_genes.png"), width = 1000, height = 1200)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = paste0("Top 50 DE Genes: ", comparison_name))
dev.off()

# Dispersion plot
cat("Creating dispersion plot...\n")
png(file.path(output_dir, "dispersion_plot.png"), width = 800, height = 600)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Summary report
cat("\n=== Creating Summary Report ===\n")
sink(file.path(output_dir, "analysis_summary.txt"))

cat("DESeq2 Differential Expression Analysis Summary\n")
cat("================================================\n")
cat(paste0("Dataset: ", dataset_id, "\n\n"))

cat("Dataset Information:\n")
cat(paste0("  Total genes analyzed: ", nrow(count_matrix), "\n"))
cat(paste0("  Total samples: ", ncol(count_matrix), "\n\n"))

cat("Sample Details:\n")
print(coldata)

cat(paste0("\n\nComparison: ", comparison_name, "\n"))
cat("--------------------------------------\n")
summary(res)
cat(paste0("\nSignificant genes (padj < 0.05, |log2FC| > 1): ", nrow(sig_genes), "\n"))

if (nrow(sig_genes) > 0) {
  cat(paste0("  Upregulated in ", comparison$numerator, ": ", 
             sum(sig_genes$log2FoldChange > 0), "\n"))
  cat(paste0("  Downregulated in ", comparison$numerator, ": ", 
             sum(sig_genes$log2FoldChange < 0), "\n"))
  
  cat("\n\nTop 20 Upregulated Genes:\n")
  top_up <- head(subset(res_ordered, log2FoldChange > 0), 20)
  if (nrow(top_up) > 0) {
    print(as.data.frame(top_up)[, c("log2FoldChange", "padj")])
  }
  
  cat("\n\nTop 20 Downregulated Genes:\n")
  top_down <- head(subset(res_ordered, log2FoldChange < 0), 20)
  if (nrow(top_down) > 0) {
    print(as.data.frame(top_down)[, c("log2FoldChange", "padj")])
  }
}

cat("\n\nOutput Files:\n")
cat("  - pseudobulk_counts.csv\n")
cat("  - sample_metadata.csv\n")
cat(paste0("  - ", comparison_name, "_results.csv\n"))
cat(paste0("  - ", comparison_name, "_significant_genes.csv\n"))
cat("  - MA_plot.png\n")
cat("  - volcano_plot.png\n")
cat("  - PCA_plot.png\n")
cat("  - heatmap_top50_genes.png\n")
cat("  - dispersion_plot.png\n")

sink()

cat("\n=================================================\n")
cat("Analysis complete!\n")
cat(paste0("Results saved to: ", normalizePath(output_dir), "\n"))
cat("=================================================\n")
