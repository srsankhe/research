#!/usr/bin/env Rscript

# Single-Cell RNA-seq Analysis with Seurat - Version 2
# Simplified and more robust version

cat("=======================================================\n")
cat("Single-Cell Analysis with Seurat\n")
cat("Dataset: GSE282344\n")
cat("=======================================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# Set output directory
output_dir <- "single_cell_analysis"
plots_dir <- file.path(output_dir, "plots")
results_dir <- file.path(output_dir, "results")
data_dir <- file.path(output_dir, "data")

# Sample information (transcriptome samples only)
samples <- c("GSM8641663_NR_JK_047",
             "GSM8641664_NR_JK_048",
             "GSM8641665_NR_JK_049",
             "GSM8641666_NR_JK_053")

# Condition assignments
conditions <- c("control", "control", "treatment", "treatment")

cat("\n=== STEP 1: Loading Single-Cell Data ===\n")

# Load data for each sample
seurat_list <- list()

for (i in 1:length(samples)) {
  sample_name <- samples[i]
  cat(paste0("\nLoading sample: ", sample_name, "\n"))

  # Read files
  genes <- read.table(gzfile(paste0(sample_name, "_genes.tsv.gz")),
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(gzfile(paste0(sample_name, "_barcodes.tsv.gz")),
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  mtx <- readMM(gzfile(paste0(sample_name, "_matrix.mtx.gz")))

  # Set names
  rownames(mtx) <- genes$V2
  colnames(mtx) <- barcodes$V1

  cat(paste0("  Raw: ", nrow(mtx), " genes x ", ncol(mtx), " cells\n"))

  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = mtx,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )

  # Add metadata
  obj$sample <- sample_name
  obj$condition <- conditions[i]
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

  seurat_list[[i]] <- obj

  cat(paste0("  Filtered: ", ncol(obj), " cells\n"))
}

cat("\n=== STEP 2: Merging Samples ===\n")

# Merge
merged <- merge(x = seurat_list[[1]],
               y = seurat_list[2:length(seurat_list)],
               add.cell.ids = samples)

cat(paste0("Total cells: ", ncol(merged), "\n"))
cat(paste0("Total genes: ", nrow(merged), "\n"))

cat("\n=== STEP 3: Quality Control ===\n")

# QC plots before filtering
pdf(file.path(plots_dir, "qc_before_filter.pdf"), width = 14, height = 6)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample")
dev.off()

# Apply QC filters
merged <- subset(merged,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 6000 &
                         percent.mt < 15)

cat(paste0("After QC: ", ncol(merged), " cells\n"))

# QC plots after filtering
pdf(file.path(plots_dir, "qc_after_filter.pdf"), width = 14, height = 6)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample")
dev.off()

cat("\n=== STEP 4: Normalization ===\n")

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)

top10 <- head(VariableFeatures(merged), 10)
cat("Top 10 variable genes:\n")
cat(paste(top10, collapse = ", "), "\n")

cat("\n=== STEP 5: Scaling and PCA ===\n")

merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)

pdf(file.path(plots_dir, "pca_elbow.pdf"), width = 8, height = 6)
ElbowPlot(merged, ndims = 50)
dev.off()

cat("\n=== STEP 6: Clustering and UMAP ===\n")

dims <- 30
merged <- FindNeighbors(merged, dims = 1:dims)
merged <- FindClusters(merged, resolution = 0.5)
merged <- RunUMAP(merged, dims = 1:dims)

cat(paste0("Clusters: ", length(unique(Idents(merged)))), "\n")

# UMAP plots
pdf(file.path(plots_dir, "umap_clusters.pdf"), width = 10, height = 8)
DimPlot(merged, reduction = "umap", label = TRUE, label.size = 5)
dev.off()

pdf(file.path(plots_dir, "umap_by_sample.pdf"), width = 10, height = 8)
DimPlot(merged, reduction = "umap", group.by = "sample")
dev.off()

pdf(file.path(plots_dir, "umap_by_condition.pdf"), width = 10, height = 8)
DimPlot(merged, reduction = "umap", group.by = "condition")
dev.off()

cat("\n=== STEP 7: Find Cluster Markers ===\n")

# Join layers for Seurat 5
merged <- JoinLayers(merged)

all_markers <- FindAllMarkers(merged,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

if (nrow(all_markers) > 0) {
  write.csv(all_markers,
            file.path(results_dir, "cluster_markers_all.csv"),
            row.names = FALSE)

  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

  write.csv(top10_markers,
            file.path(results_dir, "cluster_markers_top10.csv"),
            row.names = FALSE)

  cat(paste0("Markers identified: ", nrow(all_markers), "\n"))
} else {
  cat("No markers identified\n")
}

cat("\n=== STEP 8: Kidney Cell Type Markers ===\n")

kidney_markers <- c(
  "MME", "ANPEP", "CD24", "CD44", "PROM1", "EPCAM",
  "AQP1", "SLC5A2", "SLC12A1", "SLC12A3", "AQP2",
  "PTPRC", "CD3D", "CD79A", "CD14", "CD68",
  "KRT18", "HAVCR1", "LCN2"
)

present <- kidney_markers[kidney_markers %in% rownames(merged)]
cat(paste0("Markers found: ", length(present), "/", length(kidney_markers), "\n"))

if (length(present) > 0) {
  pdf(file.path(plots_dir, "kidney_markers_dotplot.pdf"), width = 12, height = 8)
  print(DotPlot(merged, features = present) + RotatedAxis())
  dev.off()

  if ("MME" %in% present) {
    pdf(file.path(plots_dir, "MME_feature_plot.pdf"), width = 10, height = 8)
    print(FeaturePlot(merged, features = "MME", label = TRUE))
    dev.off()

    pdf(file.path(plots_dir, "MME_violin.pdf"), width = 8, height = 6)
    print(VlnPlot(merged, features = "MME", group.by = "condition"))
    dev.off()
  }
}

cat("\n=== STEP 9: Differential Expression (All Cells) ===\n")

Idents(merged) <- "condition"

de_all <- FindMarkers(merged,
                     ident.1 = "treatment",
                     ident.2 = "control",
                     min.pct = 0.1,
                     logfc.threshold = 0.25)

de_all$gene <- rownames(de_all)
de_all <- de_all %>% arrange(p_val_adj)

write.csv(de_all,
          file.path(results_dir, "DE_treatment_vs_control_all_cells.csv"),
          row.names = FALSE)

sig <- de_all %>% filter(p_val_adj < 0.05)
cat(paste0("Significant DEGs: ", nrow(sig), "\n"))

# Top genes
top_up <- de_all %>% filter(avg_log2FC > 0) %>% head(20)
top_down <- de_all %>% filter(avg_log2FC < 0) %>% head(20)

write.csv(top_up, file.path(results_dir, "top20_up_all_cells.csv"), row.names = FALSE)
write.csv(top_down, file.path(results_dir, "top20_down_all_cells.csv"), row.names = FALSE)

cat("\nTop 10 upregulated:\n")
print(head(top_up[, c("gene", "avg_log2FC", "p_val_adj")], 10))

cat("\nTop 10 downregulated:\n")
print(head(top_down[, c("gene", "avg_log2FC", "p_val_adj")], 10))

if ("MME" %in% rownames(de_all)) {
  cat("\nMME results:\n")
  print(de_all[de_all$gene == "MME", ])
}

cat("\n=== STEP 10: Cluster-Specific DE ===\n")

Idents(merged) <- "seurat_clusters"

for (cluster_id in levels(Idents(merged))) {
  cat(paste0("\nCluster ", cluster_id, "...\n"))

  cluster_cells <- subset(merged, idents = cluster_id)
  n_ctrl <- sum(cluster_cells$condition == "control")
  n_trt <- sum(cluster_cells$condition == "treatment")

  cat(paste0("  Control: ", n_ctrl, ", Treatment: ", n_trt, "\n"))

  if (n_ctrl >= 10 && n_trt >= 10) {
    Idents(cluster_cells) <- "condition"

    tryCatch({
      de_cluster <- FindMarkers(cluster_cells,
                               ident.1 = "treatment",
                               ident.2 = "control",
                               min.pct = 0.1,
                               logfc.threshold = 0.25)

      de_cluster$gene <- rownames(de_cluster)
      de_cluster$cluster <- cluster_id

      write.csv(de_cluster,
                file.path(results_dir, paste0("DE_cluster_", cluster_id, ".csv")),
                row.names = FALSE)

      sig_cluster <- de_cluster %>% filter(p_val_adj < 0.05)
      cat(paste0("  Significant: ", nrow(sig_cluster), "\n"))

    }, error = function(e) {
      cat(paste0("  Error: ", e$message, "\n"))
    })
  } else {
    cat("  Skipped (too few cells)\n")
  }
}

cat("\n=== STEP 11: Save Object ===\n")

saveRDS(merged, file.path(data_dir, "seurat_object.rds"))
cat("Saved to: single_cell_analysis/data/seurat_object.rds\n")

cat("\n=== STEP 12: Summary ===\n")

summary_text <- paste0(
  "=======================================================\n",
  "Single-Cell Analysis Summary\n",
  "=======================================================\n\n",
  "Samples: ", paste(samples, collapse = ", "), "\n",
  "Conditions: ", paste(conditions, collapse = ", "), "\n\n",
  "Cells analyzed: ", ncol(merged), "\n",
  "Genes analyzed: ", nrow(merged), "\n",
  "Clusters identified: ", length(unique(merged$seurat_clusters)), "\n\n",
  "Differential Expression (all cells):\n",
  "  Total significant (padj < 0.05): ", nrow(sig), "\n",
  "  Upregulated: ", sum(sig$avg_log2FC > 0), "\n",
  "  Downregulated: ", sum(sig$avg_log2FC < 0), "\n\n"
)

if ("MME" %in% de_all$gene) {
  mme <- de_all[de_all$gene == "MME", ]
  summary_text <- paste0(summary_text,
    "MME (CD10):\n",
    "  log2FC: ", round(mme$avg_log2FC, 3), "\n",
    "  p-value: ", format(mme$p_val, digits = 3), "\n",
    "  adj p-value: ", format(mme$p_val_adj, digits = 3), "\n\n"
  )
}

summary_text <- paste0(summary_text,
  "Output:\n",
  "  Plots: single_cell_analysis/plots/\n",
  "  Results: single_cell_analysis/results/\n",
  "  Data: single_cell_analysis/data/\n\n",
  "=======================================================\n"
)

writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
cat(summary_text)

cat("\nAnalysis complete!\n")
