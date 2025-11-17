#!/usr/bin/env Rscript
# Generic Single-Cell RNA-seq Analysis with Seurat
# Works with any 10X Genomics format dataset

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript single_cell_seurat_analysis.R <config.json>\n")
  cat("\nConfig JSON format:\n")
  cat("{\n")
  cat('  "dataset_id": "GSE123456",\n')
  cat('  "data_format": "standard" or "gsm_prefix",\n')
  cat('  "archive_dir": "archives" (optional),\n')
  cat('  "samples": {\n')
  cat('    "sample1": {"condition": "control", "additional_metadata": {...}},\n')
  cat('    "sample2": {"condition": "treatment", "additional_metadata": {...}}\n')
  cat('  },\n')
  cat('  "output_dir": "single_cell_analysis",\n')
  cat('  "qc_filters": {\n')
  cat('    "min_features": 200,\n')
  cat('    "max_features": 5000,\n')
  cat('    "max_mt_percent": 10\n')
  cat('  },\n')
  cat('  "analysis_params": {\n')
  cat('    "n_variable_features": 2000,\n')
  cat('    "n_pca_dims": 30,\n')
  cat('    "clustering_resolution": 0.5\n')
  cat('  },\n')
  cat('  "genes_of_interest": ["Gene1", "Gene2"] (optional)\n')
  cat("}\n")
  quit(status = 1)
}

config_file <- args[1]
cat("======================================================\n")
cat("Single-Cell RNA-seq Analysis with Seurat\n")
cat(paste0("Config: ", config_file, "\n"))
cat("======================================================\n\n")

# Load configuration
config <- fromJSON(config_file)

# Extract configuration parameters
dataset_id <- config$dataset_id
data_format <- config$data_format
archive_dir <- ifelse(is.null(config$archive_dir), NULL, config$archive_dir)
samples <- config$samples
output_dir <- config$output_dir
qc_filters <- config$qc_filters
analysis_params <- config$analysis_params
genes_of_interest <- config$genes_of_interest

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "data"), showWarnings = FALSE)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)

# Function to load 10X data
load_sample_data <- function(sample_name, sample_info, data_format, archive_dir) {
  cat(paste0("  Loading ", sample_name, "...\n"))
  
  if (data_format == "standard") {
    # Standard 10X format in sample directories
    data_dir <- ifelse(is.null(archive_dir), 
                      sample_name,
                      file.path(archive_dir, sample_name))
    data <- Read10X(data.dir = data_dir, gene.column = 2)
    
  } else if (data_format == "gsm_prefix") {
    # Files with GSM prefix naming
    mtx_file <- paste0(archive_dir, "/", sample_name, "_matrix.mtx.gz")
    genes_file <- paste0(archive_dir, "/", sample_name, "_genes.tsv.gz")
    barcodes_file <- paste0(archive_dir, "/", sample_name, "_barcodes.tsv.gz")
    
    genes <- read.table(gzfile(genes_file), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    barcodes <- read.table(gzfile(barcodes_file), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    mtx <- readMM(gzfile(mtx_file))
    
    # Set names
    if (ncol(genes) >= 2) {
      rownames(mtx) <- genes$V2
    } else {
      rownames(mtx) <- genes$V1
    }
    colnames(mtx) <- barcodes$V1
    data <- mtx
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,
    min.features = qc_filters$min_features
  )
  
  # Add metadata
  seurat_obj$sample <- sample_name
  for (meta_name in names(sample_info)) {
    seurat_obj[[meta_name]] <- sample_info[[meta_name]]
  }
  
  return(seurat_obj)
}

# Load all samples
cat("\n=== STEP 1: Loading Data ===\n")
seurat_objects <- list()
sample_names <- names(samples)

for (sample_name in sample_names) {
  sample_info <- samples[[sample_name]]
  seurat_obj <- load_sample_data(sample_name, sample_info, data_format, archive_dir)
  seurat_objects[[sample_name]] <- seurat_obj
  cat(paste0("    ", sample_name, ": ", ncol(seurat_obj), " cells, ", 
             nrow(seurat_obj), " genes\n"))
}

# Merge all samples
cat("\n=== STEP 2: Merging Samples ===\n")
if (length(seurat_objects) > 1) {
  combined <- merge(seurat_objects[[1]], 
                   y = seurat_objects[2:length(seurat_objects)],
                   add.cell.ids = sample_names)
} else {
  combined <- seurat_objects[[1]]
}

cat(paste0("Combined dataset: ", ncol(combined), " cells, ", 
           nrow(combined), " genes\n"))

# Calculate QC metrics
cat("\n=== STEP 3: Quality Control ===\n")
cat("Calculating QC metrics...\n")

# Mitochondrial percentage (try both mouse and human patterns)
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-|^Mt-")

# QC plots before filtering
pdf(file.path(output_dir, "plots", "01_QC_before_filter.pdf"), width = 12, height = 8)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1, group.by = "sample")
dev.off()

# Feature scatter
pdf(file.path(output_dir, "plots", "02_QC_feature_scatter.pdf"), width = 12, height = 5)
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# Filter cells
cat("Filtering cells...\n")
cat(paste0("  Criteria: min features = ", qc_filters$min_features, 
           ", max features = ", qc_filters$max_features,
           ", max MT% = ", qc_filters$max_mt_percent, "\n"))

combined <- subset(combined, 
                  subset = nFeature_RNA > qc_filters$min_features & 
                          nFeature_RNA < qc_filters$max_features & 
                          percent.mt < qc_filters$max_mt_percent)

cat(paste0("After filtering: ", ncol(combined), " cells\n"))

# QC plots after filtering
pdf(file.path(output_dir, "plots", "03_QC_after_filter.pdf"), width = 12, height = 8)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1, group.by = "sample")
dev.off()

# Normalize and find variable features
cat("\n=== STEP 4: Normalization and Feature Selection ===\n")
cat("Normalizing data...\n")
combined <- NormalizeData(combined)

cat("Finding variable features...\n")
combined <- FindVariableFeatures(combined, 
                                 selection.method = "vst", 
                                 nfeatures = analysis_params$n_variable_features)

# Plot variable features
pdf(file.path(output_dir, "plots", "04_variable_features.pdf"), width = 10, height = 6)
top10 <- head(VariableFeatures(combined), 10)
plot1 <- VariableFeaturePlot(combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()

# Scale data
cat("Scaling data...\n")
all_genes <- rownames(combined)
combined <- ScaleData(combined, features = all_genes)

# PCA
cat("\n=== STEP 5: Dimensionality Reduction ===\n")
cat("Running PCA...\n")
combined <- RunPCA(combined, features = VariableFeatures(object = combined))

# PCA plots
pdf(file.path(output_dir, "plots", "05_PCA_overview.pdf"), width = 12, height = 10)
print(DimPlot(combined, reduction = "pca"))
print(DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE))
dev.off()

# Elbow plot
pdf(file.path(output_dir, "plots", "06_elbow_plot.pdf"), width = 8, height = 6)
ElbowPlot(combined, ndims = 50)
dev.off()

# Clustering
cat("\n=== STEP 6: Clustering ===\n")
cat(paste0("Using ", analysis_params$n_pca_dims, " PCA dimensions\n"))
cat(paste0("Resolution: ", analysis_params$clustering_resolution, "\n"))

combined <- FindNeighbors(combined, dims = 1:analysis_params$n_pca_dims)
combined <- FindClusters(combined, resolution = analysis_params$clustering_resolution)

# UMAP
cat("Running UMAP...\n")
combined <- RunUMAP(combined, dims = 1:analysis_params$n_pca_dims)

# Save Seurat object
cat("\nSaving Seurat object...\n")
saveRDS(combined, file.path(output_dir, "data", "seurat_object.rds"))

# Visualization
cat("\n=== STEP 7: Visualization ===\n")

# UMAP by different groupings
pdf(file.path(output_dir, "plots", "07_UMAP_overview.pdf"), width = 16, height = 6)
p1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(combined, reduction = "umap", group.by = "condition")
p3 <- DimPlot(combined, reduction = "umap", group.by = "sample")
print(p1 + p2 + p3)
dev.off()

# UMAP split by condition
pdf(file.path(output_dir, "plots", "08_UMAP_by_condition.pdf"), width = 12, height = 5)
print(DimPlot(combined, reduction = "umap", split.by = "condition", label = TRUE))
dev.off()

# Find markers for all clusters
cat("\n=== STEP 8: Finding Cluster Markers ===\n")
cat("This may take several minutes...\n")
all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, file.path(output_dir, "results", "cluster_markers.csv"))

# Top markers per cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(top_markers, file.path(output_dir, "results", "top10_markers_per_cluster.csv"))

# Heatmap of top markers
pdf(file.path(output_dir, "plots", "09_marker_heatmap.pdf"), width = 14, height = 12)
top5 <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
DoHeatmap(combined, features = top5$gene) + NoLegend()
dev.off()

# Genes of interest analysis
if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
  cat("\n=== STEP 9: Genes of Interest Analysis ===\n")
  
  # Check which genes are present
  genes_present <- genes_of_interest[genes_of_interest %in% rownames(combined)]
  genes_missing <- genes_of_interest[!genes_of_interest %in% rownames(combined)]
  
  if (length(genes_missing) > 0) {
    cat("Warning: The following genes were not found in the dataset:\n")
    cat(paste("  ", genes_missing, collapse = "\n"))
    cat("\n")
  }
  
  if (length(genes_present) > 0) {
    cat(paste0("Analyzing ", length(genes_present), " genes of interest\n"))
    
    # Feature plots
    pdf(file.path(output_dir, "plots", "10_genes_of_interest_UMAP.pdf"), 
        width = 12, height = 4 * ceiling(length(genes_present) / 3))
    print(FeaturePlot(combined, features = genes_present, ncol = 3))
    dev.off()
    
    # Violin plots
    pdf(file.path(output_dir, "plots", "11_genes_of_interest_violin.pdf"), 
        width = 12, height = 4 * length(genes_present))
    for (gene in genes_present) {
      p <- VlnPlot(combined, features = gene, group.by = "condition", pt.size = 0.1)
      print(p)
    }
    dev.off()
    
    # Dot plot
    pdf(file.path(output_dir, "plots", "12_genes_of_interest_dotplot.pdf"), 
        width = 10, height = 8)
    print(DotPlot(combined, features = genes_present, group.by = "seurat_clusters") + 
          RotatedAxis())
    dev.off()
  }
}

# Summary report
cat("\n=== Creating Summary Report ===\n")
sink(file.path(output_dir, "analysis_summary.txt"))

cat("Single-Cell RNA-seq Analysis Summary\n")
cat("=====================================\n")
cat(paste0("Dataset: ", dataset_id, "\n\n"))

cat("Sample Information:\n")
for (sample_name in sample_names) {
  cat(paste0("  ", sample_name, ": "))
  sample_cells <- sum(combined$sample == sample_name)
  cat(paste0(sample_cells, " cells\n"))
}

cat(paste0("\nTotal cells after QC: ", ncol(combined), "\n"))
cat(paste0("Total genes: ", nrow(combined), "\n"))
cat(paste0("Variable features: ", length(VariableFeatures(combined)), "\n"))
cat(paste0("Number of clusters: ", length(unique(combined$seurat_clusters)), "\n"))

cat("\nAnalysis Parameters:\n")
cat(paste0("  Min features: ", qc_filters$min_features, "\n"))
cat(paste0("  Max features: ", qc_filters$max_features, "\n"))
cat(paste0("  Max MT%: ", qc_filters$max_mt_percent, "\n"))
cat(paste0("  Variable features: ", analysis_params$n_variable_features, "\n"))
cat(paste0("  PCA dimensions: ", analysis_params$n_pca_dims, "\n"))
cat(paste0("  Clustering resolution: ", analysis_params$clustering_resolution, "\n"))

cat("\nOutput Files:\n")
cat("  - seurat_object.rds: Complete Seurat object\n")
cat("  - cluster_markers.csv: Marker genes for all clusters\n")
cat("  - top10_markers_per_cluster.csv: Top 10 markers per cluster\n")
cat("  - Various visualization plots in plots/ directory\n")

sink()

cat("\n======================================================\n")
cat("Analysis complete!\n")
cat(paste0("Results saved to: ", normalizePath(output_dir), "\n"))
cat("======================================================\n")
