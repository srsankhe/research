#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)

# Read one sample
genes <- read.table(gzfile("GSM8641663_NR_JK_047_genes.tsv.gz"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(gzfile("GSM8641663_NR_JK_047_barcodes.tsv.gz"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
mtx <- readMM(gzfile("GSM8641663_NR_JK_047_matrix.mtx.gz"))

rownames(mtx) <- genes$V2
colnames(mtx) <- barcodes$V1

cat("Matrix dimensions:", nrow(mtx), "x", ncol(mtx), "\n")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = mtx, project = "test", min.cells = 3, min.features = 200)

cat("Seurat object cells:", ncol(seurat_obj), "\n")
cat("Cell names (first 5):\n")
print(head(colnames(seurat_obj)))

cat("\nMetadata before adding:\n")
print(head(seurat_obj@meta.data))

# Try adding metadata
cat("\nTrying to add metadata...\n")
seurat_obj$test_meta <- "test_value"

cat("\nMetadata after adding:\n")
print(head(seurat_obj@meta.data))

cat("\nSuccess!\n")
