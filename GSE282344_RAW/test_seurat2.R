#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)

# Read one sample
genes <- read.table(gzfile("GSM8641663_NR_JK_047_genes.tsv.gz"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(gzfile("GSM8641663_NR_JK_047_barcodes.tsv.gz"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
mtx <- readMM(gzfile("GSM8641663_NR_JK_047_matrix.mtx.gz"))

rownames(mtx) <- genes$V2
colnames(mtx) <- barcodes$V1

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = mtx, project = "test", min.cells = 3, min.features = 200)

cat("Trying PercentageFeatureSet...\n")
seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

cat("\nMetadata after PercentageFeatureSet:\n")
print(head(seurat_obj@meta.data))

cat("\nTrying to add sample metadata...\n")
seurat_obj$sample <- "GSM8641663_NR_JK_047"

cat("\nMetadata after adding sample:\n")
print(head(seurat_obj@meta.data))

cat("\nSuccess!\n")
