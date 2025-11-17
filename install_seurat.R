#!/usr/bin/env Rscript
# Install Seurat and dependencies for single-cell analysis
# Generic installation script

cat("Installing Seurat and dependencies...\n")
cat("This may take 15-20 minutes. Please be patient.\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", quiet = FALSE)
}

# Required packages
required_packages <- c(
  # Core Seurat and single-cell packages
  "Seurat",
  "Matrix",
  
  # Data manipulation
  "dplyr",
  
  # Plotting
  "ggplot2",
  "patchwork",
  "cowplot",
  "RColorBrewer",
  "viridis",
  "ggrepel",
  
  # Additional utilities
  "jsonlite"
)

cat("Installing/checking packages...\n")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, "...\n"))
    install.packages(pkg, quiet = FALSE)
  } else {
    cat(paste0(pkg, " already installed.\n"))
  }
}

cat("\n=== Verifying Installations ===\n")
all_ok <- TRUE
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("✓ ", pkg, "\n"))
  } else {
    cat(paste0("✗ ", pkg, " - FAILED\n"))
    all_ok <- FALSE
  }
}

if (all_ok) {
  cat("\n=== Installation Complete ===\n")
  cat("Testing Seurat installation...\n")
  library(Seurat)
  cat(paste0("Seurat version: ", packageVersion("Seurat"), "\n"))
  cat("\nAll packages installed successfully!\n")
} else {
  cat("\n⚠ Some packages failed to install. Please check the error messages above.\n")
  quit(status = 1)
}
