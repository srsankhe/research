#!/usr/bin/env Rscript
# Install required R packages for DESeq2 analysis

cat("Installing required R packages...\n\n")

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "http://cran.rstudio.com/")
}

# Install DESeq2 and dependencies
cat("Installing DESeq2...\n")
BiocManager::install("DESeq2", ask = FALSE, update = FALSE)

# Install other required packages
required_packages <- c("ggplot2", "pheatmap", "RColorBrewer")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, "...\n"))
    install.packages(pkg, repos = "http://cran.rstudio.com/")
  }
}

cat("\nAll packages installed successfully!\n")
