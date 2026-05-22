#!/usr/bin/env Rscript
# Create annotated plots highlighting MME gene
# Dataset: GSE282344

library(DESeq2)
library(ggplot2)
library(ggrepel)  # For gene labels

# Set working directory
setwd("/Users/ssankhe/Downloads/GSE282344_RAW")

cat("Loading results...\n")

# Read the DESeq2 results
results <- read.csv("deseq2_output/treatment_vs_control_results.csv", row.names = 1)

# Create output directory
dir.create("deseq2_output/annotated_plots", showWarnings = FALSE)

cat("Creating annotated plots with MME highlighted...\n\n")

# ================================================================
# 1. VOLCANO PLOT with MME annotation
# ================================================================
cat("1. Creating annotated volcano plot...\n")

volcano_data <- results
volcano_data$gene <- rownames(volcano_data)

# Define significance
volcano_data$significant <- ifelse(
  volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1,
  "Significant",
  "Not significant"
)

# Special highlighting for MME
volcano_data$highlight <- "Other"
volcano_data$highlight[volcano_data$gene == "MME"] <- "MME"
volcano_data$highlight[volcano_data$significant == "Significant"] <- "Significant"

# Get top genes for labeling
top_genes <- rbind(
  head(volcano_data[order(volcano_data$padj), ], 10),  # Top 10 by p-value
  volcano_data[volcano_data$gene == "MME", ]  # Always include MME
)
top_genes <- unique(top_genes)

# Create volcano plot
p1 <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(volcano_data, highlight == "Other"),
             aes(color = highlight), alpha = 0.3, size = 1) +
  geom_point(data = subset(volcano_data, highlight == "Significant"),
             aes(color = highlight), alpha = 0.6, size = 2) +
  geom_point(data = subset(volcano_data, highlight == "MME"),
             aes(color = highlight), size = 4, shape = 17) +  # Triangle for MME
  scale_color_manual(values = c("Other" = "gray70",
                                 "Significant" = "red",
                                 "MME" = "blue"),
                     name = "Gene Status") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", alpha = 0.5) +
  geom_label_repel(data = top_genes,
                   aes(label = gene),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.3,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   force = 2) +
  theme_minimal() +
  labs(title = "Volcano Plot: Treatment vs Control (GSE282344)",
       subtitle = "MME highlighted in blue (log2FC = 1.04, padj = 0.925 - NOT significant)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "blue"))

ggsave("deseq2_output/annotated_plots/volcano_plot_MME_annotated.png",
       p1, width = 12, height = 8, dpi = 300)

cat("  ✓ Saved: volcano_plot_MME_annotated.png\n")

# ================================================================
# 2. MA PLOT with MME annotation
# ================================================================
cat("2. Creating annotated MA plot...\n")

ma_data <- results
ma_data$gene <- rownames(ma_data)
ma_data$significant <- ifelse(
  ma_data$padj < 0.05 & abs(ma_data$log2FoldChange) > 1,
  "Significant",
  "Not significant"
)

# Highlight MME
ma_data$highlight <- "Other"
ma_data$highlight[ma_data$gene == "MME"] <- "MME"
ma_data$highlight[ma_data$significant == "Significant"] <- "Significant"

# Get genes to label
label_genes <- rbind(
  head(ma_data[order(ma_data$padj), ], 15),
  ma_data[ma_data$gene == "MME", ]
)
label_genes <- unique(label_genes)

p2 <- ggplot(ma_data, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(data = subset(ma_data, highlight == "Other"),
             aes(color = highlight), alpha = 0.3, size = 1) +
  geom_point(data = subset(ma_data, highlight == "Significant"),
             aes(color = highlight), alpha = 0.6, size = 2) +
  geom_point(data = subset(ma_data, highlight == "MME"),
             aes(color = highlight), size = 4, shape = 17) +
  scale_color_manual(values = c("Other" = "gray70",
                                 "Significant" = "red",
                                 "MME" = "blue"),
                     name = "Gene Status") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40", alpha = 0.5) +
  geom_label_repel(data = label_genes,
                   aes(label = gene),
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.3,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   force = 2) +
  theme_minimal() +
  labs(title = "MA Plot: Treatment vs Control (GSE282344)",
       subtitle = "MME: baseMean = 1487, log2FC = 1.04 (blue triangle)",
       x = "Log10 Mean Expression",
       y = "Log2 Fold Change") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "blue"))

ggsave("deseq2_output/annotated_plots/MA_plot_MME_annotated.png",
       p2, width = 12, height = 8, dpi = 300)

cat("  ✓ Saved: MA_plot_MME_annotated.png\n")

# ================================================================
# 3. MME-FOCUSED PLOT
# ================================================================
cat("3. Creating MME-focused comparison plot...\n")

# Get MME data
mme_data <- ma_data[ma_data$gene == "MME", ]

# Create a focused plot showing MME's position
p3 <- ggplot() +
  geom_point(data = ma_data,
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "gray80", alpha = 0.3, size = 1) +
  geom_point(data = subset(ma_data, significant == "Significant"),
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "red", alpha = 0.5, size = 2) +
  geom_point(data = mme_data,
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "blue", size = 6, shape = 17) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  annotate("text", x = mme_data$log2FoldChange + 1.5,
           y = -log10(mme_data$padj) + 0.5,
           label = paste0("MME\nlog2FC: ", round(mme_data$log2FoldChange, 2),
                         "\npadj: ", round(mme_data$padj, 3),
                         "\nNOT significant"),
           size = 4, color = "blue", fontface = "bold",
           hjust = 0) +
  annotate("rect",
           xmin = mme_data$log2FoldChange - 0.2,
           xmax = mme_data$log2FoldChange + 0.2,
           ymin = -log10(mme_data$padj) - 0.02,
           ymax = -log10(mme_data$padj) + 0.02,
           alpha = 0, color = "blue", size = 1.5) +
  theme_minimal() +
  labs(title = "MME Gene Location in Volcano Plot",
       subtitle = "Blue triangle = MME; Red dots = Significant genes (padj < 0.05, |log2FC| > 1)",
       x = "Log2 Fold Change (Treatment vs Control)",
       y = "-Log10 Adjusted P-value") +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10))

ggsave("deseq2_output/annotated_plots/MME_focused_plot.png",
       p3, width = 10, height = 8, dpi = 300)

cat("  ✓ Saved: MME_focused_plot.png\n")

# ================================================================
# 4. CREATE MME SUMMARY TABLE
# ================================================================
cat("4. Creating MME summary information...\n")

sink("deseq2_output/annotated_plots/MME_summary.txt")

cat("====================================================================\n")
cat("MME Gene Analysis Summary - GSE282344 Dataset\n")
cat("====================================================================\n\n")

cat("GENE INFORMATION:\n")
cat("-----------------\n")
cat("Gene Name: MME (Membrane Metalloendopeptidase, also known as CD10, Neprilysin)\n")
cat("Function: Cell surface peptidase that degrades various peptides including\n")
cat("          substance P, enkephalins, bradykinin, and amyloid beta.\n")
cat("          Important in neuropeptide regulation and Alzheimer's disease.\n\n")

cat("EXPRESSION RESULTS:\n")
cat("-------------------\n")
cat(sprintf("Base Mean Expression: %.2f\n", mme_data$baseMean))
cat(sprintf("Log2 Fold Change: %.2f\n", mme_data$log2FoldChange))
cat(sprintf("Actual Fold Change: %.2f-fold increase in treatment\n", 2^mme_data$log2FoldChange))
cat(sprintf("Standard Error: %.2f\n", mme_data$lfcSE))
cat(sprintf("Test Statistic: %.2f\n", mme_data$stat))
cat(sprintf("P-value (raw): %.4f\n", mme_data$pvalue))
cat(sprintf("P-value (adjusted): %.4f\n", mme_data$padj))
cat("\n")

cat("SIGNIFICANCE ASSESSMENT:\n")
cat("------------------------\n")
if (mme_data$padj < 0.05 & abs(mme_data$log2FoldChange) > 1) {
  cat("STATUS: *** SIGNIFICANT ***\n")
} else if (mme_data$padj < 0.05) {
  cat("STATUS: Significant by p-value, but fold change < 2\n")
} else if (abs(mme_data$log2FoldChange) > 1) {
  cat("STATUS: ** NOT SIGNIFICANT ** (high p-value despite fold change > 2)\n")
} else {
  cat("STATUS: NOT SIGNIFICANT\n")
}

cat(sprintf("\nCriteria for significance:\n"))
cat(sprintf("  - Adjusted p-value < 0.05: %s (actual: %.4f)\n",
           ifelse(mme_data$padj < 0.05, "PASS", "FAIL"), mme_data$padj))
cat(sprintf("  - |log2FC| > 1: %s (actual: %.2f)\n",
           ifelse(abs(mme_data$log2FoldChange) > 1, "PASS", "FAIL"),
           mme_data$log2FoldChange))

cat("\n\nINTERPRETATION:\n")
cat("---------------\n")
cat("MME shows a 2-fold upregulation in treatment samples compared to controls.\n")
cat("However, this change is NOT statistically significant after multiple testing\n")
cat("correction (adjusted p-value = 0.925).\n\n")

cat("Possible reasons for non-significance:\n")
cat("1. High variability between replicates\n")
cat("2. Low statistical power (only 2 replicates per condition)\n")
cat("3. The observed change may be due to random variation\n\n")

cat("RECOMMENDATIONS:\n")
cat("----------------\n")
cat("1. If MME is a gene of interest, consider:\n")
cat("   - Validating with qRT-PCR in more samples\n")
cat("   - Checking expression in individual cells (single-cell analysis)\n")
cat("   - Looking at MME in related datasets\n\n")

cat("2. For publication:\n")
cat("   - Report as 'trending toward upregulation but not significant'\n")
cat("   - DO NOT report as a differentially expressed gene\n")
cat("   - Can mention in discussion as potential candidate for follow-up\n\n")

cat("BIOLOGICAL CONTEXT:\n")
cat("-------------------\n")
cat("MME/CD10/Neprilysin is particularly relevant in:\n")
cat("- Cancer biology (tumor marker in various cancers)\n")
cat("- Neurodegeneration (cleaves amyloid-beta in Alzheimer's)\n")
cat("- Immune regulation (expressed on B-cells and neutrophils)\n")
cat("- Pain modulation (degrades enkephalins)\n\n")

# Compare to other genes
cat("COMPARISON TO SIGNIFICANT GENES:\n")
cat("--------------------------------\n")
sig_genes <- ma_data[ma_data$significant == "Significant", ]
cat(sprintf("Number of significant genes: %d\n", nrow(sig_genes)))
if (nrow(sig_genes) > 0) {
  cat(sprintf("Range of significant log2FC: %.2f to %.2f\n",
             min(sig_genes$log2FoldChange), max(sig_genes$log2FoldChange)))
  cat(sprintf("Range of significant padj: %.2e to %.2e\n",
             min(sig_genes$padj, na.rm = TRUE), max(sig_genes$padj, na.rm = TRUE)))
  cat(sprintf("\nMME log2FC (%.2f) is %s than median significant gene (%.2f)\n",
             mme_data$log2FoldChange,
             ifelse(mme_data$log2FoldChange > median(sig_genes$log2FoldChange), "higher", "lower"),
             median(sig_genes$log2FoldChange)))
}

cat("\n====================================================================\n")

sink()

cat("  ✓ Saved: MME_summary.txt\n")

# ================================================================
# FINAL SUMMARY
# ================================================================
cat("\n====================================================================\n")
cat("ANNOTATED PLOTS CREATED SUCCESSFULLY!\n")
cat("====================================================================\n\n")

cat("Output files saved to: deseq2_output/annotated_plots/\n\n")

cat("Files created:\n")
cat("  1. volcano_plot_MME_annotated.png - Full volcano plot with MME highlighted\n")
cat("  2. MA_plot_MME_annotated.png - MA plot with MME highlighted\n")
cat("  3. MME_focused_plot.png - Focused view of MME position\n")
cat("  4. MME_summary.txt - Detailed MME analysis report\n\n")

cat("MME QUICK SUMMARY:\n")
cat(sprintf("  Expression: %.0f (moderate)\n", mme_data$baseMean))
cat(sprintf("  Fold Change: %.2f-fold increase\n", 2^mme_data$log2FoldChange))
cat(sprintf("  Significance: %s (padj = %.3f)\n",
           ifelse(mme_data$padj < 0.05, "YES", "NO"), mme_data$padj))

cat("\n====================================================================\n")
