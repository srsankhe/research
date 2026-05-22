# Single-Cell RNA-seq Analysis Guide for Beginners

## Complete Guide: From Raw Data to Differential Expression Results

This guide will walk you through analyzing single-cell RNA-seq data from the GSE208526 dataset, comparing Ischemia-Reperfusion Injury (IRI) samples with sham controls in rat kidney tissue.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Understanding Your Data](#understanding-your-data)
3. [Setting Up Your Environment](#setting-up-your-environment)
4. [Step 1: Converting Raw Data to CSV](#step-1-converting-raw-data-to-csv)
5. [Step 2: Running DESeq2 Analysis](#step-2-running-deseq2-analysis)
6. [Understanding the Results](#understanding-the-results)
7. [Troubleshooting](#troubleshooting)
8. [Next Steps](#next-steps)

---

## Prerequisites

### What You Need

- **A Mac, Linux, or Windows computer** with at least 8GB RAM
- **Basic familiarity with the command line/terminal** (don't worry, we'll guide you!)
- **Internet connection** for downloading software packages
- **The GSE208526_RAW data** (which you already have!)

### No Prior Experience Needed In:
- R programming
- Python programming
- Bioinformatics
- Statistics (though it helps!)

---

## Understanding Your Data

### What is Single-Cell RNA-seq?

Single-cell RNA-sequencing (scRNA-seq) measures gene expression in individual cells. Traditional RNA-seq averages gene expression across millions of cells, but scRNA-seq lets us see what's happening in each cell individually.

### Your Dataset: GSE208526

This dataset contains rat kidney samples from:
- **Sham controls** (normal/healthy tissue)
- **IRI samples** (Ischemia-Reperfusion Injury - tissue that was deprived of blood flow and then restored)

Each sample is further divided by cell type:
- **CD45+** cells (immune cells)
- **CD45-** cells (non-immune cells, likely epithelial kidney cells)

### Files in Your Folder

**Directory Structure:**
```
archives/                            (All compressed data files)
├── GSM6347353_sham.tar.gz          (Original compressed archive)
├── GSM6347354_sham_CD45-.tar.gz
├── GSM6347355_IRI.tar.gz
├── GSM6347356_IRI_CD45-.tar.gz
├── sham/                            (Extracted data for sham sample)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── sham_CD45-/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── IRI/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── IRI_CD45-/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

Each sample's data files are in 10X Genomics format:
- **matrix.mtx.gz** - Gene expression counts (sparse matrix format)
- **features.tsv.gz** - Gene information (gene IDs and names)
- **barcodes.tsv.gz** - Cell identifiers (unique barcode for each cell)

**Note:** All compressed files (.gz and .tar.gz) have been organized into the `archives/` folder to keep the workspace clean. The conversion script automatically reads from this location.

---

## Setting Up Your Environment

### Step 1: Install Python 3

**Check if you already have Python:**
```bash
python3 --version
```

If you see a version number (like `Python 3.9.6`), you're good! If not:

**On Mac:**
```bash
# Install Homebrew first (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python
brew install python3
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

**On Windows:**
Download from [python.org](https://www.python.org/downloads/) and run the installer.

### Step 2: Install R

**Check if you have R:**
```bash
R --version
```

If not installed:

**On Mac:**
```bash
brew install r
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install r-base
```

**On Windows:**
Download from [r-project.org](https://www.r-project.org/) and run the installer.

### Step 3: Navigate to Your Data Folder

Open Terminal (Mac/Linux) or Command Prompt (Windows) and navigate to your data:

```bash
cd /path/to/GSE208526_RAW
```

**Tip:** You can type `cd` followed by a space, then drag-and-drop the folder into the terminal window to auto-fill the path!

---

## Step 1: Converting Raw Data to CSV

The raw data is in 10X Genomics format (sparse matrix), which is space-efficient but hard to work with. We'll convert it to CSV (comma-separated values) files that are easier to read and analyze.

### 1.1: Data Organization

All compressed data files (`.gz` and `.tar.gz`) have been organized into the `archives/` directory:
- Original `.tar.gz` archives
- Extracted sample folders with the actual data files

The conversion script automatically reads from `archives/` - no manual extraction needed!

### 1.2: Create a Python Virtual Environment

This keeps your project dependencies isolated:

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate  # Mac/Linux
# OR
venv\Scripts\activate     # Windows

# Your prompt should now show (venv) at the beginning
```

### 1.3: Install Required Python Packages

```bash
pip install pandas scipy numpy
```

**What these do:**
- **pandas** - Data manipulation (like Excel for Python)
- **scipy** - Scientific computing (reads sparse matrices)
- **numpy** - Numerical operations (math with arrays)

### 1.4: Run the Conversion Script

```bash
python3 convert_to_csv.py
```

**What's happening:**
1. Script reads data from `archives/` subdirectories
2. Reads the sparse matrix files for each sample
3. Converts them to dense CSV format
4. Saves everything to `csv_output/` folder

**Expected output:**
```
Processing sample: sham
  Reading features...
  ✓ Saved: .../sham_features.csv
  Reading barcodes...
  ✓ Saved: .../sham_barcodes.csv
  Reading matrix (this may take a while)...
  ✓ Saved: .../sham_matrix.csv
  Matrix shape: 25339 genes × 8785 cells
```

This will take **5-10 minutes** depending on your computer speed.

### 1.5: Check Your Output

```bash
ls -lh csv_output/
```

You should see 12 files (3 per sample):
- `*_matrix.csv` - Expression data
- `*_features.csv` - Gene information
- `*_barcodes.csv` - Cell identifiers

---

## Step 2: Running DESeq2 Analysis

DESeq2 is a popular R package for differential expression analysis. We'll use it to find genes that are significantly different between IRI and sham samples.

### 2.1: Install R Packages

First, install the required packages:

```bash
Rscript install_packages.R
```

**What's being installed:**
- **DESeq2** - Differential expression analysis
- **ggplot2** - Publication-quality plots
- **pheatmap** - Heatmaps
- **RColorBrewer** - Color palettes

This will take **10-15 minutes**. You'll see lots of output - this is normal!

### 2.2: Run DESeq2 Analysis

```bash
Rscript deseq2_analysis.R
```

**What's happening:**

1. **Loading data** - Reads all CSV matrices
2. **Creating pseudo-bulk** - Combines single cells into bulk samples (DESeq2 works on bulk RNA-seq)
3. **Filtering** - Removes genes with zero counts
4. **Normalization** - Accounts for different library sizes
5. **Statistical testing** - Identifies differentially expressed genes
6. **Generating plots** - Creates visualizations

**Expected output:**
```
Loading data...
  Loading sham...
  Loading sham_CD45-...
  Loading IRI...
  Loading IRI_CD45-...

Creating pseudo-bulk count matrix...
  Matrix dimensions: 25339 genes x 4 samples
  After filtering zero counts: 18724 genes

=== Running DESeq2 Analysis ===

Significant genes (padj < 0.05, |log2FC| > 1): 356
  Upregulated in IRI: 239
  Downregulated in IRI: 117
```

This will take **3-5 minutes**.

### 2.3: Check Your Results

```bash
ls -lh deseq2_output/
```

You should see:
- CSV files with results
- PNG images with plots
- A text summary file

---

## Understanding the Results

### Result Files Explained

#### 1. `analysis_summary.txt`
A human-readable summary of the analysis.

**View it:**
```bash
cat deseq2_output/analysis_summary.txt
```

**Key information:**
- Total genes analyzed
- Number of significant genes
- Top upregulated/downregulated genes

#### 2. `IRI_vs_sham_results.csv`
**Complete results for all genes.**

**Columns explained:**
- `baseMean` - Average expression across all samples (higher = more expressed)
- `log2FoldChange` - How much the gene changes in IRI vs sham
  - Positive = upregulated in IRI
  - Negative = downregulated in IRI
  - log2FC of 1 = 2-fold increase
  - log2FC of -1 = 2-fold decrease
- `lfcSE` - Standard error of the fold change
- `stat` - Test statistic
- `pvalue` - Raw p-value (probability this happened by chance)
- `padj` - Adjusted p-value (accounts for testing thousands of genes)

**Opening it:**
- **In Excel/Numbers:** Double-click the file
- **In R:**
  ```r
  results <- read.csv("deseq2_output/IRI_vs_sham_results.csv", row.names=1)
  head(results)
  ```
- **In Python:**
  ```python
  import pandas as pd
  results = pd.read_csv("deseq2_output/IRI_vs_sham_results.csv", index_col=0)
  results.head()
  ```

#### 3. `IRI_vs_sham_significant_genes.csv`
**Only genes that pass significance thresholds:**
- `padj < 0.05` (less than 5% chance of false positive)
- `|log2FoldChange| > 1` (at least 2-fold change)

These are your **high-confidence candidates** for further investigation!

#### 4. `pseudobulk_counts.csv`
The aggregated count matrix (genes × samples) used for analysis.

### Visualization Files

#### 1. `volcano_plot.png`
**Shows:** All genes plotted by fold change (x-axis) vs significance (y-axis)

**How to read:**
- **X-axis:** log2 fold change (left = downregulated, right = upregulated)
- **Y-axis:** -log10(adjusted p-value) (higher = more significant)
- **Red dots:** Significant genes
- **Gray dots:** Non-significant genes

**Look for:** Genes in the upper-left and upper-right corners (big changes + significant)

#### 2. `MA_plot.png`
**Shows:** Expression level vs fold change

**How to read:**
- **X-axis:** Average expression (baseMean)
- **Y-axis:** log2 fold change
- **Blue dots:** Significant genes
- **Gray dots:** Non-significant genes

**Look for:** Blue dots far from zero on the y-axis

#### 3. `PCA_plot.png`
**Shows:** How similar/different your samples are

**How to read:**
- Each dot = one sample
- Samples close together = similar gene expression
- Samples far apart = different gene expression

**Good result:** IRI samples cluster together, sham samples cluster together, clearly separated

#### 4. `heatmap_top50_genes.png`
**Shows:** Top 50 most significantly changed genes across all samples

**How to read:**
- **Rows:** Genes
- **Columns:** Samples
- **Red:** High expression
- **Blue:** Low expression

**Look for:** Patterns where genes are red in IRI samples and blue in sham (or vice versa)

#### 5. `dispersion_plot.png`
**Shows:** Quality control for the statistical model

**How to read:**
- Black dots = individual genes
- Red line = fitted relationship
- Blue dots = genes with high variability

**Good result:** Points follow the red curve smoothly

---

## Interpreting Your Results

### What Did We Find?

From the analysis, we identified **356 significantly changed genes** between IRI and sham:

**Top Upregulated in IRI (Injury):**
- **Ccl12, Ccl7, Ccl2** - Chemokines (immune cell recruitment)
- **Cxcl1, Cxcl10, Cxcl11** - More chemokines (inflammation)
- **Timp1** - Tissue inhibitor of metalloproteinases (wound repair)

**Biology:** These results make sense! After ischemia-reperfusion injury, the kidney:
1. Recruits immune cells (chemokines)
2. Activates inflammatory pathways
3. Initiates tissue repair

**Top Downregulated in IRI:**
- **Egf** - Epidermal growth factor (cell growth)
- **Igf1** - Insulin-like growth factor (tissue maintenance)
- **Umod** - Uromodulin (kidney-specific protein)

**Biology:** Loss of normal kidney function markers suggests tissue damage.

### What is "Significant"?

A gene is considered significant if:
1. **padj < 0.05** - Less than 5% chance this is a false positive
2. **|log2FoldChange| > 1** - At least 2-fold change (doubled or halved)

**Example:**
- **Ccl12**: log2FC = 12.6, padj = 4e-11
  - This means ~6,000-fold increase (2^12.6)
  - Basically zero chance this is random
  - **VERY significant finding!**

### Searching for Specific Genes

Want to know about a specific gene? Use `grep`:

```bash
# Find the gene Mme (for example)
grep "Mme" deseq2_output/IRI_vs_sham_results.csv
```

Or open the CSV in Excel and use Ctrl+F (Cmd+F on Mac) to search.

---

## Troubleshooting

### Problem: "ModuleNotFoundError: No module named 'pandas'"

**Solution:**
```bash
# Make sure virtual environment is activated
source venv/bin/activate  # You should see (venv) in your prompt

# Install packages
pip install pandas scipy numpy
```

### Problem: "Error in library(DESeq2) : there is no package called 'DESeq2'"

**Solution:**
```bash
# Run the package installation script
Rscript install_packages.R

# If that fails, try installing manually in R:
R
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("DESeq2")
> quit()
```

### Problem: "duplicate rownames" error

**Solution:** This is handled in the updated script. Make sure you're using the latest `deseq2_analysis.R` file.

### Problem: Out of memory errors

**Solution:** Your computer may not have enough RAM. Try:
1. Close other applications
2. Work with fewer samples at a time
3. Use a computer with more RAM (8GB minimum, 16GB recommended)

### Problem: Plots won't open

**Solution:** The plots are PNG files. Try:
```bash
# On Mac
open deseq2_output/volcano_plot.png

# On Linux
xdg-open deseq2_output/volcano_plot.png

# On Windows
start deseq2_output/volcano_plot.png
```

Or navigate to the folder in Finder/File Explorer and double-click.

---

## Next Steps

### 1. Pathway Analysis

Want to know what biological processes these genes are involved in?

**Try:**
- **DAVID** (https://david.ncifcrf.gov/) - Free web tool
- **Enrichr** (https://maayanlab.cloud/Enrichr/) - Easy interface
- **Gene Ontology** (http://geneontology.org/) - Comprehensive database

**How:**
1. Copy the list of significant gene names from `IRI_vs_sham_significant_genes.csv`
2. Paste into one of the tools above
3. Get pathways, GO terms, and more!

### 2. Network Analysis

Visualize how these genes interact:

**Try:**
- **STRING** (https://string-db.org/) - Protein-protein interactions
- **Cytoscape** (https://cytoscape.org/) - Network visualization software

### 3. Single-Cell Analysis

We converted to pseudo-bulk for DESeq2, but you could also:
- Use **Seurat** (R package) for single-cell specific analysis
- Cluster cells by type
- Find cell-type specific markers
- Trajectory analysis (how cells change over time)

### 4. Validation

Key genes should be validated:
- **qRT-PCR** - Measure specific genes in lab
- **Immunohistochemistry** - Visualize protein localization
- **Western blot** - Confirm protein levels

### 5. Literature Search

Research your top genes:
- **PubMed** (https://pubmed.ncbi.nlm.nih.gov/)
- **Google Scholar** (https://scholar.google.com/)
- Search: "[gene name] kidney injury" or "[gene name] ischemia reperfusion"

---

## Key Concepts Explained

### Pseudo-bulk Aggregation

**What:** Summing counts across all cells in a sample
**Why:** DESeq2 was designed for bulk RNA-seq, so we convert single-cell to bulk
**How:** For each gene, add up expression across all cells in that sample

**Note:** For true single-cell analysis, use specialized tools like Seurat or Scanpy

### Log2 Fold Change

**Why log2?**
- Makes fold changes symmetric
- 2-fold increase = log2FC of +1
- 2-fold decrease = log2FC of -1
- 4-fold increase = log2FC of +2
- 4-fold decrease = log2FC of -2

**Conversion:**
```
Actual fold change = 2^(log2FoldChange)

Example:
log2FC = 3 → 2^3 = 8-fold increase
log2FC = -2 → 2^-2 = 0.25 (4-fold decrease)
```

### P-value vs Adjusted P-value

**P-value:** Probability that the result happened by chance
- p = 0.05 means 5% chance it's random
- p = 0.001 means 0.1% chance it's random

**Adjusted P-value (padj):**
- Accounts for testing thousands of genes simultaneously
- Without adjustment, testing 20,000 genes at p < 0.05 would give ~1,000 false positives!
- Methods: Benjamini-Hochberg (used here), Bonferroni
- **Always use padj, not p, for significance**

### Base Mean

Average expression across all samples
- **High baseMean (>1000):** Highly expressed gene
- **Medium baseMean (100-1000):** Moderately expressed
- **Low baseMean (<100):** Lowly expressed

**Note:** Changes in highly expressed genes are often more reliable

---

## File Structure Reference

```
GSE208526_RAW/
├── readme.md                          (This file!)
├── convert_to_csv.py                  (Python script for conversion)
├── consolidate_data.py                (Data consolidation script)
├── GSE208526_sample_metadata.csv      (Sample metadata with file paths)
├── install_packages.R                 (R package installer)
├── deseq2_analysis.R                  (Main analysis script)
│
├── archives/                          (All compressed data files)
│   ├── GSM6347353_sham.tar.gz        (Original compressed archives)
│   ├── GSM6347354_sham_CD45-.tar.gz
│   ├── GSM6347355_IRI.tar.gz
│   ├── GSM6347356_IRI_CD45-.tar.gz
│   ├── sham/                          (Extracted sample data)
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   ├── sham_CD45-/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   ├── IRI/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── IRI_CD45-/
│       ├── matrix.mtx.gz
│       ├── features.tsv.gz
│       └── barcodes.tsv.gz
│
├── csv_output/                        (Converted CSV files)
│   ├── sham_matrix.csv
│   ├── sham_features.csv
│   ├── sham_barcodes.csv
│   ├── sham_CD45-_matrix.csv
│   ├── sham_CD45-_features.csv
│   ├── sham_CD45-_barcodes.csv
│   ├── IRI_matrix.csv
│   ├── IRI_features.csv
│   ├── IRI_barcodes.csv
│   ├── IRI_CD45-_matrix.csv
│   ├── IRI_CD45-_features.csv
│   └── IRI_CD45-_barcodes.csv
│
├── deseq2_output/                     (Analysis results)
│   ├── analysis_summary.txt           (Human-readable summary)
│   ├── IRI_vs_sham_results.csv        (All genes)
│   ├── IRI_vs_sham_significant_genes.csv (Significant genes only)
│   ├── pseudobulk_counts.csv          (Aggregated count matrix)
│   ├── sample_metadata.csv            (Sample information)
│   ├── volcano_plot.png               (Volcano plot)
│   ├── MA_plot.png                    (MA plot)
│   ├── PCA_plot.png                   (PCA plot)
│   ├── heatmap_top50_genes.png        (Heatmap)
│   └── dispersion_plot.png            (Dispersion estimates)
│
└── venv/                              (Python virtual environment)
```

---

## Quick Command Reference

### Activate Python environment
```bash
source venv/bin/activate
```

### Convert data to CSV
```bash
python3 convert_to_csv.py
```

### Install R packages
```bash
Rscript install_packages.R
```

### Run DESeq2 analysis
```bash
Rscript deseq2_analysis.R
```

### View summary
```bash
cat deseq2_output/analysis_summary.txt
```

### Open a plot
```bash
open deseq2_output/volcano_plot.png     # Mac
xdg-open deseq2_output/volcano_plot.png # Linux
start deseq2_output/volcano_plot.png    # Windows
```

### Search for a gene
```bash
grep "GeneNameHere" deseq2_output/IRI_vs_sham_results.csv
```

---

## Getting Help

### Online Resources

1. **DESeq2 Manual:** https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
2. **Bioconductor Support:** https://support.bioconductor.org/
3. **Seurat Tutorials:** https://satijalab.org/seurat/
4. **Python for Data Analysis:** https://pandas.pydata.org/docs/
5. **R for Data Science:** https://r4ds.had.co.nz/

### Key Terms to Search

- "DESeq2 tutorial"
- "differential expression analysis RNA-seq"
- "single cell RNA-seq analysis"
- "10X Genomics data analysis"
- "interpreting volcano plot"

### Communities

- **Biostars:** https://www.biostars.org/ (Q&A for bioinformatics)
- **Stack Overflow:** https://stackoverflow.com/ (Programming questions)
- **Reddit r/bioinformatics:** https://reddit.com/r/bioinformatics/

---

## Citation

If you use these results in a publication, please cite:

**Dataset:**
- GEO Accession: GSE208526
- Citation: [Include the paper citation when published]

**Software:**
- **DESeq2:** Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550.
- **pandas:** McKinney, W. (2010). Data structures for statistical computing in python.
- **R Core Team** (2024). R: A language and environment for statistical computing.

---

## License & Acknowledgments

This analysis pipeline was created for educational purposes.

**Data source:** NCBI GEO (GSE208526)

**Author:** Analysis pipeline created on 2025-11-16

**Questions?** Check the troubleshooting section or search online using the resources above!

---

**Good luck with your analysis! 🧬📊**
