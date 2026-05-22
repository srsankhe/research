# Single-Cell RNA-seq Analysis Guide: GSE282344 Dataset

## Complete Guide: From Raw Data to Differential Expression Results

This guide will walk you through analyzing single-cell RNA-seq data from the GSE282344 dataset. This is a human single-cell transcriptomics study.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Understanding Your Data](#understanding-your-data)
3. [Setting Up Your Environment](#setting-up-your-environment)
4. [Step 1: Converting Raw Data to CSV](#step-1-converting-raw-data-to-csv)
5. [Step 2: Running DESeq2 Analysis](#step-2-running-deseq2-analysis)
6. [Understanding the Results](#understanding-the-results)
7. [Important Notes](#important-notes)
8. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### What You Need

- **Python 3** (for data conversion)
- **R** (for DESeq2 analysis)
- **8GB+ RAM** recommended
- **Basic command line knowledge**

See the detailed setup instructions in [Setting Up Your Environment](#setting-up-your-environment).

---

## Understanding Your Data

### Dataset: GSE282344

This dataset contains **human** single-cell RNA-seq data with:
- **8 samples total** (GSM8641663-GSM8641670)
- **Two types of data:**
  - **4 transcriptome samples** (047, 048, 049, 053) - ~33,538 genes each
  - **4 protein marker samples** (050, 051, 052, 054) - only 6 CD markers

### Files in Your Folder

**Directory Structure:**
```
archives/                            (All compressed data files)
├── GSM8641663_NR_JK_047_barcodes.tsv.gz
├── GSM8641663_NR_JK_047_genes.tsv.gz
├── GSM8641663_NR_JK_047_matrix.mtx.gz
├── GSM8641664_NR_JK_048_*.gz       (Sample 048)
├── GSM8641665_NR_JK_049_*.gz       (Sample 049)
├── GSM8641666_NR_JK_053_*.gz       (Sample 053)
├── GSM8641667_NR_JK_050_*.gz       (Sample 050 - protein markers only)
├── GSM8641668_NR_JK_051_*.gz       (Sample 051 - protein markers only)
├── GSM8641669_NR_JK_052_*.gz       (Sample 052 - protein markers only)
└── GSM8641670_NR_JK_054_*.gz       (Sample 054 - protein markers only)
```

### File Formats

Each transcriptome sample has three files in the archives folder:
- **genes.tsv.gz** - Gene names
- **barcodes.tsv.gz** - Cell barcodes
- **matrix.mtx.gz** - Expression counts (sparse matrix format)

**Note:** All `.gz` files have been organized into the `archives/` folder. The conversion script automatically reads from this location.

---

## Setting Up Your Environment

### Step 1: Install Python 3

**Check if you have Python:**
```bash
python3 --version
```

**On Mac:**
```bash
brew install python3
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

### Step 2: Install R

**Check if you have R:**
```bash
R --version
```

**On Mac:**
```bash
brew install r
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt install r-base
```

### Step 3: Navigate to Your Data

```bash
cd /path/to/GSE282344_RAW
```

**Tip:** You can drag-and-drop the folder into Terminal to auto-fill the path!

---

## Step 1: Converting Raw Data to CSV

### 1.1: Data Organization

All compressed `.gz` files have been organized into the `archives/` directory to keep the workspace clean. The conversion script automatically reads from this location - no manual extraction needed!

### 1.2: Create Python Virtual Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate  # Mac/Linux
# OR
venv\Scripts\activate     # Windows
```

### 1.3: Install Required Packages

```bash
pip install pandas scipy numpy
```

**What these do:**
- **pandas** - Data manipulation
- **scipy** - Reads sparse matrices
- **numpy** - Numerical operations

### 1.4: Run the Conversion Script

```bash
python3 convert_to_csv.py
```

**What's happening:**
1. Reads compressed `.gz` files from `archives/` folder
2. Converts sparse matrices to CSV
3. Saves to `csv_output/` folder
4. **Automatically skips** protein-only samples (050, 051, 052, 054)

**Expected output:**
```
Processing sample: GSM8641663_NR_JK_047
  ✓ Saved: GSM8641663_NR_JK_047_genes.csv
  ✓ Saved: GSM8641663_NR_JK_047_barcodes.csv
  ✓ Saved: GSM8641663_NR_JK_047_matrix.csv
  Matrix shape: 33538 genes × 2772 cells

... (similar for 048, 049, 053)
```

Takes **3-5 minutes**.

### 1.5: Check Your Output

```bash
ls -lh csv_output/
```

You should see 12 CSV files (3 per transcriptome sample).

---

## Step 2: Running DESeq2 Analysis

### 2.1: Install R Packages

```bash
Rscript install_packages.R
```

**What's being installed:**
- **DESeq2** - Differential expression analysis
- **ggplot2** - Plots
- **pheatmap** - Heatmaps
- **RColorBrewer** - Color palettes

Takes **10-15 minutes**. Lots of output is normal!

### 2.2: Run DESeq2 Analysis

```bash
Rscript deseq2_analysis.R
```

**What's happening:**

1. **Loading data** - Reads CSV matrices
2. **Creating pseudo-bulk** - Sums cells into bulk samples
3. **Sample assignment** - **IMPORTANT!** Assigns conditions:
   - Samples 047, 048 → **control**
   - Samples 049, 053 → **treatment**
4. **DESeq2 analysis** - Statistical testing
5. **Generating plots** - Visualizations

**⚠️ WARNING:** Sample condition assignments are **ASSUMED**! Verify they match your experiment!

**Expected output:**
```
Sample assignments:
                       sample condition
GSM8641663_NR_JK_047 ...      control
GSM8641664_NR_JK_048 ...      control
GSM8641665_NR_JK_049 ...      treatment
GSM8641666_NR_JK_053 ...      treatment

Significant genes (padj < 0.05, |log2FC| > 1): 79
  Upregulated in treatment: 69
  Downregulated in treatment: 10
```

Takes **3-5 minutes**.

### 2.3: Check Your Results

```bash
ls -lh deseq2_output/
```

You should see:
- CSV files with results
- PNG images with plots
- `analysis_summary.txt`

---

## Understanding the Results

### Key Findings (Based on ASSUMED Condition Assignments)

- **Total genes analyzed:** 26,703
- **Significant genes:** 79 (padj < 0.05, |log2FC| > 1)
  - **Upregulated in treatment:** 69 genes
  - **Downregulated in treatment:** 10 genes

### Top Upregulated Genes

1. **CYP2A7** (log2FC = 23.1) - Cytochrome P450 enzyme
2. **TPSB2** (log2FC = 6.7) - Tryptase beta 2
3. **MUC5AC** (log2FC = 5.2) - Mucin (airway secretion)
4. **AC020656.1** (log2FC = 4.6) - Long non-coding RNA
5. **FLT3** (log2FC = 4.2) - FMS-like tyrosine kinase

### Top Downregulated Genes

1. **IL22** (log2FC = -23.8) - Interleukin 22
2. **CXCL9** (log2FC = -10.3) - Chemokine (immune response)
3. **CXCL10** (log2FC = -9.2) - Chemokine (inflammation)
4. **CXCL11** (log2FC = -6.2) - Chemokine
5. **IDO1** (log2FC = -5.8) - Indoleamine 2,3-dioxygenase

### Result Files

#### 1. `treatment_vs_control_results.csv`
All genes with statistics.

**Columns:**
- `baseMean` - Average expression
- `log2FoldChange` - Fold change (treatment vs control)
- `lfcSE` - Standard error
- `stat` - Test statistic
- `pvalue` - Raw p-value
- `padj` - Adjusted p-value (use this!)

#### 2. `treatment_vs_control_significant_genes.csv`
Only significant genes (padj < 0.05, |log2FC| > 1).

#### 3. Plots

**volcano_plot.png:**
- Shows all genes
- Red = significant
- Genes far from center = big changes

**MA_plot.png:**
- Expression vs fold change
- Blue = significant

**PCA_plot.png:**
- Shows sample clustering
- Control samples should cluster together
- Treatment samples should cluster together

**heatmap_top50_genes.png:**
- Top 50 most changed genes
- Red = high expression
- Blue = low expression

---

## Important Notes

### ⚠️ Sample Condition Assignments

**CRITICAL:** The analysis assumes:
- **Samples 047, 048 = control**
- **Samples 049, 053 = treatment**

**How to verify/change:**

1. Check GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282344
2. Look at the sample descriptions
3. If assignments are wrong, edit `deseq2_analysis.R`:

```r
# Find this section (around line 35):
coldata <- data.frame(
  sample = samples,
  condition = c("control", "control", "treatment", "treatment"),  # MODIFY THIS
  row.names = samples
)
```

4. Re-run: `Rscript deseq2_analysis.R`

### Why Some Samples Were Excluded

Samples 050, 051, 052, 054 contain **only 6 genes**:
- CD24.1, CD44.1, CD13, CD10, CD326, CD133

These are **protein markers** (likely CyTOF or CITE-seq data), not transcriptomics. DESeq2 requires full transcriptome data, so they were automatically excluded.

---

## Troubleshooting

### Problem: "ModuleNotFoundError: No module named 'pandas'"

**Solution:**
```bash
# Make sure venv is activated
source venv/bin/activate  # You should see (venv) in prompt

# Install packages
pip install pandas scipy numpy
```

### Problem: "Error in library(DESeq2)"

**Solution:**
```bash
# Run package installer
Rscript install_packages.R
```

### Problem: Wrong sample assignments

**Solution:**
1. Open `deseq2_analysis.R` in a text editor
2. Find line ~35 with condition assignments
3. Modify the condition vector
4. Save and re-run

### Problem: Out of memory

**Solution:**
- Close other applications
- Use a computer with more RAM (8GB minimum)

---

## File Structure

```
GSE282344_RAW/
├── readme.md                          (This file!)
├── convert_to_csv.py                  (Python conversion script)
├── consolidate_data.py                (Data consolidation script)
├── GSE282344_sample_metadata.csv      (Sample metadata with file paths)
├── GSE282344_feature_reference.csv    (Feature reference data)
├── install_packages.R                 (R package installer)
├── deseq2_analysis.R                  (Main analysis script)
│
├── archives/                          (All compressed data files)
│   ├── GSM8641663_NR_JK_047_*.gz     (Sample 047 - control)
│   ├── GSM8641664_NR_JK_048_*.gz     (Sample 048 - control)
│   ├── GSM8641665_NR_JK_049_*.gz     (Sample 049 - treatment)
│   ├── GSM8641666_NR_JK_053_*.gz     (Sample 053 - treatment)
│   ├── GSM8641667_NR_JK_050_*.gz     (Excluded - protein markers only)
│   ├── GSM8641668_NR_JK_051_*.gz     (Excluded - protein markers only)
│   ├── GSM8641669_NR_JK_052_*.gz     (Excluded - protein markers only)
│   └── GSM8641670_NR_JK_054_*.gz     (Excluded - protein markers only)
│
├── csv_output/                        (Converted CSV files)
│   ├── GSM8641663_NR_JK_047_matrix.csv
│   ├── GSM8641663_NR_JK_047_genes.csv
│   ├── GSM8641663_NR_JK_047_barcodes.csv
│   └── ... (similar for 048, 049, 053)
│
├── deseq2_output/                     (Analysis results)
│   ├── analysis_summary.txt           (Human-readable summary)
│   ├── treatment_vs_control_results.csv (All genes)
│   ├── treatment_vs_control_significant_genes.csv (Significant only)
│   ├── pseudobulk_counts.csv          (Count matrix)
│   ├── sample_metadata.csv            (Sample info - VERIFY THIS!)
│   ├── volcano_plot.png
│   ├── MA_plot.png
│   ├── PCA_plot.png
│   ├── heatmap_top50_genes.png
│   └── dispersion_plot.png
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
grep "GENE_NAME" deseq2_output/treatment_vs_control_results.csv
```

---

## Next Steps

### 1. Verify Sample Assignments!
**MOST IMPORTANT:** Check that control/treatment assignments are correct!

1. Visit GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282344
2. Check sample descriptions
3. Modify `deseq2_analysis.R` if needed
4. Re-run analysis

### 2. Pathway Analysis

Submit significant genes to:
- **DAVID:** https://david.ncifcrf.gov/
- **Enrichr:** https://maayanlab.cloud/Enrichr/
- **Reactome:** https://reactome.org/

### 3. Literature Search

Research top genes on:
- **PubMed:** https://pubmed.ncbi.nlm.nih.gov/
- **GeneCards:** https://www.genecards.org/

### 4. Validation

For publication, validate with:
- qRT-PCR
- Western blot
- Immunohistochemistry

---

## Dataset Information

**GEO Accession:** GSE282344

**Organism:** *Homo sapiens* (Human)

**Data Type:** Single-cell RNA-seq

**Samples Analyzed:**
- GSM8641663_NR_JK_047 - 33,538 genes × 2,772 cells
- GSM8641664_NR_JK_048 - 33,538 genes × 2,741 cells
- GSM8641665_NR_JK_049 - 33,538 genes × 3,136 cells
- GSM8641666_NR_JK_053 - 33,538 genes × 1,914 cells

**Samples Excluded:**
- GSM8641667-70 (samples 050-054) - protein markers only

---

## Citation

If you use this analysis in a publication, cite:

**Dataset:**
- GEO Accession: GSE282344

**Software:**
- **DESeq2:** Love, M.I., Huber, W., Anders, S. (2014). Genome Biology, 15:550.
- **Python packages:** pandas, scipy, numpy
- **R:** R Core Team (2024)

---

## Getting Help

### Online Resources

1. **GEO Database:** https://www.ncbi.nlm.nih.gov/geo/
2. **DESeq2 Vignette:** https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
3. **Bioconductor Support:** https://support.bioconductor.org/
4. **Biostars:** https://www.biostars.org/

### Common Search Terms

- "DESeq2 tutorial"
- "single cell RNA-seq to pseudo-bulk"
- "interpreting volcano plot RNA-seq"
- "GSE282344" (for dataset-specific info)

---

## Summary

This dataset contains human single-cell RNA-seq data with 4 transcriptome samples. The analysis pipeline:

1. ✅ Converted sparse matrices to CSV
2. ✅ Created pseudo-bulk samples
3. ✅ Ran DESeq2 differential expression
4. ✅ Generated plots and results

**⚠️ IMPORTANT:** Verify sample condition assignments before trusting results!

**Good luck with your analysis! 🧬📊**
