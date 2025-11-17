# Generalized Single-Cell Analysis Scripts

This directory contains generalized scripts for single-cell RNA-seq analysis that work with any GSE dataset. All scripts are parameterized using JSON configuration files.

## Available Scripts

### 1. `consolidate_data.py`
Creates sample metadata CSV files from dataset information.

**Usage:**
```bash
python3 consolidate_data.py --config config_consolidate.json --dataset-id GSE123456
```

**Config Format:**
```json
{
  "dataset_id": "GSE123456",
  "samples": {
    "sample1": {
      "condition": "control",
      "cell_type": "CD45+",
      "any_other_metadata": "value"
    },
    "sample2": {
      "condition": "treatment",
      "cell_type": "CD45-"
    }
  },
  "archive_mode": true,
  "file_pattern": "standard"
}
```

### 2. `convert_to_csv.py`
Converts 10X Genomics format (MTX, TSV) to CSV files for downstream analysis.

**Usage:**
```bash
python3 convert_to_csv.py --config config_convert.json [--output-dir csv_output]
```

**Config Format:**
```json
{
  "dataset_id": "GSE123456",
  "samples": {
    "sample1": {},
    "sample2": {}
  },
  "archive_dir": "archives",
  "file_pattern": "standard",
  "gene_column": "gene_name"
}
```

**File Patterns:**
- `"standard"`: Standard 10X format in sample directories (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
- `"gsm_prefix"`: Files with GSM prefix (GSM123_sample_barcodes.tsv.gz, etc.)

### 3. `deseq2_analysis.R`
DESeq2 differential expression analysis on pseudo-bulk data.

**Usage:**
```bash
Rscript deseq2_analysis.R config_deseq2.json
```

**Config Format:**
```json
{
  "dataset_id": "GSE123456",
  "csv_dir": "csv_output",
  "output_dir": "deseq2_output",
  "samples": ["sample1", "sample2", "sample3", "sample4"],
  "conditions": ["control", "control", "treatment", "treatment"],
  "comparison": {
    "numerator": "treatment",
    "denominator": "control"
  },
  "additional_metadata": {
    "cell_type": ["CD45+", "CD45-", "CD45+", "CD45-"]
  }
}
```

### 4. `single_cell_seurat_analysis.R`
Complete single-cell analysis pipeline using Seurat.

**Usage:**
```bash
Rscript single_cell_seurat_analysis.R config_seurat.json
```

**Config Format:**
```json
{
  "dataset_id": "GSE123456",
  "data_format": "standard",
  "archive_dir": "archives",
  "samples": {
    "sample1": {
      "condition": "control",
      "batch": "1"
    },
    "sample2": {
      "condition": "treatment",
      "batch": "2"
    }
  },
  "output_dir": "single_cell_analysis",
  "qc_filters": {
    "min_features": 200,
    "max_features": 5000,
    "max_mt_percent": 10
  },
  "analysis_params": {
    "n_variable_features": 2000,
    "n_pca_dims": 30,
    "clustering_resolution": 0.5
  },
  "genes_of_interest": ["Gene1", "Gene2"]
}
```

### 5. `install_packages.R` and `install_seurat.R`
Install required R packages for analysis.

**Usage:**
```bash
Rscript install_packages.R
Rscript install_seurat.R
```

## Directory Structure

```
research/
├── consolidate_data.py          # Generalized scripts
├── convert_to_csv.py
├── deseq2_analysis.R
├── install_packages.R
├── install_seurat.R
├── single_cell_seurat_analysis.R
├── GSE208526_RAW/
│   ├── consolidate_data.py      # Symlink to parent
│   ├── convert_to_csv.py        # Symlink to parent
│   ├── deseq2_analysis.R        # Symlink to parent
│   ├── install_packages.R       # Symlink to parent
│   ├── install_seurat.R         # Symlink to parent
│   ├── single_cell_seurat_analysis.R  # Symlink to parent
│   ├── config_consolidate.json  # Dataset-specific config
│   ├── config_convert.json      # Dataset-specific config
│   ├── config_deseq2.json       # Dataset-specific config
│   └── config_seurat.json       # Dataset-specific config
└── GSE282344_RAW/
    ├── (same symlinks and configs)
    └── ...
```

## Example Workflows

### Complete Analysis Pipeline

#### For GSE208526:
```bash
cd GSE208526_RAW

# 1. Create sample metadata
python3 consolidate_data.py --config config_consolidate.json --dataset-id GSE208526

# 2. Convert to CSV (if needed for DESeq2)
python3 convert_to_csv.py --config config_convert.json

# 3. Run DESeq2 analysis
Rscript deseq2_analysis.R config_deseq2.json

# 4. Run Seurat single-cell analysis
Rscript single_cell_seurat_analysis.R config_seurat.json
```

#### For GSE282344:
```bash
cd GSE282344_RAW

# Same commands, different configs
python3 consolidate_data.py --config config_consolidate.json --dataset-id GSE282344
python3 convert_to_csv.py --config config_convert.json
Rscript deseq2_analysis.R config_deseq2.json
Rscript single_cell_seurat_analysis.R config_seurat.json
```

## Benefits of This Approach

1. **Single Source of Truth**: All scripts are in the parent directory. Updates benefit all datasets.
2. **Easy to Maintain**: Fix bugs or add features in one place.
3. **Dataset-Specific Configuration**: Each dataset has its own config files for customization.
4. **Reusable**: Easy to add new datasets - just create new config files.
5. **Version Control Friendly**: Symlinks ensure all datasets use the same script version.

## Adding a New Dataset

1. Create a new GSE folder
2. Create symlinks to the generalized scripts:
   ```bash
   cd GSE_NEW_DATASET_RAW
   ln -s ../consolidate_data.py .
   ln -s ../convert_to_csv.py .
   ln -s ../deseq2_analysis.R .
   ln -s ../install_packages.R .
   ln -s ../install_seurat.R .
   ln -s ../single_cell_seurat_analysis.R .
   ```
3. Create config files specific to your dataset
4. Run the analysis pipeline

## Notes

- **Python Requirements**: scipy, pandas, numpy
- **R Requirements**: DESeq2, Seurat, ggplot2, pheatmap, jsonlite
- Install R packages first: `Rscript install_packages.R && Rscript install_seurat.R`
- Modify QC filters and analysis parameters in config files based on your data quality
- Review condition assignments in DESeq2 config to match your experimental design
