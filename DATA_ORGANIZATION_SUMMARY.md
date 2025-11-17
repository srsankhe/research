# Data Organization Summary

This document describes the reorganized structure of the GSE datasets in this repository.

## Changes Made

### 1. Archive Organization
All compressed data files (`.gz` and `.tar.gz`) have been moved to `archives/` subdirectories within each GSE folder to keep the workspace clean and organized.

### 2. Sample Metadata Files
Created consolidated metadata files for easy data access:
- **GSE208526_sample_metadata.csv** - Sample information with file paths for GSE208526
- **GSE282344_sample_metadata.csv** - Sample information with file paths for GSE282344

### 3. Updated Scripts
Conversion scripts have been updated to read data from the new `archives/` locations.

---

## GSE208526_RAW Structure

```
GSE208526_RAW/
├── archives/                          # All compressed files
│   ├── GSM6347353_sham.tar.gz        # Original tar.gz archives
│   ├── GSM6347354_sham_CD45-.tar.gz
│   ├── GSM6347355_IRI.tar.gz
│   ├── GSM6347356_IRI_CD45-.tar.gz
│   ├── sham/                          # Extracted sample data
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── sham_CD45-/
│   ├── IRI/
│   └── IRI_CD45-/
├── GSE208526_sample_metadata.csv     # Sample metadata file
├── consolidate_data.py                # Data consolidation script
├── convert_to_csv.py                  # Conversion script (updated)
├── deseq2_analysis.R                  # R analysis script
├── install_packages.R
└── readme.md                          # Comprehensive guide (updated)
```

**Dataset Details:**
- **Organism:** Rat (Rattus norvegicus)
- **Tissue:** Kidney
- **Conditions:** IRI (Ischemia-Reperfusion Injury) vs Sham
- **Cell Types:** CD45+ (immune) and CD45- (non-immune)
- **4 samples:** sham, sham_CD45-, IRI, IRI_CD45-
- **Format:** 10X Genomics single-cell RNA-seq

---

## GSE282344_RAW Structure

```
GSE282344_RAW/
├── archives/                          # All compressed files
│   ├── GSM8641663_NR_JK_047_*.gz     # Full transcriptome samples
│   ├── GSM8641664_NR_JK_048_*.gz
│   ├── GSM8641665_NR_JK_049_*.gz
│   ├── GSM8641666_NR_JK_053_*.gz
│   ├── GSM8641667_NR_JK_050_*.gz     # Protein markers only
│   ├── GSM8641668_NR_JK_051_*.gz
│   ├── GSM8641669_NR_JK_052_*.gz
│   └── GSM8641670_NR_JK_054_*.gz
├── GSE282344_sample_metadata.csv     # Sample metadata file
├── GSE282344_feature_reference.csv   # Feature reference
├── consolidate_data.py                # Data consolidation script
├── convert_to_csv.py                  # Conversion script (updated)
├── deseq2_analysis.R                  # R analysis script
├── install_packages.R
└── readme.md                          # Comprehensive guide (updated)
```

**Dataset Details:**
- **Organism:** Human (Homo sapiens)
- **8 samples total:**
  - 4 full transcriptome samples (047, 048, 049, 053) - ~33,538 genes
  - 4 protein marker samples (050, 051, 052, 054) - 6 CD markers only
- **Format:** 10X Genomics single-cell RNA-seq

---

## Sample Metadata Files

### GSE208526_sample_metadata.csv
Contains sample information for the rat kidney dataset:
- `sample_id`: Sample name (sham, sham_CD45-, IRI, IRI_CD45-)
- `condition`: Experimental condition (sham or IRI)
- `cell_type`: Cell type (CD45+ or CD45-)
- `barcodes_file`: Path to barcodes file in archives
- `features_file`: Path to features file in archives
- `matrix_file`: Path to matrix file in archives

### GSE282344_sample_metadata.csv
Contains sample information for the human dataset:
- `gsm_id`: GEO sample accession (GSM8641663_NR_JK_047, etc.)
- `sample_id`: Short sample ID (NR_JK_047, etc.)
- `barcode_file`: Path to barcodes file in archives
- `genes_file`: Path to genes file in archives
- `matrix_file`: Path to matrix file in archives

---

## Workflow

### Step 1: Data is organized in archives/
All raw compressed files are stored in the `archives/` folder.

### Step 2: Convert to CSV (optional)
Run the conversion script to create CSV files from the compressed data:
```bash
cd GSE208526_RAW  # or GSE282344_RAW
source ../venv/bin/activate
python convert_to_csv.py
```

This creates a `csv_output/` folder with readable CSV files.

### Step 3: Run Analysis
Use the provided R scripts or your own analysis tools:
```bash
Rscript install_packages.R  # First time only
Rscript deseq2_analysis.R
```

---

## Important Notes

1. **Archives are read-only:** The `archives/` folder contains the original compressed data. The analysis scripts read from here without modifying the originals.

2. **CSV output is optional:** CSV files are created for convenience but take more disk space. You can analyze the data directly from the compressed formats using Seurat, Scanpy, or other tools.

3. **Sample assignments:** For GSE282344, verify the control/treatment assignments in the analysis scripts match your experimental design.

4. **Virtual environment:** A Python virtual environment (`venv/`) in the repository root is used for all Python scripts.

---

## File Sizes

**GSE208526 archives/:** ~154MB (tar.gz archives) + sample data
**GSE282344 archives/:** ~54MB (all .gz files)

---

## Questions or Issues?

Refer to the `readme.md` files in each GSE folder for detailed analysis guides and troubleshooting.
