# Script Generalization Summary

## What Was Done

Successfully generalized all common scripts from GSE208526_RAW and GSE282344_RAW folders and extracted them to the parent directory.

## Scripts Generalized (6 total)

1. **consolidate_data.py** - Creates sample metadata CSV
2. **convert_to_csv.py** - Converts 10X format to CSV
3. **deseq2_analysis.R** - DESeq2 differential expression analysis
4. **install_packages.R** - Installs R packages (already identical)
5. **install_seurat.R** - Installs Seurat and dependencies
6. **single_cell_seurat_analysis.R** - Complete Seurat single-cell pipeline

## Key Improvements

### Before
- Each dataset had its own hardcoded scripts
- Hardcoded paths (e.g., `/Users/ssankhe/Downloads/GSE208526_RAW/`)
- Hardcoded sample names and conditions
- Different implementations for similar tasks
- Difficult to maintain and update

### After
- **Single source of truth**: All scripts in parent directory
- **Configuration-driven**: JSON configs for each dataset
- **Parameterized**: No hardcoded paths or sample names
- **Reusable**: Works with any GSE dataset
- **Easy to maintain**: Update once, benefits all datasets
- **Symlinks**: Each dataset folder links to parent scripts

## Directory Structure

```
research/
├── consolidate_data.py               # Generalized scripts
├── convert_to_csv.py
├── deseq2_analysis.R
├── install_packages.R
├── install_seurat.R
├── single_cell_seurat_analysis.R
├── GENERALIZED_SCRIPTS_README.md     # Complete documentation
│
├── GSE208526_RAW/
│   ├── consolidate_data.py           # → symlink to ../consolidate_data.py
│   ├── convert_to_csv.py             # → symlink to ../convert_to_csv.py
│   ├── deseq2_analysis.R             # → symlink to ../deseq2_analysis.R
│   ├── install_packages.R            # → symlink to ../install_packages.R
│   ├── install_seurat.R              # → symlink to ../install_seurat.R
│   ├── single_cell_seurat_analysis.R # → symlink
│   ├── config_consolidate.json       # Dataset-specific config
│   ├── config_convert.json
│   ├── config_deseq2.json
│   └── config_seurat.json
│
└── GSE282344_RAW/
    ├── (same structure as GSE208526_RAW)
    └── ...
```

## Usage Examples

### GSE208526 Analysis
```bash
cd GSE208526_RAW

# Run with dataset-specific config
python3 consolidate_data.py --config config_consolidate.json --dataset-id GSE208526
python3 convert_to_csv.py --config config_convert.json
Rscript deseq2_analysis.R config_deseq2.json
Rscript single_cell_seurat_analysis.R config_seurat.json
```

### GSE282344 Analysis
```bash
cd GSE282344_RAW

# Same commands, different configs
python3 consolidate_data.py --config config_consolidate.json --dataset-id GSE282344
python3 convert_to_csv.py --config config_convert.json
Rscript deseq2_analysis.R config_deseq2.json
Rscript single_cell_seurat_analysis.R config_seurat.json
```

## Configuration Files

Each dataset has 4 config files:
- `config_consolidate.json` - For consolidate_data.py
- `config_convert.json` - For convert_to_csv.py
- `config_deseq2.json` - For deseq2_analysis.R
- `config_seurat.json` - For single_cell_seurat_analysis.R

## Benefits

1. ✅ **Maintainability**: Fix bugs in one place
2. ✅ **Consistency**: All datasets use same analysis pipeline
3. ✅ **Scalability**: Easy to add new datasets
4. ✅ **Flexibility**: JSON configs allow customization per dataset
5. ✅ **Documentation**: Self-documenting through config structure
6. ✅ **Version Control**: Symlinks ensure consistency

## Adding New Datasets

1. Create new GSE folder
2. Create symlinks to parent scripts (or copy the pattern)
3. Create dataset-specific config JSON files
4. Run the pipeline

See `GENERALIZED_SCRIPTS_README.md` for complete documentation.
