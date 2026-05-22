# MME/Mme Pathway Analysis

Comprehensive pathway analysis identifying all pathways containing the target gene MME (human) / Mme (mouse/rat) using MSigDB Hallmark and PathDIP5 Curated databases.

## Quick Summary

**Objective**: Find unique pathways in which MME/Mme exists

**Key Results**:
- ❌ **MSigDB Hallmark**: MME/Mme NOT found (0 pathways)
- ✅ **PathDIP5 Curated**: 
  - Human MME: **17 pathways**
  - Mouse/Rat Mme: **15 pathways**
  - **15 conserved** pathways (88% conservation)

## Primary Functions of MME/Mme

1. **Immune System** (5 pathways)
   - Innate immune response
   - Neutrophil degranulation
   - Hematopoietic cell lineage

2. **Protein Metabolism** (3 pathways)
   - Peptide hormone metabolism
   - Angiotensinogen metabolism
   - Protein digestion

3. **Neurodegenerative Disease** (3 pathways)
   - Alzheimer's disease (degrades amyloid-beta)

4. **Cardiovascular/Muscular** (3 pathways)
   - Cardiac conduction
   - Muscle contraction

5. **Kidney Disease** (1 pathway)
   - Focal segmental glomerulosclerosis (FSGS)

## Datasets Analyzed

- **Human**: GSE282344 (treatment vs control)
- **Rat**: GSE208526 (IRI vs sham, using mouse pathways as proxy)

## Files Generated

### Main Reports
- `PATHWAY_ANALYSIS_SUMMARY.txt` - Quick reference summary
- `MME_PATHWAY_ANALYSIS_REPORT.md` - Comprehensive markdown report

### Data Files
- `human_MME_pathdip5_curated_pathways.csv` - All human MME pathways
- `mouse_Mme_pathdip5_curated_pathways.csv` - All mouse Mme pathways
- `human_MME_pathdip5_enrichment.csv` - Human DE gene enrichment
- `rat_Mme_pathdip5_enrichment.csv` - Rat DE gene enrichment
- `pathway_species_comparison.csv` - Human vs Mouse comparison
- `pathway_analysis_summary.json` - Summary statistics

## Scripts

### 1. `mme_pathway_analysis.py`
Main analysis script that:
- Loads MSigDB Hallmark and PathDIP5 Curated databases
- Identifies pathways containing MME/Mme
- Performs enrichment analysis with DESeq2 results
- Calculates overlap with differentially expressed genes

**Usage**:
```bash
python3 mme_pathway_analysis.py
```

### 2. `create_pathway_summary.py`
Generates comprehensive reports and comparisons:
- Creates markdown report
- Compares human vs mouse pathways
- Summarizes pathway categories
- Identifies species-specific pathways

**Usage**:
```bash
python3 create_pathway_summary.py
```

## Key Findings

### 1. High Conservation
88% of pathways containing MME/Mme are conserved between human and mouse/rat

### 2. Human-Specific Pathways (2)
- Alzheimer's disease and miRNA effects (WikiPathways)
- MARKERS_NEUTROPHIL (ACSN2)

### 3. Strongest Enrichment
**Human (GSE282344)**:
- MARKERS_NEUTROPHIL: 12.5% pathway hit
- Hematopoietic cell lineage: 3.2% hit
- Neutrophil degranulation: 1.9% hit

**Rat (GSE208526)**:
- Physiological factors: 14.3% pathway hit
- Protein digestion and absorption: 9.3% hit
- Muscle contraction: 3.1% hit

### 4. Pathway Sources
- REACTOME: 9 pathways (most comprehensive)
- KEGG: 4 pathways
- WikiPathways: 3 pathways
- ACSN2: 1 pathway

## Biological Significance

**MME/Mme (Neprilysin)** is a zinc-dependent metalloproteinase that:

1. **Degrades peptides** including:
   - Amyloid-beta (Alzheimer's disease)
   - Angiotensins (blood pressure regulation)
   - Peptide hormones

2. **Critical in immune function**:
   - Neutrophil activity
   - Innate immunity
   - Hematopoiesis

3. **Disease associations**:
   - Alzheimer's disease (neuroprotective)
   - Kidney disease (FSGS)
   - Cardiovascular disease

## Databases Used

- **MSigDB Hallmark**: v2025.1.Hs (Human), v2025.1.Mm (Mouse)
- **PathDIP5 Curated**: Integrates KEGG, Reactome, WikiPathways, ACSN2
- **DESeq2 Results**: From GSE282344 and GSE208526

## Citation

If you use this analysis, please cite:
- MSigDB: Liberzon et al. (2011, 2015)
- PathDIP: Rahmati et al. (2020)
- DESeq2: Love et al. (2014)

## Contact

Analysis performed: 2024-11-16

## Visualization Plots

### Generated Visualizations (9 plots total)

All plots are available in the `plots/` directory:

#### 1. Overview Plots
- **01_pathway_categories.png** - Pathway distribution by biological category
- **02_pathway_sources.png** - Distribution of pathways across databases (KEGG, Reactome, WikiPathways)

#### 2. Enrichment Bar Charts
- **03_human_enrichment_top15.png** - Top 15 most enriched pathways in human data
- **04_rat_enrichment_top15.png** - Top 15 most enriched pathways in rat data

#### 3. Advanced Visualizations
- **05_size_vs_enrichment_scatter.png** - Relationship between pathway size and enrichment level
- **06_enrichment_heatmap_by_source.png** - Heatmaps showing enrichment patterns by database source
- **07_all_pathways_dotplot.png** - Comprehensive dot plot of all pathway enrichments
- **08_species_comparison.png** - Direct comparison of conserved pathway enrichment between species

#### 4. Summary
- **09_summary_statistics.png** - Multi-panel overview with key statistics

### How to View Plots

```bash
# Navigate to plots directory
cd pathway_analysis_results/plots/

# View on Mac
open *.png

# View on Linux with default image viewer
xdg-open *.png

# Or open individual plots
open 03_human_enrichment_top15.png
```

### Plot Color Coding

- 🔴 **Red** (>5%): High enrichment - pathway strongly affected
- 🔵 **Blue** (2-5%): Moderate enrichment - notable pathway involvement
- ⚪ **Gray** (<2%): Low enrichment - minimal pathway involvement

### Understanding the Plots

**Dot/Scatter Plots:**
- **Dot Size**: Proportional to number of overlapping DE genes
- **Color Intensity**: Represents pathway size (total genes in pathway)
- **Position**: X or Y axis shows enrichment percentage

**Heatmaps:**
- Darker colors indicate higher enrichment
- Values show exact percentage of pathway hit rate
- Organized by pathway source database

