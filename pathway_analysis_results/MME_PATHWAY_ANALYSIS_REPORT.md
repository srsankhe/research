# MME/Mme Pathway Analysis Report

Analysis Date: 2024-11-16

## Executive Summary

- **Human MME**: Found in 17 unique pathways (PathDIP5 Curated)
- **Mouse/Rat Mme**: Found in 15 unique pathways (PathDIP5 Curated)
- **MSigDB Hallmark**: MME/Mme not found in Hallmark gene sets

## 1. Human MME Pathways (PathDIP5 Curated)

### All Pathways Containing MME

#### ACSN2 (1 pathways)

- MARKERS_NEUTROPHIL

#### KEGG (4 pathways)

- Alzheimer disease
- Hematopoietic cell lineage
- Protein digestion and absorption
- Renin-angiotensin system

#### REACTOME (9 pathways)

- Cardiac conduction
- Immune System
- Innate Immune System
- Metabolism of Angiotensinogen to Angiotensins
- Metabolism of proteins
- Muscle contraction
- Neutrophil degranulation
- Peptide hormone metabolism
- Physiological factors

#### WikiPathways (3 pathways)

- Alzheimer's disease
- Alzheimer's disease and miRNA effects
- Primary focal segmental glomerulosclerosis (FSGS)

### Enrichment Analysis (GSE282344 - Human)

Pathways with highest overlap of differentially expressed genes:

| Pathway | Source | Size | DE Genes | % Hit | Top Overlapping Genes |
|---------|--------|------|----------|-------|-----------------------|
| MARKERS_NEUTROPHIL | ACSN2 | 8 | 1 | 12.5% | PECAM1 |
| Hematopoietic cell lineage | KEGG | 94 | 3 | 3.2% | FLT3, CD1D, CSF3R |
| Neutrophil degranulation | REACTOME | 478 | 9 | 1.9% | FCAR, CFD, CRISPLD2, ATP8B4, SLC11A1, C5AR1, PECAM... |
| Innate Immune System | REACTOME | 1151 | 14 | 1.2% | FCAR, C5AR2, MUC5AC, CFD, IGHG1, CRISPLD2, ATP8B4,... |
| Immune System | REACTOME | 2114 | 24 | 1.1% | FCAR, FLT3, CD1D, IL10, C5AR2, IL22, MUC5AC, CFD, ... |
| Cardiac conduction | REACTOME | 130 | 1 | 0.8% | KCNE1 |
| Muscle contraction | REACTOME | 204 | 1 | 0.5% | KCNE1 |
| Alzheimer's disease | WikiPathways | 264 | 1 | 0.4% | LRP1 |
| Metabolism of proteins | REACTOME | 1937 | 6 | 0.3% | MUC5AC, VNN2, ADAMTS2, NUP214, LTF, FN1 |
| Alzheimer's disease and miRNA effects | WikiPathways | 329 | 1 | 0.3% | LRP1 |

## 2. Mouse/Rat Mme Pathways (PathDIP5 Curated)

### All Pathways Containing Mme

#### KEGG (4 pathways)

- Alzheimer disease
- Hematopoietic cell lineage
- Protein digestion and absorption
- Renin-angiotensin system

#### REACTOME (9 pathways)

- Cardiac conduction
- Immune System
- Innate Immune System
- Metabolism of Angiotensinogen to Angiotensins
- Metabolism of proteins
- Muscle contraction
- Neutrophil degranulation
- Peptide hormone metabolism
- Physiological factors

#### WikiPathways (2 pathways)

- Alzheimer's disease
- Primary focal segmental glomerulosclerosis (FSGS)

### Enrichment Analysis (GSE208526 - Rat IRI model)

Pathways with highest overlap of differentially expressed genes:

| Pathway | Source | Size | DE Genes | % Hit | Top Overlapping Genes |
|---------|--------|------|----------|-------|-----------------------|
| Physiological factors | REACTOME | 7 | 1 | 14.3% | NPR1 |
| Protein digestion and absorption | KEGG | 108 | 10 | 9.3% | COL6A2, COL6A1, COL12A1, COL5A3, COL11A2, MEP1B, C... |
| Peptide hormone metabolism | REACTOME | 64 | 2 | 3.1% | IGF1, INHBA |
| Muscle contraction | REACTOME | 160 | 5 | 3.1% | TPM2, DES, NPR1, MYL9, FXYD3 |
| Alzheimer disease | KEGG | 371 | 9 | 2.4% | FZD9, TUBB6, TUBB3, LPL, SLC39A2, WNT7A, TUBB4B, D... |
| Metabolism of proteins | REACTOME | 1635 | 37 | 2.3% | APOL9A, TIMP1, TUBB6, TNC, TUBB3, SPP1, MSLN, FBN1... |
| Cardiac conduction | REACTOME | 98 | 2 | 2.0% | NPR1, FXYD3 |
| Immune System | REACTOME | 1617 | 32 | 2.0% | IFITM3, CXCL1, TUBB6, CLU, TUBB3, C7, MX2, SERPING... |
| Innate Immune System | REACTOME | 959 | 14 | 1.5% | CXCL1, CLU, C7, SERPING1, LCN2, CXCL3, FABP5, DYNL... |
| Alzheimer's disease | WikiPathways | 74 | 1 | 1.4% | LPL |

## 3. Pathway Comparison: Human vs Mouse/Rat

### Shared Pathways (15)

Pathways containing MME/Mme in both species:

- Alzheimer disease
- Alzheimer's disease
- Cardiac conduction
- Hematopoietic cell lineage
- Immune System
- Innate Immune System
- Metabolism of Angiotensinogen to Angiotensins
- Metabolism of proteins
- Muscle contraction
- Neutrophil degranulation
- Peptide hormone metabolism
- Physiological factors
- Primary focal segmental glomerulosclerosis (FSGS)
- Protein digestion and absorption
- Renin-angiotensin system

### Human-Only Pathways (2)

Pathways containing MME only in human:

- Alzheimer's disease and miRNA effects
- MARKERS_NEUTROPHIL

### Mouse-Only Pathways (0)

Pathways containing Mme only in mouse/rat:


## 4. Key Findings

### Major Categories

#### Human MME Pathway Categories:

- **Immune system**: 5 pathways
- **Neuropsychiatric and neurodevelopmental disorders**: 3 pathways
- **Folding, sorting and degradation**: 3 pathways
- **Muscular and bone system**: 3 pathways
- **Endocrine system**: 1 pathways
- **Digestive system**: 1 pathways
- **Immune diseases**: 1 pathways

#### Mouse/Rat Mme Pathway Categories:

- **Immune system**: 4 pathways
- **Folding, sorting and degradation**: 3 pathways
- **Muscular and bone system**: 3 pathways
- **Neuropsychiatric and neurodevelopmental disorders**: 2 pathways
- **Digestive system**: 1 pathways
- **Endocrine system**: 1 pathways
- **Immune diseases**: 1 pathways

### Notable Pathways

1. **Immune System**: MME/Mme is found in multiple immune-related pathways including:
   - Innate Immune System
   - Neutrophil degranulation
   - Hematopoietic cell lineage

2. **Alzheimer's Disease**: MME/Mme is involved in Alzheimer's disease pathways, consistent with its role as a neprilysin enzyme that degrades amyloid-beta.

3. **Renin-Angiotensin System**: Present in both species, important for cardiovascular and kidney function.

4. **Protein Metabolism**: Multiple pathways related to peptide hormone metabolism and protein digestion.

## 5. Conclusions

- MME/Mme is primarily involved in **immune system** and **protein metabolism** pathways
- Strong conservation between human and mouse/rat (most pathways are shared)
- Enrichment analysis shows overlap with DE genes, particularly in immune pathways
- Not found in MSigDB Hallmark gene sets, but well-represented in PathDIP5 curated pathways

## Data Sources

- **PathDIP5 Curated**: Curated pathway database integrating KEGG, Reactome, WikiPathways, ACSN2
- **MSigDB Hallmark**: Hallmark gene sets from MSigDB (v2025.1)
- **Human Dataset**: GSE282344
- **Rat Dataset**: GSE208526 (using mouse pathways as proxy)
