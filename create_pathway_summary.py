#!/usr/bin/env python3
"""
Create summary report and visualizations for MME/Mme pathway analysis
"""
import pandas as pd
import os

def create_summary_report():
    results_dir = '/Users/ssankhe/Documents/Github/research/pathway_analysis_results'
    
    # Load data
    human_pathways = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_curated_pathways.csv'))
    human_enrichment = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_enrichment.csv'))
    mouse_pathways = pd.read_csv(os.path.join(results_dir, 'mouse_Mme_pathdip5_curated_pathways.csv'))
    rat_enrichment = pd.read_csv(os.path.join(results_dir, 'rat_Mme_pathdip5_enrichment.csv'))
    
    report_file = os.path.join(results_dir, 'MME_PATHWAY_ANALYSIS_REPORT.md')
    
    with open(report_file, 'w') as f:
        f.write("# MME/Mme Pathway Analysis Report\n\n")
        f.write("Analysis Date: 2024-11-16\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"- **Human MME**: Found in {len(human_pathways)} unique pathways (PathDIP5 Curated)\n")
        f.write(f"- **Mouse/Rat Mme**: Found in {len(mouse_pathways)} unique pathways (PathDIP5 Curated)\n")
        f.write("- **MSigDB Hallmark**: MME/Mme not found in Hallmark gene sets\n\n")
        
        # Human pathways
        f.write("## 1. Human MME Pathways (PathDIP5 Curated)\n\n")
        f.write("### All Pathways Containing MME\n\n")
        
        # Group by source
        human_by_source = human_pathways.groupby('source')['pathway_name'].apply(list).to_dict()
        
        for source, pathways in sorted(human_by_source.items()):
            f.write(f"#### {source} ({len(pathways)} pathways)\n\n")
            for pathway in sorted(pathways):
                f.write(f"- {pathway}\n")
            f.write("\n")
        
        # Enrichment analysis
        f.write("### Enrichment Analysis (GSE282344 - Human)\n\n")
        f.write("Pathways with highest overlap of differentially expressed genes:\n\n")
        
        top_enriched = human_enrichment.nlargest(10, 'percent_pathway_hit')[
            ['pathway_name', 'source', 'pathway_size', 'num_overlap', 'percent_pathway_hit', 'overlap_genes']
        ]
        
        f.write("| Pathway | Source | Size | DE Genes | % Hit | Top Overlapping Genes |\n")
        f.write("|---------|--------|------|----------|-------|-----------------------|\n")
        for _, row in top_enriched.iterrows():
            genes = row['overlap_genes'] if row['overlap_genes'] != 'None' else '-'
            if len(genes) > 50:
                genes = genes[:50] + '...'
            f.write(f"| {row['pathway_name']} | {row['source']} | {row['pathway_size']} | "
                   f"{row['num_overlap']} | {row['percent_pathway_hit']:.1f}% | {genes} |\n")
        
        f.write("\n")
        
        # Mouse/Rat pathways
        f.write("## 2. Mouse/Rat Mme Pathways (PathDIP5 Curated)\n\n")
        f.write("### All Pathways Containing Mme\n\n")
        
        # Group by source
        mouse_by_source = mouse_pathways.groupby('source')['pathway_name'].apply(list).to_dict()
        
        for source, pathways in sorted(mouse_by_source.items()):
            f.write(f"#### {source} ({len(pathways)} pathways)\n\n")
            for pathway in sorted(pathways):
                f.write(f"- {pathway}\n")
            f.write("\n")
        
        # Enrichment analysis
        f.write("### Enrichment Analysis (GSE208526 - Rat IRI model)\n\n")
        f.write("Pathways with highest overlap of differentially expressed genes:\n\n")
        
        top_enriched_rat = rat_enrichment.nlargest(10, 'percent_pathway_hit')[
            ['pathway_name', 'source', 'pathway_size', 'num_overlap', 'percent_pathway_hit', 'overlap_genes']
        ]
        
        f.write("| Pathway | Source | Size | DE Genes | % Hit | Top Overlapping Genes |\n")
        f.write("|---------|--------|------|----------|-------|-----------------------|\n")
        for _, row in top_enriched_rat.iterrows():
            genes = row['overlap_genes'] if row['overlap_genes'] != 'None' else '-'
            if len(genes) > 50:
                genes = genes[:50] + '...'
            f.write(f"| {row['pathway_name']} | {row['source']} | {row['pathway_size']} | "
                   f"{row['num_overlap']} | {row['percent_pathway_hit']:.1f}% | {genes} |\n")
        
        f.write("\n")
        
        # Unique and shared pathways
        f.write("## 3. Pathway Comparison: Human vs Mouse/Rat\n\n")
        
        human_pathway_names = set(human_pathways['pathway_name'].unique())
        mouse_pathway_names = set(mouse_pathways['pathway_name'].unique())
        
        shared = human_pathway_names & mouse_pathway_names
        human_only = human_pathway_names - mouse_pathway_names
        mouse_only = mouse_pathway_names - human_pathway_names
        
        f.write(f"### Shared Pathways ({len(shared)})\n\n")
        f.write("Pathways containing MME/Mme in both species:\n\n")
        for pathway in sorted(shared):
            f.write(f"- {pathway}\n")
        f.write("\n")
        
        f.write(f"### Human-Only Pathways ({len(human_only)})\n\n")
        f.write("Pathways containing MME only in human:\n\n")
        for pathway in sorted(human_only):
            f.write(f"- {pathway}\n")
        f.write("\n")
        
        f.write(f"### Mouse-Only Pathways ({len(mouse_only)})\n\n")
        f.write("Pathways containing Mme only in mouse/rat:\n\n")
        for pathway in sorted(mouse_only):
            f.write(f"- {pathway}\n")
        f.write("\n")
        
        # Key findings
        f.write("## 4. Key Findings\n\n")
        
        f.write("### Major Categories\n\n")
        
        # Human categories
        human_categories = human_pathways['category'].value_counts()
        f.write("#### Human MME Pathway Categories:\n\n")
        for cat, count in human_categories.items():
            f.write(f"- **{cat}**: {count} pathways\n")
        f.write("\n")
        
        # Mouse categories
        mouse_categories = mouse_pathways['category'].value_counts()
        f.write("#### Mouse/Rat Mme Pathway Categories:\n\n")
        for cat, count in mouse_categories.items():
            f.write(f"- **{cat}**: {count} pathways\n")
        f.write("\n")
        
        f.write("### Notable Pathways\n\n")
        f.write("1. **Immune System**: MME/Mme is found in multiple immune-related pathways including:\n")
        f.write("   - Innate Immune System\n")
        f.write("   - Neutrophil degranulation\n")
        f.write("   - Hematopoietic cell lineage\n\n")
        
        f.write("2. **Alzheimer's Disease**: MME/Mme is involved in Alzheimer's disease pathways, "
               "consistent with its role as a neprilysin enzyme that degrades amyloid-beta.\n\n")
        
        f.write("3. **Renin-Angiotensin System**: Present in both species, important for cardiovascular "
               "and kidney function.\n\n")
        
        f.write("4. **Protein Metabolism**: Multiple pathways related to peptide hormone metabolism and "
               "protein digestion.\n\n")
        
        f.write("## 5. Conclusions\n\n")
        f.write("- MME/Mme is primarily involved in **immune system** and **protein metabolism** pathways\n")
        f.write("- Strong conservation between human and mouse/rat (most pathways are shared)\n")
        f.write("- Enrichment analysis shows overlap with DE genes, particularly in immune pathways\n")
        f.write("- Not found in MSigDB Hallmark gene sets, but well-represented in PathDIP5 curated pathways\n\n")
        
        f.write("## Data Sources\n\n")
        f.write("- **PathDIP5 Curated**: Curated pathway database integrating KEGG, Reactome, WikiPathways, ACSN2\n")
        f.write("- **MSigDB Hallmark**: Hallmark gene sets from MSigDB (v2025.1)\n")
        f.write("- **Human Dataset**: GSE282344\n")
        f.write("- **Rat Dataset**: GSE208526 (using mouse pathways as proxy)\n")
    
    print(f"✓ Created comprehensive report: {report_file}")
    
    # Create a pathway comparison CSV
    comparison_data = []
    
    all_pathways = human_pathway_names | mouse_pathway_names
    for pathway in sorted(all_pathways):
        in_human = pathway in human_pathway_names
        in_mouse = pathway in mouse_pathway_names
        
        # Get source and category
        if in_human:
            row = human_pathways[human_pathways['pathway_name'] == pathway].iloc[0]
            source = row['source']
            category = row['category']
        else:
            row = mouse_pathways[mouse_pathways['pathway_name'] == pathway].iloc[0]
            source = row['source']
            category = row['category']
        
        comparison_data.append({
            'pathway_name': pathway,
            'source': source,
            'category': category,
            'in_human': in_human,
            'in_mouse_rat': in_mouse,
            'conserved': in_human and in_mouse
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_file = os.path.join(results_dir, 'pathway_species_comparison.csv')
    comparison_df.to_csv(comparison_file, index=False)
    print(f"✓ Created species comparison: {comparison_file}")

if __name__ == '__main__':
    create_summary_report()
    print("\n✓ Summary report generation complete!")
