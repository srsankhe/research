#!/usr/bin/env python3
"""
Create enrichment visualization plots for MME/Mme pathway analysis
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def create_enrichment_plots():
    results_dir = '/Users/ssankhe/Documents/Github/research/pathway_analysis_results'
    plots_dir = os.path.join(results_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 10
    
    print("Creating enrichment visualization plots...\n")
    
    # Load data
    human_enrichment = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_enrichment.csv'))
    rat_enrichment = pd.read_csv(os.path.join(results_dir, 'rat_Mme_pathdip5_enrichment.csv'))
    human_pathways = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_curated_pathways.csv'))
    mouse_pathways = pd.read_csv(os.path.join(results_dir, 'mouse_Mme_pathdip5_curated_pathways.csv'))
    
    # ========== PLOT 1: Pathway Category Distribution ==========
    print("1. Creating pathway category distribution plots...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Human categories
    human_cat = human_pathways['category'].value_counts()
    colors1 = sns.color_palette("husl", len(human_cat))
    ax1.barh(range(len(human_cat)), human_cat.values, color=colors1)
    ax1.set_yticks(range(len(human_cat)))
    ax1.set_yticklabels(human_cat.index, fontsize=9)
    ax1.set_xlabel('Number of Pathways', fontsize=11)
    ax1.set_title('Human MME Pathway Categories', fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Mouse categories
    mouse_cat = mouse_pathways['category'].value_counts()
    colors2 = sns.color_palette("husl", len(mouse_cat))
    ax2.barh(range(len(mouse_cat)), mouse_cat.values, color=colors2)
    ax2.set_yticks(range(len(mouse_cat)))
    ax2.set_yticklabels(mouse_cat.index, fontsize=9)
    ax2.set_xlabel('Number of Pathways', fontsize=11)
    ax2.set_title('Mouse/Rat Mme Pathway Categories', fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '01_pathway_categories.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 01_pathway_categories.png")
    
    # ========== PLOT 2: Pathway Source Distribution ==========
    print("2. Creating pathway source distribution plots...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Human sources
    human_src = human_pathways['source'].value_counts()
    ax1.pie(human_src.values, labels=human_src.index, autopct='%1.1f%%',
            colors=sns.color_palette("Set2"), startangle=90)
    ax1.set_title('Human MME - Pathway Sources', fontsize=12, fontweight='bold')
    
    # Mouse sources
    mouse_src = mouse_pathways['source'].value_counts()
    ax2.pie(mouse_src.values, labels=mouse_src.index, autopct='%1.1f%%',
            colors=sns.color_palette("Set2"), startangle=90)
    ax2.set_title('Mouse/Rat Mme - Pathway Sources', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '02_pathway_sources.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 02_pathway_sources.png")
    
    # ========== PLOT 3: Enrichment Bar Plots (Top 15) ==========
    print("3. Creating enrichment bar plots (top 15 pathways)...")
    
    # Human enrichment
    fig, ax = plt.subplots(figsize=(12, 8))
    top_human = human_enrichment.nlargest(15, 'percent_pathway_hit')
    
    colors = ['#e74c3c' if x > 5 else '#3498db' if x > 2 else '#95a5a6' 
              for x in top_human['percent_pathway_hit']]
    
    y_pos = np.arange(len(top_human))
    ax.barh(y_pos, top_human['percent_pathway_hit'], color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{p[:50]}..." if len(p) > 50 else p 
                        for p in top_human['pathway_name']], fontsize=9)
    ax.set_xlabel('Pathway Hit Rate (%)', fontsize=11)
    ax.set_title('Human MME - Top 15 Enriched Pathways (GSE282344)', 
                 fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    # Add number of overlapping genes as text
    for i, (idx, row) in enumerate(top_human.iterrows()):
        ax.text(row['percent_pathway_hit'] + 0.2, i, 
                f"{row['num_overlap']} genes", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '03_human_enrichment_top15.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 03_human_enrichment_top15.png")
    
    # Rat enrichment
    fig, ax = plt.subplots(figsize=(12, 8))
    top_rat = rat_enrichment.nlargest(15, 'percent_pathway_hit')
    
    colors = ['#e74c3c' if x > 5 else '#3498db' if x > 2 else '#95a5a6' 
              for x in top_rat['percent_pathway_hit']]
    
    y_pos = np.arange(len(top_rat))
    ax.barh(y_pos, top_rat['percent_pathway_hit'], color=colors, alpha=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{p[:50]}..." if len(p) > 50 else p 
                        for p in top_rat['pathway_name']], fontsize=9)
    ax.set_xlabel('Pathway Hit Rate (%)', fontsize=11)
    ax.set_title('Rat Mme - Top 15 Enriched Pathways (GSE208526)', 
                 fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    for i, (idx, row) in enumerate(top_rat.iterrows()):
        ax.text(row['percent_pathway_hit'] + 0.2, i, 
                f"{row['num_overlap']} genes", 
                va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '04_rat_enrichment_top15.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 04_rat_enrichment_top15.png")
    
    # ========== PLOT 4: Scatter Plot - Pathway Size vs Enrichment ==========
    print("4. Creating scatter plots (pathway size vs enrichment)...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Human scatter
    scatter1 = ax1.scatter(human_enrichment['pathway_size'], 
                          human_enrichment['percent_pathway_hit'],
                          c=human_enrichment['num_overlap'],
                          s=100, alpha=0.6, cmap='viridis', edgecolors='black', linewidth=0.5)
    ax1.set_xlabel('Pathway Size (number of genes)', fontsize=11)
    ax1.set_ylabel('Pathway Hit Rate (%)', fontsize=11)
    ax1.set_title('Human MME - Pathway Size vs Enrichment', fontsize=12, fontweight='bold')
    ax1.set_xscale('log')
    ax1.grid(alpha=0.3)
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label('DE Genes Overlap', fontsize=10)
    
    # Annotate top 3
    top3_human = human_enrichment.nlargest(3, 'percent_pathway_hit')
    for _, row in top3_human.iterrows():
        name = row['pathway_name'][:20] + '...' if len(row['pathway_name']) > 20 else row['pathway_name']
        ax1.annotate(name, 
                    xy=(row['pathway_size'], row['percent_pathway_hit']),
                    xytext=(10, 10), textcoords='offset points',
                    fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    # Rat scatter
    scatter2 = ax2.scatter(rat_enrichment['pathway_size'], 
                          rat_enrichment['percent_pathway_hit'],
                          c=rat_enrichment['num_overlap'],
                          s=100, alpha=0.6, cmap='viridis', edgecolors='black', linewidth=0.5)
    ax2.set_xlabel('Pathway Size (number of genes)', fontsize=11)
    ax2.set_ylabel('Pathway Hit Rate (%)', fontsize=11)
    ax2.set_title('Rat Mme - Pathway Size vs Enrichment', fontsize=12, fontweight='bold')
    ax2.set_xscale('log')
    ax2.grid(alpha=0.3)
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label('DE Genes Overlap', fontsize=10)
    
    # Annotate top 3
    top3_rat = rat_enrichment.nlargest(3, 'percent_pathway_hit')
    for _, row in top3_rat.iterrows():
        name = row['pathway_name'][:20] + '...' if len(row['pathway_name']) > 20 else row['pathway_name']
        ax2.annotate(name,
                    xy=(row['pathway_size'], row['percent_pathway_hit']),
                    xytext=(10, 10), textcoords='offset points',
                    fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '05_size_vs_enrichment_scatter.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 05_size_vs_enrichment_scatter.png")
    
    # ========== PLOT 5: Heatmap of Top Pathways by Category ==========
    print("5. Creating category-based enrichment heatmap...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Human heatmap
    human_top = human_enrichment.nlargest(20, 'percent_pathway_hit')
    heatmap_data_h = human_top.pivot_table(
        values='percent_pathway_hit',
        index='pathway_name',
        columns='source',
        fill_value=0
    )
    
    sns.heatmap(heatmap_data_h, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Pathway Hit Rate (%)'}, ax=ax1,
                linewidths=0.5, linecolor='gray')
    ax1.set_title('Human MME - Top 20 Pathways by Source', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Pathway Source', fontsize=11)
    ax1.set_ylabel('Pathway Name', fontsize=11)
    ax1.set_yticklabels([label.get_text()[:40] + '...' if len(label.get_text()) > 40 
                         else label.get_text() for label in ax1.get_yticklabels()], 
                        fontsize=8)
    
    # Rat heatmap
    rat_top = rat_enrichment.nlargest(20, 'percent_pathway_hit')
    heatmap_data_r = rat_top.pivot_table(
        values='percent_pathway_hit',
        index='pathway_name',
        columns='source',
        fill_value=0
    )
    
    sns.heatmap(heatmap_data_r, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Pathway Hit Rate (%)'}, ax=ax2,
                linewidths=0.5, linecolor='gray')
    ax2.set_title('Rat Mme - Top 20 Pathways by Source', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Pathway Source', fontsize=11)
    ax2.set_ylabel('Pathway Name', fontsize=11)
    ax2.set_yticklabels([label.get_text()[:40] + '...' if len(label.get_text()) > 40 
                         else label.get_text() for label in ax2.get_yticklabels()], 
                        fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '06_enrichment_heatmap_by_source.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 06_enrichment_heatmap_by_source.png")
    
    # ========== PLOT 6: Dot Plot - All Pathways ==========
    print("6. Creating comprehensive dot plot...")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
    
    # Human dot plot
    human_sorted = human_enrichment.sort_values('percent_pathway_hit', ascending=True)
    
    scatter1 = ax1.scatter(human_sorted['percent_pathway_hit'],
                          range(len(human_sorted)),
                          s=human_sorted['num_overlap'] * 20,
                          c=human_sorted['pathway_size'],
                          cmap='coolwarm', alpha=0.7, edgecolors='black', linewidth=0.5)
    
    ax1.set_yticks(range(len(human_sorted)))
    ax1.set_yticklabels([f"{p[:45]}..." if len(p) > 45 else p 
                         for p in human_sorted['pathway_name']], fontsize=8)
    ax1.set_xlabel('Pathway Hit Rate (%)', fontsize=11)
    ax1.set_title('Human MME - All Pathway Enrichment', fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label('Pathway Size', fontsize=10)
    
    # Rat dot plot
    rat_sorted = rat_enrichment.sort_values('percent_pathway_hit', ascending=True)
    
    scatter2 = ax2.scatter(rat_sorted['percent_pathway_hit'],
                          range(len(rat_sorted)),
                          s=rat_sorted['num_overlap'] * 20,
                          c=rat_sorted['pathway_size'],
                          cmap='coolwarm', alpha=0.7, edgecolors='black', linewidth=0.5)
    
    ax2.set_yticks(range(len(rat_sorted)))
    ax2.set_yticklabels([f"{p[:45]}..." if len(p) > 45 else p 
                         for p in rat_sorted['pathway_name']], fontsize=8)
    ax2.set_xlabel('Pathway Hit Rate (%)', fontsize=11)
    ax2.set_title('Rat Mme - All Pathway Enrichment', fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label('Pathway Size', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '07_all_pathways_dotplot.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 07_all_pathways_dotplot.png")
    
    # ========== PLOT 7: Comparison Between Species ==========
    print("7. Creating species comparison plot...")
    
    # Match pathways between species
    comparison_df = pd.read_csv(os.path.join(results_dir, 'pathway_species_comparison.csv'))
    conserved = comparison_df[comparison_df['conserved'] == True]['pathway_name'].tolist()
    
    # Get enrichment values for conserved pathways
    human_conserved = human_enrichment[human_enrichment['pathway_name'].isin(conserved)].copy()
    rat_conserved = rat_enrichment[rat_enrichment['pathway_name'].isin(conserved)].copy()
    
    # Merge
    comparison = human_conserved[['pathway_name', 'percent_pathway_hit', 'num_overlap']].merge(
        rat_conserved[['pathway_name', 'percent_pathway_hit', 'num_overlap']],
        on='pathway_name', suffixes=('_human', '_rat')
    )
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    scatter = ax.scatter(comparison['percent_pathway_hit_human'],
                        comparison['percent_pathway_hit_rat'],
                        s=100, alpha=0.6, c=range(len(comparison)),
                        cmap='viridis', edgecolors='black', linewidth=0.5)
    
    # Add diagonal line
    max_val = max(comparison['percent_pathway_hit_human'].max(), 
                  comparison['percent_pathway_hit_rat'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('Human Enrichment (%)', fontsize=11)
    ax.set_ylabel('Rat Enrichment (%)', fontsize=11)
    ax.set_title('Species Comparison - Conserved Pathway Enrichment', 
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Annotate interesting points
    for _, row in comparison.iterrows():
        if abs(row['percent_pathway_hit_human'] - row['percent_pathway_hit_rat']) > 3:
            name = row['pathway_name'][:20] + '...' if len(row['pathway_name']) > 20 else row['pathway_name']
            ax.annotate(name,
                       xy=(row['percent_pathway_hit_human'], row['percent_pathway_hit_rat']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=7, bbox=dict(boxstyle='round,pad=0.3', 
                                            facecolor='yellow', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, '08_species_comparison.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 08_species_comparison.png")
    
    # ========== PLOT 8: Summary Statistics ==========
    print("8. Creating summary statistics panel...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Panel 1: Total pathways
    ax1 = fig.add_subplot(gs[0, 0])
    pathways_data = [len(human_pathways), len(mouse_pathways)]
    ax1.bar(['Human MME', 'Mouse Mme'], pathways_data, color=['#3498db', '#e74c3c'])
    ax1.set_ylabel('Number of Pathways', fontsize=10)
    ax1.set_title('Total Pathways Containing MME/Mme', fontsize=11, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    for i, v in enumerate(pathways_data):
        ax1.text(i, v + 0.5, str(v), ha='center', fontweight='bold')
    
    # Panel 2: Conservation
    ax2 = fig.add_subplot(gs[0, 1])
    conservation = [len(comparison_df[comparison_df['conserved'] == True]),
                   len(comparison_df[comparison_df['conserved'] == False])]
    colors = ['#2ecc71', '#95a5a6']
    ax2.pie(conservation, labels=['Conserved', 'Species-specific'], autopct='%1.1f%%',
           colors=colors, startangle=90)
    ax2.set_title('Pathway Conservation', fontsize=11, fontweight='bold')
    
    # Panel 3: Average enrichment
    ax3 = fig.add_subplot(gs[0, 2])
    avg_enrich = [human_enrichment['percent_pathway_hit'].mean(),
                  rat_enrichment['percent_pathway_hit'].mean()]
    ax3.bar(['Human', 'Rat'], avg_enrich, color=['#3498db', '#e74c3c'])
    ax3.set_ylabel('Average Hit Rate (%)', fontsize=10)
    ax3.set_title('Average Pathway Enrichment', fontsize=11, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    for i, v in enumerate(avg_enrich):
        ax3.text(i, v + 0.1, f'{v:.2f}%', ha='center')
    
    # Panel 4: Top categories human
    ax4 = fig.add_subplot(gs[1, :])
    top_cat_human = human_enrichment.groupby('category')['percent_pathway_hit'].mean().nlargest(7)
    ax4.barh(range(len(top_cat_human)), top_cat_human.values, 
            color=sns.color_palette("husl", len(top_cat_human)))
    ax4.set_yticks(range(len(top_cat_human)))
    ax4.set_yticklabels(top_cat_human.index, fontsize=9)
    ax4.set_xlabel('Average Pathway Hit Rate (%)', fontsize=10)
    ax4.set_title('Human - Average Enrichment by Category', fontsize=11, fontweight='bold')
    ax4.grid(axis='x', alpha=0.3)
    
    # Panel 5: Top categories rat
    ax5 = fig.add_subplot(gs[2, :])
    top_cat_rat = rat_enrichment.groupby('category')['percent_pathway_hit'].mean().nlargest(7)
    ax5.barh(range(len(top_cat_rat)), top_cat_rat.values,
            color=sns.color_palette("husl", len(top_cat_rat)))
    ax5.set_yticks(range(len(top_cat_rat)))
    ax5.set_yticklabels(top_cat_rat.index, fontsize=9)
    ax5.set_xlabel('Average Pathway Hit Rate (%)', fontsize=10)
    ax5.set_title('Rat - Average Enrichment by Category', fontsize=11, fontweight='bold')
    ax5.grid(axis='x', alpha=0.3)
    
    plt.suptitle('MME/Mme Pathway Analysis - Summary Statistics', 
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.savefig(os.path.join(plots_dir, '09_summary_statistics.png'), bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 09_summary_statistics.png")
    
    print(f"\n✓ All plots saved to: {plots_dir}/")
    print(f"✓ Total plots created: 9")
    
    # Create plot index
    plot_index = """# Enrichment Plots Index

## Generated Plots

1. **01_pathway_categories.png**
   - Bar charts showing pathway distribution by category for human and mouse

2. **02_pathway_sources.png**
   - Pie charts showing pathway source distribution (KEGG, Reactome, etc.)

3. **03_human_enrichment_top15.png**
   - Top 15 enriched pathways in human (GSE282344) with hit rates

4. **04_rat_enrichment_top15.png**
   - Top 15 enriched pathways in rat (GSE208526) with hit rates

5. **05_size_vs_enrichment_scatter.png**
   - Scatter plots showing relationship between pathway size and enrichment

6. **06_enrichment_heatmap_by_source.png**
   - Heatmaps of top 20 pathways organized by source database

7. **07_all_pathways_dotplot.png**
   - Comprehensive dot plots showing all pathway enrichments

8. **08_species_comparison.png**
   - Scatter plot comparing enrichment between human and rat for conserved pathways

9. **09_summary_statistics.png**
   - Multi-panel summary with key statistics and comparisons

## Color Coding

- Red (>5%): High enrichment
- Blue (2-5%): Moderate enrichment  
- Gray (<2%): Low enrichment

## Dot/Scatter Plot Legend

- **Size**: Number of overlapping DE genes
- **Color**: Pathway size (number of genes in pathway)
"""
    
    with open(os.path.join(plots_dir, 'PLOTS_README.md'), 'w') as f:
        f.write(plot_index)
    
    print("✓ Created: PLOTS_README.md")

if __name__ == '__main__':
    create_enrichment_plots()
    print("\n✓ Enrichment plot generation complete!")
