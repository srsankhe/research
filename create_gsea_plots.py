#!/usr/bin/env python3
"""
Create GSEA-style enrichment plots for MME/Mme pathway analysis
Includes running enrichment score plots for each pathway
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from matplotlib.patches import Rectangle

def calculate_enrichment_score(pathway_genes, gene_scores):
    """
    Calculate GSEA-style enrichment score
    Similar to GSEA running enrichment score
    """
    # Rank genes by score (absolute value for differential expression)
    ranked_genes = sorted(gene_scores.items(), key=lambda x: abs(x[1]), reverse=True)
    ranked_gene_names = [g[0].upper() for g in ranked_genes]
    
    # Mark pathway gene positions
    pathway_genes_upper = [g.upper() for g in pathway_genes]
    n_pathway = len(pathway_genes_upper)
    n_total = len(ranked_genes)
    
    # Calculate running enrichment score
    es_score = 0
    max_es = 0
    min_es = 0
    es_profile = []
    hit_indices = []
    
    # Weight for hits vs misses
    p = 1  # Weight for correlation
    
    # Calculate weighted sum of pathway genes
    pathway_scores = [abs(gene_scores.get(g, 0)) for g in pathway_genes if g in gene_scores]
    nr = sum([s**p for s in pathway_scores]) if pathway_scores else 1
    
    for i, gene in enumerate(ranked_gene_names):
        if gene in pathway_genes_upper:
            # Hit: gene is in pathway
            score = abs(gene_scores.get(gene, 0))
            es_score += (score**p) / nr
            hit_indices.append(i)
        else:
            # Miss: gene not in pathway
            es_score -= 1 / (n_total - n_pathway)
        
        es_profile.append(es_score)
        max_es = max(max_es, es_score)
        min_es = min(min_es, es_score)
    
    # Use the ES with the largest absolute value
    es = max_es if abs(max_es) > abs(min_es) else min_es
    
    return es, es_profile, hit_indices, ranked_gene_names

def create_gsea_enrichment_plot(pathway_name, pathway_genes, de_results, output_file,
                                gene_col='log2FoldChange', title_suffix=''):
    """
    Create GSEA-style enrichment plot for a single pathway
    """
    # Prepare data
    gene_scores = de_results[gene_col].to_dict()
    
    # Calculate enrichment
    es, es_profile, hit_indices, ranked_genes = calculate_enrichment_score(
        pathway_genes, gene_scores
    )
    
    # Create figure with 3 subplots
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(4, 1, hspace=0.3, height_ratios=[1.5, 0.5, 0.5, 1])
    
    # Plot 1: Running Enrichment Score
    ax1 = fig.add_subplot(gs[0])
    x = np.arange(len(es_profile))
    
    # Color based on positive/negative
    colors = ['green' if es > 0 else 'red' for es in es_profile]
    ax1.plot(x, es_profile, linewidth=2, color='green' if es > 0 else 'red')
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    
    # Mark maximum/minimum
    max_idx = np.argmax(np.abs(es_profile))
    ax1.plot(max_idx, es_profile[max_idx], 'o', markersize=10, 
             color='darkred' if es_profile[max_idx] < 0 else 'darkgreen')
    
    ax1.set_ylabel('Enrichment Score (ES)', fontsize=11, fontweight='bold')
    ax1.set_title(f'{pathway_name}\nES = {es:.3f}', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0, len(es_profile))
    
    # Plot 2: Hit marks (vertical lines where pathway genes appear)
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    for hit_idx in hit_indices:
        ax2.axvline(x=hit_idx, color='black', linewidth=0.5, alpha=0.8)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_ylabel('Hits', fontsize=10, fontweight='bold', rotation=0, labelpad=30)
    ax2.set_xlim(0, len(es_profile))
    
    # Plot 3: Ranking metric (gene scores)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    gene_score_values = [gene_scores.get(g, 0) for g in ranked_genes]
    
    # Color bar based on positive/negative
    colors_bar = ['red' if s < 0 else 'blue' for s in gene_score_values]
    ax3.bar(x, gene_score_values, color=colors_bar, width=1, alpha=0.7, edgecolor='none')
    ax3.axhline(y=0, color='black', linewidth=0.5)
    ax3.set_ylabel('Ranked\nMetric', fontsize=10, fontweight='bold', rotation=0, labelpad=30)
    ax3.set_xlim(0, len(es_profile))
    ax3.grid(axis='y', alpha=0.3)
    
    # Plot 4: Gene rank
    ax4 = fig.add_subplot(gs[3], sharex=ax1)
    ax4.plot(x, x, color='gray', linewidth=1)
    ax4.fill_between(x, 0, x, alpha=0.2, color='gray')
    ax4.set_xlabel('Gene Rank', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Rank', fontsize=10, fontweight='bold')
    ax4.set_xlim(0, len(es_profile))
    ax4.grid(alpha=0.3)
    
    # Add text annotations
    text_str = f'Pathway Size: {len(pathway_genes)}\nHits: {len(hit_indices)}\nHit Rate: {len(hit_indices)/len(pathway_genes)*100:.1f}%'
    ax1.text(0.98, 0.02, text_str, transform=ax1.transAxes,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=9)
    
    plt.suptitle(f'GSEA Enrichment Plot{title_suffix}', fontsize=14, fontweight='bold', y=0.995)
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_all_gsea_plots():
    """Create GSEA-style plots for all pathways"""
    
    results_dir = '/Users/ssankhe/Documents/Github/research/pathway_analysis_results'
    plots_dir = os.path.join(results_dir, 'plots', 'gsea_style')
    os.makedirs(plots_dir, exist_ok=True)
    
    print("Creating GSEA-style enrichment plots...\n")
    
    # Load pathway databases
    pathway_db_dir = '/Users/ssankhe/Documents/Github/research/pathway_db'
    base_dir = '/Users/ssankhe/Documents/Github/research'
    
    # Load pathway lists
    human_pathways = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_curated_pathways.csv'))
    mouse_pathways = pd.read_csv(os.path.join(results_dir, 'mouse_Mme_pathdip5_curated_pathways.csv'))
    
    # Load enrichment results
    human_enrichment = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_enrichment.csv'))
    rat_enrichment = pd.read_csv(os.path.join(results_dir, 'rat_Mme_pathdip5_enrichment.csv'))
    
    # Load PathDIP5 databases
    pathdip_human_file = os.path.join(pathway_db_dir, 'pathDip5/pathDip5CuratedHuman.csv')
    pathdip_mouse_file = os.path.join(pathway_db_dir, 'pathDip5/pathDip5CuratedMouse.csv')
    
    pathdip_human_db = pd.read_csv(pathdip_human_file)
    pathdip_mouse_db = pd.read_csv(pathdip_mouse_file)
    
    # Load DE results
    de_human = pd.read_csv(os.path.join(base_dir, 'GSE282344_RAW/deseq2_output/treatment_vs_control_results.csv'),
                           index_col=0)
    de_rat = pd.read_csv(os.path.join(base_dir, 'GSE208526_RAW/deseq2_output/IRI_vs_sham_results.csv'),
                         index_col=0)
    
    # Get top 10 enriched pathways for each species
    top_human = human_enrichment.nlargest(10, 'percent_pathway_hit')
    top_rat = rat_enrichment.nlargest(10, 'percent_pathway_hit')
    
    # Create GSEA plots for top human pathways
    print("Creating plots for top 10 human pathways...")
    for idx, (_, row) in enumerate(top_human.iterrows(), 1):
        pathway_name = row['pathway_name']
        print(f"  {idx}/10: {pathway_name[:60]}...")
        
        # Get all genes in this pathway
        pathway_genes = pathdip_human_db[
            pathdip_human_db['PATHWAY'] == pathway_name
        ]['GENESYMBOL'].tolist()
        
        if len(pathway_genes) > 0:
            output_file = os.path.join(plots_dir, f'human_{idx:02d}_{pathway_name[:50].replace("/", "_")}.png')
            
            try:
                create_gsea_enrichment_plot(
                    pathway_name, pathway_genes, de_human, output_file,
                    title_suffix=' - Human (GSE282344)'
                )
                print(f"     ✓ Saved")
            except Exception as e:
                print(f"     ✗ Error: {e}")
    
    # Create GSEA plots for top rat pathways
    print("\nCreating plots for top 10 rat pathways...")
    for idx, (_, row) in enumerate(top_rat.iterrows(), 1):
        pathway_name = row['pathway_name']
        print(f"  {idx}/10: {pathway_name[:60]}...")
        
        # Get all genes in this pathway
        pathway_genes = pathdip_mouse_db[
            pathdip_mouse_db['PATHWAY'] == pathway_name
        ]['GENESYMBOL'].tolist()
        
        if len(pathway_genes) > 0:
            output_file = os.path.join(plots_dir, f'rat_{idx:02d}_{pathway_name[:50].replace("/", "_")}.png')
            
            try:
                create_gsea_enrichment_plot(
                    pathway_name, pathway_genes, de_rat, output_file,
                    title_suffix=' - Rat (GSE208526)'
                )
                print(f"     ✓ Saved")
            except Exception as e:
                print(f"     ✗ Error: {e}")
    
    print(f"\n✓ All GSEA-style plots saved to: {plots_dir}/")
    print(f"✓ Total GSEA plots created: 20")

if __name__ == '__main__':
    create_all_gsea_plots()
    print("\n✓ GSEA-style enrichment plot generation complete!")
