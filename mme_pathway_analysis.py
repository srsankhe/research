#!/usr/bin/env python3
"""
Pathway Analysis for MME/Mme Gene
Identifies pathways containing MME/Mme from MSigDB Hallmark and PathDIP5 curated databases
"""
import pandas as pd
import os
import json
from collections import defaultdict

def load_msigdb_hallmark(filepath, target_gene):
    """Load MSigDB Hallmark pathways and filter for target gene"""
    df = pd.read_csv(filepath)
    
    # Filter for target gene
    gene_pathways = df[df['gene_symbol'].str.upper() == target_gene.upper()].copy()
    
    results = []
    for _, row in gene_pathways.iterrows():
        results.append({
            'database': 'MSigDB Hallmark',
            'pathway_id': row['gs_id'],
            'pathway_name': row['gs_name'],
            'description': row['gs_description'],
            'gene_symbol': row['gene_symbol'],
            'ncbi_gene': row['ncbi_gene'],
            'ensembl_gene': row['ensembl_gene']
        })
    
    return results

def load_pathdip5_curated(filepath, target_gene, is_human=True):
    """Load PathDIP5 curated pathways and filter for target gene"""
    df = pd.read_csv(filepath)
    
    # Filter for target gene
    gene_pathways = df[df['GENESYMBOL'].str.upper() == target_gene.upper()].copy()
    
    results = []
    for _, row in gene_pathways.iterrows():
        pathway_dict = {
            'database': 'PathDIP5 Curated',
            'pathway_name': row['PATHWAY'],
            'source': row['SOURCE'],
            'gene_symbol': row['GENESYMBOL'],
            'category': row.get('CATEGORY', 'N/A'),
            'type': row.get('TYPE', 'N/A')
        }
        
        if is_human:
            pathway_dict.update({
                'uniprot_id': row.get('UPID', 'N/A'),
                'entrez_id': row.get('EGID', 'N/A'),
                'gene_name': row.get('GENENAME', 'N/A')
            })
        else:
            pathway_dict.update({
                'uniprot_id': row.get('UPID', 'N/A'),
                'entrez_id': row.get('EGID', 'N/A'),
                'gene_name': row.get('GENENAME', 'N/A'),
                'gene_symbol_alt': row.get('GENESYMBOL_A', 'N/A')
            })
        
        results.append(pathway_dict)
    
    return results

def get_all_genes_in_pathway(filepath, pathway_name, is_msigdb=True):
    """Get all genes in a specific pathway"""
    df = pd.read_csv(filepath)
    
    if is_msigdb:
        pathway_genes = df[df['gs_name'] == pathway_name]['gene_symbol'].tolist()
    else:
        pathway_genes = df[df['PATHWAY'] == pathway_name]['GENESYMBOL'].tolist()
    
    return pathway_genes

def analyze_pathway_enrichment(de_results_file, pathways_list, pathway_db_file, 
                                is_msigdb=True, padj_cutoff=0.05, log2fc_cutoff=1.0):
    """
    Analyze if pathways containing MME/Mme are enriched in DE results
    """
    # Load DE results
    de_results = pd.read_csv(de_results_file, index_col=0)
    
    # Get significant DE genes
    sig_genes = de_results[
        (de_results['padj'] < padj_cutoff) & 
        (abs(de_results['log2FoldChange']) > log2fc_cutoff)
    ].index.tolist()
    sig_genes_upper = [g.upper() for g in sig_genes]
    
    enrichment_results = []
    
    for pathway_info in pathways_list:
        pathway_name = pathway_info['pathway_name']
        
        # Get all genes in this pathway
        pathway_genes = get_all_genes_in_pathway(pathway_db_file, pathway_name, is_msigdb)
        pathway_genes_upper = [g.upper() for g in pathway_genes]
        
        # Find overlap
        overlap = [g for g in sig_genes_upper if g in pathway_genes_upper]
        
        enrichment_results.append({
            **pathway_info,
            'pathway_size': len(pathway_genes),
            'num_sig_genes': len(sig_genes),
            'num_overlap': len(overlap),
            'overlap_genes': ', '.join(overlap[:20]) if overlap else 'None',
            'percent_pathway_hit': (len(overlap) / len(pathway_genes) * 100) if pathway_genes else 0
        })
    
    return pd.DataFrame(enrichment_results)

def main():
    # Configuration
    pathway_db_dir = '/Users/ssankhe/Documents/Github/research/pathway_db'
    base_dir = '/Users/ssankhe/Documents/Github/research'
    output_dir = os.path.join(base_dir, 'pathway_analysis_results')
    os.makedirs(output_dir, exist_ok=True)
    
    # Target genes
    human_gene = 'MME'
    mouse_gene = 'Mme'  # Will use for rat data as proxy
    
    print("="*60)
    print("MME/Mme Pathway Analysis")
    print("="*60)
    
    # ========== HUMAN ANALYSIS (GSE282344) ==========
    print("\n=== HUMAN ANALYSIS (GSE282344) ===\n")
    
    # 1. MSigDB Hallmark Human
    print("1. Analyzing MSigDB Hallmark Human...")
    hallmark_human_file = os.path.join(pathway_db_dir, 'msigdb/msigdbHallmarkHuman.csv')
    mme_hallmark_human = load_msigdb_hallmark(hallmark_human_file, human_gene)
    print(f"   Found {len(mme_hallmark_human)} MME pathways in MSigDB Hallmark")
    
    # 2. PathDIP5 Curated Human
    print("2. Analyzing PathDIP5 Curated Human...")
    pathdip_human_file = os.path.join(pathway_db_dir, 'pathDip5/pathDip5CuratedHuman.csv')
    mme_pathdip_human = load_pathdip5_curated(pathdip_human_file, human_gene, is_human=True)
    print(f"   Found {len(mme_pathdip_human)} MME pathways in PathDIP5 Curated")
    
    # Save human pathway lists
    if mme_hallmark_human:
        df_h_hall = pd.DataFrame(mme_hallmark_human)
        df_h_hall.to_csv(os.path.join(output_dir, 'human_MME_msigdb_hallmark_pathways.csv'), index=False)
        print(f"   ✓ Saved: human_MME_msigdb_hallmark_pathways.csv")
    
    if mme_pathdip_human:
        df_h_path = pd.DataFrame(mme_pathdip_human)
        df_h_path.to_csv(os.path.join(output_dir, 'human_MME_pathdip5_curated_pathways.csv'), index=False)
        print(f"   ✓ Saved: human_MME_pathdip5_curated_pathways.csv")
    
    # 3. Enrichment analysis with DE results (GSE282344)
    print("\n3. Performing enrichment analysis with GSE282344 DE results...")
    de_human_file = os.path.join(base_dir, 'GSE282344_RAW/deseq2_output/treatment_vs_control_results.csv')
    
    if os.path.exists(de_human_file):
        if mme_hallmark_human:
            enrich_h_hall = analyze_pathway_enrichment(
                de_human_file, mme_hallmark_human, hallmark_human_file, 
                is_msigdb=True
            )
            enrich_h_hall.to_csv(os.path.join(output_dir, 'human_MME_hallmark_enrichment.csv'), index=False)
            print(f"   ✓ Saved: human_MME_hallmark_enrichment.csv")
        
        if mme_pathdip_human:
            enrich_h_path = analyze_pathway_enrichment(
                de_human_file, mme_pathdip_human, pathdip_human_file, 
                is_msigdb=False
            )
            enrich_h_path.to_csv(os.path.join(output_dir, 'human_MME_pathdip5_enrichment.csv'), index=False)
            print(f"   ✓ Saved: human_MME_pathdip5_enrichment.csv")
    
    # ========== MOUSE ANALYSIS (for Rat GSE208526) ==========
    print("\n=== MOUSE ANALYSIS (for Rat GSE208526 as proxy) ===\n")
    
    # 1. MSigDB Hallmark Mouse
    print("1. Analyzing MSigDB Hallmark Mouse...")
    hallmark_mouse_file = os.path.join(pathway_db_dir, 'msigdb/msigdbHallmarkMouse.csv')
    mme_hallmark_mouse = load_msigdb_hallmark(hallmark_mouse_file, mouse_gene)
    print(f"   Found {len(mme_hallmark_mouse)} Mme pathways in MSigDB Hallmark")
    
    # 2. PathDIP5 Curated Mouse
    print("2. Analyzing PathDIP5 Curated Mouse...")
    pathdip_mouse_file = os.path.join(pathway_db_dir, 'pathDip5/pathDip5CuratedMouse.csv')
    mme_pathdip_mouse = load_pathdip5_curated(pathdip_mouse_file, mouse_gene, is_human=False)
    print(f"   Found {len(mme_pathdip_mouse)} Mme pathways in PathDIP5 Curated")
    
    # Save mouse pathway lists
    if mme_hallmark_mouse:
        df_m_hall = pd.DataFrame(mme_hallmark_mouse)
        df_m_hall.to_csv(os.path.join(output_dir, 'mouse_Mme_msigdb_hallmark_pathways.csv'), index=False)
        print(f"   ✓ Saved: mouse_Mme_msigdb_hallmark_pathways.csv")
    
    if mme_pathdip_mouse:
        df_m_path = pd.DataFrame(mme_pathdip_mouse)
        df_m_path.to_csv(os.path.join(output_dir, 'mouse_Mme_pathdip5_curated_pathways.csv'), index=False)
        print(f"   ✓ Saved: mouse_Mme_pathdip5_curated_pathways.csv")
    
    # 3. Enrichment analysis with DE results (GSE208526 - Rat)
    print("\n3. Performing enrichment analysis with GSE208526 (Rat) DE results...")
    de_mouse_file = os.path.join(base_dir, 'GSE208526_RAW/deseq2_output/IRI_vs_sham_results.csv')
    
    if os.path.exists(de_mouse_file):
        if mme_hallmark_mouse:
            enrich_m_hall = analyze_pathway_enrichment(
                de_mouse_file, mme_hallmark_mouse, hallmark_mouse_file, 
                is_msigdb=True
            )
            enrich_m_hall.to_csv(os.path.join(output_dir, 'rat_Mme_hallmark_enrichment.csv'), index=False)
            print(f"   ✓ Saved: rat_Mme_hallmark_enrichment.csv")
        
        if mme_pathdip_mouse:
            enrich_m_path = analyze_pathway_enrichment(
                de_mouse_file, mme_pathdip_mouse, pathdip_mouse_file, 
                is_msigdb=False
            )
            enrich_m_path.to_csv(os.path.join(output_dir, 'rat_Mme_pathdip5_enrichment.csv'), index=False)
            print(f"   ✓ Saved: rat_Mme_pathdip5_enrichment.csv")
    
    # ========== SUMMARY ==========
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    summary = {
        'Human (GSE282344)': {
            'MSigDB Hallmark pathways with MME': len(mme_hallmark_human),
            'PathDIP5 Curated pathways with MME': len(mme_pathdip_human)
        },
        'Mouse/Rat (GSE208526)': {
            'MSigDB Hallmark pathways with Mme': len(mme_hallmark_mouse),
            'PathDIP5 Curated pathways with Mme': len(mme_pathdip_mouse)
        }
    }
    
    for species, counts in summary.items():
        print(f"\n{species}:")
        for db, count in counts.items():
            print(f"  {db}: {count}")
    
    # Save summary
    with open(os.path.join(output_dir, 'pathway_analysis_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ All results saved to: {output_dir}")
    print("="*60)

if __name__ == '__main__':
    main()
