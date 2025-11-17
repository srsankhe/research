#!/usr/bin/env python3
"""
Convert 10X Genomics format (MTX, TSV) to CSV files for DE analysis
Generic script that works with any dataset
"""
import os
import gzip
import pandas as pd
from scipy.io import mmread
import numpy as np
import argparse
import json

def read_mtx_file(filepath):
    """Read MTX file (compressed or uncompressed)"""
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rb') as f:
            matrix = mmread(f)
    else:
        matrix = mmread(filepath)
    return matrix

def read_tsv_file(filepath):
    """Read TSV file (compressed or uncompressed)"""
    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', header=None, compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t', header=None)
    return df

def find_file(base_dir, sample_name, file_type, prefixes):
    """Find file with different possible naming patterns"""
    for prefix in prefixes:
        for ext in ['', '.gz']:
            filepath = os.path.join(base_dir, f"{sample_name}_{prefix}.{file_type}{ext}")
            if os.path.exists(filepath):
                return filepath
            filepath = os.path.join(base_dir, f"{prefix}.{file_type}{ext}")
            if os.path.exists(filepath):
                return filepath
    return None

def convert_sample_to_csv(sample_name, base_dir, archive_dir, output_dir, file_pattern, gene_col_name='gene_name'):
    """Convert one sample to CSV files"""
    print(f"\nProcessing sample: {sample_name}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine data directory
    if archive_dir:
        data_dir = os.path.join(base_dir, archive_dir, sample_name)
    else:
        data_dir = os.path.join(base_dir, sample_name)
    
    # Find files based on pattern
    if file_pattern == 'gsm_prefix':
        # Files like GSM8641663_NR_JK_047_matrix.mtx.gz
        matrix_file = find_file(os.path.join(base_dir, archive_dir) if archive_dir else base_dir, 
                               sample_name, 'mtx', ['matrix'])
        genes_file = find_file(os.path.join(base_dir, archive_dir) if archive_dir else base_dir, 
                              sample_name, 'tsv', ['genes', 'features'])
        barcodes_file = find_file(os.path.join(base_dir, archive_dir) if archive_dir else base_dir, 
                                  sample_name, 'tsv', ['barcodes'])
    else:
        # Standard 10X format in sample directories
        matrix_file = find_file(data_dir, '', 'mtx', ['matrix'])
        genes_file = find_file(data_dir, '', 'tsv', ['features', 'genes'])
        barcodes_file = find_file(data_dir, '', 'tsv', ['barcodes'])
    
    if not all([matrix_file, genes_file, barcodes_file]):
        print(f"  ✗ Could not find all required files for {sample_name}")
        return None
    
    # Read features/genes
    print(f"  Reading features...")
    features_df = read_tsv_file(genes_file)
    
    # Handle different TSV formats
    if features_df.shape[1] >= 3:
        features_df.columns = ['gene_id', 'gene_name', 'feature_type']
    elif features_df.shape[1] == 2:
        features_df.columns = ['gene_name', 'gene_name_dup']
        features_df = features_df[['gene_name']]
    else:
        features_df.columns = ['gene_name']
    
    # Check if this is full transcriptome or just markers
    if len(features_df) < 100:
        print(f"  ⚠ Only {len(features_df)} genes found - likely protein markers, skipping")
        return None
    
    # Save features as CSV
    features_output = os.path.join(output_dir, f'{sample_name}_features.csv')
    features_df.to_csv(features_output, index=False)
    print(f"  ✓ Saved: {features_output}")
    
    # Read barcodes
    print(f"  Reading barcodes...")
    barcodes_df = read_tsv_file(barcodes_file)
    barcodes_df.columns = ['barcode']
    
    # Save barcodes as CSV
    barcodes_output = os.path.join(output_dir, f'{sample_name}_barcodes.csv')
    barcodes_df.to_csv(barcodes_output, index=False)
    print(f"  ✓ Saved: {barcodes_output}")
    
    # Read matrix
    print(f"  Reading matrix (this may take a while)...")
    matrix = read_mtx_file(matrix_file)
    
    # Convert sparse matrix to dense
    print(f"  Converting sparse matrix to dense format...")
    print(f"  Matrix dimensions: {matrix.shape[0]} genes × {matrix.shape[1]} cells")
    
    dense_matrix = matrix.toarray()
    
    # Create DataFrame with gene names as rows and barcodes as columns
    matrix_df = pd.DataFrame(
        dense_matrix,
        index=features_df[gene_col_name if gene_col_name in features_df.columns else 'gene_name'],
        columns=barcodes_df['barcode']
    )
    
    # Save matrix as CSV
    matrix_output = os.path.join(output_dir, f'{sample_name}_matrix.csv')
    print(f"  Saving matrix to CSV...")
    matrix_df.to_csv(matrix_output)
    print(f"  ✓ Saved: {matrix_output}")
    print(f"  Matrix shape: {matrix_df.shape[0]} genes × {matrix_df.shape[1]} cells")
    
    return {
        'sample': sample_name,
        'features_file': features_output,
        'barcodes_file': barcodes_output,
        'matrix_file': matrix_output,
        'n_genes': matrix_df.shape[0],
        'n_cells': matrix_df.shape[1]
    }

def main():
    parser = argparse.ArgumentParser(description='Convert 10X format to CSV')
    parser.add_argument('--config', required=True, help='Path to JSON config file')
    parser.add_argument('--output-dir', default='csv_output', help='Output directory')
    
    args = parser.parse_args()
    
    # Use current directory or config directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(base_dir, args.output_dir)
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    samples = list(config.get('samples', {}).keys())
    archive_dir = config.get('archive_dir', 'archives')
    file_pattern = config.get('file_pattern', 'standard')
    gene_col = config.get('gene_column', 'gene_name')
    
    results = []
    for sample in samples:
        result = convert_sample_to_csv(sample, base_dir, archive_dir, output_dir, 
                                      file_pattern, gene_col)
        if result:
            results.append(result)
    
    # Print summary
    print("\n" + "="*60)
    print("CONVERSION SUMMARY")
    print("="*60)
    
    if not results:
        print("No samples were successfully converted!")
        return
    
    for r in results:
        print(f"\nSample: {r['sample']}")
        print(f"  Genes: {r['n_genes']}")
        print(f"  Cells: {r['n_cells']}")
        print(f"  Output files:")
        print(f"    - {os.path.basename(r['features_file'])}")
        print(f"    - {os.path.basename(r['barcodes_file'])}")
        print(f"    - {os.path.basename(r['matrix_file'])}")
    
    print(f"\nAll CSV files saved to: {output_dir}")

if __name__ == '__main__':
    main()
