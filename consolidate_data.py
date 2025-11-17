#!/usr/bin/env python3
"""
Consolidate single-cell data with sample information
Generic script that works with any GSE dataset
"""
import pandas as pd
import os
import sys
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description='Consolidate single-cell sample metadata')
    parser.add_argument('--config', required=True, help='Path to JSON config file with sample information')
    parser.add_argument('--dataset-id', required=True, help='Dataset ID (e.g., GSE208526)')
    parser.add_argument('--output', help='Output CSV file name (default: <dataset_id>_sample_metadata.csv)')
    
    args = parser.parse_args()
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    samples = config.get('samples', {})
    archive_mode = config.get('archive_mode', False)
    file_pattern = config.get('file_pattern', 'standard')  # 'standard' or 'gsm_prefix'
    
    # Create sample metadata file with file locations
    sample_list = []
    for sample_name, metadata in samples.items():
        if archive_mode:
            base_path = f"archives/{sample_name}"
        else:
            base_path = sample_name
        
        if os.path.exists(base_path) or archive_mode:
            sample_entry = {
                'sample_id': sample_name,
                **metadata  # Include all metadata fields
            }
            
            # Add file paths based on pattern
            if file_pattern == 'gsm_prefix':
                sample_entry.update({
                    'barcode_file': f"archives/{sample_name}_barcodes.tsv.gz",
                    'genes_file': f"archives/{sample_name}_genes.tsv.gz",
                    'matrix_file': f"archives/{sample_name}_matrix.mtx.gz"
                })
            else:  # standard 10X format
                sample_entry.update({
                    'barcodes_file': f"{base_path}/barcodes.tsv.gz",
                    'features_file': f"{base_path}/features.tsv.gz",
                    'matrix_file': f"{base_path}/matrix.mtx.gz"
                })
            
            sample_list.append(sample_entry)
    
    # Create DataFrame and save
    sample_metadata_df = pd.DataFrame(sample_list)
    output_file = args.output or f'{args.dataset_id}_sample_metadata.csv'
    sample_metadata_df.to_csv(output_file, index=False)
    
    print(f"Sample metadata saved to {output_file}")
    print(f"Total samples: {len(sample_list)}")
    print("\nNote: This is single-cell RNA-seq data in 10X format (MTX matrices)")
    print("Matrix data is stored in compressed format")
    print("Use Seurat, Scanpy, or similar tools to load the matrix files for analysis")

if __name__ == '__main__':
    main()
