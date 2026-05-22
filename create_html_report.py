#!/usr/bin/env python3
"""
Generate comprehensive HTML report for MME/Mme pathway analysis
Includes all plots, data tables, and GSEA-style enrichment plots
"""
import pandas as pd
import os
import base64
from pathlib import Path

def encode_image(image_path):
    """Encode image to base64 for embedding in HTML"""
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode()

def create_html_report():
    results_dir = '/Users/ssankhe/Documents/Github/research/pathway_analysis_results'
    plots_dir = os.path.join(results_dir, 'plots')
    gsea_plots_dir = os.path.join(plots_dir, 'gsea_style')
    
    # Load data
    human_pathways = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_curated_pathways.csv'))
    human_enrichment = pd.read_csv(os.path.join(results_dir, 'human_MME_pathdip5_enrichment.csv'))
    mouse_pathways = pd.read_csv(os.path.join(results_dir, 'mouse_Mme_pathdip5_curated_pathways.csv'))
    rat_enrichment = pd.read_csv(os.path.join(results_dir, 'rat_Mme_pathdip5_enrichment.csv'))
    comparison = pd.read_csv(os.path.join(results_dir, 'pathway_species_comparison.csv'))
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MME/Mme Pathway Analysis - Complete Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }}
        .header p {{
            margin: 5px 0;
            font-size: 1.1em;
            opacity: 0.9;
        }}
        .summary-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .card {{
            background: white;
            padding: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            border-left: 4px solid #667eea;
        }}
        .card h3 {{
            margin-top: 0;
            color: #667eea;
            font-size: 1.1em;
        }}
        .card .number {{
            font-size: 2.5em;
            font-weight: bold;
            color: #333;
            margin: 10px 0;
        }}
        .section {{
            background: white;
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .section h3 {{
            color: #764ba2;
            margin-top: 25px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.9em;
        }}
        th {{
            background-color: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ddd;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }}
        .plot-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .plot-container h4 {{
            margin-top: 0;
            color: #764ba2;
            font-size: 1em;
        }}
        .plot-container img {{
            width: 100%;
            height: auto;
            border-radius: 4px;
        }}
        .gsea-plot {{
            margin: 20px 0;
            padding: 15px;
            background: #f9f9f9;
            border-radius: 8px;
        }}
        .gsea-plot img {{
            width: 100%;
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
            margin-right: 5px;
        }}
        .badge-high {{ background-color: #e74c3c; color: white; }}
        .badge-medium {{ background-color: #3498db; color: white; }}
        .badge-low {{ background-color: #95a5a6; color: white; }}
        .badge-conserved {{ background-color: #2ecc71; color: white; }}
        .toc {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        .toc h2 {{
            color: #667eea;
            margin-top: 0;
        }}
        .toc ul {{
            list-style: none;
            padding-left: 0;
        }}
        .toc li {{
            padding: 8px 0;
            border-bottom: 1px solid #eee;
        }}
        .toc a {{
            color: #667eea;
            text-decoration: none;
            font-weight: 500;
        }}
        .toc a:hover {{
            text-decoration: underline;
        }}
        .highlight {{
            background-color: #fff3cd;
            padding: 15px;
            border-left: 4px solid #ffc107;
            border-radius: 4px;
            margin: 15px 0;
        }}
        .footer {{
            text-align: center;
            padding: 30px;
            color: #666;
            border-top: 2px solid #ddd;
            margin-top: 50px;
        }}
    </style>
</head>
<body>
"""

    # Header
    html_content += """
    <div class="header">
        <h1>🧬 MME/Mme Pathway Analysis</h1>
        <p><strong>Complete Report with GSEA-Style Enrichment Plots</strong></p>
        <p>Analysis Date: November 17, 2024</p>
        <p>Datasets: GSE282344 (Human) | GSE208526 (Rat/Mouse)</p>
    </div>
"""

    # Summary Cards
    html_content += f"""
    <div class="summary-cards">
        <div class="card">
            <h3>Total Pathways</h3>
            <div class="number">{len(human_pathways) + len(mouse_pathways)}</div>
            <p>{len(human_pathways)} Human | {len(mouse_pathways)} Mouse/Rat</p>
        </div>
        <div class="card">
            <h3>Conserved Pathways</h3>
            <div class="number">{len(comparison[comparison['conserved'] == True])}</div>
            <p>{len(comparison[comparison['conserved'] == True])/len(comparison)*100:.0f}% Conservation Rate</p>
        </div>
        <div class="card">
            <h3>Databases Searched</h3>
            <div class="number">2</div>
            <p>MSigDB Hallmark (0 hits)<br>PathDIP5 Curated ({len(human_pathways)} hits)</p>
        </div>
        <div class="card">
            <h3>Visualizations</h3>
            <div class="number">29</div>
            <p>9 Summary plots<br>20 GSEA-style plots</p>
        </div>
    </div>
"""

    # Table of Contents
    html_content += """
    <div class="toc">
        <h2>📑 Table of Contents</h2>
        <ul>
            <li><a href="#executive-summary">Executive Summary</a></li>
            <li><a href="#pathways">Pathway Lists</a></li>
            <li><a href="#enrichment">Enrichment Analysis</a></li>
            <li><a href="#visualizations">Standard Visualizations</a></li>
            <li><a href="#gsea-plots">GSEA-Style Enrichment Plots</a></li>
            <li><a href="#species-comparison">Species Comparison</a></li>
            <li><a href="#conclusions">Conclusions</a></li>
        </ul>
    </div>
"""

    # Executive Summary
    html_content += f"""
    <div class="section" id="executive-summary">
        <h2>📊 Executive Summary</h2>
        
        <div class="highlight">
            <strong>Key Finding:</strong> MME/Mme (Neprilysin/CD10) is a highly conserved zinc-dependent 
            metalloproteinase primarily involved in <strong>immune system function</strong> and 
            <strong>protein metabolism</strong>.
        </div>
        
        <h3>Pathway Discovery</h3>
        <ul>
            <li><strong>MSigDB Hallmark:</strong> ❌ MME/Mme not found (0 pathways)</li>
            <li><strong>PathDIP5 Curated:</strong> ✅ Found in {len(human_pathways)} human and {len(mouse_pathways)} mouse pathways</li>
            <li><strong>Conservation:</strong> {len(comparison[comparison['conserved'] == True])} pathways conserved ({len(comparison[comparison['conserved'] == True])/len(comparison)*100:.1f}%)</li>
        </ul>
        
        <h3>Top Functional Categories</h3>
        <ol>
            <li><strong>Immune System</strong> - 5 pathways (neutrophils, innate immunity)</li>
            <li><strong>Protein Metabolism</strong> - 3 pathways (peptide degradation)</li>
            <li><strong>Alzheimer's Disease</strong> - 3 pathways (amyloid-β degradation)</li>
            <li><strong>Cardiovascular</strong> - 3 pathways (cardiac conduction, muscle)</li>
            <li><strong>Kidney Disease</strong> - 1 pathway (FSGS)</li>
        </ol>
    </div>
"""

    # Pathway Lists
    html_content += f"""
    <div class="section" id="pathways">
        <h2>🔬 Pathway Lists</h2>
        
        <h3>Human MME Pathways ({len(human_pathways)})</h3>
        <table>
            <thead>
                <tr>
                    <th>Pathway Name</th>
                    <th>Source</th>
                    <th>Category</th>
                    <th>Conservation</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for _, row in human_pathways.iterrows():
        pathway_name = row['pathway_name']
        is_conserved = pathway_name in comparison[comparison['conserved'] == True]['pathway_name'].values
        conserved_badge = '<span class="badge badge-conserved">Conserved</span>' if is_conserved else ''
        html_content += f"""
                <tr>
                    <td><strong>{row['pathway_name']}</strong></td>
                    <td>{row['source']}</td>
                    <td>{row['category']}</td>
                    <td>{conserved_badge}</td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
    </div>
"""

    # Enrichment Results
    html_content += f"""
    <div class="section" id="enrichment">
        <h2>📈 Enrichment Analysis</h2>
        
        <h3>Human (GSE282344) - Top 10 Enriched Pathways</h3>
        <table>
            <thead>
                <tr>
                    <th>Pathway</th>
                    <th>Source</th>
                    <th>Pathway Size</th>
                    <th>DE Genes</th>
                    <th>Hit Rate</th>
                    <th>Top Genes</th>
                </tr>
            </thead>
            <tbody>
"""
    
    top_human = human_enrichment.nlargest(10, 'percent_pathway_hit')
    for _, row in top_human.iterrows():
        hit_rate = row['percent_pathway_hit']
        badge = 'badge-high' if hit_rate > 5 else 'badge-medium' if hit_rate > 2 else 'badge-low'
        genes = row['overlap_genes'][:100] + '...' if len(str(row['overlap_genes'])) > 100 else row['overlap_genes']
        html_content += f"""
                <tr>
                    <td><strong>{row['pathway_name']}</strong></td>
                    <td>{row['source']}</td>
                    <td>{row['pathway_size']}</td>
                    <td>{row['num_overlap']}</td>
                    <td><span class="badge {badge}">{hit_rate:.1f}%</span></td>
                    <td><small>{genes}</small></td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
        
        <h3>Rat (GSE208526) - Top 10 Enriched Pathways</h3>
        <table>
            <thead>
                <tr>
                    <th>Pathway</th>
                    <th>Source</th>
                    <th>Pathway Size</th>
                    <th>DE Genes</th>
                    <th>Hit Rate</th>
                    <th>Top Genes</th>
                </tr>
            </thead>
            <tbody>
"""
    
    top_rat = rat_enrichment.nlargest(10, 'percent_pathway_hit')
    for _, row in top_rat.iterrows():
        hit_rate = row['percent_pathway_hit']
        badge = 'badge-high' if hit_rate > 5 else 'badge-medium' if hit_rate > 2 else 'badge-low'
        genes = row['overlap_genes'][:100] + '...' if len(str(row['overlap_genes'])) > 100 else row['overlap_genes']
        html_content += f"""
                <tr>
                    <td><strong>{row['pathway_name']}</strong></td>
                    <td>{row['source']}</td>
                    <td>{row['pathway_size']}</td>
                    <td>{row['num_overlap']}</td>
                    <td><span class="badge {badge}">{hit_rate:.1f}%</span></td>
                    <td><small>{genes}</small></td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
    </div>
"""

    # Standard Visualizations
    html_content += """
    <div class="section" id="visualizations">
        <h2>📊 Standard Visualizations</h2>
        <div class="plot-grid">
"""
    
    standard_plots = [
        ('01_pathway_categories.png', 'Pathway Categories'),
        ('02_pathway_sources.png', 'Pathway Sources'),
        ('03_human_enrichment_top15.png', 'Human Top 15 Enriched'),
        ('04_rat_enrichment_top15.png', 'Rat Top 15 Enriched'),
        ('05_size_vs_enrichment_scatter.png', 'Size vs Enrichment'),
        ('06_enrichment_heatmap_by_source.png', 'Enrichment Heatmap'),
        ('07_all_pathways_dotplot.png', 'All Pathways Dotplot'),
        ('08_species_comparison.png', 'Species Comparison'),
        ('09_summary_statistics.png', 'Summary Statistics')
    ]
    
    for plot_file, plot_title in standard_plots:
        plot_path = os.path.join(plots_dir, plot_file)
        if os.path.exists(plot_path):
            img_data = encode_image(plot_path)
            html_content += f"""
            <div class="plot-container">
                <h4>{plot_title}</h4>
                <img src="data:image/png;base64,{img_data}" alt="{plot_title}">
            </div>
"""
    
    html_content += """
        </div>
    </div>
"""

    # GSEA-Style Plots
    html_content += """
    <div class="section" id="gsea-plots">
        <h2>🎯 GSEA-Style Enrichment Plots</h2>
        <p>These plots show running enrichment scores similar to Gene Set Enrichment Analysis (GSEA), 
        displaying how pathway genes are distributed across the ranked gene list.</p>
        
        <h3>Human (GSE282344) - Top 10 Pathways</h3>
"""
    
    # Human GSEA plots
    human_gsea_files = sorted([f for f in os.listdir(gsea_plots_dir) if f.startswith('human_')])
    for gsea_file in human_gsea_files:
        plot_path = os.path.join(gsea_plots_dir, gsea_file)
        img_data = encode_image(plot_path)
        pathway_name = gsea_file.replace('human_', '').replace('.png', '').split('_', 1)[1].replace('_', ' ')
        html_content += f"""
        <div class="gsea-plot">
            <h4>{pathway_name}</h4>
            <img src="data:image/png;base64,{img_data}" alt="GSEA {pathway_name}">
        </div>
"""
    
    html_content += """
        <h3>Rat (GSE208526) - Top 10 Pathways</h3>
"""
    
    # Rat GSEA plots
    rat_gsea_files = sorted([f for f in os.listdir(gsea_plots_dir) if f.startswith('rat_')])
    for gsea_file in rat_gsea_files:
        plot_path = os.path.join(gsea_plots_dir, gsea_file)
        img_data = encode_image(plot_path)
        pathway_name = gsea_file.replace('rat_', '').replace('.png', '').split('_', 1)[1].replace('_', ' ')
        html_content += f"""
        <div class="gsea-plot">
            <h4>{pathway_name}</h4>
            <img src="data:image/png;base64,{img_data}" alt="GSEA {pathway_name}">
        </div>
"""
    
    html_content += """
    </div>
"""

    # Conclusions
    html_content += """
    <div class="section" id="conclusions">
        <h2>💡 Conclusions</h2>
        
        <h3>Key Biological Insights</h3>
        <ol>
            <li><strong>Highly Conserved Function:</strong> 88% of pathways are shared between human and mouse/rat,
            indicating strong evolutionary conservation of MME/Mme function.</li>
            
            <li><strong>Primary Role in Immunity:</strong> Strongest enrichment in immune-related pathways, 
            particularly neutrophil degranulation and innate immune response.</li>
            
            <li><strong>Neuroprotective Function:</strong> Presence in multiple Alzheimer's disease pathways 
            supports role in amyloid-β degradation.</li>
            
            <li><strong>Not a Hallmark Gene:</strong> Absence from MSigDB Hallmark indicates context-specific 
            rather than broadly dysregulated function.</li>
            
            <li><strong>Multi-System Involvement:</strong> Found in pathways spanning immune, nervous, 
            cardiovascular, and renal systems.</li>
        </ol>
        
        <h3>Clinical Relevance</h3>
        <ul>
            <li><strong>Alzheimer's Disease:</strong> Therapeutic target for increasing amyloid-β clearance</li>
            <li><strong>Kidney Disease:</strong> Involved in FSGS and acute kidney injury</li>
            <li><strong>Hypertension:</strong> Regulates angiotensin peptides in renin-angiotensin system</li>
            <li><strong>Cancer:</strong> CD10 marker expression in various cancers</li>
        </ul>
        
        <h3>Methodological Notes</h3>
        <ul>
            <li>Analysis based on PathDIP5 curated pathways (KEGG, Reactome, WikiPathways, ACSN2)</li>
            <li>Enrichment calculated using differential expression results from DESeq2</li>
            <li>GSEA-style plots show running enrichment scores for pathway gene distributions</li>
            <li>Rat data analyzed using mouse pathway annotations (high homology)</li>
        </ul>
    </div>
"""

    # Footer
    html_content += """
    <div class="footer">
        <p><strong>MME/Mme Pathway Analysis Report</strong></p>
        <p>Generated: November 17, 2024</p>
        <p>Databases: MSigDB v2025.1, PathDIP5 Curated</p>
        <p>Analysis Tools: Python 3, pandas, matplotlib, seaborn</p>
        <p>Location: /Users/ssankhe/Documents/Github/research/pathway_analysis_results/</p>
    </div>
</body>
</html>
"""

    # Save HTML file
    output_file = os.path.join(results_dir, 'COMPLETE_PATHWAY_ANALYSIS.html')
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✓ HTML report created: {output_file}")
    print(f"✓ File size: {os.path.getsize(output_file) / 1024 / 1024:.2f} MB")
    print(f"\nTo view the report:")
    print(f"  open {output_file}")

if __name__ == '__main__':
    print("Creating comprehensive HTML report with all plots...\n")
    create_html_report()
    print("\n✓ HTML report generation complete!")
