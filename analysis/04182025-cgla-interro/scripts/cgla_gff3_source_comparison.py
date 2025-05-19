#!/usr/bin/env python3
"""
cgla_gff3_source_comparison.py - Compare NCBI and CGD GFF3 files for C. glabrata
focusing on exon vs CDS feature relationships
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
import os

def parse_gff3(file_path):
    """Parse GFF3 file and extract gene, mRNA, exon, and CDS features"""
    genes = {}
    exons = {}
    cds_regions = {}
    parent_map = {}
    curr_gene = None
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # Store features by type
            if feature_type == 'gene':
                gene_id = attr_dict.get('ID', '')
                if gene_id:
                    curr_gene = gene_id
                    genes[gene_id] = {
                        'seqid': seqid,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'transcripts': []
                    }
            
            elif feature_type == 'mRNA':
                parent = attr_dict.get('Parent', '')
                transcript_id = attr_dict.get('ID', '')
                if parent and transcript_id:
                    if parent in genes:
                        genes[parent]['transcripts'].append(transcript_id)
                    parent_map[transcript_id] = parent
            
            elif feature_type == 'exon':
                parent = attr_dict.get('Parent', '')
                exon_id = attr_dict.get('ID', '')
                if parent:
                    if parent not in exons:
                        exons[parent] = []
                    exons[parent].append({
                        'seqid': seqid,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'id': exon_id
                    })
            
            elif feature_type == 'CDS':
                parent = attr_dict.get('Parent', '')
                cds_id = attr_dict.get('ID', '')
                if parent:
                    if parent not in cds_regions:
                        cds_regions[parent] = []
                    cds_regions[parent].append({
                        'seqid': seqid,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'id': cds_id
                    })
    
    return genes, exons, cds_regions, parent_map

def calculate_feature_metrics(genes, exons, cds_regions, parent_map):
    """Calculate metrics for feature relationships"""
    metrics = {
        'gene_count': len(genes),
        'transcript_count': sum(len(g['transcripts']) for g in genes.values()),
        'multiexon_transcript_count': 0,
        'exon_count': sum(len(exons[t]) for t in exons),
        'cds_count': sum(len(cds_regions[t]) for t in cds_regions if t in cds_regions),
        'exon_lengths': [],
        'cds_lengths': [],
        'utr_percentages': [],
        'exon_cds_match': 0,
        'exon_cds_diff': 0,
        'exon_cds_diff_sizes': []
    }
    
    # Transcript metrics
    for transcript_id in exons:
        if len(exons[transcript_id]) > 1:
            metrics['multiexon_transcript_count'] += 1
    
    # Compare exon vs CDS
    transcript_data = []
    
    for transcript_id in exons:
        if transcript_id not in cds_regions:
            continue
        
        transcript_exons = exons[transcript_id]
        transcript_cds = cds_regions[transcript_id]
        
        # Handle single exon cases
        if len(transcript_exons) == 1 and len(transcript_cds) == 1:
            exon = transcript_exons[0]
            cds = transcript_cds[0]
            
            exon_length = exon['end'] - exon['start'] + 1
            cds_length = cds['end'] - cds['start'] + 1
            
            metrics['exon_lengths'].append(exon_length)
            metrics['cds_lengths'].append(cds_length)
            
            # Calculate UTR percentage
            utr_size = exon_length - cds_length
            utr_percentage = (utr_size / exon_length) * 100 if exon_length > 0 else 0
            metrics['utr_percentages'].append(utr_percentage)
            
            if exon['start'] == cds['start'] and exon['end'] == cds['end']:
                metrics['exon_cds_match'] += 1
            else:
                metrics['exon_cds_diff'] += 1
                metrics['exon_cds_diff_sizes'].append(utr_size)
            
            gene_id = parent_map.get(transcript_id, '')
            transcript_data.append({
                'gene_id': gene_id,
                'transcript_id': transcript_id,
                'exon_length': exon_length,
                'cds_length': cds_length,
                'utr_size': utr_size,
                'utr_percentage': utr_percentage,
                'match': exon['start'] == cds['start'] and exon['end'] == cds['end']
            })
        
        # Handle multi-exon cases (simplified for this analysis)
        else:
            total_exon_length = sum(e['end'] - e['start'] + 1 for e in transcript_exons)
            total_cds_length = sum(c['end'] - c['start'] + 1 for c in transcript_cds)
            
            metrics['exon_lengths'].append(total_exon_length)
            metrics['cds_lengths'].append(total_cds_length)
            
            # Calculate UTR percentage
            utr_size = total_exon_length - total_cds_length
            utr_percentage = (utr_size / total_exon_length) * 100 if total_exon_length > 0 else 0
            metrics['utr_percentages'].append(utr_percentage)
            
            # For multi-exon, we'll consider it a match if total lengths are the same
            if total_exon_length == total_cds_length:
                metrics['exon_cds_match'] += 1
            else:
                metrics['exon_cds_diff'] += 1
                metrics['exon_cds_diff_sizes'].append(utr_size)
            
            gene_id = parent_map.get(transcript_id, '')
            transcript_data.append({
                'gene_id': gene_id,
                'transcript_id': transcript_id,
                'exon_length': total_exon_length,
                'cds_length': total_cds_length,
                'utr_size': utr_size,
                'utr_percentage': utr_percentage,
                'match': total_exon_length == total_cds_length
            })
    
    # Create dataframe for more detailed analysis
    transcript_df = pd.DataFrame(transcript_data)
    
    # Calculate summary statistics
    if metrics['utr_percentages']:
        metrics['mean_utr_percentage'] = np.mean(metrics['utr_percentages'])
        metrics['median_utr_percentage'] = np.median(metrics['utr_percentages'])
        metrics['max_utr_percentage'] = np.max(metrics['utr_percentages'])
    else:
        metrics['mean_utr_percentage'] = 0
        metrics['median_utr_percentage'] = 0
        metrics['max_utr_percentage'] = 0
    
    if metrics['exon_cds_diff_sizes']:
        metrics['mean_diff_size'] = np.mean(metrics['exon_cds_diff_sizes'])
        metrics['median_diff_size'] = np.median(metrics['exon_cds_diff_sizes'])
        metrics['max_diff_size'] = np.max(metrics['exon_cds_diff_sizes'])
    else:
        metrics['mean_diff_size'] = 0
        metrics['median_diff_size'] = 0
        metrics['max_diff_size'] = 0
    
    return metrics, transcript_df

def generate_visualizations(cgla_cgd_metrics, cgla_ncbi_metrics, cgla_cgd_df, cgla_ncbi_df, output_dir):
    """Generate visualizations comparing CGD and NCBI annotations for C. glabrata"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Set the aesthetics for the plots
    sns.set(style="whitegrid")
    plt.rcParams.update({'font.size': 12})
    
    # 1. UTR Percentage Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(cgla_cgd_df['utr_percentage'], kde=True, color='blue', label='CGD', alpha=0.5)
    sns.histplot(cgla_ncbi_df['utr_percentage'], kde=True, color='orange', label='NCBI', alpha=0.5)
    plt.axvline(x=cgla_cgd_metrics['mean_utr_percentage'], color='blue', linestyle='--', label=f'CGD mean: {cgla_cgd_metrics["mean_utr_percentage"]:.2f}%')
    plt.axvline(x=cgla_ncbi_metrics['mean_utr_percentage'], color='orange', linestyle='--', label=f'NCBI mean: {cgla_ncbi_metrics["mean_utr_percentage"]:.2f}%')
    plt.title('UTR Percentage Distribution (% of exon that is UTR)')
    plt.xlabel('UTR Percentage')
    plt.ylabel('Count')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'utr_percentage_distribution.png'), dpi=300)
    plt.close()
    
    # 2. Exon-CDS Match vs Diff
    match_diff_data = [
        ['CGD', cgla_cgd_metrics['exon_cds_match'], cgla_cgd_metrics['exon_cds_diff']],
        ['NCBI', cgla_ncbi_metrics['exon_cds_match'], cgla_ncbi_metrics['exon_cds_diff']]
    ]
    match_diff_df = pd.DataFrame(match_diff_data, columns=['Source', 'Matching', 'Different'])
    
    # Calculate percentages
    match_diff_df['Total'] = match_diff_df['Matching'] + match_diff_df['Different']
    match_diff_df['Match_Pct'] = (match_diff_df['Matching'] / match_diff_df['Total']) * 100
    match_diff_df['Diff_Pct'] = (match_diff_df['Different'] / match_diff_df['Total']) * 100
    
    plt.figure(figsize=(10, 6))
    x = np.arange(2)
    width = 0.35
    
    plt.bar(x - width/2, match_diff_df['Match_Pct'], width, label='Exon-CDS Match')
    plt.bar(x + width/2, match_diff_df['Diff_Pct'], width, label='Exon-CDS Different')
    
    plt.xticks(x, match_diff_df['Source'])
    plt.ylabel('Percentage')
    plt.title('Exon-CDS Coordinate Matching vs. Different')
    
    for i, source in enumerate(match_diff_df['Source']):
        plt.text(i - width/2, match_diff_df['Match_Pct'][i] + 1, f"{match_diff_df['Match_Pct'][i]:.1f}%", ha='center')
        plt.text(i + width/2, match_diff_df['Diff_Pct'][i] + 1, f"{match_diff_df['Diff_Pct'][i]:.1f}%", ha='center')
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'exon_cds_match_vs_diff.png'), dpi=300)
    plt.close()
    
    # 3. Length differences scatter plot
    plt.figure(figsize=(10, 8))
    
    # CGD
    plt.scatter(cgla_cgd_df['cds_length'], cgla_cgd_df['exon_length'], 
                alpha=0.5, label='CGD', color='blue')
    
    # NCBI
    plt.scatter(cgla_ncbi_df['cds_length'], cgla_ncbi_df['exon_length'], 
                alpha=0.5, label='NCBI', color='orange')
    
    # Add diagonal line for perfect match
    max_val = max(cgla_cgd_df['exon_length'].max(), cgla_ncbi_df['exon_length'].max())
    plt.plot([0, max_val], [0, max_val], 'k--', label='Perfect Match Line')
    
    plt.title('CDS Length vs. Exon Length')
    plt.xlabel('CDS Length (bp)')
    plt.ylabel('Exon Length (bp)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cds_vs_exon_length.png'), dpi=300)
    plt.close()
    
    # 4. UTR size boxplot
    plt.figure(figsize=(10, 6))
    
    # Prepare data for plotting
    data = {
        'CGD': cgla_cgd_df['utr_size'],
        'NCBI': cgla_ncbi_df['utr_size']
    }
    df = pd.DataFrame({k: pd.Series(v) for k, v in data.items()})
    
    # Create boxplot
    sns.boxplot(data=df)
    plt.title('UTR Size Distribution')
    plt.ylabel('UTR Size (bp)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'utr_size_boxplot.png'), dpi=300)
    plt.close()
    
    # 5. Summary metrics table
    summary_data = {
        'Metric': [
            'Gene Count',
            'Transcript Count',
            'Multi-exon Transcript Count',
            'Exon Count',
            'CDS Count',
            'Exon-CDS Match Count',
            'Exon-CDS Different Count',
            'Mean UTR Percentage',
            'Median UTR Percentage',
            'Max UTR Percentage',
            'Mean UTR Size (bp)',
            'Median UTR Size (bp)',
            'Max UTR Size (bp)'
        ],
        'CGD': [
            cgla_cgd_metrics['gene_count'],
            cgla_cgd_metrics['transcript_count'],
            cgla_cgd_metrics['multiexon_transcript_count'],
            cgla_cgd_metrics['exon_count'],
            cgla_cgd_metrics['cds_count'],
            cgla_cgd_metrics['exon_cds_match'],
            cgla_cgd_metrics['exon_cds_diff'],
            f"{cgla_cgd_metrics['mean_utr_percentage']:.2f}%",
            f"{cgla_cgd_metrics['median_utr_percentage']:.2f}%",
            f"{cgla_cgd_metrics['max_utr_percentage']:.2f}%",
            f"{cgla_cgd_metrics['mean_diff_size']:.2f}",
            f"{cgla_cgd_metrics['median_diff_size']:.2f}",
            f"{cgla_cgd_metrics['max_diff_size']}"
        ],
        'NCBI': [
            cgla_ncbi_metrics['gene_count'],
            cgla_ncbi_metrics['transcript_count'],
            cgla_ncbi_metrics['multiexon_transcript_count'],
            cgla_ncbi_metrics['exon_count'],
            cgla_ncbi_metrics['cds_count'],
            cgla_ncbi_metrics['exon_cds_match'],
            cgla_ncbi_metrics['exon_cds_diff'],
            f"{cgla_ncbi_metrics['mean_utr_percentage']:.2f}%",
            f"{cgla_ncbi_metrics['median_utr_percentage']:.2f}%",
            f"{cgla_ncbi_metrics['max_utr_percentage']:.2f}%",
            f"{cgla_ncbi_metrics['mean_diff_size']:.2f}",
            f"{cgla_ncbi_metrics['median_diff_size']:.2f}",
            f"{cgla_ncbi_metrics['max_diff_size']}"
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'comparison_summary.csv'), index=False)
    
    # Create a second CSV with raw exon/CDS/UTR data for further analysis
    combined_df = pd.concat([
        cgla_cgd_df.assign(source='CGD'),
        cgla_ncbi_df.assign(source='NCBI')
    ])
    combined_df.to_csv(os.path.join(output_dir, 'transcript_data.csv'), index=False)
    
    # 6. Create a figure showing examples of matching vs non-matching genes
    # Select a few examples from each source
    match_examples = {
        'CGD': cgla_cgd_df[cgla_cgd_df['match'] == True].head(3) if not cgla_cgd_df[cgla_cgd_df['match'] == True].empty else pd.DataFrame(),
        'NCBI': cgla_ncbi_df[cgla_ncbi_df['match'] == True].head(3) if not cgla_ncbi_df[cgla_ncbi_df['match'] == True].empty else pd.DataFrame()
    }
    
    diff_examples = {
        'CGD': cgla_cgd_df[cgla_cgd_df['match'] == False].head(3) if not cgla_cgd_df[cgla_cgd_df['match'] == False].empty else pd.DataFrame(),
        'NCBI': cgla_ncbi_df[cgla_ncbi_df['match'] == False].head(3) if not cgla_ncbi_df[cgla_ncbi_df['match'] == False].empty else pd.DataFrame()
    }
    
    # Create summary of examples
    example_data = []
    
    for source in ['CGD', 'NCBI']:
        # Matching examples
        if not match_examples[source].empty:
            for _, row in match_examples[source].iterrows():
                example_data.append({
                    'Source': source,
                    'Gene ID': row['gene_id'],
                    'Transcript ID': row['transcript_id'],
                    'Type': 'Match',
                    'Exon Length': row['exon_length'],
                    'CDS Length': row['cds_length'],
                    'UTR Size': row['utr_size'],
                    'UTR Percentage': f"{row['utr_percentage']:.2f}%"
                })
        
        # Different examples
        if not diff_examples[source].empty:
            for _, row in diff_examples[source].iterrows():
                example_data.append({
                    'Source': source,
                    'Gene ID': row['gene_id'],
                    'Transcript ID': row['transcript_id'],
                    'Type': 'Different',
                    'Exon Length': row['exon_length'],
                    'CDS Length': row['cds_length'],
                    'UTR Size': row['utr_size'],
                    'UTR Percentage': f"{row['utr_percentage']:.2f}%"
                })
    
    if example_data:
        example_df = pd.DataFrame(example_data)
        example_df.to_csv(os.path.join(output_dir, 'example_genes.csv'), index=False)
    
    # 7. GFF3 structure comparison visualization
    plt.figure(figsize=(12, 8))
    
    # Set up the comparison
    plt.subplot(2, 1, 1)
    
    # CGD Example
    exon_length = 1000
    cds_length = 750
    utr_size = exon_length - cds_length
    
    # Draw CGD gene structure
    plt.barh(0, exon_length, height=0.5, color='lightgrey', label='Exon')
    plt.barh(0, cds_length, height=0.5, color='blue', alpha=0.7, label='CDS')
    
    # Highlight UTR regions
    plt.barh(0, utr_size/2, height=0.5, color='red', alpha=0.3, label='5\' UTR')
    plt.barh(0, utr_size/2, height=0.5, left=cds_length, color='red', alpha=0.3, label='3\' UTR')
    
    plt.yticks([0], ['CGD'])
    plt.title('CGD Annotation Structure')
    plt.xlabel('Position (bp)')
    plt.legend(loc='upper right')
    
    # NCBI Example
    plt.subplot(2, 1, 2)
    
    # Get typical values from NCBI data
    exon_length = np.mean(cgla_ncbi_metrics['exon_lengths']) if cgla_ncbi_metrics['exon_lengths'] else 1000
    cds_length = np.mean(cgla_ncbi_metrics['cds_lengths']) if cgla_ncbi_metrics['cds_lengths'] else 1000
    utr_size = exon_length - cds_length
    
    # Draw NCBI gene structure
    plt.barh(0, exon_length, height=0.5, color='lightgrey', label='Exon')
    plt.barh(0, cds_length, height=0.5, color='orange', alpha=0.7, label='CDS')
    
    # Highlight UTR regions
    if utr_size > 0:
        plt.barh(0, utr_size/2, height=0.5, color='red', alpha=0.3, label='5\' UTR')
        plt.barh(0, utr_size/2, height=0.5, left=cds_length, color='red', alpha=0.3, label='3\' UTR')
    
    plt.yticks([0], ['NCBI'])
    plt.title('NCBI Annotation Structure')
    plt.xlabel('Position (bp)')
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'annotation_structure_comparison.png'), dpi=300)
    plt.close()
    
def main():
    parser = argparse.ArgumentParser(description='Compare CGD and NCBI GFF3 files for C. glabrata')
    parser.add_argument('--cgla_cgd', required=True, help='Path to CGD C. glabrata GFF3 file')
    parser.add_argument('--cgla_ncbi', required=True, help='Path to NCBI C. glabrata GFF3 file')
    parser.add_argument('--output', default='gff3_comparison_output', help='Output directory for results')
    
    args = parser.parse_args()
    
    print(f"Parsing CGD C. glabrata GFF3 file: {args.cgla_cgd}")
    cgla_cgd_genes, cgla_cgd_exons, cgla_cgd_cds, cgla_cgd_parent_map = parse_gff3(args.cgla_cgd)
    
    print(f"Parsing NCBI C. glabrata GFF3 file: {args.cgla_ncbi}")
    cgla_ncbi_genes, cgla_ncbi_exons, cgla_ncbi_cds, cgla_ncbi_parent_map = parse_gff3(args.cgla_ncbi)
    
    print("Calculating metrics for CGD C. glabrata")
    cgla_cgd_metrics, cgla_cgd_df = calculate_feature_metrics(
        cgla_cgd_genes, cgla_cgd_exons, cgla_cgd_cds, cgla_cgd_parent_map
    )
    
    print("Calculating metrics for NCBI C. glabrata")
    cgla_ncbi_metrics, cgla_ncbi_df = calculate_feature_metrics(
        cgla_ncbi_genes, cgla_ncbi_exons, cgla_ncbi_cds, cgla_ncbi_parent_map
    )
    
    print("Generating visualizations")
    generate_visualizations(
        cgla_cgd_metrics, cgla_ncbi_metrics, 
        cgla_cgd_df, cgla_ncbi_df, 
        args.output
    )
    
    print(f"Analysis complete. Results saved to: {args.output}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"{'Metric':<30} {'CGD':<15} {'NCBI':<15}")
    print("-" * 60)
    print(f"{'Gene Count':<30} {cgla_cgd_metrics['gene_count']:<15} {cgla_ncbi_metrics['gene_count']:<15}")
    print(f"{'Transcript Count':<30} {cgla_cgd_metrics['transcript_count']:<15} {cgla_ncbi_metrics['transcript_count']:<15}")
    print(f"{'Exon-CDS Match Count':<30} {cgla_cgd_metrics['exon_cds_match']:<15} {cgla_ncbi_metrics['exon_cds_match']:<15}")
    print(f"{'Exon-CDS Different Count':<30} {cgla_cgd_metrics['exon_cds_diff']:<15} {cgla_ncbi_metrics['exon_cds_diff']:<15}")
    
    # Calculate match percentages
    cgla_cgd_match_pct = (cgla_cgd_metrics['exon_cds_match'] / 
                      (cgla_cgd_metrics['exon_cds_match'] + cgla_cgd_metrics['exon_cds_diff'])) * 100 if (cgla_cgd_metrics['exon_cds_match'] + cgla_cgd_metrics['exon_cds_diff']) > 0 else 0
    cgla_ncbi_match_pct = (cgla_ncbi_metrics['exon_cds_match'] / 
                      (cgla_ncbi_metrics['exon_cds_match'] + cgla_ncbi_metrics['exon_cds_diff'])) * 100 if (cgla_ncbi_metrics['exon_cds_match'] + cgla_ncbi_metrics['exon_cds_diff']) > 0 else 0
    
    print(f"{'Exon-CDS Match Percentage':<30} {cgla_cgd_match_pct:.2f}%{'':<9} {cgla_ncbi_match_pct:.2f}%")
    print(f"{'Mean UTR Percentage':<30} {cgla_cgd_metrics['mean_utr_percentage']:.2f}%{'':<9} {cgla_ncbi_metrics['mean_utr_percentage']:.2f}%")
    
    print("\nConclusion: This analysis compares the CGD and NCBI annotations for C. glabrata")
    print(f"CGD has {cgla_cgd_match_pct:.2f}% of genes with matching exon-CDS coordinates")
    print(f"NCBI has {cgla_ncbi_match_pct:.2f}% of genes with matching exon-CDS coordinates")
    if abs(cgla_cgd_match_pct - cgla_ncbi_match_pct) < 10:
        print("Both annotation sources show similar patterns, confirming this is a consistent property of C. glabrata annotation.")
    else:
        print("There are significant differences between annotation sources, suggesting annotation methodology affects this property.")

if __name__ == "__main__":
    main()
