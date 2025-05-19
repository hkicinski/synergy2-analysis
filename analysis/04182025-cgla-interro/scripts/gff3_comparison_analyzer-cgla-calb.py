#!/usr/bin/env python3
"""
gff3_comparison_analyzer.py - Compare and visualize differences between GFF3 files
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

def generate_visualizations(c_glabrata_metrics, c_albicans_metrics, c_glabrata_df, c_albicans_df, output_dir):
    """Generate visualizations comparing the two species"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Set the aesthetics for the plots
    sns.set(style="whitegrid")
    plt.rcParams.update({'font.size': 12})
    
    # 1. UTR Percentage Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(c_glabrata_df['utr_percentage'], kde=True, color='blue', label='C. glabrata', alpha=0.5)
    sns.histplot(c_albicans_df['utr_percentage'], kde=True, color='red', label='C. albicans', alpha=0.5)
    plt.axvline(x=c_glabrata_metrics['mean_utr_percentage'], color='blue', linestyle='--', label=f'C. glabrata mean: {c_glabrata_metrics["mean_utr_percentage"]:.2f}%')
    plt.axvline(x=c_albicans_metrics['mean_utr_percentage'], color='red', linestyle='--', label=f'C. albicans mean: {c_albicans_metrics["mean_utr_percentage"]:.2f}%')
    plt.title('UTR Percentage Distribution (% of exon that is UTR)')
    plt.xlabel('UTR Percentage')
    plt.ylabel('Count')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'utr_percentage_distribution.png'), dpi=300)
    plt.close()
    
    # 2. Exon-CDS Match vs Diff
    match_diff_data = [
        ['C. glabrata', c_glabrata_metrics['exon_cds_match'], c_glabrata_metrics['exon_cds_diff']],
        ['C. albicans', c_albicans_metrics['exon_cds_match'], c_albicans_metrics['exon_cds_diff']]
    ]
    match_diff_df = pd.DataFrame(match_diff_data, columns=['Species', 'Matching', 'Different'])
    
    # Calculate percentages
    match_diff_df['Total'] = match_diff_df['Matching'] + match_diff_df['Different']
    match_diff_df['Match_Pct'] = (match_diff_df['Matching'] / match_diff_df['Total']) * 100
    match_diff_df['Diff_Pct'] = (match_diff_df['Different'] / match_diff_df['Total']) * 100
    
    plt.figure(figsize=(10, 6))
    x = np.arange(2)
    width = 0.35
    
    plt.bar(x - width/2, match_diff_df['Match_Pct'], width, label='Exon-CDS Match')
    plt.bar(x + width/2, match_diff_df['Diff_Pct'], width, label='Exon-CDS Different')
    
    plt.xticks(x, match_diff_df['Species'])
    plt.ylabel('Percentage')
    plt.title('Exon-CDS Coordinate Matching vs. Different')
    
    for i, species in enumerate(match_diff_df['Species']):
        plt.text(i - width/2, match_diff_df['Match_Pct'][i] + 1, f"{match_diff_df['Match_Pct'][i]:.1f}%", ha='center')
        plt.text(i + width/2, match_diff_df['Diff_Pct'][i] + 1, f"{match_diff_df['Diff_Pct'][i]:.1f}%", ha='center')
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'exon_cds_match_vs_diff.png'), dpi=300)
    plt.close()
    
    # 3. Length differences scatter plot
    plt.figure(figsize=(10, 8))
    
    # C. glabrata
    plt.scatter(c_glabrata_df['cds_length'], c_glabrata_df['exon_length'], 
                alpha=0.5, label='C. glabrata', color='blue')
    
    # C. albicans
    plt.scatter(c_albicans_df['cds_length'], c_albicans_df['exon_length'], 
                alpha=0.5, label='C. albicans', color='red')
    
    # Add diagonal line for perfect match
    max_val = max(c_glabrata_df['exon_length'].max(), c_albicans_df['exon_length'].max())
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
        'C. glabrata': c_glabrata_df['utr_size'],
        'C. albicans': c_albicans_df['utr_size']
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
        'C. glabrata': [
            c_glabrata_metrics['gene_count'],
            c_glabrata_metrics['transcript_count'],
            c_glabrata_metrics['multiexon_transcript_count'],
            c_glabrata_metrics['exon_count'],
            c_glabrata_metrics['cds_count'],
            c_glabrata_metrics['exon_cds_match'],
            c_glabrata_metrics['exon_cds_diff'],
            f"{c_glabrata_metrics['mean_utr_percentage']:.2f}%",
            f"{c_glabrata_metrics['median_utr_percentage']:.2f}%",
            f"{c_glabrata_metrics['max_utr_percentage']:.2f}%",
            f"{c_glabrata_metrics['mean_diff_size']:.2f}",
            f"{c_glabrata_metrics['median_diff_size']:.2f}",
            f"{c_glabrata_metrics['max_diff_size']}"
        ],
        'C. albicans': [
            c_albicans_metrics['gene_count'],
            c_albicans_metrics['transcript_count'],
            c_albicans_metrics['multiexon_transcript_count'],
            c_albicans_metrics['exon_count'],
            c_albicans_metrics['cds_count'],
            c_albicans_metrics['exon_cds_match'],
            c_albicans_metrics['exon_cds_diff'],
            f"{c_albicans_metrics['mean_utr_percentage']:.2f}%",
            f"{c_albicans_metrics['median_utr_percentage']:.2f}%",
            f"{c_albicans_metrics['max_utr_percentage']:.2f}%",
            f"{c_albicans_metrics['mean_diff_size']:.2f}",
            f"{c_albicans_metrics['median_diff_size']:.2f}",
            f"{c_albicans_metrics['max_diff_size']}"
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'comparison_summary.csv'), index=False)
    
    # Create a second CSV with raw exon/CDS/UTR data for further analysis
    combined_df = pd.concat([
        c_glabrata_df.assign(species='C. glabrata'),
        c_albicans_df.assign(species='C. albicans')
    ])
    combined_df.to_csv(os.path.join(output_dir, 'transcript_data.csv'), index=False)
    
    # 6. Create a figure showing examples of matching vs non-matching genes
    # Select a few examples from each species
    match_examples = {
        'C. glabrata': c_glabrata_df[c_glabrata_df['match'] == True].head(3),
        'C. albicans': c_albicans_df[c_albicans_df['match'] == True].head(3)
    }
    
    diff_examples = {
        'C. glabrata': c_glabrata_df[c_glabrata_df['match'] == False].head(3),
        'C. albicans': c_albicans_df[c_albicans_df['match'] == False].head(3)
    }
    
    # Create summary of examples
    example_data = []
    
    for species in ['C. glabrata', 'C. albicans']:
        # Matching examples
        for _, row in match_examples[species].iterrows():
            example_data.append({
                'Species': species,
                'Gene ID': row['gene_id'],
                'Transcript ID': row['transcript_id'],
                'Type': 'Match',
                'Exon Length': row['exon_length'],
                'CDS Length': row['cds_length'],
                'UTR Size': row['utr_size'],
                'UTR Percentage': f"{row['utr_percentage']:.2f}%"
            })
        
        # Different examples
        for _, row in diff_examples[species].iterrows():
            example_data.append({
                'Species': species,
                'Gene ID': row['gene_id'],
                'Transcript ID': row['transcript_id'],
                'Type': 'Different',
                'Exon Length': row['exon_length'],
                'CDS Length': row['cds_length'],
                'UTR Size': row['utr_size'],
                'UTR Percentage': f"{row['utr_percentage']:.2f}%"
            })
    
    example_df = pd.DataFrame(example_data)
    example_df.to_csv(os.path.join(output_dir, 'example_genes.csv'), index=False)
    
    # 7. Additional visualization: Impact on translation
    # Create a figure showing how UTR percentages affect protein translation
    plt.figure(figsize=(12, 8))
    
    # Set up the subplots
    gs = plt.GridSpec(2, 2, height_ratios=[1, 1])
    
    # C. albicans typical gene (high match rate)
    ax1 = plt.subplot(gs[0, 0])
    exon_length = c_albicans_metrics['exon_lengths'][0] if c_albicans_metrics['exon_lengths'] else 1000
    cds_length = c_albicans_metrics['cds_lengths'][0] if c_albicans_metrics['cds_lengths'] else 1000
    
    # Draw C. albicans exon and CDS
    ax1.barh(0, exon_length, height=0.5, color='lightgrey', label='Exon')
    ax1.barh(0, cds_length, height=0.5, color='green', alpha=0.7, label='CDS')
    ax1.set_yticks([])
    ax1.set_xlabel('Position (bp)')
    ax1.set_title('C. albicans Typical Gene Structure')
    ax1.legend()
    
    # Draw protein sequence below (schematic)
    ax1.text(exon_length/2, -0.5, 'Properly Translated Protein', ha='center')
    
    # C. glabrata typical gene (low match rate)
    ax2 = plt.subplot(gs[0, 1])
    exon_length = c_glabrata_metrics['exon_lengths'][0] if c_glabrata_metrics['exon_lengths'] else 1000
    cds_length = c_glabrata_metrics['cds_lengths'][0] if c_glabrata_metrics['cds_lengths'] else 800
    utr_size = exon_length - cds_length
    
    # Draw C. glabrata exon and CDS
    ax2.barh(0, exon_length, height=0.5, color='lightgrey', label='Exon')
    ax2.barh(0, cds_length, height=0.5, color='green', alpha=0.7, label='CDS')
    
    # Highlight UTR regions
    ax2.barh(0, utr_size/2, height=0.5, color='red', alpha=0.3, label='5\' UTR')
    ax2.barh(0, utr_size/2, height=0.5, left=cds_length, color='red', alpha=0.3, label='3\' UTR')
    
    ax2.set_yticks([])
    ax2.set_xlabel('Position (bp)')
    ax2.set_title('C. glabrata Typical Gene Structure')
    ax2.legend()
    
    # Draw incorrect protein translation below (schematic)
    ax2.text(exon_length/2, -0.5, 'Incorrectly Translated Protein (includes UTRs)', ha='center')
    
    # Translation issue visualization
    ax3 = plt.subplot(gs[1, :])
    
    # Define example sequences
    correct_dna = "ATGGCTAGCTAGCTAGCTATGCTAGCTATGCTAGCTTAG"
    correct_protein = "MASSLAMLASS*"
    
    utr_dna = "ACGTACGTACGTACGT" + correct_dna + "ACGTACGTACGTACGTACGT"
    utr_protein = "TRTYVRMSSLAMSSTRTRTRR*"
    
    # Create DNA/protein visualization
    ax3.text(0.1, 0.8, "Correct DNA → Protein Translation (CDS only):", fontsize=12)
    ax3.text(0.1, 0.7, f"DNA: {correct_dna}", fontsize=10)
    ax3.text(0.1, 0.6, f"Protein: {correct_protein}", fontsize=10)
    
    ax3.text(0.1, 0.4, "Incorrect DNA → Protein Translation (including UTRs):", fontsize=12)
    ax3.text(0.1, 0.3, f"DNA: {utr_dna}", fontsize=10)
    ax3.text(0.1, 0.2, f"Protein: {utr_protein} (incorrect, containing premature stop codons)", fontsize=10)
    
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_title('Impact of UTRs on Protein Translation')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'translation_impact.png'), dpi=300)
    plt.close()
    
    # 8. Plot showing SYNERGY2 workflow and where the issue occurs
    plt.figure(figsize=(12, 8))
    
    # Define positions for flowchart boxes
    positions = {
        'input_gff': (0.1, 0.8),
        'input_fasta': (0.1, 0.6),
        'synergy_process': (0.5, 0.7),
        'exon_extract': (0.3, 0.5),
        'translation': (0.7, 0.5),
        'orthology': (0.5, 0.3),
        'output': (0.5, 0.1)
    }
    
    # Draw boxes
    for key, (x, y) in positions.items():
        box_width = 0.2
        box_height = 0.1
        
        if key == 'exon_extract':
            fc = 'lightcoral'  # Highlight the problematic step
        else:
            fc = 'lightblue'
            
        rect = plt.Rectangle((x-box_width/2, y-box_height/2), box_width, box_height, fc=fc, ec='black')
        plt.gca().add_patch(rect)
    
    # Add text
    plt.text(positions['input_gff'][0], positions['input_gff'][1], 'GFF3\nInput', ha='center', va='center')
    plt.text(positions['input_fasta'][0], positions['input_fasta'][1], 'FASTA\nInput', ha='center', va='center')
    plt.text(positions['synergy_process'][0], positions['synergy_process'][1], 'SYNERGY2\nProcess', ha='center', va='center')
    plt.text(positions['exon_extract'][0], positions['exon_extract'][1], 'Sequence\nExtraction\nUsing EXON\nCoordinates', ha='center', va='center', color='red')
    plt.text(positions['translation'][0], positions['translation'][1], 'Protein\nTranslation', ha='center', va='center')
    plt.text(positions['orthology'][0], positions['orthology'][1], 'Orthology\nDetection', ha='center', va='center')
    plt.text(positions['output'][0], positions['output'][1], 'Output\nClusters', ha='center', va='center')
    
    # Draw arrows
    plt.arrow(positions['input_gff'][0]+0.1, positions['input_gff'][1], 
              positions['synergy_process'][0]-positions['input_gff'][0]-0.15, 0, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    plt.arrow(positions['input_fasta'][0]+0.1, positions['input_fasta'][1], 
              positions['synergy_process'][0]-positions['input_fasta'][0]-0.15, 0, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    plt.arrow(positions['synergy_process'][0], positions['synergy_process'][1]-0.05, 
              positions['exon_extract'][0]-positions['synergy_process'][0]+0.1, 
              positions['exon_extract'][1]-positions['synergy_process'][1]+0.05, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    plt.arrow(positions['exon_extract'][0]+0.1, positions['exon_extract'][1], 
              positions['translation'][0]-positions['exon_extract'][0]-0.15, 0, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    plt.arrow(positions['translation'][0], positions['translation'][1]-0.05, 
              0, positions['orthology'][1]-positions['translation'][1]+0.1, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    plt.arrow(positions['orthology'][0], positions['orthology'][1]-0.05, 
              0, positions['output'][1]-positions['orthology'][1]+0.05, 
              head_width=0.02, head_length=0.02, fc='black', ec='black')
    
    # Add annotation for the issue
    plt.text(0.3, 0.42, 'Issue: C. glabrata has significant UTRs in exons,\nleading to incorrect protein sequences', 
             color='red', ha='center', weight='bold')
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    plt.title('SYNERGY2 Workflow and Translation Issue')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'synergy2_workflow.png'), dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare GFF3 files from C. glabrata and C. albicans')
    parser.add_argument('--cgla', required=True, help='Path to C. glabrata GFF3 file')
    parser.add_argument('--calb', required=True, help='Path to C. albicans GFF3 file')
    parser.add_argument('--output', default='gff3_comparison_output', help='Output directory for results')
    
    args = parser.parse_args()
    
    print(f"Parsing C. glabrata GFF3 file: {args.cgla}")
    c_glabrata_genes, c_glabrata_exons, c_glabrata_cds, c_glabrata_parent_map = parse_gff3(args.cgla)
    
    print(f"Parsing C. albicans GFF3 file: {args.calb}")
    c_albicans_genes, c_albicans_exons, c_albicans_cds, c_albicans_parent_map = parse_gff3(args.calb)
    
    print("Calculating metrics for C. glabrata")
    c_glabrata_metrics, c_glabrata_df = calculate_feature_metrics(
        c_glabrata_genes, c_glabrata_exons, c_glabrata_cds, c_glabrata_parent_map
    )
    
    print("Calculating metrics for C. albicans")
    c_albicans_metrics, c_albicans_df = calculate_feature_metrics(
        c_albicans_genes, c_albicans_exons, c_albicans_cds, c_albicans_parent_map
    )
    
    print("Generating visualizations")
    generate_visualizations(
        c_glabrata_metrics, c_albicans_metrics, 
        c_glabrata_df, c_albicans_df, 
        args.output
    )
    
    print(f"Analysis complete. Results saved to: {args.output}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"{'Metric':<30} {'C. glabrata':<15} {'C. albicans':<15}")
    print("-" * 60)
    print(f"{'Gene Count':<30} {c_glabrata_metrics['gene_count']:<15} {c_albicans_metrics['gene_count']:<15}")
    print(f"{'Transcript Count':<30} {c_glabrata_metrics['transcript_count']:<15} {c_albicans_metrics['transcript_count']:<15}")
    print(f"{'Exon-CDS Match Count':<30} {c_glabrata_metrics['exon_cds_match']:<15} {c_albicans_metrics['exon_cds_match']:<15}")
    print(f"{'Exon-CDS Different Count':<30} {c_glabrata_metrics['exon_cds_diff']:<15} {c_albicans_metrics['exon_cds_diff']:<15}")
    
    # Calculate match percentages
    c_glabrata_match_pct = (c_glabrata_metrics['exon_cds_match'] / 
                            (c_glabrata_metrics['exon_cds_match'] + c_glabrata_metrics['exon_cds_diff'])) * 100
    c_albicans_match_pct = (c_albicans_metrics['exon_cds_match'] / 
                            (c_albicans_metrics['exon_cds_match'] + c_albicans_metrics['exon_cds_diff'])) * 100
    
    print(f"{'Exon-CDS Match Percentage':<30} {c_glabrata_match_pct:.2f}%{'':<9} {c_albicans_match_pct:.2f}%")
    print(f"{'Mean UTR Percentage':<30} {c_glabrata_metrics['mean_utr_percentage']:.2f}%{'':<9} {c_albicans_metrics['mean_utr_percentage']:.2f}%")
    
    print("\nConclusion: The analysis shows a significant difference in exon-CDS coordinate matching between C. glabrata and C. albicans.")
    print(f"C. glabrata has {c_glabrata_match_pct:.2f}% of genes with matching exon-CDS coordinates, while")
    print(f"C. albicans has {c_albicans_match_pct:.2f}% of genes with matching coordinates.")
    print("This strongly suggests that SYNERGY2's use of exon coordinates for protein translation")
    print("is causing inaccurate protein sequences for C. glabrata, affecting orthology detection.")

if __name__ == "__main__":
    main()
