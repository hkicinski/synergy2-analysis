#!/usr/bin/env python3
"""
compare_protein_sequences.py - Compare protein sequences between two FASTA files
"""

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import argparse
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare protein sequences between two FASTA files')
    parser.add_argument('--file1', required=True, help='First protein FASTA file (gffread)')
    parser.add_argument('--file2', required=True, help='Second protein FASTA file (root.pep)')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    return parser.parse_args()

def run_blastp(query, subject, output_file):
    """Run BLASTp between query and subject"""
    cmd = [
        'blastp',
        '-query', query,
        '-subject', subject,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid pident length qlen slen evalue bitscore',
        '-num_threads', '4'
    ]
    
    print(f"Running BLASTp: {' '.join(cmd)}")
    start_time = time.time()
    subprocess.run(cmd, check=True)
    end_time = time.time()
    print(f"BLASTp completed in {end_time - start_time:.2f} seconds")

def analyze_blast_results(blast_file, output_dir, direction):
    """Analyze BLASTp results and generate statistics and visualizations"""
    # Check if the blast file exists and has content
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        print(f"No BLAST results found in {blast_file}")
        return None
    
    # Load BLASTp results
    columns = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_file, sep='\t', names=columns)
    
    # Calculate query coverage
    blast_df['qcovs'] = (blast_df['length'] / blast_df['qlen']) * 100
    
    # Set up plotting directory
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # Group by query ID and get best hit for each query
    best_hits = blast_df.loc[blast_df.groupby('qseqid')['bitscore'].idxmax()]
    
    # Calculate summary statistics
    total_queries = len(set(blast_df['qseqid']))
    matched_queries = len(best_hits)
    high_quality_matches = len(best_hits[(best_hits['pident'] >= 90) & (best_hits['qcovs'] >= 90)])
    match_rate = matched_queries / total_queries * 100 if total_queries > 0 else 0
    high_quality_rate = high_quality_matches / matched_queries * 100 if matched_queries > 0 else 0
    
    # Save results to CSV
    best_hits.to_csv(os.path.join(output_dir, f"best_hits_{direction}.csv"), index=False)
    
    # Create percent identity histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(best_hits['pident'], bins=20, kde=True)
    plt.axvline(x=90, color='red', linestyle='--', label='90% Identity Threshold')
    plt.title(f'Distribution of Percent Identity ({direction})')
    plt.xlabel('Percent Identity')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, f"percent_identity_distribution_{direction}.png"))
    plt.close()
    
    # Create query coverage histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(best_hits['qcovs'], bins=20, kde=True)
    plt.axvline(x=90, color='red', linestyle='--', label='90% Coverage Threshold')
    plt.title(f'Distribution of Query Coverage ({direction})')
    plt.xlabel('Query Coverage')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, f"query_coverage_distribution_{direction}.png"))
    plt.close()
    
    # Create scatter plot of identity vs coverage
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='pident', y='qcovs', data=best_hits, alpha=0.5)
    plt.axhline(y=90, color='red', linestyle='--')
    plt.axvline(x=90, color='red', linestyle='--')
    plt.title(f'Percent Identity vs. Query Coverage ({direction})')
    plt.xlabel('Percent Identity')
    plt.ylabel('Query Coverage')
    plt.savefig(os.path.join(plot_dir, f"identity_vs_coverage_{direction}.png"))
    plt.close()
    
    # Print summary
    print(f"\nBLASTp Analysis Summary ({direction}):")
    print(f"Total queries: {total_queries}")
    print(f"Queries with hits: {matched_queries} ({match_rate:.2f}%)")
    print(f"High quality matches (≥90% identity & coverage): {high_quality_matches} ({high_quality_rate:.2f}% of matches)")
    
    return {
        'total_queries': total_queries,
        'matched_queries': matched_queries,
        'high_quality_matches': high_quality_matches,
        'match_rate': match_rate,
        'high_quality_rate': high_quality_rate,
        'direction': direction
    }

def main():
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Count sequences in each file
    file1_count = sum(1 for _ in SeqIO.parse(args.file1, "fasta"))
    file2_count = sum(1 for _ in SeqIO.parse(args.file2, "fasta"))
    
    print(f"File 1 (gffread): {file1_count} sequences")
    print(f"File 2 (root.pep): {file2_count} sequences")
    
    # Run BLASTp in both directions
    # Direction 1: gffread → root.pep
    blast_output1 = os.path.join(args.output_dir, "gffread_to_root.tsv")
    run_blastp(args.file1, args.file2, blast_output1)
    stats1 = analyze_blast_results(blast_output1, args.output_dir, "gffread_to_root")
    
    # Direction 2: root.pep → gffread
    blast_output2 = os.path.join(args.output_dir, "root_to_gffread.tsv")
    run_blastp(args.file2, args.file1, blast_output2)
    stats2 = analyze_blast_results(blast_output2, args.output_dir, "root_to_gffread")
    
    # Create summary bar chart if both analyses succeeded
    if stats1 and stats2:
        plt.figure(figsize=(12, 7))
        
        metrics = ['Match Rate', 'High Quality Match Rate']
        values1 = [stats1['match_rate'], stats1['high_quality_rate']]
        values2 = [stats2['match_rate'], stats2['high_quality_rate']]
        
        x = range(len(metrics))
        width = 0.35
        
        plt.bar([i - width/2 for i in x], values1, width, label='gffread → root.pep')
        plt.bar([i + width/2 for i in x], values2, width, label='root.pep → gffread')
        
        plt.axhline(y=90, color='red', linestyle='--', label='90% Threshold')
        plt.ylabel('Percentage (%)')
        plt.title('BLASTp Analysis Summary Metrics')
        plt.xticks(x, metrics)
        plt.legend()
        
        # Add value labels
        for i, v in enumerate(values1):
            plt.text(i - width/2, v + 1, f'{v:.1f}%', ha='center')
        for i, v in enumerate(values2):
            plt.text(i + width/2, v + 1, f'{v:.1f}%', ha='center')
        
        plt.ylim(0, 105)
        plt.savefig(os.path.join(args.output_dir, "plots", "summary_metrics.png"))
        plt.close()
        
        # Provide overall interpretation
        print("\nOverall Interpretation:")
        if stats1['high_quality_rate'] >= 90 and stats2['high_quality_rate'] >= 90:
            print("✓ EXCELLENT MATCH: The vast majority of proteins match well in both directions")
            print("  This suggests the gffread translations match very well with SYNERGY2's proteins.")
        elif stats1['high_quality_rate'] >= 70 or stats2['high_quality_rate'] >= 70:
            print("⚠ GOOD MATCH WITH SOME DIFFERENCES: Many proteins match well, but there are some differences")
            print("  This suggests some systematic differences in how proteins were translated.")
        else:
            print("✗ SIGNIFICANT DIFFERENCES: Relatively few proteins match well between the two sets")
            print("  This indicates significant differences in protein translations, which may")
            print("  explain the orthology discrepancies observed with C. glabrata.")
    
    print(f"\nAll analysis results saved to: {args.output_dir}")

if __name__ == "__main__":
    main()