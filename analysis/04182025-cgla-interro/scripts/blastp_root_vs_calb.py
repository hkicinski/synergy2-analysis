#!/usr/bin/env python3
"""
blastp_root_vs_calb.py - Perform stringent BLASTp analysis between root.pep and C. albicans proteome
"""

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import argparse
import time

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare root.pep vs C. albicans proteome using BLASTp')
    parser.add_argument('--root_pep', required=True, help='Path to root.pep file')
    parser.add_argument('--calb_ref', required=True, help='Path to C. albicans reference proteome')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    parser.add_argument('--evalue', default='1e-10', help='E-value threshold for BLASTp (default: 1e-10)')
    parser.add_argument('--identity', type=float, default=50.0, help='Minimum percent identity for reporting (default: 50.0)')
    parser.add_argument('--coverage', type=float, default=50.0, help='Minimum query coverage for reporting (default: 50.0)')
    return parser.parse_args()

def run_blast(query, subject, output, evalue='1e-10', threads=4):
    """Run BLASTp between query and subject databases"""
    cmd = [
        'blastp',
        '-query', query,
        '-subject', subject,
        '-out', output,
        '-evalue', evalue,
        '-outfmt', '6 qseqid sseqid pident qcovs evalue bitscore',
        '-num_threads', str(threads)
    ]
    
    print(f"Running BLASTp: {' '.join(cmd)}")
    start_time = time.time()
    subprocess.run(cmd, check=True)
    end_time = time.time()
    print(f"BLASTp completed in {end_time - start_time:.2f} seconds")

def run_bidirectional_blast(root_pep, calb_ref, output_dir, evalue):
    """Run BLASTp in both directions and identify reciprocal best hits"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Forward blast: root.pep → calb_ref
    forward_output = os.path.join(output_dir, "root_to_calb.tsv")
    run_blast(root_pep, calb_ref, forward_output, evalue)
    
    # Reverse blast: calb_ref → root.pep
    reverse_output = os.path.join(output_dir, "calb_to_root.tsv")
    run_blast(calb_ref, root_pep, reverse_output, evalue)
    
    return forward_output, reverse_output

def load_blast_results(blast_file):
    """Load BLASTp results into a DataFrame"""
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'qcovs', 'evalue', 'bitscore'])
    
    return pd.read_csv(blast_file, sep='\t', 
                      names=['qseqid', 'sseqid', 'pident', 'qcovs', 'evalue', 'bitscore'])

def identify_best_hits(blast_df):
    """Identify best hits for each query sequence based on bit score"""
    if blast_df.empty:
        return pd.DataFrame(columns=blast_df.columns)
    
    # Sort by query ID and bit score (descending)
    sorted_df = blast_df.sort_values(['qseqid', 'bitscore'], ascending=[True, False])
    
    # Keep only the best hit for each query
    best_hits = sorted_df.drop_duplicates(subset=['qseqid'], keep='first')
    
    return best_hits

def find_reciprocal_best_hits(forward_best, reverse_best):
    """Identify reciprocal best hits between the two datasets"""
    rbh_list = []
    
    # Check each forward best hit
    for _, row in forward_best.iterrows():
        query_id = row['qseqid']
        subject_id = row['sseqid']
        
        # Find if there's a corresponding reverse hit
        reverse_hit = reverse_best[reverse_best['qseqid'] == subject_id]
        
        if not reverse_hit.empty and reverse_hit.iloc[0]['sseqid'] == query_id:
            # This is a reciprocal best hit
            rbh_list.append({
                'root_id': query_id,
                'calb_id': subject_id,
                'forward_pident': row['pident'],
                'forward_qcovs': row['qcovs'],
                'forward_evalue': row['evalue'],
                'forward_bitscore': row['bitscore'],
                'reverse_pident': reverse_hit.iloc[0]['pident'],
                'reverse_qcovs': reverse_hit.iloc[0]['qcovs'],
                'reverse_evalue': reverse_hit.iloc[0]['evalue'],
                'reverse_bitscore': reverse_hit.iloc[0]['bitscore']
            })
    
    return pd.DataFrame(rbh_list)

def analyze_results(forward_file, reverse_file, root_pep, calb_ref, output_dir, min_identity, min_coverage):
    """Analyze BLASTp results and generate statistics and visualizations"""
    print("Loading BLASTp results...")
    
    # Load BLAST results
    forward_df = load_blast_results(forward_file)
    reverse_df = load_blast_results(reverse_file)
    
    # Get best hits in each direction
    forward_best = identify_best_hits(forward_df)
    reverse_best = identify_best_hits(reverse_df)
    
    # Find reciprocal best hits
    rbh_df = find_reciprocal_best_hits(forward_best, reverse_best)
    
    # Count sequences in each file
    root_count = sum(1 for _ in SeqIO.parse(root_pep, "fasta"))
    calb_count = sum(1 for _ in SeqIO.parse(calb_ref, "fasta"))
    
    # Filter by identity and coverage thresholds
    forward_filtered = forward_df[(forward_df['pident'] >= min_identity) & 
                                  (forward_df['qcovs'] >= min_coverage)]
    reverse_filtered = reverse_df[(reverse_df['pident'] >= min_identity) & 
                                  (reverse_df['qcovs'] >= min_coverage)]
    
    # Count unique queries and subjects after filtering
    forward_queries = forward_filtered['qseqid'].nunique()
    forward_subjects = forward_filtered['sseqid'].nunique()
    reverse_queries = reverse_filtered['qseqid'].nunique()
    reverse_subjects = reverse_filtered['sseqid'].nunique()
    
    # Calculate statistics
    forward_match_rate = forward_queries / root_count if root_count > 0 else 0
    reverse_match_rate = reverse_queries / calb_count if calb_count > 0 else 0
    rbh_count = len(rbh_df)
    rbh_rate = rbh_count / min(root_count, calb_count) if min(root_count, calb_count) > 0 else 0
    
    # High quality matches (≥90% identity, ≥90% coverage)
    high_quality_matches = forward_df[(forward_df['pident'] >= 90) & (forward_df['qcovs'] >= 90)]
    high_quality_count = high_quality_matches['qseqid'].nunique()
    high_quality_rate = high_quality_count / root_count if root_count > 0 else 0
    
    # Prepare summary statistics
    stats = {
        'total_root_sequences': root_count,
        'total_calb_sequences': calb_count,
        'forward_matches': forward_queries,
        'forward_match_rate': forward_match_rate,
        'reverse_matches': reverse_queries,
        'reverse_match_rate': reverse_match_rate,
        'reciprocal_best_hits': rbh_count,
        'rbh_rate': rbh_rate,
        'high_quality_matches': high_quality_count,
        'high_quality_rate': high_quality_rate
    }
    
    # Save statistics to file
    stats_df = pd.DataFrame([stats])
    stats_file = os.path.join(output_dir, "blast_stats.csv")
    stats_df.to_csv(stats_file, index=False)
    
    # Save filtered results
    forward_filtered.to_csv(os.path.join(output_dir, "root_to_calb_filtered.csv"), index=False)
    reverse_filtered.to_csv(os.path.join(output_dir, "calb_to_root_filtered.csv"), index=False)
    
    # Save RBH results
    if not rbh_df.empty:
        rbh_df.to_csv(os.path.join(output_dir, "reciprocal_best_hits.csv"), index=False)
    
    # Create visualizations
    create_visualizations(forward_df, reverse_df, rbh_df, stats, output_dir)
    
    # Print summary
    print("\nBLASTp Analysis Summary:")
    print(f"Total sequences in root.pep: {root_count}")
    print(f"Total C. albicans reference proteins: {calb_count}")
    print(f"root.pep sequences with C. albicans matches (≥{min_identity}% identity, ≥{min_coverage}% coverage): "
          f"{forward_queries} ({forward_match_rate*100:.1f}%)")
    print(f"C. albicans proteins with root.pep matches (≥{min_identity}% identity, ≥{min_coverage}% coverage): "
          f"{reverse_queries} ({reverse_match_rate*100:.1f}%)")
    print(f"Reciprocal best hits: {rbh_count} ({rbh_rate*100:.1f}%)")
    print(f"High quality matches (≥90% identity, ≥90% coverage): {high_quality_count} ({high_quality_rate*100:.1f}%)")
    
    # Provide interpretation
    print("\nInterpretation:")
    if rbh_rate >= 0.9:
        print("✓ TRANSLATION QUALITY IS GOOD: Most C. albicans proteins have high-quality matches in root.pep")
        print("  The high F1 scores in orthology detection likely reflect this good translation quality.")
    elif rbh_rate <= 0.5:
        print("✗ TRANSLATION QUALITY ISSUES DETECTED: Many C. albicans proteins don't have good matches in root.pep")
        print("  The high F1 scores are surprising given the translation quality issues and may reflect")
        print("  other factors in the orthology detection process.")
    else:
        print("⚠ PARTIAL TRANSLATION ISSUES: Some C. albicans proteins have good matches in root.pep, but others don't")
        print("  The high F1 scores suggest that SYNERGY2's algorithm may be robust to these partial translation issues")
        print("  when determining orthology for C. albicans.")
    
    return stats

def create_visualizations(forward_df, reverse_df, rbh_df, stats, output_dir):
    """Create visualizations to help interpret BLASTp results"""
    # Set up plotting directory
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # Set seaborn style
    sns.set(style="whitegrid")
    
    # 1. Percent identity distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=forward_df, x='pident', bins=20, kde=True)
    plt.title('Distribution of Percent Identity (root.pep → C. albicans)')
    plt.xlabel('Percent Identity')
    plt.ylabel('Count')
    plt.axvline(x=90, color='red', linestyle='--', label='90% Identity Threshold')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, "identity_distribution.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Query coverage distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=forward_df, x='qcovs', bins=20, kde=True)
    plt.title('Distribution of Query Coverage (root.pep → C. albicans)')
    plt.xlabel('Query Coverage')
    plt.ylabel('Count')
    plt.axvline(x=90, color='red', linestyle='--', label='90% Coverage Threshold')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, "coverage_distribution.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Scatter plot of identity vs. coverage
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=forward_df, x='pident', y='qcovs', alpha=0.5, s=10)
    plt.title('Percent Identity vs. Query Coverage (root.pep → C. albicans)')
    plt.xlabel('Percent Identity')
    plt.ylabel('Query Coverage')
    plt.axhline(y=90, color='red', linestyle='--')
    plt.axvline(x=90, color='red', linestyle='--')
    plt.savefig(os.path.join(plot_dir, "identity_vs_coverage.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Reciprocal best hits identity comparison
    if not rbh_df.empty:
        plt.figure(figsize=(10, 6))
        plt.scatter(rbh_df['forward_pident'], rbh_df['reverse_pident'], alpha=0.5, s=10)
        plt.title('Reciprocal Best Hits: Forward vs. Reverse Percent Identity')
        plt.xlabel('Forward Percent Identity (root.pep → C. albicans)')
        plt.ylabel('Reverse Percent Identity (C. albicans → root.pep)')
        plt.plot([0, 100], [0, 100], 'r--')
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.savefig(os.path.join(plot_dir, "rbh_identity_comparison.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 5. Summary bar chart
    plt.figure(figsize=(12, 7))
    metrics = ['forward_match_rate', 'reverse_match_rate', 'rbh_rate', 'high_quality_rate']
    values = [stats[m] * 100 for m in metrics]
    labels = ['root.pep → C. albicans\nMatch Rate', 
              'C. albicans → root.pep\nMatch Rate',
              'Reciprocal Best Hit\nRate',
              'High Quality Match\nRate (≥90% id/cov)']
    
    colors = ['#3498db', '#2ecc71', '#f39c12', '#9b59b6']
    threshold_color = '#e74c3c'
    
    bars = plt.bar(labels, values, color=colors)
    
    # Add a horizontal line at 90%
    plt.axhline(y=90, color=threshold_color, linestyle='--', label='90% Threshold')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                 f'{height:.1f}%', ha='center', va='bottom')
    
    plt.ylim(0, 105)
    plt.ylabel('Percentage (%)')
    plt.title('BLASTp Analysis Summary Metrics')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "summary_metrics.png"), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    args = parse_arguments()
    
    # Run bidirectional BLASTp
    forward_file, reverse_file = run_bidirectional_blast(
        args.root_pep, args.calb_ref, args.output_dir, args.evalue
    )
    
    # Analyze results
    stats = analyze_results(
        forward_file, reverse_file, args.root_pep, args.calb_ref, 
        args.output_dir, args.identity, args.coverage
    )
    
    print(f"\nAll analysis results saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
