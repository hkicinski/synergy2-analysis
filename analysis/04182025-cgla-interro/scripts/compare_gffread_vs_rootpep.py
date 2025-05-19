#!/usr/bin/env python3
"""
compare_gffread_vs_rootpep.py - Compare gffread-generated proteins with root.pep
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
    parser = argparse.ArgumentParser(description='Compare gffread-generated proteins with root.pep')
    parser.add_argument('--gffread', required=True, help='Path to gffread-generated protein FASTA')
    parser.add_argument('--root_pep', required=True, help='Path to root.pep file')
    parser.add_argument('--output_dir', required=True, help='Directory for output files')
    parser.add_argument('--species_id', default='CAGL0', help='Species identifier to filter in root.pep')
    return parser.parse_args()

def extract_species_proteins(input_file, output_file, species_id):
    """Extract species-specific proteins from root.pep"""
    records_count = 0
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(input_file, "fasta"):
            if species_id in record.id:
                out_f.write(f">{record.id}\n{record.seq}\n")
                records_count += 1
    
    print(f"Extracted {records_count} {species_id} proteins to {output_file}")
    return records_count

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

def analyze_blast_results(blast_file, output_dir):
    """Analyze BLASTp results and generate statistics and visualizations"""
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
    best_hits.to_csv(os.path.join(output_dir, "best_hits.csv"), index=False)
    
    # Create percent identity histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(best_hits['pident'], bins=20, kde=True)
    plt.axvline(x=90, color='red', linestyle='--', label='90% Identity Threshold')
    plt.title('Distribution of Percent Identity')
    plt.xlabel('Percent Identity')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, "percent_identity_distribution.png"))
    plt.close()
    
    # Create query coverage histogram
    plt.figure(figsize=(10, 6))
    sns.histplot(best_hits['qcovs'], bins=20, kde=True)
    plt.axvline(x=90, color='red', linestyle='--', label='90% Coverage Threshold')
    plt.title('Distribution of Query Coverage')
    plt.xlabel('Query Coverage')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, "query_coverage_distribution.png"))
    plt.close()
    
    # Create scatter plot of identity vs coverage
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='pident', y='qcovs', data=best_hits, alpha=0.5)
    plt.axhline(y=90, color='red', linestyle='--')
    plt.axvline(x=90, color='red', linestyle='--')
    plt.title('Percent Identity vs. Query Coverage')
    plt.xlabel('Percent Identity')
    plt.ylabel('Query Coverage')
    plt.savefig(os.path.join(plot_dir, "identity_vs_coverage.png"))
    plt.close()
    
    # Create summary bar chart
    plt.figure(figsize=(10, 6))
    metrics = ['Match Rate', 'High Quality Match\nRate (≥90% id/cov)']
    values = [match_rate, high_quality_rate]
    
    bars = plt.bar(metrics, values, color=['#3498db', '#9b59b6'])
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                 f'{height:.1f}%', ha='center', va='bottom')
    
    plt.axhline(y=90, color='red', linestyle='--', label='90% Threshold')
    plt.ylim(0, 105)
    plt.ylabel('Percentage (%)')
    plt.title('BLASTp Analysis Summary Metrics')
    plt.legend()
    plt.savefig(os.path.join(plot_dir, "summary_metrics.png"))
    plt.close()
    
    # Print summary
    print("\nBLASTp Analysis Summary:")
    print(f"Total queries: {total_queries}")
    print(f"Queries with hits: {matched_queries} ({match_rate:.2f}%)")
    print(f"High quality matches (≥90% identity & coverage): {high_quality_matches} ({high_quality_rate:.2f}% of matches)")
    
    # Provide interpretation
    print("\nInterpretation:")
    if high_quality_rate >= 90:
        print("✓ EXCELLENT MATCH: The vast majority of proteins have high-quality matches")
        print("  This suggests the gffread translations match very well with SYNERGY2's proteins.")
    elif high_quality_rate >= 50:
        print("⚠ MODERATE MATCH: Many proteins match well, but there are substantial differences")
        print("  This suggests some systematic differences in how proteins were translated.")
    else:
        print("✗ POOR MATCH: Few proteins match well between the two sets")
        print("  This indicates significant differences in protein translations, which may")
        print("  explain the orthology discrepancies observed with C. glabrata.")
    
    return {
        'total_queries': total_queries,
        'matched_queries': matched_queries,
        'high_quality_matches': high_quality_matches,
        'match_rate': match_rate,
        'high_quality_rate': high_quality_rate
    }

def main():
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Extract C. glabrata proteins from root.pep
    root_species_file = os.path.join(args.output_dir, "root_species_proteins.fasta")
    extract_species_proteins(args.root_pep, root_species_file, args.species_id)
    
    # Run BLASTp: gffread proteins vs root.pep C. glabrata proteins
    blast_output = os.path.join(args.output_dir, "gffread_vs_root_blast.tsv")
    run_blastp(args.gffread, root_species_file, blast_output)
    
    # Analyze results
    analyze_blast_results(blast_output, args.output_dir)
    
    print(f"\nAll analysis results saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
