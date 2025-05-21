#!/usr/bin/env python3
"""
validate_synergy_parsing.py - Independent validation script for SYNERGY2 outputs
Author: Hubert Kicinski
Date: May 20, 2025

This script performs independent verification of SYNERGY2 parsing results
by checking the parsed output against the original input file.
"""

import pandas as pd
import re
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Input files
    original_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis_run2/nodes/root/final_clusters.txt"
    parsed_clusters_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-parsed/parsed_clusters.csv"
    genes_by_cluster_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-parsed/genes_by_cluster.csv"

    # Output directory
    validation_dir = os.path.dirname(parsed_clusters_file) + "/validation"
    os.makedirs(validation_dir, exist_ok=True)

    #validation report as MD
    report_file = f"{validation_dir}/validation_report.md"
    with open(report_file, "w") as report:
        report.write("# SYNERGY2 Parsing Validation Report\n\n")
        report.write(f"Generated: {pd.Timestamp.now()}\n\n")

        #debug
        report.write("## File Existence Check\n\n")
        files_exist = True
        for file_path in [original_file, parsed_clusters_file, genes_by_cluster_file]:
            exists = os.path.exists(file_path)
            report.write(f"- {os.path.basename(file_path)}: {'✓' if exists else '✗'}\n")
            files_exist = files_exist and exists

        if not files_exist:
            report.write("\n**ERROR**: One or more required files are missing.\n")
            print("ERROR: One or more required files are missing. See validation report for details.")
            return

        report.write("\n## Independent Original File Analysis\n\n")

        # analyze the original file (final_clusters)
        with open(original_file, 'r') as f:
            original_lines = f.readlines()

        total_original_lines = len(original_lines)
        report.write(f"Total lines in original file: {total_original_lines}\n\n")

        #  original file structure
        original_cluster_sizes = []
        original_gene_count = 0
        sample_original_clusters = []

        for i, line in enumerate(original_lines):
            if i < 5:  # Keep first 5 lines as samples
                sample_original_clusters.append(line.strip())

            # parsed line to count genes
            parts = [p for p in re.split(r'\t+|\s{2,}', line.strip()) if p]
            if len(parts) >= 2:
                genes_text = parts[1]
                genes = genes_text.split()
                original_cluster_sizes.append(len(genes))
                original_gene_count += len(genes)

        # reported original file metrics
        avg_original_cluster_size = sum(original_cluster_sizes) / len(original_cluster_sizes) if original_cluster_sizes else 0
        report.write(f"Original file metrics:\n")
        report.write(f"- Total clusters: {len(original_cluster_sizes)}\n")
        report.write(f"- Total genes: {original_gene_count}\n")
        report.write(f"- Average genes per cluster: {avg_original_cluster_size:.2f}\n\n")

        report.write("Sample lines from original file:\n")
        for sample in sample_original_clusters:
            report.write(f"```\n{sample}\n```\n")

        # loaded parsed data
        report.write("\n## Parsed Data Analysis\n\n")

        try:
            parsed_clusters = pd.read_csv(parsed_clusters_file)
            genes_by_cluster = pd.read_csv(genes_by_cluster_file)
        except Exception as e:
            report.write(f"**ERROR**: Failed to read parsed files: {e}\n")
            print(f"ERROR: Failed to read parsed files: {e}")
            return

        # reported parsed file metrics
        total_parsed_clusters = len(parsed_clusters)
        total_parsed_genes = len(genes_by_cluster)
        avg_parsed_cluster_size = parsed_clusters['size'].mean() if 'size' in parsed_clusters.columns else 0

        report.write(f"Parsed data metrics:\n")
        report.write(f"- Total clusters: {total_parsed_clusters}\n")
        report.write(f"- Total genes: {total_parsed_genes}\n")
        report.write(f"- Average genes per cluster: {avg_parsed_cluster_size:.2f}\n\n")

        # debug: Count Validation
        report.write("## Count Validation\n\n")

        cluster_count_match = len(original_cluster_sizes) == total_parsed_clusters
        gene_count_match = original_gene_count == total_parsed_genes

        report.write(f"- Cluster count match: {'✓' if cluster_count_match else '✗'} (Original: {len(original_cluster_sizes)}, Parsed: {total_parsed_clusters})\n")
        report.write(f"- Gene count match: {'✓' if gene_count_match else '✗'} (Original: {original_gene_count}, Parsed: {total_parsed_genes})\n\n")

        if not cluster_count_match or not gene_count_match:
            report.write("**WARNING**: Count mismatch between original and parsed data.\n")

        # Cluster Size Distribution
        report.write("## Cluster Size Distribution\n\n")

        # Original size distribution (i.e. from the final_clusters as a sanity check)
        original_size_counts = {}
        for size in original_cluster_sizes:
            original_size_counts[size] = original_size_counts.get(size, 0) + 1

        # Parsed size distribution if 'size' column exists
        parsed_size_counts = {}
        if 'size' in parsed_clusters.columns:
            for size in parsed_clusters['size']:
                parsed_size_counts[size] = parsed_size_counts.get(size, 0) + 1

        # Compare distributions
        report.write("### Size Distribution Comparison\n\n")
        report.write("| Size | Original Count | Parsed Count |\n")
        report.write("|------|---------------|-------------|\n")

        all_sizes = sorted(set(list(original_size_counts.keys()) + list(parsed_size_counts.keys())))
        for size in all_sizes:
            orig_count = original_size_counts.get(size, 0)
            parsed_count = parsed_size_counts.get(size, 0)
            report.write(f"| {size} | {orig_count} | {parsed_count} |\n")

        # Check for single-gene clusters
        orig_single_gene = original_size_counts.get(1, 0)
        parsed_single_gene = parsed_size_counts.get(1, 0)

        orig_single_gene_pct = orig_single_gene / len(original_cluster_sizes) * 100 if original_cluster_sizes else 0
        parsed_single_gene_pct = parsed_single_gene / total_parsed_clusters * 100 if total_parsed_clusters > 0 else 0

        report.write(f"\n**Single-gene clusters**:\n")
        report.write(f"- Original: {orig_single_gene} ({orig_single_gene_pct:.2f}%)\n")
        report.write(f"- Parsed: {parsed_single_gene} ({parsed_single_gene_pct:.2f}%)\n\n")

        # Critical check: If all clusters are parsed as single-gene but original has multi-gene clusters
        if parsed_single_gene == total_parsed_clusters and orig_single_gene < len(original_cluster_sizes):
            report.write("**CRITICAL ERROR**: All parsed clusters contain only one gene, but original data has multi-gene clusters.\n")
            report.write("This suggests a parsing error that failed to properly split gene lists.\n\n")

        # Visualize size distributions
        plt.figure(figsize=(10, 6))

        # Plot original sizes
        sizes, counts = zip(*sorted(original_size_counts.items()))
        plt.bar([str(s) for s in sizes], counts, alpha=0.5, label='Original')

        # Plot parsed sizes
        if parsed_size_counts:
            sizes, counts = zip(*sorted(parsed_size_counts.items()))
            plt.bar([str(s) for s in sizes], counts, alpha=0.5, label='Parsed')

        plt.title('Cluster Size Distribution Comparison')
        plt.xlabel('Cluster Size (Number of Genes)')
        plt.ylabel('Count of Clusters')
        plt.legend()
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{validation_dir}/size_distribution_comparison.png")
        plt.close()

        report.write(f"![Size Distribution Comparison]({os.path.basename(validation_dir)}/size_distribution_comparison.png)\n\n")

        # Gene-to-Cluster Ratio
        report.write("## Gene-to-Cluster Ratio\n\n")

        orig_ratio = original_gene_count / len(original_cluster_sizes) if original_cluster_sizes else 0
        parsed_ratio = total_parsed_genes / total_parsed_clusters if total_parsed_clusters > 0 else 0

        ratio_match = abs(orig_ratio - parsed_ratio) < 0.1  # Allow small difference due to rounding

        report.write(f"- Original ratio: {orig_ratio:.2f} genes per cluster\n")
        report.write(f"- Parsed ratio: {parsed_ratio:.2f} genes per cluster\n")
        report.write(f"- Ratio match: {'✓' if ratio_match else '✗'}\n\n")

        if not ratio_match:
            report.write("**WARNING**: Gene-to-cluster ratio differs significantly between original and parsed data.\n")
            report.write("This suggests that parsing may have split clusters incorrectly or failed to capture all genes.\n\n")

        # Final validation verdict
        report.write("## Validation Verdict\n\n")

        critical_errors = []
        warnings = []

        # Check for critical errors
        if not files_exist:
            critical_errors.append("One or more required files are missing")

        if parsed_single_gene == total_parsed_clusters and orig_single_gene < len(original_cluster_sizes):
            critical_errors.append("All parsed clusters contain only one gene, but original data has multi-gene clusters")

        # Check for warnings
        if not cluster_count_match:
            warnings.append("Cluster count mismatch")

        if not gene_count_match:
            warnings.append("Gene count mismatch")

        if not ratio_match:
            warnings.append("Gene-to-cluster ratio mismatch")

        # Final verdict
        if critical_errors:
            report.write("**VALIDATION FAILED: CRITICAL ERRORS DETECTED**\n\n")
            report.write("Critical errors:\n")
            for error in critical_errors:
                report.write(f"- {error}\n")
        elif warnings:
            report.write("**VALIDATION PASSED WITH WARNINGS**\n\n")
            report.write("Warnings:\n")
            for warning in warnings:
                report.write(f"- {warning}\n")
        else:
            report.write("**VALIDATION PASSED SUCCESSFULLY**\n\n")
            report.write("The parsed results match the expected patterns from the original data.\n")

        report.write("\nThis validation was performed by directly comparing the parsed results against an independent analysis of the original file.\n")

    # Print summary to console
    print(f"Validation complete! Report saved to {report_file}")

    if critical_errors:
        print("VALIDATION FAILED: CRITICAL ERRORS DETECTED")
        for error in critical_errors:
            print(f"- {error}")
        sys.exit(1)
    elif warnings:
        print("VALIDATION PASSED WITH WARNINGS")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("VALIDATION PASSED SUCCESSFULLY")
        print("The parsed results match the expected patterns from the original data.")

if __name__ == "__main__":
    main()
