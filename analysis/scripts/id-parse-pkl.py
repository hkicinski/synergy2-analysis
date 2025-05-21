#!/usr/bin/env python3
"""
id-parse-pkl.py - Parse SYNERGY2 cluster output files with built-in sanity checks
Author: Hubert Kicinski
Date: May 20, 2025
"""

import pandas as pd
import pickle
import re
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# Input and output paths
clusters_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis_run2/nodes/root/final_clusters.txt"
locus_mapping_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis_run2/nodes/root/locus_mappings.pkl"
output_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-parsed"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Create a sanity check log file
sanity_log = open(f"{output_dir}/sanity_check_log.txt", "w")
def log(message):
    """Write to both console and log file"""
    print(message)
    sanity_log.write(message + "\n")

log(f"SYNERGY2 Cluster Parser with Sanity Checks")
log(f"==========================================")
log(f"Input file: {clusters_file}")
log(f"Output directory: {output_dir}")
log(f"Started at: {pd.Timestamp.now()}")
log(f"")

# SANITY CHECK 1: Verify input file exists
if not os.path.exists(clusters_file):
    log(f"ERROR: Input file not found: {clusters_file}")
    sanity_log.close()
    sys.exit(1)

# Parse the final_clusters.txt file correctly
log(f"Parsing cluster file...")
clusters = []
line_count = 0
multi_gene_clusters = 0
multi_species_clusters = 0

# Sample lines for validation
sample_lines = []

with open(clusters_file, 'r') as f:
    for line_number, line in enumerate(f, 1):
        line_count += 1
        if line_number <= 5:  # Keep first 5 lines for sample
            sample_lines.append(line.strip())
            
        # Handle different possible delimiters - tab or multiple spaces
        parts = [p for p in re.split(r'\t+|\s{2,}', line.strip()) if p]
        
        if len(parts) >= 2:
            # First part contains cluster ID and metadata
            cluster_info = parts[0]
            # Second part contains space-separated genes
            genes_text = parts[1]
            
            # Extract cluster ID
            cluster_match = re.match(r'(Cluster\d+)', cluster_info)
            if cluster_match:
                cluster_id = cluster_match.group(1)
                
                # Split genes by space
                genes = genes_text.split()
                
                # Count species in this cluster
                species_in_cluster = set()
                for gene in genes:
                    if re.search(r'CAGL0', gene):
                        species_in_cluster.add("C_glabrata")
                    elif re.search(r'orf19', gene):
                        species_in_cluster.add("C_albicans")
                    elif re.search(r'KLLA0', gene):
                        species_in_cluster.add("K_lactis")
                    elif re.search(r'gene-Y[A-Z][A-Z]', gene) or re.search(r'Y[A-Z][A-Z]\d', gene):
                        species_in_cluster.add("S_cerevisiae")
                
                # Count multiple gene and species clusters
                if len(genes) > 1:
                    multi_gene_clusters += 1
                if len(species_in_cluster) > 1:
                    multi_species_clusters += 1
                
                clusters.append({
                    'cluster_id': cluster_id,
                    'genes': genes_text,
                    'gene_list': genes,
                    'size': len(genes),
                    'species_count': len(species_in_cluster)
                })

# SANITY CHECK 2: Verify parsing produced results
if len(clusters) == 0:
    log(f"ERROR: No clusters were parsed from the input file.")
    log(f"Sample of the input file:")
    for sample in sample_lines:
        log(f"  {sample}")
    log(f"Please check the file format and update the parsing logic if needed.")
    sanity_log.close()
    sys.exit(1)

# Convert to DataFrame
clusters_df = pd.DataFrame(clusters)
log(f"Parsed {len(clusters_df)} clusters with {clusters_df['size'].sum()} total genes")

# SANITY CHECK 3: Verify multiple genes per cluster
if multi_gene_clusters == 0:
    log(f"WARNING: No clusters with multiple genes were detected. This is highly unusual!")
    log(f"Sample of the input file:")
    for sample in sample_lines:
        log(f"  {sample}")
    log(f"The current parsing logic may not be correctly identifying multiple genes per cluster.")
    
    # Don't exit since one-gene-per-cluster is technically possible but unlikely
else:
    log(f"Found {multi_gene_clusters} clusters ({multi_gene_clusters/len(clusters_df)*100:.2f}%) with multiple genes")

# SANITY CHECK 4: Verify multiple species per cluster
if multi_species_clusters == 0:
    log(f"WARNING: No clusters with multiple species were detected. This is highly unusual!")
    log(f"The current parsing logic may not be correctly identifying species from gene IDs.")
    
    # Don't exit since one-species-per-cluster is technically possible but unlikely
else:
    log(f"Found {multi_species_clusters} clusters ({multi_species_clusters/len(clusters_df)*100:.2f}%) with multiple species")

# Extract all genes and identify species
all_genes = []
species_pattern_matches = {
    "C_glabrata": 0,
    "C_albicans": 0,
    "K_lactis": 0,
    "S_cerevisiae": 0,
    "unknown": 0
}

for _, row in clusters_df.iterrows():
    cluster_id = row['cluster_id']
    
    for gene_id in row['gene_list']:
        # Identify species based on gene ID pattern
        if re.search(r'CAGL0', gene_id):
            species = "C_glabrata"
            species_pattern_matches[species] += 1
        elif re.search(r'orf19', gene_id):
            species = "C_albicans"
            species_pattern_matches[species] += 1
        elif re.search(r'KLLA0', gene_id):
            species = "K_lactis"
            species_pattern_matches[species] += 1
        elif re.search(r'gene-Y[A-Z][A-Z]', gene_id) or re.search(r'Y[A-Z][A-Z]\d', gene_id):
            species = "S_cerevisiae"
            species_pattern_matches[species] += 1
        else:
            species = "unknown"
            species_pattern_matches[species] += 1
        
        all_genes.append({
            'cluster_id': cluster_id,
            'original_id': gene_id,
            'species': species
        })

# Create genes dataframe
genes_df = pd.DataFrame(all_genes)
log(f"Extracted {len(genes_df)} genes across {len(clusters_df)} clusters")

# SANITY CHECK 5: Verify all species were detected
missing_species = []
low_species = []
expected_species = ["C_glabrata", "C_albicans", "K_lactis", "S_cerevisiae"]

for species in expected_species:
    if species_pattern_matches[species] == 0:
        missing_species.append(species)
    elif species_pattern_matches[species] < 100:  # Arbitrary low threshold
        low_species.append((species, species_pattern_matches[species]))

if missing_species:
    log(f"WARNING: The following expected species were not detected: {', '.join(missing_species)}")
    log(f"Please check the species identification patterns in the script.")

if low_species:
    log(f"WARNING: The following species have unusually low counts:")
    for species, count in low_species:
        log(f"  {species}: {count} genes")
    log(f"This might indicate issues with the species identification patterns.")

# Calculate species distribution
species_counts = genes_df['species'].value_counts()
log(f"\nSpecies distribution:")
for species, count in species_counts.items():
    percentage = count / len(genes_df) * 100
    log(f"  {species}: {count} genes ({percentage:.2f}%)")

# Count species per cluster
species_per_cluster = genes_df.groupby('cluster_id')['species'].nunique()
log(f"\nClusters by number of species:")
for species_count, count in species_per_cluster.value_counts().sort_index().items():
    percentage = count / len(clusters_df) * 100
    log(f"  {species_count} species: {count} clusters ({percentage:.2f}%)")

# SANITY CHECK 6: Verify single-copy orthologs
# Count clusters with exactly one gene from each species
single_copy_count = 0
for cluster_id, group in genes_df.groupby('cluster_id'):
    species_counts = group['species'].value_counts()
    if all(species_counts.get(species, 0) == 1 for species in expected_species):
        single_copy_count += 1

log(f"\nSingle-copy orthologs (one gene from each species): {single_copy_count} clusters")
if single_copy_count == 0:
    log(f"WARNING: No single-copy orthologs were detected. This is unusual for related species.")
    log(f"This might indicate issues with the parsing or species identification.")
else:
    log(f"Single-copy orthologs represent {single_copy_count/len(clusters_df)*100:.2f}% of all clusters")

# SANITY CHECK 7: Generate validation visualizations
log(f"\nGenerating validation visualizations...")

# 1. Cluster size distribution
plt.figure(figsize=(10, 6))
sns.histplot(clusters_df['size'], bins=20, kde=False)
plt.title('Distribution of Cluster Sizes')
plt.xlabel('Number of Genes in Cluster')
plt.ylabel('Count of Clusters')
plt.savefig(f"{output_dir}/validation_cluster_sizes.png")
plt.close()

# 2. Species distribution
plt.figure(figsize=(10, 6))
sns.barplot(x=species_counts.index, y=species_counts.values)
plt.title('Number of Genes by Species')
plt.xlabel('Species')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_dir}/validation_species_distribution.png")
plt.close()

# 3. Species count distribution
plt.figure(figsize=(10, 6))
species_count_dist = species_per_cluster.value_counts().sort_index()
sns.barplot(x=species_count_dist.index.astype(str), y=species_count_dist.values)
plt.title('Clusters by Number of Species')
plt.xlabel('Number of Species in Cluster')
plt.ylabel('Count of Clusters')
plt.savefig(f"{output_dir}/validation_species_per_cluster.png")
plt.close()

# Save output files
log(f"\nSaving output files to {output_dir}")
clusters_df.to_csv(f"{output_dir}/parsed_clusters.csv", index=False)
genes_df.to_csv(f"{output_dir}/genes_by_cluster.csv", index=False)

# Create a mapping file (simplified since we're using original IDs directly)
id_mapping = genes_df[['original_id', 'species']].drop_duplicates()
id_mapping.to_csv(f"{output_dir}/id_mapping.csv", index=False)

# SANITY CHECK 8: Final validation summary
log(f"\nSANITY CHECK SUMMARY:")
log(f"  1. Input file exists and was parsed: ✓")
log(f"  2. Parsed {len(clusters_df)} clusters: ✓")
log(f"  3. Multi-gene clusters: {multi_gene_clusters > 0 and '✓' or '✗'}")
log(f"  4. Multi-species clusters: {multi_species_clusters > 0 and '✓' or '✗'}")
log(f"  5. All expected species detected: {(not missing_species) and '✓' or '✗'}")
log(f"  6. Single-copy orthologs detected: {single_copy_count > 0 and '✓' or '✗'}")
log(f"  7. Validation visualizations generated: ✓")
log(f"  8. Output files created: ✓")

# Overall pass/fail assessment
major_issues = len(missing_species) > 0 or multi_gene_clusters == 0 or multi_species_clusters == 0
if major_issues:
    log(f"\nWARNING: Major issues were detected in the parsing results!")
    log(f"Please review the warnings above and correct any issues before proceeding.")
else:
    log(f"\nAll sanity checks PASSED. Results appear valid and ready for further analysis.")

log(f"\nCompleted at: {pd.Timestamp.now()}")
sanity_log.close()

# Create a simple report file for quick reference
with open(f"{output_dir}/parse_summary.md", "w") as f:
    f.write("# SYNERGY2 Cluster Parsing Summary\n\n")
    f.write(f"- **Total clusters:** {len(clusters_df)}\n")
    f.write(f"- **Total genes:** {len(genes_df)}\n")
    f.write(f"- **Multi-gene clusters:** {multi_gene_clusters} ({multi_gene_clusters/len(clusters_df)*100:.2f}%)\n")
    f.write(f"- **Multi-species clusters:** {multi_species_clusters} ({multi_species_clusters/len(clusters_df)*100:.2f}%)\n")
    f.write(f"- **Single-copy orthologs:** {single_copy_count} ({single_copy_count/len(clusters_df)*100:.2f}%)\n\n")
    
    f.write("## Species Distribution\n\n")
    f.write("| Species | Gene Count | Percentage |\n")
    f.write("|---------|------------|------------|\n")
    for species, count in species_counts.items():
        percentage = count / len(genes_df) * 100
        f.write(f"| {species} | {count} | {percentage:.2f}% |\n")
    
    f.write("\n## Clusters by Number of Species\n\n")
    f.write("| Species Count | Cluster Count | Percentage |\n")
    f.write("|--------------|---------------|------------|\n")
    for species_count, count in species_per_cluster.value_counts().sort_index().items():
        percentage = count / len(clusters_df) * 100
        f.write(f"| {species_count} | {count} | {percentage:.2f}% |\n")
    
    f.write("\n## Validation Status\n\n")
    f.write(f"Status: {'All sanity checks PASSED' if not major_issues else 'ISSUES DETECTED'}\n")
    f.write(f"Ready for next analysis step: {'YES' if not major_issues else 'NO'}\n")

print("\nDone! Files saved in", output_dir)
print(f"Check {output_dir}/sanity_check_log.txt for detailed validation results")
print(f"Check {output_dir}/parse_summary.md for a summary of the parsing results")
