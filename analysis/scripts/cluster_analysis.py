#!/usr/bin/env python3
"""
cluster_analysis.py - Analyze SYNERGY2 clusters using output from id-parse-pkl.py
Author: Hubert Kicinski
Date: May 20, 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import numpy as np

# Input and output paths
input_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-parsed"
output_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro"
os.makedirs(output_dir, exist_ok=True)

# Create a log file
log_file = open(f"{output_dir}/cluster_analysis_log.txt", "w")

def log(message):
    """Write to both console and log file"""
    print(message)
    log_file.write(message + "\n")

log(f"SYNERGY2 Cluster Analysis")
log(f"=======================")
log(f"Input directory: {input_dir}")
log(f"Output directory: {output_dir}")
log(f"Started at: {pd.Timestamp.now()}")
log("")

# Validate input data format
def validate_input_data(genes_df):
    """Verify the genes dataframe has expected structure and content"""
    expected_columns = ['cluster_id', 'original_id', 'species']
    missing_columns = [col for col in expected_columns if col not in genes_df.columns]
    
    if missing_columns:
        log(f"ERROR: Input data missing required columns: {', '.join(missing_columns)}")
        log_file.close()
        sys.exit(1)
    
    # Check species composition
    species_counts = genes_df['species'].value_counts()
    expected_species = ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']
    missing_species = [sp for sp in expected_species if sp not in species_counts.index]
    
    if missing_species:
        log(f"WARNING: Input data missing expected species: {', '.join(missing_species)}")
    
    # Check multi-gene clusters
    genes_per_cluster = genes_df.groupby('cluster_id').size()
    single_gene_clusters = sum(genes_per_cluster == 1)
    
    if single_gene_clusters == len(genes_per_cluster):
        log("ERROR: All clusters contain only one gene. This suggests parsing issues.")
        log_file.close()
        sys.exit(1)
    
    multi_gene_percentage = (len(genes_per_cluster) - single_gene_clusters) / len(genes_per_cluster) * 100
    log(f"Validation: {multi_gene_percentage:.2f}% of clusters contain multiple genes")
    
    # Input data looks valid
    return True

# Load the gene data from id-parse-pkl.py output
log("Loading parsed data...")
genes_file = os.path.join(input_dir, "genes_by_cluster.csv")

if not os.path.exists(genes_file):
    log(f"ERROR: Input file {genes_file} not found.")
    log("Make sure to run id-parse-pkl.py first to generate the required input files.")
    log_file.close()
    sys.exit(1)

try:
    genes_df = pd.read_csv(genes_file)
    log(f"Loaded {len(genes_df)} genes across {genes_df['cluster_id'].nunique()} clusters")
except Exception as e:
    log(f"ERROR loading input file: {e}")
    log_file.close()
    sys.exit(1)

# Validate the loaded data
validate_input_data(genes_df)

# Step 1: Calculate species composition of each cluster
log("\nAnalyzing cluster compositions...")
cluster_composition = genes_df.groupby('cluster_id')['species'].value_counts().unstack(fill_value=0)

# Add species count column
cluster_composition['species_count'] = (cluster_composition > 0).sum(axis=1)

# Save the composition data
cluster_composition.to_csv(f"{output_dir}/cluster_compositions.csv")
log(f"Saved cluster composition data to {output_dir}/cluster_compositions.csv")

# Step 2: Classify clusters by composition type
log("\nClassifying clusters by composition pattern...")

def classify_cluster(row):
    """Classify cluster based on species presence and counts"""
    species_present = set()
    for species in ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']:
        if species in row and row[species] > 0:
            species_present.add(species)
    
    # Single-copy orthologs (one gene from each species)
    if (row.get('C_albicans', 0) == 1 and row.get('C_glabrata', 0) == 1 and 
        row.get('K_lactis', 0) == 1 and row.get('S_cerevisiae', 0) == 1):
        return "single_copy_ortholog"
    
    # Multi-copy orthologs (at least one gene from each species)
    elif len(species_present) == 4:
        return "multi_copy_ortholog"
    
    # Partial orthologs (missing one or more species)
    elif len(species_present) > 1:
        return "partial_ortholog"
    
    # Species-specific clusters
    elif len(species_present) == 1:
        # Return the species name
        return list(species_present)[0] + "_specific"
    
    return "other"

# Apply classification
cluster_composition['classification'] = cluster_composition.apply(classify_cluster, axis=1)

# Count clusters by classification
classification_counts = cluster_composition['classification'].value_counts()
log("Cluster classifications:")
for category, count in classification_counts.items():
    percentage = count / len(cluster_composition) * 100
    log(f"  {category}: {count} clusters ({percentage:.2f}%)")

# Save classifications
cluster_composition.to_csv(f"{output_dir}/cluster_classifications.csv")
log(f"Saved cluster classifications to {output_dir}/cluster_classifications.csv")

# Step 3: Analyze single-copy orthologs (the most reliable relationships)
log("\nAnalyzing single-copy orthologs...")
single_copy_clusters = cluster_composition[
    cluster_composition['classification'] == 'single_copy_ortholog'
].index.tolist()

# Extract all genes in single-copy ortholog clusters
single_copy_genes = genes_df[genes_df['cluster_id'].isin(single_copy_clusters)]

# Create a pivot table for easy reference
single_copy_pivot = single_copy_genes.pivot_table(
    index='cluster_id', 
    columns='species', 
    values='original_id', 
    aggfunc='first'
)

# Save single-copy orthologs
single_copy_pivot.to_csv(f"{output_dir}/single_copy_orthologs.csv")
log(f"Saved {len(single_copy_pivot)} single-copy ortholog clusters to {output_dir}/single_copy_orthologs.csv")

# Step 4: Analyze species-specific expansions
log("\nAnalyzing species-specific expansions...")
for species in ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']:
    if species not in cluster_composition.columns:
        log(f"  Warning: Species {species} not found in data")
        continue
        
    # Find clusters with 3 or more genes from this species
    expanded_clusters = cluster_composition[cluster_composition[species] >= 3].sort_values(by=species, ascending=False)
    
    if len(expanded_clusters) > 0:
        log(f"  {species}: {len(expanded_clusters)} clusters with gene expansions")
        
        # Save the top 20 expansions
        top_expansions = expanded_clusters.head(20)
        top_expansions.to_csv(f"{output_dir}/{species}_expansions.csv")
        log(f"  Saved top expansions to {output_dir}/{species}_expansions.csv")
        
        # Extract genes in these expanded clusters
        expansion_genes = genes_df[
            (genes_df['cluster_id'].isin(top_expansions.index)) & 
            (genes_df['species'] == species)
        ]
        expansion_genes.to_csv(f"{output_dir}/{species}_expansion_genes.csv", index=False)

# Step 5: Analyze ortholog conservation patterns
log("\nAnalyzing ortholog conservation patterns...")

# Calculate the percentage of genes with orthologs in each species
species_with_orthologs = {}
for species in ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']:
    if species not in genes_df['species'].unique():
        continue
        
    # Get all genes for this species
    species_genes = genes_df[genes_df['species'] == species]
    
    # Find clusters containing these genes
    clusters_with_species = species_genes['cluster_id'].unique()
    
    # Find clusters that also contain other species
    multi_species_clusters = []
    for cluster_id in clusters_with_species:
        cluster_species = genes_df[genes_df['cluster_id'] == cluster_id]['species'].unique()
        if len(cluster_species) > 1:
            multi_species_clusters.append(cluster_id)
    
    # Calculate percentage of genes with orthologs
    genes_with_orthologs = species_genes[species_genes['cluster_id'].isin(multi_species_clusters)]
    percent_with_orthologs = len(genes_with_orthologs) / len(species_genes) * 100
    
    species_with_orthologs[species] = percent_with_orthologs
    log(f"  {species}: {percent_with_orthologs:.2f}% of genes have orthologs in other species")

# Step 6: Visualize the key findings
log("\nGenerating visualizations...")

# 1. Cluster classification distribution
plt.figure(figsize=(12, 6))
sorted_counts = classification_counts.sort_values(ascending=False)
sns.barplot(x=sorted_counts.index, y=sorted_counts.values)
plt.title('Cluster Classification Distribution')
plt.xlabel('Classification')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{output_dir}/cluster_classifications.png")
plt.close()

# 2. Species count distribution
species_counts = cluster_composition['species_count'].value_counts().sort_index()
plt.figure(figsize=(10, 6))
sns.barplot(x=species_counts.index.astype(str), y=species_counts.values)
plt.title('Number of Species per Cluster')
plt.xlabel('Number of Species')
plt.ylabel('Count of Clusters')
plt.tight_layout()
plt.savefig(f"{output_dir}/species_per_cluster.png")
plt.close()

# 3. Species distribution in genes
species_gene_counts = genes_df['species'].value_counts()
plt.figure(figsize=(10, 6))
sns.barplot(x=species_gene_counts.index, y=species_gene_counts.values)
plt.title('Gene Distribution by Species')
plt.xlabel('Species')
plt.ylabel('Number of Genes')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_dir}/species_gene_distribution.png")
plt.close()

# 4. Heatmap of species distribution in multi-copy orthologs
species_columns = [col for col in cluster_composition.columns 
                  if col in ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']]

multi_copy_clusters = cluster_composition[
    cluster_composition['classification'] == 'multi_copy_ortholog'
].sort_values(by=species_columns)

if len(multi_copy_clusters) > 0:
    # Take a subset for better visualization if there are many
    sample_size = min(50, len(multi_copy_clusters))
    sample_clusters = multi_copy_clusters.head(sample_size)
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        sample_clusters[species_columns], 
        cmap='viridis', 
        annot=True, 
        fmt='d'
    )
    plt.title(f'Species Distribution in Multi-Copy Ortholog Clusters (Sample of {sample_size})')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/multi_copy_heatmap.png")
    plt.close()

# 5. Ortholog conservation by species
if species_with_orthologs:
    plt.figure(figsize=(10, 6))
    species_list = list(species_with_orthologs.keys())
    values = [species_with_orthologs[sp] for sp in species_list]
    sns.barplot(x=species_list, y=values)
    plt.title('Percentage of Genes with Orthologs in Other Species')
    plt.xlabel('Species')
    plt.ylabel('Percentage')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/ortholog_conservation.png")
    plt.close()

# 6. Cluster size distribution
cluster_sizes = genes_df.groupby('cluster_id').size().value_counts().sort_index()
plt.figure(figsize=(12, 6))
plt.bar(cluster_sizes.index.astype(str), cluster_sizes.values)
plt.title('Cluster Size Distribution')
plt.xlabel('Number of Genes in Cluster')
plt.ylabel('Number of Clusters')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{output_dir}/cluster_size_distribution.png")
plt.close()

# Generate summary report
log("\nGenerating summary report...")
with open(f"{output_dir}/cluster_analysis_summary.md", "w") as f:
    f.write("# SYNERGY2 Cluster Analysis Summary\n\n")
    
    f.write("## Dataset Overview\n\n")
    f.write(f"- Total genes analyzed: {len(genes_df)}\n")
    f.write(f"- Total clusters: {genes_df['cluster_id'].nunique()}\n")
    f.write(f"- Average genes per cluster: {len(genes_df) / genes_df['cluster_id'].nunique():.2f}\n\n")
    
    f.write("## Species Distribution\n\n")
    f.write("| Species | Gene Count | Percentage |\n")
    f.write("|---------|------------|------------|\n")
    
    species_distribution = genes_df['species'].value_counts()
    total_genes = len(genes_df)
    
    for species, count in species_distribution.items():
        percentage = count / total_genes * 100
        f.write(f"| {species} | {count} | {percentage:.2f}% |\n")
    
    f.write("\n## Cluster Classifications\n\n")
    f.write("| Classification | Count | Percentage |\n")
    f.write("|----------------|-------|------------|\n")
    
    total_clusters = len(cluster_composition)
    sorted_classifications = classification_counts.sort_values(ascending=False)
    
    for classification, count in sorted_classifications.items():
        percentage = count / total_clusters * 100
        f.write(f"| {classification} | {count} | {percentage:.2f}% |\n")
    
    f.write("\n## Single-Copy Orthologs\n\n")
    single_copy_count = len(single_copy_pivot)
    f.write(f"Found {single_copy_count} clusters with exactly one gene from each species.\n")
    f.write(f"These represent {single_copy_count / total_clusters * 100:.2f}% of all clusters.\n\n")
    
    f.write("## Ortholog Conservation\n\n")
    f.write("| Species | Percentage with Orthologs |\n")
    f.write("|---------|-------------------------|\n")
    
    for species, percentage in species_with_orthologs.items():
        f.write(f"| {species} | {percentage:.2f}% |\n")
    
    f.write("\n## Gene Family Expansions\n\n")
    for species in ['C_albicans', 'C_glabrata', 'K_lactis', 'S_cerevisiae']:
        if species not in cluster_composition.columns:
            continue
            
        expanded_clusters = cluster_composition[cluster_composition[species] >= 3]
        if len(expanded_clusters) > 0:
            f.write(f"### {species}\n\n")
            f.write(f"- {len(expanded_clusters)} clusters with gene expansions\n")
            
            # List top 5 expansions
            top5 = expanded_clusters.sort_values(by=species, ascending=False).head(5)
            f.write("- Top 5 expansions:\n")
            
            for idx, row in top5.iterrows():
                f.write(f"  - {idx}: {row[species]} genes\n")
            
            f.write("\n")
    
    f.write("\n## Evolutionary Insights\n\n")
    
    f.write(f"- **Core genome**: {single_copy_count} single-copy orthologs ({single_copy_count / total_clusters * 100:.2f}%) represent the core conserved genome.\n")
    
    species_specific = sum(classification_counts.filter(regex='_specific'))
    f.write(f"- **Species-specific genes**: {species_specific} clusters ({species_specific / total_clusters * 100:.2f}%) contain genes unique to one species.\n")
    
    multi_copy = classification_counts.get('multi_copy_ortholog', 0)
    f.write(f"- **Gene family expansions**: {multi_copy} clusters ({multi_copy / total_clusters * 100:.2f}%) show evidence of gene duplication while maintaining representation in all species.\n")
    
    partial = classification_counts.get('partial_ortholog', 0)
    f.write(f"- **Partial orthologs**: {partial} clusters ({partial / total_clusters * 100:.2f}%) show gene loss in one or more species.\n")

log(f"\nAnalysis complete! Results saved to {output_dir}/")
log(f"Completed at: {pd.Timestamp.now()}")
log_file.close()

print(f"\nDone! Files saved in {output_dir}/")
print(f"Check {output_dir}/cluster_analysis_summary.md for a comprehensive summary of the analysis")
