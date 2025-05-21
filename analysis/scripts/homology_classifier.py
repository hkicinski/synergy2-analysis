#!/usr/bin/env python3
"""
homology_classifier.py - Extract and classify homologous relationships from SYNERGY2 output

This script analyzes SYNERGY2 output files to classify genes into different homology categories:
- Pure orthologs: 1:1 orthologous relationships across species
- Ohnologs: Paralogs derived from whole genome duplication (WGD)
- Paralogs: Duplicated genes within a species (not ohnologs)

Author: Hubert Kicinski
Date: May 20, 2025
"""

import pandas as pd
import numpy as np
import os
import sys
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

# inputs and output directories
BASE_INPUT_DIR = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro"
CLUSTER_PARSED_DIR = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-parsed"
CLUSTER_ANALYSIS_DIR = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/cluster-analysis-results"
OUTPUT_DIR = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/05202025-synergy2-output-interro/homology_classification"

# species info
WGD_SPECIES = ['S_cerevisiae', 'C_glabrata']  # Species that underwent whole genome duplication
NON_WGD_SPECIES = ['C_albicans', 'K_lactis']   # Species that did not undergo WGD
ALL_SPECIES = WGD_SPECIES + NON_WGD_SPECIES

def load_synergy2_data():
    """Load required SYNERGY2 output files."""
    print("Loading SYNERGY2 data...")

    files = {
        'genes': os.path.join(CLUSTER_PARSED_DIR, 'genes_by_cluster.csv'),
        'clusters': os.path.join(CLUSTER_ANALYSIS_DIR, 'cluster_classifications.csv'),
        'compositions': os.path.join(CLUSTER_ANALYSIS_DIR, 'cluster_compositions.csv'),
        'single_copy': os.path.join(CLUSTER_ANALYSIS_DIR, 'single_copy_orthologs.csv')
    }

    #verifies file existence (debug)
    missing_files = []
    for file_type, file_path in files.items():
        if not os.path.exists(file_path):
            missing_files.append(file_type)
            print(f"WARNING: {file_type} file not found at {file_path}")

    # alt paths if file is not OG found
    if 'single_copy' in missing_files:
        alt_path = os.path.join(BASE_INPUT_DIR, 'single_copy_orthologs.csv')
        if os.path.exists(alt_path):
            print(f"Found single_copy at alternate location: {alt_path}")
            files['single_copy'] = alt_path
            missing_files.remove('single_copy')

    if 'genes' in missing_files:
        alt_path = os.path.join(BASE_INPUT_DIR, 'genes_by_cluster.csv')
        if os.path.exists(alt_path):
            print(f"Found genes at alternate location: {alt_path}")
            files['genes'] = alt_path
            missing_files.remove('genes')

    if missing_files:
        print("Critical files are still missing. Will continue with available files.")

    # Load files
    data = {}

    # Load genes
    if 'genes' not in missing_files:
        data['genes'] = pd.read_csv(files['genes'])
        print(f"Loaded genes data with {len(data['genes'])} entries")
    else:
        print("ERROR: Cannot continue without genes data")
        sys.exit(1)

    # Load clusters
    if 'clusters' not in missing_files:
        data['clusters'] = pd.read_csv(files['clusters'])
        print(f"Loaded cluster classifications with {len(data['clusters'])} entries")
    else:
        print("ERROR: Cannot continue without cluster classifications")
        sys.exit(1)

    # load compositions
    if 'compositions' not in missing_files:
        data['compositions'] = pd.read_csv(files['compositions'])
        print(f"Loaded cluster compositions with {len(data['compositions'])} entries")
    else:
        print("WARNING: Missing cluster compositions, using cluster classifications instead")
        data['compositions'] = data['clusters']

    # load single_copy if available
    if 'single_copy' not in missing_files:
        data['single_copy'] = pd.read_csv(files['single_copy'])
        print(f"Loaded single-copy orthologs with {len(data['single_copy'])} entries")
    else:
        print("WARNING: Single-copy orthologs file not found. Will derive from cluster classifications.")
        # We'll derive single-copy orthologs from cluster classifications later

    return data

def classify_pure_orthologs(data):
    """Identify pure 1:1 orthologous relationships (exactly one gene from each species)."""
    print("Classifying pure orthologous relationships...")

    # get single-copy ortholog clusters
    single_copy_clusters = data['clusters'][data['clusters']['classification'] == 'single_copy_ortholog']

    # get genes in these clusters
    pure_ortholog_genes = data['genes'][data['genes']['cluster_id'].isin(single_copy_clusters['cluster_id'])]

    print(f"Found {len(pure_ortholog_genes)} genes in {len(single_copy_clusters)} pure ortholog clusters")

    # added ortholog partners information
    orthologs_with_partners = []

    # processed each cluster to identify ortholog relationships
    for cluster_id, group in pure_ortholog_genes.groupby('cluster_id'):
        # get genes for each species
        species_genes = {}
        for _, row in group.iterrows():
            species_genes[row['species']] = row['original_id']

        # for each gene, list its orthologs in other species
        for species, gene_id in species_genes.items():
            orthologs = {sp: sp_gene for sp, sp_gene in species_genes.items() if sp != species}

            orthologs_with_partners.append({
                'cluster_id': cluster_id,
                'gene_id': gene_id,
                'species': species,
                'homology_class': 'pure_ortholog',
                'ortholog_partners': str(orthologs)
            })

    #  DataFrame conversion
    pure_orthologs = pd.DataFrame(orthologs_with_partners)

    return pure_orthologs

def identify_ohnologs(data):
    """Identify ohnologs (WGD-derived paralogs) from SYNERGY2 output."""
    print("Identifying ohnologs...")

   # This is where assumptions are made. I assumed a 2:1 ratio where 2 post-WGD genes exist for 1 non-wgd
   # This isn't a robust assumption as this discredits gene loss history.
   #Even in the YGOB paper, only 11% of the ohnologs post-wgd were retained in S. cerevisiae. Very quick and dirty.

    ohnolog_pattern_clusters = []

    for wgd_species in WGD_SPECIES:
        # finds clusters where this WGD species has exactly 2 genes
        wgd_candidates = data['clusters'][data['clusters'][wgd_species] == 2]

        # for these clusters, check if any non-WGD species has exactly 1 gene
        for non_wgd_species in NON_WGD_SPECIES:
            if non_wgd_species in data['clusters'].columns:
                pattern_matches = wgd_candidates[wgd_candidates[non_wgd_species] == 1]
                ohnolog_pattern_clusters.extend(pattern_matches['cluster_id'].tolist())

    # removes duplicates
    ohnolog_pattern_clusters = list(set(ohnolog_pattern_clusters))

    print(f"Found {len(ohnolog_pattern_clusters)} clusters with ohnolog pattern")

    # extracts genes from these clusters
    ohnolog_genes = []

    for cluster_id in ohnolog_pattern_clusters:
        cluster_genes = data['genes'][data['genes']['cluster_id'] == cluster_id]

        # Group genes by species
        species_genes = {}
        for _, row in cluster_genes.iterrows():
            if row['species'] not in species_genes:
                species_genes[row['species']] = []
            species_genes[row['species']].append(row['original_id'])

        # For each WGD species with exactly 2 genes (ohnolog pair)
        for wgd_species in WGD_SPECIES:
            if wgd_species in species_genes and len(species_genes[wgd_species]) == 2:
                gene1, gene2 = species_genes[wgd_species]

                # Get orthologs in non-WGD species (if any)
                orthologs = {}
                for non_wgd_species in NON_WGD_SPECIES:
                    if non_wgd_species in species_genes and len(species_genes[non_wgd_species]) == 1:
                        orthologs[non_wgd_species] = species_genes[non_wgd_species][0]

                # Add each gene as an ohnolog
                for gene_id in [gene1, gene2]:
                    ohnolog_genes.append({
                        'cluster_id': cluster_id,
                        'gene_id': gene_id,
                        'species': wgd_species,
                        'homology_class': 'ohnolog',
                        'ohnolog_partner': gene2 if gene_id == gene1 else gene1,
                        'ohnolog_evidence': 'WGD_2:1_pattern',
                        'outgroup_orthologs': str(orthologs)
                    })

    # Convert to DataFrame
    ohnologs = pd.DataFrame(ohnolog_genes)

    if len(ohnologs) > 0:
        print(f"Found {len(ohnologs)} genes in {len(ohnolog_pattern_clusters)} ohnolog pairs")
    else:
        print("No ohnologs found")
        ohnologs = pd.DataFrame(columns=['cluster_id', 'gene_id', 'species', 'homology_class',
                                         'ohnolog_partner', 'ohnolog_evidence', 'outgroup_orthologs'])

    return ohnologs

def identify_paralogs(data, pure_orthologs, ohnologs):
    """Identify paralogous relationships (intra-species duplications) that are not ohnologs."""
    print("Identifying non-ohnolog paralogs...")

    # Get clusters with paralogs (multiple genes from the same species)
    # that are not already classified as ohnologs

    # First, identify all clusters with paralogs
    paralog_containing_clusters = data['clusters'][
        (data['clusters']['C_albicans'] > 1) |
        (data['clusters']['C_glabrata'] > 1) |
        (data['clusters']['K_lactis'] > 1) |
        (data['clusters']['S_cerevisiae'] > 1)
    ]

    # Extract genes from these clusters
    all_paralog_genes = []

    for _, cluster in paralog_containing_clusters.iterrows():
        cluster_id = cluster['cluster_id']

        # Get genes in this cluster
        cluster_genes = data['genes'][data['genes']['cluster_id'] == cluster_id]

        # Count genes by species
        species_counts = cluster_genes['species'].value_counts().to_dict()

        # Process species with multiple genes (paralogs)
        for species, count in species_counts.items():
            if count > 1:
                # Get all genes for this species
                species_genes = cluster_genes[cluster_genes['species'] == species]['original_id'].tolist()

                # Create entries for each paralogous gene
                for i, gene_id in enumerate(species_genes):
                    # Generate list of paralog partners (all other genes of this species)
                    paralog_partners = [g for g in species_genes if g != gene_id]

                    # Generate list of ortholog partners (one gene per other species)
                    ortholog_partners = {}
                    for other_species in [s for s in ALL_SPECIES if s != species]:
                        other_species_genes = cluster_genes[cluster_genes['species'] == other_species]['original_id'].tolist()
                        if other_species_genes:
                            # Just take the first gene as the ortholog representative
                            ortholog_partners[other_species] = other_species_genes[0]

                    # Determine duplication type
                    duplication_type = "small_scale_duplication"

                    all_paralog_genes.append({
                        'cluster_id': cluster_id,
                        'gene_id': gene_id,
                        'species': species,
                        'homology_class': 'paralog',
                        'paralog_partners': str(paralog_partners),
                        'ortholog_partners': str(ortholog_partners),
                        'duplication_type': duplication_type
                    })

    # Convert to DataFrame
    all_paralogs = pd.DataFrame(all_paralog_genes)

    # Filter out genes already classified as ohnologs
    if len(ohnologs) > 0 and len(all_paralogs) > 0:
        ohnolog_genes = set(ohnologs['gene_id'])
        paralogs = all_paralogs[~all_paralogs['gene_id'].isin(ohnolog_genes)]
    else:
        paralogs = all_paralogs

    if len(paralogs) > 0:
        print(f"Found {len(paralogs)} genes in paralogous relationships (non-ohnolog)")
    else:
        print("No non-ohnolog paralogs found")
        paralogs = pd.DataFrame(columns=['cluster_id', 'gene_id', 'species', 'homology_class',
                                         'paralog_partners', 'ortholog_partners', 'duplication_type'])

    return paralogs

def refine_paralogs(paralogs):
    """Refine paralog classification to identify tandem, segmental, and dispersed duplications."""
    print("Refining paralog classifications...")

    # Without genome annotation data, we can't directly determine duplication types
    # For now, just set all to "unclassified"
    if len(paralogs) > 0:
        paralogs['duplication_subtype'] = "unclassified"

    return paralogs

def summarize_homology_classifications(pure_orthologs, ohnologs, paralogs):
    """Create summary statistics and visualizations."""
    print("Generating summary statistics and visualizations...")

    # Combine all classifications
    homology_classes = []

    if len(pure_orthologs) > 0:
        homology_classes.append(pure_orthologs[['cluster_id', 'gene_id', 'species', 'homology_class']])

    if len(ohnologs) > 0:
        homology_classes.append(ohnologs[['cluster_id', 'gene_id', 'species', 'homology_class']])

    if len(paralogs) > 0:
        homology_classes.append(paralogs[['cluster_id', 'gene_id', 'species', 'homology_class']])

    if not homology_classes:
        print("Error: No homology classifications found")
        return None

    all_classified = pd.concat(homology_classes, ignore_index=True)

    # Create summary DataFrame
    summary_data = all_classified['homology_class'].value_counts().reset_index()
    summary_data.columns = ['category', 'count']

    # Summary by species
    species_summary = all_classified.groupby(['species', 'homology_class']).size().reset_index()
    species_summary.columns = ['species', 'category', 'count']

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Save classification data
    pure_orthologs.to_csv(os.path.join(OUTPUT_DIR, 'pure_orthologs.csv'), index=False)
    if len(ohnologs) > 0:
        ohnologs.to_csv(os.path.join(OUTPUT_DIR, 'ohnologs.csv'), index=False)
    if len(paralogs) > 0:
        paralogs.to_csv(os.path.join(OUTPUT_DIR, 'paralogs.csv'), index=False)

    # Save combined classification
    all_classified.to_csv(os.path.join(OUTPUT_DIR, 'all_homology_classifications.csv'), index=False)

    # Save summary statistics
    summary_data.to_csv(os.path.join(OUTPUT_DIR, 'homology_summary.csv'), index=False)
    species_summary.to_csv(os.path.join(OUTPUT_DIR, 'homology_by_species.csv'), index=False)

    # Generate visualizations
    plt.figure(figsize=(10, 6))
    sns.barplot(data=summary_data, x='category', y='count')
    plt.title('Gene Counts by Homology Classification')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'homology_classification.png'))
    plt.close()

    # Generate species-specific visualization
    plt.figure(figsize=(12, 8))
    sns.barplot(data=species_summary, x='species', y='count', hue='category')
    plt.title('Homology Classifications by Species')
    plt.xticks(rotation=45)
    plt.legend(title='Homology Type')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'homology_by_species.png'))
    plt.close()

    # Generate markdown summary
    with open(os.path.join(OUTPUT_DIR, 'homology_classification_summary.md'), 'w') as f:
        f.write("# SYNERGY2 Homology Classification Summary\n\n")

        f.write("## Overall Classification\n\n")
        f.write("| Category | Gene Count | Percentage |\n")
        f.write("|----------|------------|------------|\n")

        total_genes = len(all_classified)
        for _, row in summary_data.iterrows():
            percentage = (row['count'] / total_genes) * 100
            f.write(f"| {row['category']} | {row['count']} | {percentage:.2f}% |\n")

        f.write("\n## Classification by Species\n\n")
        f.write("| Species | Pure Orthologs | Ohnologs | Paralogs | Total |\n")
        f.write("|---------|---------------|----------|----------|-------|\n")

        # Create pivot table for species summary
        pivot_data = species_summary.pivot_table(
            values='count',
            index='species',
            columns='category',
            fill_value=0
        ).reset_index()

        # Ensure all columns exist
        for col in ['pure_ortholog', 'ohnolog', 'paralog']:
            if col not in pivot_data.columns:
                pivot_data[col] = 0

        for _, row in pivot_data.iterrows():
            species = row['species']
            pure_count = row.get('pure_ortholog', 0)
            ohnolog_count = row.get('ohnolog', 0)
            paralog_count = row.get('paralog', 0)
            total = pure_count + ohnolog_count + paralog_count

            f.write(f"| {species} | {pure_count} | {ohnolog_count} | {paralog_count} | {total} |\n")

        if len(ohnologs) > 0 and 'ohnolog_evidence' in ohnologs.columns:
            f.write("\n## Ohnolog Evidence\n\n")
            f.write("| Evidence Type | Count | Percentage |\n")
            f.write("|---------------|-------|------------|\n")

            for evidence, count in ohnologs['ohnolog_evidence'].value_counts().items():
                percentage = (count / len(ohnologs)) * 100
                f.write(f"| {evidence} | {count} | {percentage:.2f}% |\n")

            # Show some examples of ohnolog pairs
            f.write("\n### Sample Ohnolog Pairs\n\n")
            f.write("| WGD Species | Ohnolog Pair | Non-WGD Orthologs |\n")
            f.write("|-------------|-------------|-------------------|\n")

            # Group by cluster and species to show pairs
            sample_pairs = []
            for (cluster_id, species), group in ohnologs.groupby(['cluster_id', 'species']):
                if len(group) == 2:  # Should always be 2 for ohnologs
                    sample_pairs.append({
                        'cluster_id': cluster_id,
                        'species': species,
                        'gene1': group['gene_id'].iloc[0],
                        'gene2': group['gene_id'].iloc[1],
                        'orthologs': group['outgroup_orthologs'].iloc[0]
                    })

            # Display up to 10 examples
            for i, pair in enumerate(sample_pairs[:10]):
                f.write(f"| {pair['species']} | {pair['gene1']}, {pair['gene2']} | {pair['orthologs']} |\n")

        # Get total numbers for conclusion
        total_pure_orthologs = summary_data.loc[summary_data['category'] == 'pure_ortholog', 'count'].iloc[0] if 'pure_ortholog' in summary_data['category'].values else 0
        total_ohnologs = summary_data.loc[summary_data['category'] == 'ohnolog', 'count'].iloc[0] if 'ohnolog' in summary_data['category'].values else 0
        total_paralogs = summary_data.loc[summary_data['category'] == 'paralog', 'count'].iloc[0] if 'paralog' in summary_data['category'].values else 0

        f.write("\n## Conclusion\n\n")
        f.write(f"This analysis identified three distinct classes of homology relationships in the SYNERGY2 output:\n\n")
        f.write(f"1. **Pure Orthologs**: {total_pure_orthologs} genes ({total_pure_orthologs/total_genes*100:.2f}%) with one-to-one orthologous relationships across species\n")
        f.write(f"2. **Ohnologs (WGD-derived paralogs)**: {total_ohnologs} genes ({total_ohnologs/total_genes*100:.2f}%) derived from whole genome duplication\n")
        f.write(f"3. **Other Paralogs**: {total_paralogs} genes ({total_paralogs/total_genes*100:.2f}%) resulting from small-scale duplication events\n\n")
        f.write(f"These classifications provide insights into the evolutionary history of genes across the analyzed yeast species, highlighting both conservation (orthologs) and diversification (paralogs and ohnologs).\n")

    print(f"Summary files saved to {OUTPUT_DIR}")
    return all_classified

def main():
    print("SYNERGY2 Homology Classification")
    print("================================")
    print(f"Base input directory: {BASE_INPUT_DIR}")
    print(f"Cluster parsed directory: {CLUSTER_PARSED_DIR}")
    print(f"Cluster analysis directory: {CLUSTER_ANALYSIS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"WGD species: {', '.join(WGD_SPECIES)}")
    print(f"Non-WGD species: {', '.join(NON_WGD_SPECIES)}")
    print("")

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load SYNERGY2 data
    data = load_synergy2_data()

    # 1. Identify pure orthologous relationships
    pure_orthologs = classify_pure_orthologs(data)

    # 2. Identify ohnologs (WGD-derived paralogs)
    ohnologs = identify_ohnologs(data)

    # 3. Identify other paralogs
    paralogs = identify_paralogs(data, pure_orthologs, ohnologs)

    # 4. Refine paralog classification
    paralogs = refine_paralogs(paralogs)

    # 5. Summarize classifications
    all_classified = summarize_homology_classifications(pure_orthologs, ohnologs, paralogs)

    print("\nHomology classification complete!")
    print(f"Results saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
