import pandas as pd
import os

# path to SYNERGY2 results
synergy_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis"
clusters_file = os.path.join(synergy_dir, "nodes/root/final_clusters.txt")
mapping_file = os.path.join(synergy_dir, "nodes/root/clust_to_trans.txt")
output_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/comparison_results"

# loads SYNERGY2 data and prints column info for debugging
try:
    clusters = pd.read_csv(clusters_file, sep="\t")
except:
    clusters = pd.read_csv(clusters_file, sep="\t", header=None)
    # assigns column names based on expected format
    if len(clusters.columns) >= 3:
        clusters.columns = ["cluster_id", "genes", "type"]
    else:
        clusters.columns = ["cluster_id", "genes"]

try:
    mapping = pd.read_csv(mapping_file, sep="\t")
except:
    mapping = pd.read_csv(mapping_file, sep="\t", header=None)
    #debugging info
    print("Columns in mapping file:", mapping.columns.tolist())
    print("First few rows of mapping file:")
    print(mapping.head())
    
    if len(mapping.columns) >= 2:
        mapping.columns = ["cluster_id", "gene_id"]
    else:
        print("ERROR: Unexpected mapping file format")
        exit(1)

# prints updated column names
print("Final mapping columns:", mapping.columns.tolist())
print("First few rows after renaming:")
print(mapping.head())

# function to get species from gene IDs based on check_gene_ids.py results
def get_species_from_gene_id(gene_id):
    gene_id = str(gene_id)
    if gene_id.startswith('orf19'):
        return "C. albicans"
    elif gene_id.startswith('CAGL0'):
        return "C. glabrata"
    elif 'KLLA0' in gene_id:
        return "K. lactis"
    elif gene_id.startswith('gene-Y') or gene_id.startswith('gene-S'):
        return "S. cerevisiae"
    else:
        # prints unidentified gene IDs for debugging
        print(f"Unidentified gene ID format: {gene_id}")
        return "Unknown"

# adds species information to mapping
gene_id_column = "gene_id"  # Default column name
if gene_id_column not in mapping.columns and len(mapping.columns) >= 2:
    gene_id_column = mapping.columns[1]
    print(f"Using '{gene_id_column}' as the gene ID column")

mapping['species'] = mapping[gene_id_column].apply(get_species_from_gene_id)

# counts species distribution for verification
species_counts = mapping['species'].value_counts()
print("\nSpecies distribution in dataset:")
print(species_counts)

# creates ortholog pairs from SYNERGY2 clusters
synergy_pairs = []

# sets cluster_id_column
cluster_id_column = "cluster_id"
if cluster_id_column not in mapping.columns and len(mapping.columns) >= 1:
    cluster_id_column = mapping.columns[0]
    print(f"Using '{cluster_id_column}' as the cluster ID column")

# groups genes by cluster
for cluster_id, group in mapping.groupby(cluster_id_column):
    genes = group[gene_id_column].tolist()
    species = group['species'].tolist()
    
    # for each pair of genes from different species
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            if species[i] != species[j] and species[i] != "Unknown" and species[j] != "Unknown":  # Only include identified inter-species pairs
                synergy_pairs.append({
                    'species1': species[i],
                    'gene1': genes[i],
                    'species2': species[j],
                    'gene2': genes[j],
                    'cluster_id': cluster_id
                })

# converts to DataFrame for saving
synergy_pairs_df = pd.DataFrame(synergy_pairs)
synergy_pairs_df.to_csv(os.path.join(output_dir, "synergy2_ortholog_pairs.csv"), index=False)
print(f"Created main ortholog pairs file with {len(synergy_pairs)} pairs")

# creates species-specific pairs files to be stored in output!
species_list = ["S. cerevisiae", "C. glabrata", "K. lactis", "C. albicans"]
for sp1 in species_list:
    for sp2 in species_list:
        if sp1 < sp2:  # To avoid duplicates
            sp_pairs = synergy_pairs_df[
                ((synergy_pairs_df['species1'] == sp1) & (synergy_pairs_df['species2'] == sp2)) |
                ((synergy_pairs_df['species1'] == sp2) & (synergy_pairs_df['species2'] == sp1))
            ]
            
            if len(sp_pairs) > 0:
                file_name = f"synergy2_{sp1.replace('. ', '_')}_{sp2.replace('. ', '_')}_pairs.csv"
                sp_pairs.to_csv(os.path.join(output_dir, file_name), index=False)
                print(f"Created SYNERGY2 {sp1}-{sp2} ortholog pairs file with {len(sp_pairs)} pairs")
