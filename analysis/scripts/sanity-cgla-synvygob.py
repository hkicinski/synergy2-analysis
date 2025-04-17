import pandas as pd
import os

# Path to your SYNERGY2 output files
synergy_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis"
cluster_file = os.path.join(synergy_dir, "nodes/root/final_clusters.txt")
mapping_file = os.path.join(synergy_dir, "nodes/root/clust_to_trans.txt")

# Sample gene pairs to investigate
gene_pairs = [
    ('CAGL0C01463g', 'YPL180W'),
    ('CAGL0H08910g', 'YGL138C'),
    ('CAGL0H03091g', 'YGL091C'),
    ('CAGL0H07667g', 'YGL248W'),
    ('CAGL0D01672g', 'YPR055W')
]

# Load SYNERGY2 mapping data
try:
    mapping = pd.read_csv(mapping_file, sep="\t", header=None)
    # Assume first column is cluster_id, second is gene_id
    mapping.columns = ["cluster_id", "gene_id"]
    print(f"Loaded mapping file with {len(mapping)} entries")
except Exception as e:
    print(f"Error loading mapping file: {e}")

# For each gene pair, check if genes exist and what clusters they're in
print("\nAnalyzing gene pairs missed by SYNERGY2:")
for cglabrata_gene, cerevisiae_gene in gene_pairs:
    print(f"\nPair: {cglabrata_gene} - {cerevisiae_gene}")
    
    # Check if genes exist in SYNERGY2 data
    cg_exists = cglabrata_gene in mapping['gene_id'].values
    sc_exists = cerevisiae_gene in mapping['gene_id'].values
    sc_with_prefix = f"gene-{cerevisiae_gene}" in mapping['gene_id'].values
    
    print(f"C. glabrata gene exists in SYNERGY2: {cg_exists}")
    print(f"S. cerevisiae gene exists in SYNERGY2 (without prefix): {sc_exists}")
    print(f"S. cerevisiae gene exists in SYNERGY2 (with 'gene-' prefix): {sc_with_prefix}")
    
    # If genes exist, check what clusters they're in
    if cg_exists:
        cg_cluster = mapping[mapping['gene_id'] == cglabrata_gene]['cluster_id'].values[0]
        print(f"C. glabrata gene is in cluster: {cg_cluster}")
        
        # Find other genes in this cluster
        cluster_genes = mapping[mapping['cluster_id'] == cg_cluster]['gene_id'].tolist()
        other_genes = [g for g in cluster_genes if g != cglabrata_gene]
        if other_genes:
            print(f"Other genes in this cluster: {', '.join(other_genes[:5])}")
            if len(other_genes) > 5:
                print(f"...and {len(other_genes)-5} more")
    
    # Check S. cerevisiae gene cluster
    sc_gene_to_check = f"gene-{cerevisiae_gene}" if sc_with_prefix else cerevisiae_gene
    if sc_exists or sc_with_prefix:
        sc_cluster = mapping[mapping['gene_id'] == sc_gene_to_check]['cluster_id'].values[0]
        print(f"S. cerevisiae gene is in cluster: {sc_cluster}")
        
        # Find other genes in this cluster
        cluster_genes = mapping[mapping['cluster_id'] == sc_cluster]['gene_id'].tolist()
        other_genes = [g for g in cluster_genes if g != sc_gene_to_check]
        if other_genes:
            print(f"Other genes in this cluster: {', '.join(other_genes[:5])}")
            if len(other_genes) > 5:
                print(f"...and {len(other_genes)-5} more")
