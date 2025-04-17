import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# data directories (inputs)
ygob_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/Homology_data"
synergy_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/comparison_results"
output_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/comparison_results"

# sets the pairs to compare
species_pairs = [
    ("C_albicans", "S_cerevisiae"),
    ("C_glabrata", "K_lactis"),
    ("C_glabrata", "S_cerevisiae"),
    ("K_lactis", "S_cerevisiae"),
    ("C_albicans", "C_glabrata"),
    ("C_albicans", "K_lactis")
]

def standardize_gene_id(gene_id):
    gene_id = str(gene_id)
    
    # converts S. cerevisiae gene IDs from SYNERGY2 format to standard SGD IDs
    if gene_id.startswith('gene-Y'):
        # removes 'gene-' prefix from S. cerevisiae genes
        return gene_id.replace('gene-', '')
    
    # handles K. lactis genes - remove gene: prefix AND underscores
    if 'KLLA' in gene_id:
        #  removes "gene:" prefix if present
        gene_id = gene_id.replace('gene:', '')
        # removes any underscores
        gene_id = gene_id.replace('_', '')
        if "KLLA0Cluster" in gene_id:
            gene_id = gene_id.replace("KLLA0Cluster", "KLLA Cluster")
        return gene_id
    
    #Some genes won't need changes, like for c. albicans and c. glabrata
    return gene_id

# function to compare ortholog datasets
def compare_ortholog_sets(ygob_pairs, synergy_pairs, species1, species2):
    # uses standardized gene IDs
    ygob_pairs['gene1_std'] = ygob_pairs['gene1'].apply(standardize_gene_id)
    ygob_pairs['gene2_std'] = ygob_pairs['gene2'].apply(standardize_gene_id)
    synergy_pairs['gene1_std'] = synergy_pairs['gene1'].apply(standardize_gene_id)
    synergy_pairs['gene2_std'] = synergy_pairs['gene2'].apply(standardize_gene_id)
    
    # creates a gene pair tuple for comparison
    ygob_pairs['pair'] = ygob_pairs.apply(
        lambda row: tuple(sorted([row['gene1_std'], row['gene2_std']])), axis=1
    )
    synergy_pairs['pair'] = synergy_pairs.apply(
        lambda row: tuple(sorted([row['gene1_std'], row['gene2_std']])), axis=1
    )
    
    # convert to sets for comparison
    ygob_set = set(ygob_pairs['pair'])
    synergy_set = set(synergy_pairs['pair'])
    
    # calculates metrics
    true_positives = len(ygob_set.intersection(synergy_set))
    false_positives = len(synergy_set - ygob_set)
    false_negatives = len(ygob_set - synergy_set)
    
    # calculates evaluation metrics
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # creates results dictionary
    results = {
        'species_pair': f"{species1.replace('_', '.')} - {species2.replace('_', '.')}",
        'ygob_pairs': len(ygob_set),
        'synergy_pairs': len(synergy_set),
        'true_positives': true_positives,
        'false_positives': false_positives,
        'false_negatives': false_negatives,
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score
    }
    
    # saves the comparison details
    comparison_df = pd.DataFrame({
        'gene_pair': list(ygob_set.union(synergy_set)),
        'in_ygob': [pair in ygob_set for pair in ygob_set.union(synergy_set)],
        'in_synergy': [pair in synergy_set for pair in ygob_set.union(synergy_set)]
    })
    comparison_df.to_csv(os.path.join(output_dir, f"{species1}_{species2}_comparison_details.csv"), index=False)
    
    return results

# compares each species pair
comparison_results = []

for sp1, sp2 in species_pairs:
    print(f"Comparing {sp1} - {sp2}...")
    
    ygob_file = os.path.join(ygob_dir, f"{sp1}_{sp2}_pairs.csv")
    if not os.path.exists(ygob_file):
        print(f"  YGOB file not found: {ygob_file}")
        continue
        
    ygob_pairs = pd.read_csv(ygob_file)
    
    synergy_file = os.path.join(synergy_dir, f"synergy2_{sp1}_{sp2}_pairs.csv")
    if not os.path.exists(synergy_file):
        print(f"  SYNERGY2 file not found: {synergy_file}")
        continue
        
    synergy_pairs = pd.read_csv(synergy_file)
    
    # compares the sets
    results = compare_ortholog_sets(ygob_pairs, synergy_pairs, sp1, sp2)
    comparison_results.append(results)
    
    #results
    print(f"  YGOB pairs: {results['ygob_pairs']}")
    print(f"  SYNERGY2 pairs: {results['synergy_pairs']}")
    print(f"  True positives: {results['true_positives']}")
    print(f"  Precision: {results['precision']:.4f}")
    print(f"  Recall: {results['recall']:.4f}")
    print(f"  F1 Score: {results['f1_score']:.4f}")

# combines results and saves them as .csv
results_df = pd.DataFrame(comparison_results)
results_df.to_csv(os.path.join(output_dir, "comparison_summary.csv"), index=False)

# creates the graphs shared
if len(comparison_results) > 0:
    # bar chart of precision, recall, and F1 score by species pair
    plt.figure(figsize=(12, 6))
    
    # melt the dataframe for easier plotting
    plot_df = pd.melt(
        results_df, 
        id_vars=['species_pair'], 
        value_vars=['precision', 'recall', 'f1_score'],
        var_name='Metric', 
        value_name='Value'
    )
    
    # create the plot
    sns.barplot(x='species_pair', y='Value', hue='Metric', data=plot_df)
    plt.title('Comparison of SYNERGY2 vs YGOB/CGOB Orthology Predictions')
    plt.xlabel('Species Pair')
    plt.ylabel('Score')
    plt.ylim(0, 1)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'orthology_comparison_metrics.png'), dpi=300)
    
    # bar chart of counts by species pair
    plt.figure(figsize=(12, 6))
    
    # melt the dataframe for easier plotting (wide to long format)
    plot_df = pd.melt(
        results_df, 
        id_vars=['species_pair'], 
        value_vars=['true_positives', 'false_positives', 'false_negatives'],
        var_name='Category', 
        value_name='Count'
    )
    
    # creatse the plot
    sns.barplot(x='species_pair', y='Count', hue='Category', data=plot_df)
    plt.title('Orthology Prediction Comparison Counts')
    plt.xlabel('Species Pair')
    plt.ylabel('Number of Gene Pairs')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'orthology_comparison_counts.png'), dpi=300)
    
    print(f"Visualization saved to {output_dir}")

# creates a summary report
with open(os.path.join(output_dir, "orthology_comparison_report.md"), "w") as f:
    f.write("# SYNERGY2 vs YGOB/CGOB Orthology Comparison\n\n")
    f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d')}\n\n")
    
    f.write("## Summary\n\n")
    f.write("This report compares orthology predictions between SYNERGY2 and YGOB/CGOB databases.\n\n")
    
    f.write("## Comparison Results\n\n")
    f.write("| Species Pair | YGOB Pairs | SYNERGY2 Pairs | True Positives | False Positives | False Negatives | Precision | Recall | F1 Score |\n")
    f.write("|-------------|------------|---------------|---------------|----------------|-----------------|-----------|--------|----------|\n")
    
    for _, row in results_df.iterrows():
        f.write(f"| {row['species_pair']} | {row['ygob_pairs']} | {row['synergy_pairs']} | ")
        f.write(f"{row['true_positives']} | {row['false_positives']} | {row['false_negatives']} | ")
        f.write(f"{row['precision']:.4f} | {row['recall']:.4f} | {row['f1_score']:.4f} |\n")
    
    f.write("\n## Interpretation\n\n")
    f.write("### Precision\n")
    f.write("Measures the fraction of SYNERGY2 ortholog predictions that are also found in YGOB/CGOB.\n")
    f.write("Higher precision indicates fewer false positive ortholog predictions.\n\n")
    
    f.write("### Recall\n")
    f.write("Measures the fraction of YGOB/CGOB orthologs that were correctly identified by SYNERGY2.\n")
    f.write("Higher recall indicates fewer false negative ortholog predictions.\n\n")
    
    f.write("### F1 Score\n")
    f.write("The harmonic mean of precision and recall, providing a single metric for overall performance.\n\n")
    
    f.write("## Conclusions\n\n")
    f.write("Determining whether SYNERGY2 results are reliable depends on the specific research question and threshold of accuracy required.\n")
    f.write("Generally, an F1 score above 0.7 is considered good, and above 0.8 is excellent for orthology prediction.\n\n")
    
    f.write("### Visualizations\n\n")
    f.write("Please refer to the generated PNG files for visual representation of these metrics.\n")

print(f"Report generated at {os.path.join(output_dir, 'orthology_comparison_report.md')}")
