import os

# Paths
cluster_file = "/space/hkicinski/Documents/E009-hk-Homology-SYNERGY_git/output/04102025-yeast-analysis1/clust_to_trans.txt"
blast_results = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis/nodes/root/blast.m8"

# Check the blast results for these genes
print("Searching blast results for our genes...")
cgla_blast_hits = []
scer_blast_hits = []

with open(blast_results, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 12:  # Standard BLAST output has 12 columns
            query = parts[0]
            subject = parts[1]
            percent_identity = float(parts[2])
            alignment_length = int(parts[3])
            e_value = float(parts[10])
            bit_score = float(parts[11])
            
            if "CAGL0M09405g" in query or "CAGL0M09405g" in subject:
                cgla_blast_hits.append({
                    'query': query,
                    'subject': subject,
                    'percent_identity': percent_identity,
                    'alignment_length': alignment_length,
                    'e_value': e_value,
                    'bit_score': bit_score
                })
            
            if "YBL103C" in query or "YBL103C" in subject:
                scer_blast_hits.append({
                    'query': query,
                    'subject': subject,
                    'percent_identity': percent_identity,
                    'alignment_length': alignment_length,
                    'e_value': e_value,
                    'bit_score': bit_score
                })

# Print the results
print(f"\nFound {len(cgla_blast_hits)} BLAST hits for CAGL0M09405g")
for hit in cgla_blast_hits[:5]:  # Show first 5 hits
    print(f"  {hit['query']} - {hit['subject']}: {hit['percent_identity']}% identity, E-value: {hit['e_value']}")

print(f"\nFound {len(scer_blast_hits)} BLAST hits for YBL103C")
for hit in scer_blast_hits[:5]:  # Show first 5 hits
    print(f"  {hit['query']} - {hit['subject']}: {hit['percent_identity']}% identity, E-value: {hit['e_value']}")

# Most importantly: check if these two genes hit each other
print("\nChecking if these genes hit each other in BLAST results...")
for hit in cgla_blast_hits:
    if "YBL103C" in hit['query'] or "YBL103C" in hit['subject']:
        print(f"  Match found! {hit['query']} - {hit['subject']}: {hit['percent_identity']}% identity, E-value: {hit['e_value']}")

# If no direct hit, check for next best hits
print("\nTop hits for each gene:")
if cgla_blast_hits:
    # Sort by bit score (higher is better)
    cgla_blast_hits.sort(key=lambda x: x['bit_score'], reverse=True)
    print("  Top hits for CAGL0M09405g:")
    for hit in cgla_blast_hits[:3]:  # Top 3 hits
        print(f"    {hit['query']} - {hit['subject']}: {hit['percent_identity']}% identity, E-value: {hit['e_value']}")

if scer_blast_hits:
    # Sort by bit score (higher is better)
    scer_blast_hits.sort(key=lambda x: x['bit_score'], reverse=True)
    print("  Top hits for YBL103C:")
    for hit in scer_blast_hits[:3]:  # Top 3 hits
        print(f"    {hit['query']} - {hit['subject']}: {hit['percent_identity']}% identity, E-value: {hit['e_value']}")
