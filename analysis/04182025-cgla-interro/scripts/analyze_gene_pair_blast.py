import os

# Paths
jft_blast_fa = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis/nodes/root/JFT.blast.fa"
jft_blast_headers = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis/nodes/root/JFT.blast_headers.txt"
blast_headers = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis/nodes/root/blast_headers.txt"

print("Searching blast headers for our genes...")

# Function to search a file for patterns
def search_file(file_path, patterns):
    results = {pattern: [] for pattern in patterns}
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                for pattern in patterns:
                    if pattern in line:
                        results[pattern].append((line_num, line.strip()))
        return results
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return results

# Search for our genes in the headers files
patterns = ["CAGL0M09405g", "YBL103C"]
jft_headers_results = search_file(jft_blast_headers, patterns)
blast_headers_results = search_file(blast_headers, patterns)

# Print results
print("Results from JFT blast headers:")
for pattern, matches in jft_headers_results.items():
    print(f"  {pattern}: {len(matches)} matches")
    for line_num, line in matches[:3]:  # Show up to 3 matches
        print(f"    Line {line_num}: {line}")

print("\nResults from blast headers:")
for pattern, matches in blast_headers_results.items():
    print(f"  {pattern}: {len(matches)} matches")
    for line_num, line in matches[:3]:  # Show up to 3 matches
        print(f"    Line {line_num}: {line}")

# If we found our genes in the headers, look for their sequences in the FASTA file
if any(len(matches) > 0 for matches in jft_headers_results.values()):
    print("\nSearching JFT.blast.fa for sequences...")
    # We would extract the sequences here based on the headers found
    # This is complex to do in a script without knowing the exact format
    print("Please check the headers file to identify the correct format for extraction")
