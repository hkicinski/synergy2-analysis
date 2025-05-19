import pickle

# Load the mapping file
with open('/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/synergy_runs/yeast_analysis/nodes/root/locus_mappings.pkl', 'rb') as f:
    mappings = pickle.load(f)

# Print information about the structure of the mappings
print(f"Type of mappings: {type(mappings)}")
if isinstance(mappings, dict):
    print(f"Number of keys: {len(mappings)}")
    print("Sample of keys and values:")
    for k, v in list(mappings.items())[:3]:  # Show first 3 entries
        print(f"  {k}: {v}")

# Try direct lookup and reverse lookup for our genes
cgla_gene = "CAGL0M09405g"
scer_gene = "gene-YBL103C"  # Try with prefix first

# Direct lookup (if gene IDs are keys)
if cgla_gene in mappings:
    print(f"\nDirect mapping for {cgla_gene}: {mappings[cgla_gene]}")
else:
    print(f"\n{cgla_gene} not found as a key")

if scer_gene in mappings:
    print(f"Direct mapping for {scer_gene}: {mappings[scer_gene]}")
else:
    print(f"{scer_gene} not found as a key")

# Reverse lookup (if gene IDs are values)
cgla_found = False
scer_found = False

# Look through all mappings for our genes
for k, v in mappings.items():
    # Check different data structures depending on what's in the mappings
    if isinstance(v, list) or isinstance(v, tuple):
        for item in v:
            if cgla_gene == item:
                print(f"\nReverse mapping for {cgla_gene}: {k}")
                cgla_found = True
            if scer_gene == item:
                print(f"Reverse mapping for {scer_gene}: {k}")
                scer_found = True
    elif isinstance(v, str):
        if cgla_gene == v:
            print(f"\nReverse mapping for {cgla_gene}: {k}")
            cgla_found = True
        if scer_gene == v:
            print(f"Reverse mapping for {scer_gene}: {k}")
            scer_found = True

# Try with alternative format for S. cerevisiae gene ID
if not scer_found:
    scer_gene_alt = "YBL103C"  # Without prefix
    for k, v in mappings.items():
        if isinstance(v, list) or isinstance(v, tuple):
            for item in v:
                if scer_gene_alt == item:
                    print(f"Reverse mapping for {scer_gene_alt}: {k}")
                    scer_found = True
        elif isinstance(v, str):
            if scer_gene_alt == v:
                print(f"Reverse mapping for {scer_gene_alt}: {k}")
                scer_found = True

if not cgla_found:
    print(f"{cgla_gene} not found in values")
if not scer_found:
    print(f"Neither {scer_gene} nor {scer_gene_alt} found in values")