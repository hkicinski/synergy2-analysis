import pandas as pd
import os

# Pillars files
cgob_pillars_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/Homology_data/pillar-homology/E009-hk-04142025_Pillars_CGOB.tab"
ygob_pillars_file = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/Homology_data/pillar-homology/E009-hk-04142025_Pillars_YGOB.tab"
output_dir = "/space/hkicinski/Documents/E009-hk-Homology_SYNERGY/analysis/04102025-analysis-run1/Homology_data"

# reads the pillars files
cgob_df = pd.read_csv(cgob_pillars_file, sep='\t', na_values=["---", ""], keep_default_na=True)
ygob_df = pd.read_csv(ygob_pillars_file, sep='\t', na_values=["---", ""], keep_default_na=True)

# CGOB indices 
cgob_ca_col_idx = 0    
cgob_sc_col_idx = 17   

# YGOB indices
ygob_sc_col_name = 'YAL068C'   
ygob_cg_col_name = '---.7'    
ygob_kl_col_name = '---.14'   

print(f"Using CGOB columns: C. albicans={cgob_df.columns[cgob_ca_col_idx]}, S. cerevisiae={cgob_df.columns[cgob_sc_col_idx]}")
print(f"Using YGOB columns: S. cerevisiae={ygob_sc_col_name}, C. glabrata={ygob_cg_col_name}, K. lactis={ygob_kl_col_name}")

# extracts direct ortholog pairs from CGOB
c_albicans_s_cerevisiae_pairs = []
for idx, row in cgob_df.iterrows():
    # checks if both cells have values
    if pd.notna(row.iloc[cgob_ca_col_idx]) and pd.notna(row.iloc[cgob_sc_col_idx]):
        ca_gene = row.iloc[cgob_ca_col_idx]
        sc_gene = row.iloc[cgob_sc_col_idx]
        
        # creates species pairs A -> B
        c_albicans_s_cerevisiae_pairs.append({
            'species1': 'C. albicans',
            'gene1': ca_gene,
            'species2': 'S. cerevisiae',
            'gene2': sc_gene,
            'source': 'CGOB',
            'pillar_id': idx,
            'gene_pair': f"('{ca_gene}', '{sc_gene}')"
        })
        
        # creates pair B -> A
        c_albicans_s_cerevisiae_pairs.append({
            'species1': 'S. cerevisiae',
            'gene1': sc_gene, 
            'species2': 'C. albicans',
            'gene2': ca_gene,
            'source': 'CGOB',
            'pillar_id': idx,
            'gene_pair': f"('{sc_gene}', '{ca_gene}')"
        })

# extracts direct ortholog pairs from YGOB
c_glabrata_s_cerevisiae_pairs = []
k_lactis_s_cerevisiae_pairs = []
c_glabrata_k_lactis_pairs = []

for idx, row in ygob_df.iterrows():
    # C. glabrata - S. cerevisiae pairs
    if pd.notna(row[ygob_cg_col_name]) and pd.notna(row[ygob_sc_col_name]):
        cg_gene = row[ygob_cg_col_name]
        sc_gene = row[ygob_sc_col_name]
        
        c_glabrata_s_cerevisiae_pairs.append({
            'species1': 'C. glabrata',
            'gene1': cg_gene,
            'species2': 'S. cerevisiae',
            'gene2': sc_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{cg_gene}', '{sc_gene}')"
        })
        
        c_glabrata_s_cerevisiae_pairs.append({
            'species1': 'S. cerevisiae',
            'gene1': sc_gene,
            'species2': 'C. glabrata',
            'gene2': cg_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{sc_gene}', '{cg_gene}')"
        })
    
    # K. lactis - S. cerevisiae pairs
    if pd.notna(row[ygob_kl_col_name]) and pd.notna(row[ygob_sc_col_name]):
        kl_gene = row[ygob_kl_col_name]
        sc_gene = row[ygob_sc_col_name]
        
        k_lactis_s_cerevisiae_pairs.append({
            'species1': 'K. lactis',
            'gene1': kl_gene,
            'species2': 'S. cerevisiae',
            'gene2': sc_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{kl_gene}', '{sc_gene}')"
        })
        
        k_lactis_s_cerevisiae_pairs.append({
            'species1': 'S. cerevisiae',
            'gene1': sc_gene,
            'species2': 'K. lactis',
            'gene2': kl_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{sc_gene}', '{kl_gene}')"
        })
    
    # C. glabrata - K. lactis pairs
    if pd.notna(row[ygob_cg_col_name]) and pd.notna(row[ygob_kl_col_name]):
        cg_gene = row[ygob_cg_col_name]
        kl_gene = row[ygob_kl_col_name]
        
        c_glabrata_k_lactis_pairs.append({
            'species1': 'C. glabrata',
            'gene1': cg_gene,
            'species2': 'K. lactis',
            'gene2': kl_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{cg_gene}', '{kl_gene}')"
        })
        
        c_glabrata_k_lactis_pairs.append({
            'species1': 'K. lactis',
            'gene1': kl_gene,
            'species2': 'C. glabrata',
            'gene2': cg_gene,
            'source': 'YGOB',
            'pillar_id': idx,
            'gene_pair': f"('{kl_gene}', '{cg_gene}')"
        })

# creates mappings for bridge species (S. cerevisiae)
sc_to_ca = {}  # S. cerevisiae to C. albicans
sc_to_cg = {}  # S. cerevisiae to C. glabrata
sc_to_kl = {}  # S. cerevisiae to K. lactis

# from CGOB: S. cerevisiae to C. albicans
for pair in c_albicans_s_cerevisiae_pairs:
    if pair['species1'] == 'S. cerevisiae':
        sc_to_ca[pair['gene1']] = pair['gene2']

# from YGOB: S. cerevisiae to C. glabrata and K. lactis
for pair in c_glabrata_s_cerevisiae_pairs:
    if pair['species1'] == 'S. cerevisiae':
        sc_to_cg[pair['gene1']] = pair['gene2']

for pair in k_lactis_s_cerevisiae_pairs:
    if pair['species1'] == 'S. cerevisiae':
        sc_to_kl[pair['gene1']] = pair['gene2']

# 4. this is the bridge of C. albicans to C. glabrata and K. lactis through S. cerevisiae
c_albicans_c_glabrata_pairs = []
c_albicans_k_lactis_pairs = []

# C. albicans - C. glabrata bridge
print(f"Number of S. cerevisiae genes that map to both C. albicans and C. glabrata: {len(set(sc_to_ca.keys()) & set(sc_to_cg.keys()))}")
for sc_gene in set(sc_to_ca.keys()) & set(sc_to_cg.keys()):
    ca_gene = sc_to_ca[sc_gene]
    cg_gene = sc_to_cg[sc_gene]
    
    c_albicans_c_glabrata_pairs.append({
        'species1': 'C. albicans',
        'gene1': ca_gene,
        'species2': 'C. glabrata',
        'gene2': cg_gene,
        'source': 'CGOB/YGOB_bridge',
        'pillar_id': -1,
        'gene_pair': f"('{ca_gene}', '{cg_gene}')"
    })
    
    c_albicans_c_glabrata_pairs.append({
        'species1': 'C. glabrata',
        'gene1': cg_gene,
        'species2': 'C. albicans',
        'gene2': ca_gene,
        'source': 'CGOB/YGOB_bridge',
        'pillar_id': -1,
        'gene_pair': f"('{cg_gene}', '{ca_gene}')"
    })

# C. albicans - K. lactis bridge
print(f"Number of S. cerevisiae genes that map to both C. albicans and K. lactis: {len(set(sc_to_ca.keys()) & set(sc_to_kl.keys()))}")
for sc_gene in set(sc_to_ca.keys()) & set(sc_to_kl.keys()):
    ca_gene = sc_to_ca[sc_gene]
    kl_gene = sc_to_kl[sc_gene]
    
    c_albicans_k_lactis_pairs.append({
        'species1': 'C. albicans',
        'gene1': ca_gene,
        'species2': 'K. lactis',
        'gene2': kl_gene,
        'source': 'CGOB/YGOB_bridge',
        'pillar_id': -1,
        'gene_pair': f"('{ca_gene}', '{kl_gene}')"
    })
    
    c_albicans_k_lactis_pairs.append({
        'species1': 'K. lactis',
        'gene1': kl_gene,
        'species2': 'C. albicans',
        'gene2': ca_gene,
        'source': 'CGOB/YGOB_bridge',
        'pillar_id': -1,
        'gene_pair': f"('{kl_gene}', '{ca_gene}')"
    })

# C. albicans - S. cerevisiae
if c_albicans_s_cerevisiae_pairs:
    df = pd.DataFrame(c_albicans_s_cerevisiae_pairs)
    df.to_csv(os.path.join(output_dir, "C_albicans_S_cerevisiae_pairs.csv"), index=False)
    print(f"Created C_albicans_S_cerevisiae_pairs.csv with {len(c_albicans_s_cerevisiae_pairs)//2} pairs")

# C. glabrata - S. cerevisiae
if c_glabrata_s_cerevisiae_pairs:
    df = pd.DataFrame(c_glabrata_s_cerevisiae_pairs)
    df.to_csv(os.path.join(output_dir, "C_glabrata_S_cerevisiae_pairs.csv"), index=False)
    print(f"Created C_glabrata_S_cerevisiae_pairs.csv with {len(c_glabrata_s_cerevisiae_pairs)//2} pairs")

# K. lactis - S. cerevisiae
if k_lactis_s_cerevisiae_pairs:
    df = pd.DataFrame(k_lactis_s_cerevisiae_pairs)
    df.to_csv(os.path.join(output_dir, "K_lactis_S_cerevisiae_pairs.csv"), index=False)
    print(f"Created K_lactis_S_cerevisiae_pairs.csv with {len(k_lactis_s_cerevisiae_pairs)//2} pairs")

# C. glabrata - K. lactis
if c_glabrata_k_lactis_pairs:
    df = pd.DataFrame(c_glabrata_k_lactis_pairs)
    df.to_csv(os.path.join(output_dir, "C_glabrata_K_lactis_pairs.csv"), index=False)
    print(f"Created C_glabrata_K_lactis_pairs.csv with {len(c_glabrata_k_lactis_pairs)//2} pairs")

# C. albicans - C. glabrata (bridged)
if c_albicans_c_glabrata_pairs:
    df = pd.DataFrame(c_albicans_c_glabrata_pairs)
    df.to_csv(os.path.join(output_dir, "C_albicans_C_glabrata_pairs.csv"), index=False)
    print(f"Created C_albicans_C_glabrata_pairs.csv with {len(c_albicans_c_glabrata_pairs)//2} pairs")

# C. albicans - K. lactis (bridged)
if c_albicans_k_lactis_pairs:
    df = pd.DataFrame(c_albicans_k_lactis_pairs)
    df.to_csv(os.path.join(output_dir, "C_albicans_K_lactis_pairs.csv"), index=False)
    print(f"Created C_albicans_K_lactis_pairs.csv with {len(c_albicans_k_lactis_pairs)//2} pairs")

# summary stats
print("\nSummary:")
print(f"S. cerevisiae - C. albicans mappings: {len(sc_to_ca)}")
print(f"S. cerevisiae - C. glabrata mappings: {len(sc_to_cg)}")
print(f"S. cerevisiae - K. lactis mappings: {len(sc_to_kl)}")
print(f"C. albicans - C. glabrata bridges: {len(c_albicans_c_glabrata_pairs)//2}")
print(f"C. albicans - K. lactis bridges: {len(c_albicans_k_lactis_pairs)//2}")
print(f"Total bridges created: {len(c_albicans_c_glabrata_pairs)//2 + len(c_albicans_k_lactis_pairs)//2}")
