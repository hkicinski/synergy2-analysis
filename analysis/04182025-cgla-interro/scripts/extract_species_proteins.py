import os
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract species-specific proteins')
    parser.add_argument('--input', required=True, help='Input protein FASTA file (root.pep)')
    parser.add_argument('--output_dir', required=True, help='Output directory for species-specific files')
    return parser.parse_args()

def extract_species_proteins(input_file, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize dictionaries for each species
    species_records = {
        'cgla': [],  # C. glabrata
        'scer': [],  # S. cerevisiae
        'klac': [],  # K. lactis
        'calb': []   # C. albicans
    }
    
    # Read input file and categorize by species
    for record in SeqIO.parse(input_file, "fasta"):
        gene_id = record.id
        
        # Identify species based on gene ID patterns
        if 'CAGL0' in gene_id:
            species_records['cgla'].append(record)
        elif gene_id.startswith('gene-Y') or (gene_id.startswith('Y') and len(gene_id) == 7):
            species_records['scer'].append(record)
        elif 'KLLA' in gene_id:
            species_records['klac'].append(record)
        elif 'orf19' in gene_id:
            species_records['calb'].append(record)
    
    # Write species-specific files
    for species, records in species_records.items():
        if records:
            output_file = os.path.join(output_dir, f"{species}_proteins.fasta")
            SeqIO.write(records, output_file, "fasta")
            print(f"Wrote {len(records)} {species} proteins to {output_file}")

def main():
    args = parse_arguments()
    extract_species_proteins(args.input, args.output_dir)

if __name__ == "__main__":
    main()
