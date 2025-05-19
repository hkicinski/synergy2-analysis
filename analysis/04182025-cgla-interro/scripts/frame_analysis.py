#!/usr/bin/env python3
"""
frame_analysis.py - Compare translations in all reading frames
"""

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import gffutils
import tempfile

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analyze reading frame issues')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA')
    parser.add_argument('--gff', required=True, help='Path to GFF3 annotations')
    parser.add_argument('--root_pep', required=True, help='Path to root.pep file')
    parser.add_argument('--cgd_proteins', required=True, help='Path to CGD proteins')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--identity_range', default='40-60', help='Identity range to focus on (e.g., 40-60)')
    return parser.parse_args()

def extract_cds_sequence(gene_id, db, genome_dict):
    """Extract coding sequence for a gene"""
    try:
        # Try to find the gene feature
        gene_features = list(db.features_of_type('gene', id=gene_id))
        if not gene_features:
            return None, f"Gene {gene_id} not found"
        
        gene = gene_features[0]  # Take the first match if multiple exist
        
        # Get CDS features
        cds_features = list(db.children(gene, featuretype='CDS', order_by='start'))
        
        if not cds_features:
            return None, f"No CDS features found for {gene_id}"
        
        strand = cds_features[0].strand
        chrom = cds_features[0].seqid
        
        # Check if chromosome exists in genome
        if chrom not in genome_dict:
            return None, f"Chromosome {chrom} not found"
        
        # Sort CDS features
        if strand == '+':
            cds_features.sort(key=lambda x: x.start)
        else:
            cds_features.sort(key=lambda x: x.start, reverse=True)
        
        # Extract and combine sequences
        cds_seq = ""
        for cds in cds_features:
            seq = str(genome_dict[chrom].seq[cds.start-1:cds.end])
            if strand == '-':
                seq = str(Seq(seq).reverse_complement())
            
            cds_seq += seq
        
        return cds_seq, "Success"
    
    except Exception as e:
        return None, f"Error: {str(e)}"

def translate_in_all_frames(sequence):
    """Translate sequence in all three forward frames"""
    seq_obj = Seq(sequence)
    translations = []
    
    for i in range(3):
        frame_seq = seq_obj[i:]
        # Adjust length to be divisible by 3
        adj_length = len(frame_seq) - (len(frame_seq) % 3)
        frame_seq = frame_seq[:adj_length]
        
        # Translate
        protein = frame_seq.translate()
        translations.append(str(protein))
    
    return translations

def find_protein_in_fasta(gene_id, fasta_file):
    """Find protein sequence in FASTA file by gene ID"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if gene_id in record.id:
            return str(record.seq)
    return None

def calculate_identity(seq1, seq2):
    """Calculate simple percent identity between two sequences"""
    # Trim terminal stop codons
    if seq1 and seq1.endswith('*'):
        seq1 = seq1[:-1]
    if seq2 and seq2.endswith('*'):
        seq2 = seq2[:-1]
    
    if not seq1 or not seq2:
        return 0.0
    
    # Count matches
    matches = sum(1 for a, b in zip(seq1[:min(len(seq1), len(seq2))], 
                                   seq2[:min(len(seq1), len(seq2))]) if a == b)
    
    # Calculate identity
    identity = (matches / min(len(seq1), len(seq2))) * 100
    
    return identity

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse identity range
    min_id, max_id = map(float, args.identity_range.split('-'))
    
    # Create temporary GFF database with a strategy to handle duplicates
    db_file = os.path.join(args.output_dir, "gff.db")
    
    try:
        # Use a merge strategy that won't fail on duplicates
        db = gffutils.create_db(
            args.gff, 
            db_file, 
            force=True,
            merge_strategy='create_unique',
            id_spec={'gene': 'ID'},
            verbose=True
        )
    except Exception as e:
        print(f"Error creating GFF database: {str(e)}")
        print("Trying alternative approach...")
        
        # Try a different approach - create a filtered GFF file
        temp_gff = os.path.join(args.output_dir, "filtered.gff")
        seen_ids = set()
        with open(args.gff, 'r') as f_in, open(temp_gff, 'w') as f_out:
            for line in f_in:
                if line.startswith('#') or line.strip() == '':
                    f_out.write(line)
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 9:
                    f_out.write(line)
                    continue
                    
                # Extract ID from 9th column
                attrs = parts[8]
                if 'ID=' not in attrs:
                    f_out.write(line)
                    continue
                
                # Get the ID
                id_part = attrs.split('ID=')[1].split(';')[0]
                
                # Only write lines with unique IDs
                if id_part not in seen_ids:
                    seen_ids.add(id_part)
                    f_out.write(line)
        
        # Create database from filtered GFF
        db = gffutils.create_db(
            temp_gff, 
            db_file, 
            force=True,
            merge_strategy='error',
            id_spec={'gene': 'ID'}
        )
    
    # Load genome
    print("Loading genome...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    
    # Get C. glabrata genes
    print("Identifying C. glabrata genes...")
    sample_genes = []
    for gene in db.features_of_type('gene'):
        if gene.id.startswith('CAGL'):
            sample_genes.append(gene.id)
            if len(sample_genes) >= 20:  # Analyze 20 genes
                break
    
    if not sample_genes:
        print("Error: No C. glabrata genes (CAGL*) found in the GFF file")
        return
    
    print(f"Found {len(sample_genes)} C. glabrata genes to analyze")
    
    results = []
    for gene_id in sample_genes:
        print(f"Processing {gene_id}...")
        
        # Extract CDS sequence
        cds_seq, cds_status = extract_cds_sequence(gene_id, db, genome_dict)
        if not cds_seq:
            print(f"  Error: {cds_status}")
            continue
        
        # Translate in all frames
        frame_translations = translate_in_all_frames(cds_seq)
        
        # Get reference protein sequences
        print(f"  Looking for {gene_id} in reference files...")
        synergy_protein = find_protein_in_fasta(gene_id, args.root_pep)
        cgd_protein = find_protein_in_fasta(gene_id, args.cgd_proteins)
        
        if not synergy_protein:
            print(f"  Warning: {gene_id} not found in SYNERGY2 root.pep")
            continue
            
        if not cgd_protein:
            print(f"  Warning: {gene_id} not found in CGD proteins")
            continue
        
        # Calculate identity to reference for each frame
        frame_identities_synergy = [calculate_identity(trans, synergy_protein) 
                                  for trans in frame_translations]
        
        frame_identities_cgd = [calculate_identity(trans, cgd_protein) 
                              for trans in frame_translations]
        
        # Determine best matching frames
        best_frame_synergy = frame_identities_synergy.index(max(frame_identities_synergy)) + 1
        best_frame_cgd = frame_identities_cgd.index(max(frame_identities_cgd)) + 1
        
        # Save results
        result = {
            'gene_id': gene_id,
            'cds_length': len(cds_seq),
            'frame1_synergy_identity': frame_identities_synergy[0],
            'frame2_synergy_identity': frame_identities_synergy[1],
            'frame3_synergy_identity': frame_identities_synergy[2],
            'frame1_cgd_identity': frame_identities_cgd[0],
            'frame2_cgd_identity': frame_identities_cgd[1],
            'frame3_cgd_identity': frame_identities_cgd[2],
            'best_frame_synergy': best_frame_synergy,
            'best_frame_cgd': best_frame_cgd,
            'synergy_cgd_identity': calculate_identity(synergy_protein, cgd_protein)
        }
        
        results.append(result)
        
        # Save the sequence files for inspection
        with open(os.path.join(args.output_dir, f"{gene_id}_sequences.txt"), 'w') as f:
            f.write(f"CDS: {cds_seq[:100]}...(truncated, total length: {len(cds_seq)})\n\n")
            for i, trans in enumerate(frame_translations):
                f.write(f"Frame {i+1}: {trans[:100]}...(truncated, total length: {len(trans)})\n\n")
            if synergy_protein:
                f.write(f"SYNERGY2: {synergy_protein[:100]}...(truncated, total length: {len(synergy_protein)})\n\n")
            if cgd_protein:
                f.write(f"CGD: {cgd_protein[:100]}...(truncated, total length: {len(cgd_protein)})\n\n")
    
    # Compile results
    if results:
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(args.output_dir, "frame_analysis_results.csv"), index=False)
        
        # Analyze frame shifts
        frame_match = sum(df['best_frame_synergy'] == df['best_frame_cgd'])
        frame_mismatch = len(df) - frame_match
        
        print("\nFrame Analysis Results:")
        print(f"Total genes analyzed: {len(df)}")
        print(f"Genes with matching best frames: {frame_match} ({frame_match/len(df)*100:.1f}%)")
        print(f"Genes with different best frames: {frame_mismatch} ({frame_mismatch/len(df)*100:.1f}%)")
        
        if frame_mismatch > 0:
            print("\nFrame shift pattern:")
            for synergy_frame in range(1, 4):
                for cgd_frame in range(1, 4):
                    count = sum((df['best_frame_synergy'] == synergy_frame) & 
                              (df['best_frame_cgd'] == cgd_frame))
                    if count > 0:
                        print(f"  SYNERGY frame {synergy_frame} vs CGD frame {cgd_frame}: {count} genes")
    else:
        print("No results to analyze. Check if the gene IDs in GFF match those in protein files.")
    
    # Clean up
    if os.path.exists(db_file):
        os.remove(db_file)

if __name__ == "__main__":
    main()
