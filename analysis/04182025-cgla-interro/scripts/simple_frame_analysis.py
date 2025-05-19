#!/usr/bin/env python3
"""
simple_frame_analysis.py - Compare translations in all reading frames
"""

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import gffutils

def parse_arguments():
    parser = argparse.ArgumentParser(description='Analyze reading frame issues')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA')
    parser.add_argument('--gff', required=True, help='Path to GFF3 annotations')
    parser.add_argument('--root_pep', required=True, help='Path to root.pep file')
    parser.add_argument('--cgd_proteins', required=True, help='Path to CGD proteins')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--num_genes', type=int, default=20, help='Number of genes to analyze')
    return parser.parse_args()

def extract_cds_sequence(gene_id, db, genome_dict):
    """Extract coding sequence for a gene"""
    try:
        # Find the gene feature
        gene = db[gene_id]
        
        # Get CDS features
        cds_features = list(db.children(gene, featuretype='CDS', order_by='start'))
        
        if not cds_features:
            return None, "No CDS features found"
        
        strand = cds_features[0].strand
        chrom = cds_features[0].seqid
        
        # Sort CDS features
        if strand == '+':
            cds_features.sort(key=lambda x: x.start)
        else:
            cds_features.sort(key=lambda x: x.start, reverse=True)
        
        # Extract and combine sequences
        cds_seq = ""
        for cds in cds_features:
            if chrom not in genome_dict:
                return None, "Chromosome not found"
            
            seq = str(genome_dict[chrom].seq[cds.start-1:cds.end])
            if strand == '-':
                seq = str(Seq(seq).reverse_complement())
            
            cds_seq += seq
        
        return cds_seq, "Success"
    
    except Exception as e:
        return None, str(e)

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
    
    # Create GFF database with a strategy to handle duplicates
    db_file = os.path.join(args.output_dir, "gff.db")
    db = gffutils.create_db(
        args.gff, 
        db_file, 
        force=True,
        merge_strategy='create_unique',
        verbose=True
    )
    
    # Load genome
    print("Loading genome...")
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    
    # Select specific C. glabrata genes to analyze
    sample_genes = [
        'CAGL0A00187g', 'CAGL0A00429g', 'CAGL0A01089g', 
        'CAGL0A01320g', 'CAGL0A01870g', 'CAGL0A02387g',
        'CAGL0A02904g', 'CAGL0A03344g', 'CAGL0A03608g',
        'CAGL0A04378g', 'CAGL0B00132g', 'CAGL0B01089g',
        'CAGL0B02563g', 'CAGL0B03157g', 'CAGL0B03641g',
        'CAGL0B04125g', 'CAGL0C00275g', 'CAGL0C01111g',
        'CAGL0C01705g', 'CAGL0C02959g'
    ]
    
    print(f"Analyzing {len(sample_genes)} C. glabrata genes")
    
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
        seq_file = os.path.join(args.output_dir, f"{gene_id}_sequences.txt")
        with open(seq_file, 'w') as f:
            f.write(f"CDS: {cds_seq[:100]}...\n\n")
            for i, trans in enumerate(frame_translations):
                f.write(f"Frame {i+1}: {trans[:100]}...\n\n")
            if synergy_protein:
                f.write(f"SYNERGY2: {synergy_protein[:100]}...\n\n")
            if cgd_protein:
                f.write(f"CGD: {cgd_protein[:100]}...\n\n")
    
    # Compile results
    if results:
        df = pd.DataFrame(results)
        csv_file = os.path.join(args.output_dir, "frame_analysis_results.csv")
        df.to_csv(csv_file, index=False)
        
        # Analyze frame shifts
        frame_match = sum(df['best_frame_synergy'] == df['best_frame_cgd'])
        frame_mismatch = len(df) - frame_match
        
        print("\nFrame Analysis Results:")
        print(f"Total genes analyzed: {len(df)}")
        print(f"Genes with matching best frames: {frame_match} ({frame_match/len(df)*100:.1f}%)")
        print(f"Genes with different best frames: {frame_mismatch} ({frame_mismatch/len(df)*100:.1f}%)")
        
        if frame_mismatch > 0:
            print("\nFrame shift pattern:")
            for frame in range(1, 4):
                count = sum((df['best_frame_synergy'] != df['best_frame_cgd']) & 
                           (df['best_frame_cgd'] == frame))
                print(f"  CGD frame {frame} translated in different frame: {count}")
            
            # Count specific frame transitions
            print("\nDetailed frame shifts:")
            for synergy_frame in range(1, 4):
                for cgd_frame in range(1, 4):
                    count = sum((df['best_frame_synergy'] == synergy_frame) & 
                              (df['best_frame_cgd'] == cgd_frame))
                    if count > 0:
                        print(f"  SYNERGY frame {synergy_frame} vs CGD frame {cgd_frame}: {count} genes")
    else:
        print("No results to analyze.")
    
    # Clean up
    if os.path.exists(db_file):
        os.remove(db_file)

if __name__ == "__main__":
    main()
