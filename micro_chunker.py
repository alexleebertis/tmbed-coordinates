#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import sys

def split_fasta_to_micro_chunks(input_fasta, output_dir, chunk_size=5):
    """Split a fasta file into micro-chunks of specified size."""
    input_path = Path(input_fasta).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Read all sequences using fasta-blast to handle comment lines
    sequences = list(SeqIO.parse(input_path, "fasta-blast"))
    total = len(sequences)
    
    if total == 0:
        print(f"❌ No sequences found in {input_path}")
        return
    
    # Create micro-chunks
    chunk_count = 0
    for i in range(0, total, chunk_size):
        chunk_seqs = sequences[i:i+chunk_size]
        chunk_name = f"{input_path.stem}_micro_{chunk_count:03d}.fasta"
        chunk_file = output_path / chunk_name
        
        SeqIO.write(chunk_seqs, chunk_file, "fasta")
        print(f"✅ Created {chunk_name} ({len(chunk_seqs)} proteins)")
        chunk_count += 1
    
    print(f"\n📊 Split {total} proteins into {chunk_count} micro-chunks of max {chunk_size} each")
    print(f"📁 Location: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTA into micro-chunks for WSL stability")
    parser.add_argument("--input", required=True, help="Input FASTA file (e.g., chunk_032b_0.fasta)")
    parser.add_argument("--output-dir", required=True, help="Output directory for micro-chunks")
    parser.add_argument("--chunk-size", type=int, default=5, help="Proteins per micro-chunk (default: 5)")
    
    args = parser.parse_args()
    split_fasta_to_micro_chunks(args.input, args.output_dir, args.chunk_size)