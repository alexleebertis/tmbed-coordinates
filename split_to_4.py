#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import sys

def split_to_four(input_fasta, output_dir):
    """Split a FASTA file into 4 equal chunks (032b_0 through _3)."""
    input_path = Path(input_fasta).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Read with fasta-blast to handle comments
    sequences = list(SeqIO.parse(input_path, "fasta-blast"))
    total = len(sequences)
    
    if total == 0:
        print(f"❌ No sequences found in {input_path}")
        return
    
    # Calculate chunk sizes (roughly equal)
    chunk_size = total // 4
    remainder = total % 4
    
    # Create 4 chunks
    chunks = []
    start = 0
    for i in range(4):
        # Add 1 extra to early chunks if there's remainder
        size = chunk_size + (1 if i < remainder else 0)
        end = start + size
        chunks.append(sequences[start:end])
        start = end
    
    # Write chunks
    base_name = input_path.stem  # chunk_032b
    for i, chunk_seqs in enumerate(chunks):
        chunk_name = f"{base_name}_{i}.fasta"
        chunk_file = output_path / chunk_name
        SeqIO.write(chunk_seqs, chunk_file, "fasta")
        print(f"✅ Created {chunk_name} ({len(chunk_seqs)} proteins)")
    
    print(f"\n📊 Split {total} proteins into 4 chunks")
    print(f"📁 Location: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTA into 4 equal parts")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    
    args = parser.parse_args()
    split_to_four(args.input, args.output_dir)