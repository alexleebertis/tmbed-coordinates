#!/usr/bin/env python3
"""
Split large FASTA files into smaller chunks for batch processing.

Usage:
    python split_fasta.py -i large.fasta -o chunks/ -n 150

Example:
    python split_fasta.py -i proteins.fasta -o data/chunks/ -n 150
"""

import argparse
from pathlib import Path
from Bio import SeqIO
import math


def split_fasta(input_file, output_dir, chunk_size=150):
    """Split FASTA file into chunks."""
    input_path = Path(input_file)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Read all sequences
    print(f"Reading {input_path}...")
    sequences = list(SeqIO.parse(input_path, "fasta"))
    total = len(sequences)

    if total == 0:
        print("Error: No sequences found in input file")
        return

    print(f"  Total sequences: {total}")
    print(f"  Chunk size: {chunk_size}")

    # Calculate number of chunks
    num_chunks = math.ceil(total / chunk_size)
    print(f"  Number of chunks: {num_chunks}")
    print()

    # Write chunks
    created_files = []
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, total)
        chunk = sequences[start_idx:end_idx]

        chunk_file = output_path / f"chunk_{i:03d}.fasta"
        SeqIO.write(chunk, chunk_file, "fasta")
        created_files.append(chunk_file)
        print(f"  chunk_{i:03d}.fasta: {len(chunk)} sequences")

    print()
    print(f"✓ Created {len(created_files)} chunk files in {output_path}")
    print()
    print("Next steps:")
    print(f"  1. Process chunks: ./batch_process.sh {output_path} results/")
    print(f"  2. Aggregate: python aggregate_results.py results/")


def main():
    parser = argparse.ArgumentParser(
        description='Split FASTA file into smaller chunks',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python split_fasta.py -i proteins.fasta -o chunks/
  python split_fasta.py -i large.fasta -o data/chunks/ -n 100
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for chunks')
    parser.add_argument('-n', '--chunk-size', type=int, default=150,
                        help='Number of sequences per chunk (default: 150)')

    args = parser.parse_args()

    split_fasta(args.input, args.output, args.chunk_size)


if __name__ == "__main__":
    main()
