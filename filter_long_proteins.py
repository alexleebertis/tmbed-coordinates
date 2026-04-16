#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from pathlib import Path

def filter_by_length(input_fasta, output_normal, output_giants, max_length=2000):
    """Separate normal proteins from giants."""
    input_path = Path(input_fasta)
    normal_path = Path(output_normal)
    giant_path = Path(output_giants)
    
    sequences = list(SeqIO.parse(input_path, "fasta"))
    
    normal = [s for s in sequences if len(s.seq) <= max_length]
    giants = [s for s in sequences if len(s.seq) > max_length]
    
    SeqIO.write(normal, normal_path, "fasta")
    SeqIO.write(giants, giant_path, "fasta")
    
    print(f"📊 Total: {len(sequences)} proteins")
    print(f"✅ Normal (≤{max_length} aa): {len(normal)} -> {normal_path}")
    print(f"🦕 Giants (>{max_length} aa): {len(giants)} -> {giant_path}")
    for g in giants:
        print(f"   - {g.id}: {len(g.seq)} aa")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-normal", required=True)
    parser.add_argument("--output-giants", required=True)
    parser.add_argument("--max-length", type=int, default=2000)
    args = parser.parse_args()
    filter_by_length(args.input, args.output_normal, args.output_giants, args.max_length)