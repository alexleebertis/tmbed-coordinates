#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO
import sys

def process_chunk(input_fasta, results_dir, max_length=2000, batch_size=5):
    """Filter giants and run TMbed in one go."""
    input_path = Path(input_fasta).resolve()
    chunk_name = input_path.stem  # e.g., chunk_033
    output_dir = Path(results_dir) / chunk_name
    
    print(f"\n{'='*60}")
    print(f"🔬 Processing {chunk_name}")
    print(f"{'='*60}")
    
    # Read and filter
    sequences = list(SeqIO.parse(input_path, "fasta-blast"))
    total = len(sequences)
    
    normal = [s for s in sequences if len(s.seq) <= max_length]
    giants = [s for s in sequences if len(s.seq) > max_length]
    
    print(f"📊 Total proteins: {total}")
    print(f"✅ Normal (≤{max_length} aa): {len(normal)}")
    
    if giants:
        print(f"🦕 Giants (>{max_length} aa): {len(giants)} - SKIPPED:")
        for g in giants:
            print(f"   - {g.id}: {len(g.seq)} aa")
    
    if len(normal) == 0:
        print(f"⚠️  No processable proteins in {chunk_name}!")
        # Create empty marker
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "ALL_PROTEINS_WERE_GIANTS.txt").write_text(
            f"All {len(giants)} proteins exceeded {max_length} aa\n"
        )
        return
    
    # Create filtered temp file
    temp_filtered = output_dir / f"{chunk_name}_filtered.fasta"
    output_dir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(normal, temp_filtered, "fasta")
    
    # Run TMbed
    print(f"\n🚀 Running TMbed on {len(normal)} proteins...")
    cmd = [
        "python", "tmbed_coords.py",
        "--fasta", str(temp_filtered),
        "-o", str(output_dir / chunk_name),
        "--batch-size", str(batch_size)
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"✅ {chunk_name} complete!")
        
        # Clean up temp file after success (optional - comment out if you want to keep)
        temp_filtered.unlink()
        print(f"🧹 Cleaned up temp file")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ TMbed failed: {e}")
        print(f"   Stderr: {e.stderr}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter giants and run TMbed cleanly")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--max-length", type=int, default=2000, help="Max protein length")
    parser.add_argument("--batch-size", type=int, default=5, help="TMbed batch size")
    
    args = parser.parse_args()
    process_chunk(args.input, args.results_dir, args.max_length, args.batch_size)