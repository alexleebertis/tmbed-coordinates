#!/usr/bin/env python3
import subprocess
import time
from pathlib import Path
import argparse
import sys

def run_with_recovery(micro_chunk_dir, tmbed_script_path, results_dir):
    """Run TMbed on micro-chunks sequentially with auto-recovery."""
    micro_chunk_dir = Path(micro_chunk_dir).resolve()
    results_dir = Path(results_dir).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Get all micro-chunk fasta files
    fasta_files = sorted(micro_chunk_dir.glob("*_micro_*.fasta"))
    total = len(fasta_files)
    
    if total == 0:
        print(f"❌ No micro-chunk files found in {micro_chunk_dir}")
        return
    
    print(f"🔬 Found {total} micro-chunks to process")
    print(f"💾 Results will save to: {results_dir}")
    print(f"⚡ Processing sequentially (WSL-safe mode)\n")
    
    completed = 0
    failed = 0
    skipped = 0
    
    for i, fasta_file in enumerate(fasta_files, 1):
        chunk_name = fasta_file.stem
        result_file = results_dir / f"{chunk_name}.tsv"
        
        # Count proteins in this chunk (count '>' characters minus 1 for empty split)
        protein_count = len(fasta_file.read_text().split('>')) - 1
        
        # Skip if already done (crash recovery)
        if result_file.exists():
            print(f"[{i}/{total}] ⏭️  SKIP: {chunk_name} ({protein_count} proteins, already done)")
            skipped += 1
            continue
        
        print(f"[{i}/{total}] 🚀 START: {chunk_name} ({protein_count} proteins)")
        
        try:
            # Run TMbed with 5-minute timeout per micro-chunk
            result = subprocess.run(
                ["python", tmbed_script_path, "--fasta", str(fasta_file), "-o", str(result_file), "--batch-size", "1"],
                capture_output=True,
                text=True,
                timeout=300  # 5 minutes max per micro-chunk
            )
            
            if result.returncode == 0:
                print(f"[{i}/{total}] ✅ DONE: {chunk_name}")
                completed += 1
            else:
                print(f"[{i}/{total}] ❌ FAIL: {chunk_name} (exit code {result.returncode})")
                failed += 1
                # Save error log
                (results_dir / f"{chunk_name}_ERROR.txt").write_text(result.stderr)
                
        except subprocess.TimeoutExpired:
            print(f"[{i}/{total}] ⏰ TIMEOUT: {chunk_name} (took >5 minutes)")
            failed += 1
        except Exception as e:
            print(f"[{i}/{total}] 💥 CRASH: {chunk_name} ({str(e)})")
            failed += 1
        
        # WSL breathing room
        time.sleep(2)
        
        # Progress summary every 5 chunks
        if i % 5 == 0:
            print(f"\n📈 Progress: {completed} done, {failed} failed, {skipped} skipped, {total - i} remaining\n")
    
    # Final report
    print(f"\n{'='*50}")
    print(f"🏁 FINAL REPORT:")
    print(f"   ✅ Completed: {completed}")
    print(f"   ❌ Failed:    {failed}")
    print(f"   ⏭️  Skipped:   {skipped}")
    print(f"   📊 Total:     {total}")
    print(f"{'='*50}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WSL-safe sequential TMbed runner")
    parser.add_argument("--micro-chunk-dir", required=True, help="Directory containing micro-chunks")
    parser.add_argument("--tmbed-script", required=True, help="Path to tmbed_coords.py")
    parser.add_argument("--results-dir", required=True, help="Where to save .tsv results")
    
    args = parser.parse_args()
    run_with_recovery(args.micro_chunk_dir, args.tmbed_script, args.results_dir)