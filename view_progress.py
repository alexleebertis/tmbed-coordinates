#!/usr/bin/env python3
from pathlib import Path
import argparse

def show_progress(micro_chunk_dir, results_dir):
    """Show easy-to-read progress without grep."""
    micro_chunk_dir = Path(micro_chunk_dir)
    results_dir = Path(results_dir)
    
    expected = list(micro_chunk_dir.glob("*.fasta"))
    done = list(results_dir.glob("*.tsv"))
    pending = [e for e in expected if not (results_dir / f"{e.stem}.tsv").exists()]
    
    print(f"\n📊 TMbed Progress Dashboard")
    print(f"{'='*40}")
    print(f"✅ Completed: {len(done)}/{len(expected)} ({100*len(done)//len(expected)}%)")
    print(f"⏳ Pending:   {len(pending)}")
    print(f"{'='*40}")
    
    if pending:
        print(f"\n🎯 Next 3 to process:")
        for p in pending[:3]:
            print(f"   - {p.name}")
    
    if len(done) > 0 and len(pending) == 0:
        print("\n🎉 ALL DONE! Ready to merge results.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--micro-chunk-dir", required=True)
    parser.add_argument("--results-dir", required=True)
    args = parser.parse_args()
    show_progress(args.micro_chunk_dir, args.results_dir)