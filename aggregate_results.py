#!/usr/bin/env python3
"""
Aggregate TMbed results from multiple chunk directories into combined TSV files.

Usage:
    python aggregate_results.py <results_dir>
    
Example:
    python aggregate_results.py results/
"""

import pandas as pd
import sys
from pathlib import Path


def aggregate_results(results_dir):
    results_path = Path(results_dir)
    
    if not results_path.exists():
        print(f"Error: Directory '{results_dir}' not found")
        sys.exit(1)
    
    summaries = []
    regions = []
    
    print(f"Scanning {results_dir}...")
    
    # Find all chunk directories
    chunk_dirs = sorted(results_path.glob("chunk_*"))
    if not chunk_dirs:
        chunk_dirs = sorted([d for d in results_path.iterdir() if d.is_dir()])
    
    print(f"Found {len(chunk_dirs)} subdirectories")
    
    for chunk_dir in chunk_dirs:
        summary_file = chunk_dir / "tmbed_summary.tsv"
        regions_file = chunk_dir / "tmbed_tm_regions.tsv"
        
        if summary_file.exists():
            summaries.append(pd.read_csv(summary_file, sep="\t"))
        if regions_file.exists():
            regions.append(pd.read_csv(regions_file, sep="\t"))
    
    # Combine and save
    if summaries:
        combined_summary = pd.concat(summaries, ignore_index=True)
        output_summary = results_path / "all_proteins_summary.tsv"
        combined_summary.to_csv(output_summary, sep="\t", index=False)
        
        print(f"\n✓ Protein summary: {len(combined_summary)} proteins")
        print(f"  Saved to: {output_summary}")
        print(f"\nStatistics:")
        print(f"  - With TM domains: {combined_summary['has_tm'].sum()}")
        print(f"  - Multi-pass (≥2 TMs): {combined_summary['is_multipass'].sum()}")
        print(f"  - Total TM domains: {combined_summary['tm_count'].sum()}")
    
    if regions:
        combined_regions = pd.concat(regions, ignore_index=True)
        output_regions = results_path / "all_tm_regions.tsv"
        combined_regions.to_csv(output_regions, sep="\t", index=False)
        
        print(f"\n✓ TM regions: {len(combined_regions)} regions")
        print(f"  Saved to: {output_regions}")
        
        # Type breakdown
        print(f"\nTM type distribution:")
        type_counts = combined_regions["type"].value_counts()
        for tm_type, count in type_counts.items():
            print(f"  - {tm_type}: {count}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python aggregate_results.py <results_dir>")
        sys.exit(1)
    
    aggregate_results(sys.argv[1])