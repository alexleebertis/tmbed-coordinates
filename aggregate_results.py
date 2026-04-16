#!/usr/bin/env python3
"""
Aggregate TMbed results from multiple chunk directories OR flat files into combined TSV files.

Usage:
    # Auto-detect (subdirs or flat)
    python aggregate_results.py <results_dir>
    
    # Force flat mode (files directly in dir)
    python aggregate_results.py <results_dir> --flat
    
    # Force recursive mode (subdirectories)
    python aggregate_results.py <results_dir> --recursive
    
    # Custom output names
    python aggregate_results.py <results_dir> --output-prefix rep2_final

Examples:
    python aggregate_results.py results/mem_global_rep2
    python aggregate_results.py results/mem_global_rep2 --flat --output-prefix rep2
"""

import pandas as pd
import sys
import argparse
from pathlib import Path


def find_files_flat(directory, pattern):
    """Find files matching pattern directly in directory."""
    return sorted(directory.glob(pattern))


def find_files_recursive(directory, pattern):
    """Find files matching pattern in all subdirectories."""
    return sorted(directory.rglob(pattern))


def aggregate_results(results_dir, mode="auto", output_prefix=None, summary_name="tmbed_summary.tsv", regions_name="tmbed_tm_regions.tsv"):
    results_path = Path(results_dir).resolve()
    
    if not results_path.exists():
        print(f"❌ Error: Directory '{results_dir}' not found")
        sys.exit(1)
    
    # Determine file discovery mode
    if mode == "auto":
        # Check if there are subdirectories with the target files
        subdirs = [d for d in results_path.iterdir() if d.is_dir()]
        has_subdir_files = any((d / summary_name).exists() for d in subdirs[:5])  # Check first 5
        
        if has_subdir_files:
            mode = "recursive"
            print(f"📁 Auto-detected: Recursive mode (subdirectories)")
        else:
            mode = "flat"
            print(f"📄 Auto-detected: Flat mode (files in directory)")
    
    # Set up file finder
    if mode == "recursive":
        find_files = lambda p: find_files_recursive(results_path, f"*/{p}")
    else:
        find_files = lambda p: find_files_flat(results_path, p)
    
    print(f"\n🔍 Scanning {results_path}...")
    
    # Find files
    summary_files = find_files(summary_name)
    regions_files = find_files(regions_name)
    
    print(f"   Found {len(summary_files)} summary files")
    print(f"   Found {len(regions_files)} regions files")
    
    if not summary_files and not regions_files:
        print("❌ No TMbed result files found")
        print(f"   (Looked for {summary_name} and {regions_name})")
        sys.exit(1)
    
    # Process summaries
    if summary_files:
        summaries = []
        for f in summary_files:
            try:
                df = pd.read_csv(f, sep="\t")
                summaries.append(df)
                print(f"   ✓ {f.name}: {len(df)} proteins")
            except Exception as e:
                print(f"   ⚠ Error reading {f}: {e}")
        
        if summaries:
            combined_summary = pd.concat(summaries, ignore_index=True)
            
            # Generate output name
            if output_prefix:
                output_summary = results_path / f"{output_prefix}_summary.tsv"
            else:
                output_summary = results_path / "all_proteins_summary.tsv"
            
            combined_summary.to_csv(output_summary, sep="\t", index=False)
            
            print(f"\n📊 Protein Summary:")
            print(f"   Total proteins: {len(combined_summary)}")
            print(f"   With TM domains: {combined_summary['has_tm'].sum()}")
            print(f"   Multi-pass (≥2 TMs): {combined_summary['is_multipass'].sum()}")
            print(f"   Total TM domains: {combined_summary['tm_count'].sum()}")
            print(f"   Saved: {output_summary}")
    
    # Process regions
    if regions_files:
        regions = []
        for f in regions_files:
            try:
                df = pd.read_csv(f, sep="\t")
                regions.append(df)
            except Exception as e:
                print(f"   ⚠ Error reading {f}: {e}")
        
        if regions:
            combined_regions = pd.concat(regions, ignore_index=True)
            
            if output_prefix:
                output_regions = results_path / f"{output_prefix}_tm_regions.tsv"
            else:
                output_regions = results_path / "all_tm_regions.tsv"
            
            combined_regions.to_csv(output_regions, sep="\t", index=False)
            
            print(f"\n🧬 TM Regions:")
            print(f"   Total regions: {len(combined_regions)}")
            
            if 'type' in combined_regions.columns:
                type_counts = combined_regions["type"].value_counts()
                for tm_type, count in type_counts.items():
                    print(f"   - {tm_type}: {count}")
            print(f"   Saved: {output_regions}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate TMbed results")
    parser.add_argument("results_dir", help="Directory containing TMbed results")
    parser.add_argument("--mode", choices=["auto", "flat", "recursive"], default="auto",
                       help="File discovery mode: auto (default), flat (files in dir), or recursive (subdirs)")
    parser.add_argument("--output-prefix", help="Prefix for output files (e.g., 'rep2' → rep2_summary.tsv)")
    parser.add_argument("--summary-name", default="tmbed_summary.tsv", 
                       help="Filename pattern for summary files (default: tmbed_summary.tsv)")
    parser.add_argument("--regions-name", default="tmbed_tm_regions.tsv",
                       help="Filename pattern for regions files (default: tmbed_tm_regions.tsv)")
    
    args = parser.parse_args()
    aggregate_results(args.results_dir, args.mode, args.output_prefix, args.summary_name, args.regions_name)