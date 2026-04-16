#!/usr/bin/env python3
"""
Convert TSV files to CSV for Excel opening.

Usage:
    # Single file
    python tsv_to_csv.py input.tsv output.csv
    
    # All TSVs in directory (auto-names .csv)
    python tsv_to_csv.py --dir /path/to/tsv/files/
"""

import pandas as pd
import argparse
from pathlib import Path
import sys

def convert_single(input_file, output_file):
    """Convert one TSV to CSV."""
    df = pd.read_csv(input_file, sep='\t')
    df.to_csv(output_file, index=False)
    print(f"✅ Converted: {input_file} → {output_file}")
    print(f"   Rows: {len(df)}, Columns: {len(df.columns)}")

def convert_directory(input_dir):
    """Convert all TSV files in directory to CSV."""
    dir_path = Path(input_dir)
    tsv_files = list(dir_path.glob("*.tsv"))
    
    if not tsv_files:
        print(f"❌ No .tsv files found in {input_dir}")
        return
    
    print(f"📁 Found {len(tsv_files)} TSV files")
    
    for tsv_file in tsv_files:
        csv_file = tsv_file.with_suffix('.csv')
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            df.to_csv(csv_file, index=False)
            print(f"   ✅ {tsv_file.name} → {csv_file.name} ({len(df)} rows)")
        except Exception as e:
            print(f"   ❌ Error with {tsv_file.name}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TSV to CSV for Excel")
    parser.add_argument("input", help="Input TSV file or directory (use with --dir flag)")
    parser.add_argument("output", nargs="?", help="Output CSV file (optional, for single file mode)")
    parser.add_argument("--dir", action="store_true", help="Process all TSVs in directory")
    
    args = parser.parse_args()
    
    if args.dir:
        convert_directory(args.input)
    else:
        if not args.output:
            # Auto-generate output name
            input_path = Path(args.input)
            args.output = str(input_path.with_suffix('.csv'))
        convert_single(args.input, args.output)