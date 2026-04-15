#!/usr/bin/env python3
"""
Convert TSV (with Protein ID column) to FASTA by fetching from UniProt
"""

import pandas as pd
import requests
import time
import sys
from Bio import SeqIO
from io import StringIO
from pathlib import Path

def fetch_uniprot_sequences(accessions, batch_size=100, delay=0.3):
    """Fetch sequences from UniProt"""
    sequences = {}
    
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        query = " OR ".join([f"accession:{acc}" for acc in batch])
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {"query": query, "format": "fasta"}
        
        try:
            time.sleep(delay)
            response = requests.get(url, params=params, timeout=60)
            if response.status_code == 200:
                for record in SeqIO.parse(StringIO(response.text), "fasta"):
                    # Extract accession from UniProt header sp|O00762|UBE2C_HUMAN
                    parts = record.id.split('|')
                    if len(parts) >= 2:
                        acc = parts[1]
                    else:
                        acc = record.id.split('-')[0]
                    sequences[acc] = str(record.seq)
                print(f"  Fetched {i+len(batch)}/{len(accessions)}...")
            else:
                print(f"  HTTP {response.status_code} for batch {i}")
        except Exception as e:
            print(f"  Error: {e}")
    
    return sequences

def main():
    import argparse
    parser = argparse.ArgumentParser(description='TSV to FASTA converter for TMbed')
    parser.add_argument('--tsv', '-i', required=True, help='Input TSV file')
    parser.add_argument('--col', '-c', default='Protein ID', 
                       help='Column name containing UniProt IDs (default: "Protein ID")')
    parser.add_argument('--output', '-o', default='output.fasta', help='Output FASTA file')
    parser.add_argument('--skip-contam', action='store_true', 
                       help='Skip rows where Protein column starts with "contam_"')
    
    args = parser.parse_args()
    
    # Load TSV
    print(f"Loading {args.tsv}...")
    df = pd.read_csv(args.tsv, sep='\t')
    print(f"  Loaded {len(df)} rows")
    
    # Extract accessions
    if args.col not in df.columns:
        print(f"ERROR: Column '{args.col}' not found. Available columns:")
        print(list(df.columns))
        sys.exit(1)
    
    accessions = df[args.col].dropna().unique().tolist()
    
    # Optional: filter out contaminants
    if args.skip_contam and 'Protein' in df.columns:
        mask = ~df['Protein'].str.startswith('contam_', na=False)
        accessions = df[mask][args.col].dropna().unique().tolist()
        print(f"  Filtered contaminants: {len(accessions)} proteins remaining")
    
    print(f"  Unique accessions to fetch: {len(accessions)}")
    
    # Fetch sequences
    print(f"\nFetching from UniProt...")
    sequences = fetch_uniprot_sequences(accessions, batch_size=100, delay=0.3)
    print(f"  Retrieved {len(sequences)}/{len(accessions)} sequences")
    
    # Write FASTA
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        for acc, seq in sequences.items():
            f.write(f">{acc}\n{seq}\n")
    
    print(f"\n✓ Saved to {output_path}")

if __name__ == "__main__":
    main()