#!/usr/bin/env python3
"""
Convert TSV to FASTA by fetching from UniProt.
Handles isoforms (e.g., A0FGR8-2) by querying base accession.
"""

import pandas as pd
import requests
import time
import sys
import re
from Bio import SeqIO
from io import StringIO

def fetch_uniprot_sequences(accessions, batch_size=100, delay=0.3):
    """Fetch sequences from UniProt - returns dict of all sequences including isoforms"""
    sequences = {}

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        query = " OR ".join([f"accession:{acc}" for acc in batch])
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {"query": query, "format": "fasta", "includeIsoform": "true"}  # <-- Added isoform flag

        try:
            time.sleep(delay)
            response = requests.get(url, params=params, timeout=60)
            if response.status_code == 200:
                for record in SeqIO.parse(StringIO(response.text), "fasta"):
                    # Extract full ID including isoform (e.g., A0FGR8-2)
                    parts = record.id.split('|')
                    if len(parts) >= 2:
                        full_id = parts[1]  # This includes the -2 for isoforms
                    else:
                        full_id = record.id.split('-')[0]
                    sequences[full_id] = str(record.seq)
                print(f"  Fetched batch {i//batch_size + 1}: {len(sequences)} total sequences")
            else:
                print(f"  HTTP {response.status_code} for batch {i}")
        except Exception as e:
            print(f"  Error in batch {i}: {e}")

    return sequences

def extract_id(raw_id):
    """Extract ID from sp|ID|NAME or contam_sp|ID|NAME format, keep isoform suffix"""
    raw_str = str(raw_id).strip()
    if '|' in raw_str:
        parts = raw_str.split('|')
        if len(parts) >= 2:
            return parts[1]  # Returns ID with isoform suffix if present (e.g., A0FGR8-2)
    return raw_str

def get_base_accession(isoform_id):
    """Strip isoform suffix for API query (A0FGR8-2 -> A0FGR8)"""
    # Remove dash and number at end
    return re.sub(r'-\d+$', '', isoform_id)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', '-i', required=True, help='Input TSV')
    parser.add_argument('--col', '-c', type=int, default=1, help='Column index (0-based)')
    parser.add_argument('--header', action='store_true', help='TSV has header row')
    parser.add_argument('--output', '-o', default='output.fasta', help='Output FASTA')
    args = parser.parse_args()

    # Load TSV
    print(f"Loading {args.tsv}...")
    if args.header:
        df = pd.read_csv(args.tsv, sep='\t', header=0, dtype=str)
    else:
        df = pd.read_csv(args.tsv, sep='\t', header=None, dtype=str)

    print(f"  Loaded {len(df)} rows")

    # Get full IDs (including isoform suffixes like -2, -3)
    raw_ids = df.iloc[:, args.col].dropna().tolist()

    # Skip header text if present
    full_ids = []
    for rid in raw_ids:
        if str(rid) == "Protein ID":
            continue
        full_ids.append(extract_id(rid))

    # Remove duplicates
    seen = set()
    unique_full_ids = []
    for fid in full_ids:
        if fid and fid not in seen:
            seen.add(fid)
            unique_full_ids.append(fid)

    print(f"  Unique IDs to fetch: {len(unique_full_ids)}")
    print(f"  Sample: {unique_full_ids[:5]}")

    # Extract base accessions for API query (strip -2, -3, etc.)
    base_accessions = []
    for fid in unique_full_ids:
        base = get_base_accession(fid)
        if base not in base_accessions:
            base_accessions.append(base)

    print(f"  Base accessions to query: {len(base_accessions)}")

    # Fetch from UniProt (this gets canonical + all isoforms)
    print(f"\nFetching from UniProt (including isoforms)...")
    all_sequences = fetch_uniprot_sequences(base_accessions, batch_size=100, delay=0.5)
    print(f"  Total sequences retrieved: {len(all_sequences)}")

    # Match requested IDs to retrieved sequences
    found_ids = []
    missing_ids = []

    for req_id in unique_full_ids:
        if req_id in all_sequences:
            found_ids.append(req_id)
        else:
            # Check if base accession exists (fallback to canonical)
            base = get_base_accession(req_id)
            if base in all_sequences:
                # Use canonical sequence if specific isoform not found
                all_sequences[req_id] = all_sequences[base]
                found_ids.append(req_id)
                print(f"  Note: {req_id} not found, using canonical {base}")
            else:
                missing_ids.append(req_id)

    print(f"\n  Found: {len(found_ids)}/{len(unique_full_ids)}")
    if missing_ids:
        print(f"  Missing: {len(missing_ids)} (saved to missing_ids.txt)")
        with open('missing_ids.txt', 'w') as f:
            f.write('\n'.join(missing_ids))

    # Write FASTA
    written = 0
    with open(args.output, 'w') as f:
        for uid in unique_full_ids:
            if uid in all_sequences:
                f.write(f">{uid}\n{all_sequences[uid]}\n")
                written += 1

    print(f"\n✓ Wrote {written} sequences to {args.output}")

if __name__ == "__main__":
    main()