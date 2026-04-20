#!/usr/bin/env python3
import os

rep1 = "results/mem_glyco_rep2/mem_glyco_rep2_all_summary.tsv"
rep2 = "results/mem_glyco_rep3/mem_glyco_rep3_all_summary.tsv"

print("DIAGNOSTIC: Are rep1 and rep2 identical?")
print("="*50)

# Check if files exist
if not os.path.exists(rep1):
    print(f"Missing: {rep1}")
if not os.path.exists(rep2):
    print(f"Missing: {rep2}")
    
if os.path.exists(rep1) and os.path.exists(rep2):
    # File size
    size1 = os.path.getsize(rep1)
    size2 = os.path.getsize(rep2)
    print(f"Rep1: {size1} bytes | Rep2: {size2} bytes")
    
    # Line count
    with open(rep1) as f:
        lines1 = sum(1 for _ in f)
    with open(rep2) as f:
        lines2 = sum(1 for _ in f)
    print(f"Rep1: {lines1} lines | Rep2: {lines2} lines")
    
    # Check first 5 protein IDs
    import pandas as pd
    df1 = pd.read_csv(rep1, sep='\t')
    df2 = pd.read_csv(rep2, sep='\t')
    
    print(f"\nRep1 proteins: {df1['protein_id'].nunique()}")
    print(f"Rep2 proteins: {df2['protein_id'].nunique()}")
    print(f"\nFirst 5 rep1: {df1['protein_id'].head().tolist()}")
    print(f"First 5 rep2: {df2['protein_id'].head().tolist()}")
    
    # Check if identical
    same_ids = set(df1['protein_id']) == set(df2['protein_id'])
    print(f"\nSame protein sets: {same_ids}")
    
    if same_ids:
        print("\n🔴 WARNING: Same proteins found in both!")
        print("Checking input TSVs...")
        # Check if input TSVs are different
        tsv1 = "data/mem_glyco_rep1/mem_glyco_rep1_protein.tsv"
        tsv2 = "data/mem_glyco_rep2/mem_glyco_rep2_protein.tsv"
        if os.path.exists(tsv1) and os.path.exists(tsv2):
            df_tsv1 = pd.read_csv(tsv1, sep='\t', header=None)
            df_tsv2 = pd.read_csv(tsv2, sep='\t', header=None)
            print(f"TSV1 rows: {len(df_tsv1)} | TSV2 rows: {len(df_tsv2)}")
            print(f"TSV1 col1 sample: {df_tsv1.iloc[:3, 1].tolist()}")
            print(f"TSV2 col1 sample: {df_tsv2.iloc[:3, 1].tolist()}")