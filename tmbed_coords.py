#!/usr/bin/env python3
"""
TMbed Coordinate Parser - Extracts exact TM region positions
Input: FASTA file with protein sequences
Output: TSV with TM domain coordinates (start, end, type)
"""

import pandas as pd
import numpy as np
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
import re
import argparse

TMBED_REPO = Path("/home/alexlee/work/onboarding/surfaceome_topology/TMbed")

def run_tmbed_embed_predict(fasta_path, output_dir, batch_size=2):
    """Run TMbed embed + predict pipeline"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # FIXED: Unique filenames per chunk
    input_basename = Path(fasta_path).stem
    embeddings_h5 = output_dir / f"{input_basename}_embeddings.h5"
    predictions_txt = output_dir / f"{input_basename}_predictions.txt"
    
    original_dir = os.getcwd()
    
    try:
        os.chdir(TMBED_REPO)
        
        # Step 1: Generate embeddings with ProtT5
        print(f"[1/2] Generating embeddings for {fasta_path}...")
        embed_cmd = [
            "python", "-m", "tmbed", "embed",
            "-f", str(fasta_path),
            "-e", str(embeddings_h5),
            "--batch-size", str(batch_size)
        ]
        result = subprocess.run(embed_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Embedding stderr: {result.stderr}")
        
        if not embeddings_h5.exists():
            raise RuntimeError(f"Embeddings failed: {embeddings_h5} not created")
        
        # Step 2: Predict topology
        print(f"[2/2] Predicting topology...")
        predict_cmd = [
            "python", "-m", "tmbed", "predict",
            "-f", str(fasta_path),
            "-e", str(embeddings_h5),
            "-p", str(predictions_txt),
            "--out-format", "1"  # 3-line format
        ]
        result = subprocess.run(predict_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Prediction stderr: {result.stderr}")
        
        if not predictions_txt.exists():
            raise RuntimeError(f"Prediction failed: {predictions_txt} not created")
        
        print(f"✓ Predictions saved to {predictions_txt}")
        return predictions_txt
        
    finally:
        os.chdir(original_dir)

def parse_tmbed_coordinates(pred_file, fasta_file):
    """
    Parse TMbed format 1 and extract exact coordinates of TM regions
    FIXED: Uses FASTA headers instead of corrupted TMbed prediction headers
    """
    
    results = []
    
    # Load protein IDs and sequences from FASTA (correct source)
    fasta_records = list(SeqIO.parse(fasta_file, "fasta"))
    protein_ids = []
    sequences = {}
    
    for record in fasta_records:
        header = record.id
        if '|' in header:
            parts = header.split('|')
            acc = parts[1] if len(parts) >= 2 else header
        else:
            acc = header.split()[0]
        protein_ids.append(acc)
        sequences[acc] = str(record.seq)
    
    # Parse TMbed predictions (topology only, ignore corrupted headers)
    with open(pred_file) as f:
        lines = f.readlines()
    
    # Extract topology lines (every 3rd line starting from line 3)
    # Format: header, sequence, topology, header, sequence, topology...
    topology_entries = []
    i = 0
    entry_idx = 0
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            i += 1
            continue
        
        # Header line (corrupted, ignore it)
        if line.startswith('>'):
            # Get sequence line
            i += 1
            if i >= len(lines):
                break
            seq_line = lines[i].strip()
            
            # Get topology line (this is what we need)
            i += 1
            if i >= len(lines):
                break
            topology = lines[i].strip()
            
            # Match with correct protein ID from FASTA
            if entry_idx < len(protein_ids):
                protein_id = protein_ids[entry_idx]
                sequence = sequences[protein_id]
                
                # Verify sequence length matches topology length
                if len(topology) != len(sequence):
                    print(f"Warning: Length mismatch for {protein_id}: seq={len(sequence)}, topo={len(topology)}")
                
                topology_entries.append({
                    'protein_id': protein_id,
                    'sequence': sequence,
                    'topology': topology
                })
                entry_idx += 1
        
        i += 1
    
    # Process each entry to find TM regions
    for entry in topology_entries:
        protein_id = entry['protein_id']
        sequence = entry['sequence']
        topology = entry['topology']
        
        # Extract TM regions with coordinates
        tm_regions = []
        in_tm = False
        tm_start = 0
        tm_type = None
        
        for idx, char in enumerate(topology):
            is_tm_char = char in ['H', 'h', 'B', 'b']
            
            if is_tm_char and not in_tm:
                # Start of new TM region
                in_tm = True
                tm_start = idx
                tm_type = char
            
            elif not is_tm_char and in_tm:
                # End of TM region
                in_tm = False
                tm_end = idx
                
                type_map = {
                    'H': 'TMH',
                    'h': 'TMH_re',
                    'B': 'TMB',
                    'b': 'TMB_re'
                }
                
                tm_regions.append({
                    'protein_id': protein_id,
                    'tm_index': len(tm_regions) + 1,
                    'start': tm_start + 1,  # 1-indexed
                    'end': tm_end,
                    'length': tm_end - tm_start,
                    'type': type_map.get(tm_type, 'unknown'),
                    'topology_code': tm_type,
                    'sequence': sequence[tm_start:tm_end]
                })
        
        # Handle case where TM extends to end of sequence
        if in_tm:
            tm_end = len(topology)
            type_map = {'H': 'TMH', 'h': 'TMH_re', 'B': 'TMB', 'b': 'TMB_re'}
            tm_regions.append({
                'protein_id': protein_id,
                'tm_index': len(tm_regions) + 1,
                'start': tm_start + 1,
                'end': tm_end,
                'length': tm_end - tm_start,
                'type': type_map.get(tm_type, 'unknown'),
                'topology_code': tm_type,
                'sequence': sequence[tm_start:tm_end]
            })
        
        # Calculate summary stats
        tm_count = len(tm_regions)
        alpha_helices = sum(1 for r in tm_regions if r['type'] == 'TMH')
        beta_barrels = sum(1 for r in tm_regions if r['type'] == 'TMB')
        
        results.append({
            'protein_id': protein_id,
            'sequence_length': len(sequence),
            'tm_count': tm_count,
            'alpha_helices': alpha_helices,
            'beta_barrels': beta_barrels,
            'reentrant_helices': sum(1 for r in tm_regions if r['type'] == 'TMH_re'),
            'has_tm': tm_count > 0,
            'is_multipass': tm_count >= 2,
            'topology_string': topology,
            'tm_regions': tm_regions
        })
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Run TMbed and extract exact TM region coordinates'
    )
    parser.add_argument('--fasta', '-f', type=str, required=True,
                       help='Input FASTA file with protein sequences')
    parser.add_argument('--output-dir', '-o', type=str, default='tmbed_coordinates',
                       help='Output directory')
    parser.add_argument('--batch-size', '-b', type=int, default=2,
                       help='TMbed batch size (default: 2 for RTX 5070)')
    parser.add_argument('--skip-tmbed', action='store_true',
                       help='Skip TMbed run, parse existing predictions.txt')
    parser.add_argument('--predictions', type=str, default=None,
                       help='Path to existing TMbed predictions.txt')
    
    args = parser.parse_args()
    
    fasta_path = Path(args.fasta).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(exist_ok=True, parents=True)
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    
    # Step 1: Run TMbed or use existing
    if args.skip_tmbed and args.predictions:
        pred_file = Path(args.predictions)
        print(f"Using existing predictions: {pred_file}")
    else:
        pred_file = run_tmbed_embed_predict(fasta_path, output_dir, args.batch_size)
    
    # Get basename for output files
    input_basename = fasta_path.stem
    
    # Step 2: Parse coordinates
    print(f"\nParsing TM coordinates from {pred_file}...")
    results = parse_tmbed_coordinates(pred_file, fasta_path)
    
    # Step 3: Export summary - FIXED unique filename
    summary_df = pd.DataFrame([
        {k: v for k, v in r.items() if k != 'tm_regions'} 
        for r in results
    ])
    summary_tsv = output_dir / f"{input_basename}_summary.tsv"
    summary_df.to_csv(summary_tsv, sep='\t', index=False)
    print(f"✓ Summary saved: {summary_tsv} ({len(summary_df)} proteins)")
    
    # Step 4: Export detailed TM regions (flattened) - FIXED unique filename
    all_regions = []
    for r in results:
        for region in r['tm_regions']:
            all_regions.append(region)
    
    if all_regions:
        regions_df = pd.DataFrame(all_regions)
        regions_tsv = output_dir / f"{input_basename}_tm_regions.tsv"
        regions_df.to_csv(regions_tsv, sep='\t', index=False)
        print(f"✓ TM regions saved: {regions_tsv} ({len(all_regions)} regions)")
        
        # Print sample
        print(f"\nSample TM regions (first 5):")
        print(regions_df.head().to_string())
        
        # Stats
        print(f"\n{'='*60}")
        print(f"SUMMARY")
        print(f"{'='*60}")
        print(f"Total proteins analyzed: {len(summary_df)}")
        print(f"Proteins with TM domains: {summary_df['has_tm'].sum()}")
        print(f"Multi-pass proteins (≥2 TM): {summary_df['is_multipass'].sum()}")
        
        print(f"\nTM type distribution:")
        print(f"  Alpha-helices (TMH): {regions_df[regions_df['type']=='TMH'].shape[0]}")
        print(f"  Beta-barrels (TMB): {regions_df[regions_df['type']=='TMB'].shape[0]}")
        print(f"  Re-entrant (TMH_re/TMB_re): {regions_df[regions_df['type'].str.contains('_re')].shape[0]}")
        
        print(f"\nTop multi-pass proteins:")
        multipass = summary_df[summary_df['is_multipass']].nlargest(5, 'tm_count')
        for _, row in multipass.iterrows():
            print(f"  {row['protein_id']}: {row['tm_count']} TM domains")
    else:
        print("No TM regions found in any protein")
    
    print(f"\n{'='*60}")
    print(f"All outputs in: {output_dir}")
    print(f"  - Summary: {input_basename}_summary.tsv")
    print(f"  - Regions: {input_basename}_tm_regions.tsv")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()