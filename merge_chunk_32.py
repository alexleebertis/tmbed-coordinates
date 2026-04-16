#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import shutil

def merge_chunk_32(results_base_dir):
    """Merge chunk_032a and chunk_032b_complete into chunk_032."""
    base_path = Path(results_base_dir)
    
    # Paths
    chunk_032a = base_path / "chunk_032a"
    chunk_032b = base_path / "chunk_032b_complete"
    chunk_032 = base_path / "chunk_032"
    
    # Ensure output exists
    chunk_032.mkdir(parents=True, exist_ok=True)
    
    print("🔬 Merging chunk_032a + chunk_032b_complete → chunk_032")
    print("="*60)
    
    # Track all proteins
    all_proteins = []
    sources = []
    
    # Process tmbed_summary.tsv
    summary_files = []
    if (chunk_032a / "tmbed_summary.tsv").exists():
        summary_files.append(("032a", chunk_032a / "tmbed_summary.tsv"))
    if (chunk_032b / "tmbed_summary.tsv").exists():
        summary_files.append(("032b", chunk_032b / "tmbed_summary.tsv"))
    
    all_summaries = []
    for source, file_path in summary_files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            df['source'] = source  # Track where it came from
            all_summaries.append(df)
            
            # Track proteins
            id_col = 'protein_id' if 'protein_id' in df.columns else 'ID'
            for pid in df[id_col].tolist():
                all_proteins.append(pid)
                sources.append(source)
                
            print(f"✅ {source}: {len(df)} proteins from tmbed_summary.tsv")
        except Exception as e:
            print(f"⚠️  Error reading {file_path}: {e}")
    
    if all_summaries:
        master_summary = pd.concat(all_summaries, ignore_index=True)
        # Reorder columns to put source first
        cols = ['source'] + [c for c in master_summary.columns if c != 'source']
        master_summary = master_summary[cols]
        master_summary.to_csv(chunk_032 / "tmbed_summary.tsv", sep='\t', index=False)
        print(f"\n📊 Master summary: {len(master_summary)} total entries → {chunk_032 / 'tmbed_summary.tsv'}")
    
    # Process tmbed_tm_regions.tsv
    tm_files = []
    if (chunk_032a / "tmbed_tm_regions.tsv").exists():
        tm_files.append(("032a", chunk_032a / "tmbed_tm_regions.tsv"))
    if (chunk_032b / "tmbed_tm_regions.tsv").exists():
        tm_files.append(("032b", chunk_032b / "tmbed_tm_regions.tsv"))
    
    all_tm = []
    for source, file_path in tm_files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            df['source_chunk'] = source
            all_tm.append(df)
            print(f"✅ {source}: {len(df)} TM regions")
        except Exception as e:
            print(f"⚠️  Error reading {file_path}: {e}")
    
    if all_tm:
        master_tm = pd.concat(all_tm, ignore_index=True)
        master_tm.to_csv(chunk_032 / "tmbed_tm_regions.tsv", sep='\t', index=False)
        print(f"📊 Master TM regions: {len(master_tm)} total entries → {chunk_032 / 'tmbed_tm_regions.tsv'}")
    
    # Copy other files
    if (chunk_032a / "predictions.txt").exists():
        shutil.copy(chunk_032a / "predictions.txt", chunk_032 / "predictions_032a.txt")
        print(f"📄 Copied predictions.txt → predictions_032a.txt")
    
    if (chunk_032b / "aggregation_report.txt").exists():
        shutil.copy(chunk_032b / "aggregation_report.txt", chunk_032 / "aggregation_report_032b.txt")
        print(f"📄 Copied aggregation_report.txt → aggregation_report_032b.txt")
    
    # Create master report
    report_file = chunk_032 / "MASTER_REPORT.txt"
    with open(report_file, 'w') as f:
        f.write("CHUNK 032 MASTER REPORT\n")
        f.write("="*60 + "\n")
        f.write(f"Combined: chunk_032a + chunk_032b\n")
        f.write(f"Total protein entries: {len(all_proteins)}\n")
        f.write(f"Unique proteins: {len(set(all_proteins))}\n\n")
        f.write("Source breakdown:\n")
        f.write(f"  - 032a: {sources.count('032a')} entries\n")
        f.write(f"  - 032b: {sources.count('032b')} entries\n\n")
        f.write("Files in this folder:\n")
        for file in sorted(chunk_032.iterdir()):
            if file.is_file():
                size = file.stat().st_size
                f.write(f"  - {file.name} ({size:,} bytes)\n")
    
    print(f"\n📄 MASTER_REPORT.txt created")
    print(f"🎉 chunk_032 is complete with all merged data!")
    print(f"\nFinal location: {chunk_032}")

if __name__ == "__main__":
    import sys
    base_dir = "/home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates/results"
    merge_chunk_32(base_dir)