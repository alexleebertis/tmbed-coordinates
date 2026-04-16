#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import argparse

def aggregate_chunk_results(chunk_results_dir, output_dir):
    """Aggregate all TMbed results from chunk_032b subdirectories."""
    chunk_path = Path(chunk_results_dir).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find ALL subdirectories (recursive) that contain tmbed_summary.tsv
    result_dirs = []
    for subdir in chunk_path.rglob("*"):
        if subdir.is_dir() and (subdir / "tmbed_summary.tsv").exists():
            result_dirs.append(subdir)
    
    if not result_dirs:
        print(f"❌ No result directories with tmbed_summary.tsv found in {chunk_path}")
        return
    
    print(f"🔬 Found {len(result_dirs)} result folders to aggregate")
    
    # Aggregate all summaries and TM regions
    all_summaries = []
    all_tm_regions = []
    processed_proteins = []
    
    for result_dir in sorted(result_dirs):
        summary_file = result_dir / "tmbed_summary.tsv"
        tm_file = result_dir / "tmbed_tm_regions.tsv"
        
        try:
            # Read summary
            df = pd.read_csv(summary_file, sep='\t')
            all_summaries.append(df)
            
            if 'protein_id' in df.columns:
                processed_proteins.extend(df['protein_id'].tolist())
            elif 'ID' in df.columns:
                processed_proteins.extend(df['ID'].tolist())
                
            # Read TM regions if exists
            if tm_file.exists():
                tm_df = pd.read_csv(tm_file, sep='\t')
                all_tm_regions.append(tm_df)
                
            print(f"✅ {result_dir.name}: {len(df)} proteins")
            
        except Exception as e:
            print(f"⚠️  Error reading {result_dir}: {e}")
    
    # Save aggregated results
    if all_summaries:
        master_summary = pd.concat(all_summaries, ignore_index=True)
        master_summary.to_csv(output_path / "tmbed_summary.tsv", sep='\t', index=False)
        print(f"\n📊 Master summary: {len(master_summary)} total entries")
        print(f"   Saved to: {output_path / 'tmbed_summary.tsv'}")
    
    if all_tm_regions:
        master_tm = pd.concat(all_tm_regions, ignore_index=True)
        master_tm.to_csv(output_path / "tmbed_tm_regions.tsv", sep='\t', index=False)
        print(f"📊 Master TM regions: {len(master_tm)} total entries")
        print(f"   Saved to: {output_path / 'tmbed_tm_regions.tsv'}")
    
    # Create completion report
    unique_proteins = sorted(set(processed_proteins))
    report = output_path / "aggregation_report.txt"
    with open(report, 'w') as f:
        f.write(f"TMbed Aggregation Report for chunk_032b\n")
        f.write("="*50 + "\n")
        f.write(f"Total result folders processed: {len(result_dirs)}\n")
        f.write(f"Total protein entries: {len(processed_proteins)}\n")
        f.write(f"Unique proteins: {len(unique_proteins)}\n\n")
        f.write("Unique proteins processed:\n")
        for pid in unique_proteins:
            f.write(f"  - {pid}\n")
    
    print(f"\n📄 Report saved: {report}")
    print(f"🎉 Aggregation complete! All chunk_032b results unified in {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate TMbed results from chunk_032b")
    parser.add_argument("--input-dir", 
                        default="/home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates/results/chunk_032b",
                        help="Directory containing all chunk_032b subresults")
    parser.add_argument("--output-dir", 
                        default="/home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates/results/chunk_032b_complete",
                        help="Output directory for aggregated results")
    
    args = parser.parse_args()
    aggregate_chunk_results(args.input_dir, args.output_dir)