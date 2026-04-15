#!/bin/bash
# batch_process.sh - Process all FASTA chunks through TMbed
# Usage: ./batch_process.sh <chunks_dir> <results_dir>

CHUNKS_DIR=${1:-"data/chunks"}
RESULTS_DIR=${2:-"results"}
BATCH_SIZE=${3:-2}

# Check if directories exist
if [ ! -d "$CHUNKS_DIR" ]; then
    echo "Error: Chunks directory '$CHUNKS_DIR' not found"
    exit 1
fi

mkdir -p "$RESULTS_DIR"

# Count total chunks
TOTAL=$(ls -1 "$CHUNKS_DIR"/*.fasta 2>/dev/null | wc -l)
if [ "$TOTAL" -eq 0 ]; then
    echo "Error: No .fasta files found in '$CHUNKS_DIR'"
    exit 1
fi

echo "========================================"
echo "Batch Processing $TOTAL chunks"
echo "Input: $CHUNKS_DIR"
echo "Output: $RESULTS_DIR"
echo "Batch size: $BATCH_SIZE"
echo "========================================"

CURRENT=0
for FASTA_FILE in "$CHUNKS_DIR"/*.fasta; do
    CURRENT=$((CURRENT + 1))
    CHUNK_NAME=$(basename "$FASTA_FILE" .fasta)
    OUTPUT_SUBDIR="$RESULTS_DIR/$CHUNK_NAME"
    
    echo ""
    echo "[$CURRENT/$TOTAL] Processing $CHUNK_NAME..."
    
    if [ -f "$OUTPUT_SUBDIR/tmbed_tm_regions.tsv" ]; then
        echo "  ✓ Already processed (skipping)"
        continue
    fi
    
    python tmbed_coords.py \\
        --fasta "$FASTA_FILE" \\
        --output-dir "$OUTPUT_SUBDIR" \\
        --batch-size "$BATCH_SIZE"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Success"
    else
        echo "  ✗ Failed (check logs)"
    fi
done

echo ""
echo "========================================"
echo "Batch processing complete"
echo "========================================"

# Aggregate results if all successful
echo ""
echo "Aggregating results..."
python -c "
import pandas as pd
import os
from pathlib import Path

results_dir = Path(''$RESULTS_DIR'')
summaries = []
regions = []

for chunk_dir in sorted(results_dir.glob(''chunk_*'')):
    summary_file = chunk_dir / ''tmbed_summary.tsv''
    regions_file = chunk_dir / ''tmbed_tm_regions.tsv''
    
    if summary_file.exists():
        summaries.append(pd.read_csv(summary_file, sep=''\t''))
    if regions_file.exists():
        regions.append(pd.read_csv(regions_file, sep=''\t''))

if summaries:
    combined_summary = pd.concat(summaries, ignore_index=True)
    combined_summary.to_csv(results_dir / ''all_proteins_summary.tsv'', sep=''\t'', index=False)
    print(f''✓ Aggregated {len(combined_summary)} proteins to all_proteins_summary.tsv'')
    
    # Print stats
    with_tm = combined_summary[''has_tm''].sum()
    multipass = combined_summary[''is_multipass''].sum()
    print(f''  - {with_tm} proteins with TM domains'')
    print(f''  - {multipass} multi-pass proteins'')

if regions:
    combined_regions = pd.concat(regions, ignore_index=True)
    combined_regions.to_csv(results_dir / ''all_tm_regions.tsv'', sep=''\t'', index=False)
    print(f''✓ Aggregated {len(combined_regions)} TM regions to all_tm_regions.tsv'')
"

echo ""
echo "Done! Check $RESULTS_DIR for aggregated results."