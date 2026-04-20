#!/bin/bash
cd /home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates
source /home/alexlee/work/onboarding/surfaceome_topology/.venv/bin/activate
for i in {033..038}; do
    if [ ! -f results/chunk_/tmbed_tm_regions.tsv ]; then
        nohup python tmbed_coords.py --fasta data/chunk_.fasta --output-dir results/chunk_ --batch-size 1 > chunk_.log 2>&1 &
        echo 
