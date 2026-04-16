# TMbed Coordinate Extractor

Extract exact transmembrane domain coordinates from protein sequences using [TMbed](https://github.com/HUBioDataLab/TMbed) and ProtT5 embeddings.

## Overview

- Input: TSV with UniProt IDs → FASTA → TMbed predictions
- Output: TSV files with exact 1-indexed TM region coordinates
- Handles large datasets via **chunked processing** with giant protein filtering
- Supports both **subdirectory** (rep1-style) and **flat file** (rep2-style) output structures
- No API keys required (local processing)

## Installation

**Prerequisites:**
- Python 3.8+
- CUDA GPU (8GB+ VRAM recommended, tested on RTX 5070)
- TMbed installed locally

**Setup:**
```bash
# Install TMbed
git clone https://github.com/HUBioDataLab/TMbed.git
cd TMbed && pip install -e .

# Install pipeline dependencies
pip install pandas numpy biopython requests
```

**Configure path:** Edit `tmbed_coords.py` and set `TMBED_REPO`:
```python
TMBED_REPO = Path("/path/to/your/TMbed")
```

## Usage

### 1. Convert TSV to FASTA
```bash
python tsv_to_fasta.py -i proteins.tsv -c "Protein ID" -o output.fasta --skip-contam
```

**Arguments:**
- `--tsv`: Input TSV file
- `--col`: Column with UniProt IDs (default: "Protein ID")
- `--output`: Output FASTA file
- `--skip-contam`: Skip contaminant rows

### 2. Filter Giant Proteins (Optional but Recommended)
For datasets with proteins >2000 amino acids (causes TMbed OOM errors):
```bash
python filter_long_proteins.py \
    --input proteins.fasta \
    --output-normal proteins_NORMAL.fasta \
    --output-giants proteins_GIANTS.fasta \
    --max-length 2000
```

**Arguments:**
- `--input`: Input FASTA file
- `--output-normal`: Output file for proteins ≤ max-length (process these)
- `--output-giants`: Output file for proteins > max-length (excluded from TMbed)
- `--max-length`: Length threshold (default: 2000)

### 3. Split Large FASTA (Optional)
For large datasets, split into manageable chunks:
```bash
# Standard chunks (150 proteins each)
python split_fasta.py -i proteins_NORMAL.fasta -o chunks/ -n 150

# Micro-chunks for problematic sections (50 proteins each)
python micro_chunker.py \
    --input proteins_NORMAL.fasta \
    --output-dir micro_chunks/ \
    --chunk-size 50
```

**Note:** Use `micro_chunker.py` when encountering OOM errors or for proteins 1000-2000aa range.

### 4. Run TMbed

**Single file:**
```bash
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 2
```

**Batch processing (multiple chunks):**
```bash
./batch_process.sh chunks/ results/
```

**Sequential loop for safety:**
```bash
for i in $(seq -w 000 158); do
    python tmbed_coords.py \
        --fasta chunks/chunk_$i.fasta \
        --output-dir results/
done
```

**Arguments:**
- `--fasta`: Input FASTA
- `--output-dir`: Output directory (default: tmbed_coordinates)
- `--batch-size`: Batch size (default: 2, for RTX 5070)
- `--skip-tmbed`: Use existing predictions.txt
- `--predictions`: Path to existing predictions.txt

### 5. Aggregate Results

**Auto-detect structure (subdirectories vs flat files):**
```bash
python aggregate_results.py results/
```

**Force specific modes:**
```bash
# For flat file structure (rep2-style: all .tsv files in one directory)
python aggregate_results.py results/ --mode flat --output-prefix rep2

# For subdirectory structure (rep1-style: chunk_000/, chunk_001/, etc.)
python aggregate_results.py results/ --mode recursive

# Custom file patterns
python aggregate_results.py results/ --summary-name "*_summary.tsv" --regions-name "*_regions.tsv"
```

**Arguments:**
- `results_dir`: Directory containing TMbed results
- `--mode`: `auto` (default), `flat`, or `recursive`
- `--output-prefix`: Prefix for output files (e.g., `rep2` → `rep2_summary.tsv`)
- `--summary-name`: Pattern for summary files (default: `tmbed_summary.tsv`)
- `--regions-name`: Pattern for regions files (default: `tmbed_tm_regions.tsv`)

### 6. Convert to CSV for Excel (Optional)
For viewing/filtering in Microsoft Excel:
```bash
# Single file
python tsv_to_csv.py input.tsv output.csv

# Convert all TSVs in directory
python tsv_to_csv.py results/ --dir
```

**Note:** CSV files are disposable copies. Original TSVs remain the master data.

## Directory Structure Examples

### Rep1-style (subdirectories):
```
results/
├── chunk_000/
│   ├── tmbed_summary.tsv
│   └── tmbed_tm_regions.tsv
├── chunk_001/
│   ├── tmbed_summary.tsv
│   └── tmbed_tm_regions.tsv
└── all_proteins_summary.tsv  (aggregated)
```

### Rep2-style (flat files):
```
results/
├── rep2_NORMAL_micro_000_summary.tsv
├── rep2_NORMAL_micro_000_tm_regions.tsv
├── rep2_NORMAL_micro_001_summary.tsv
├── rep2_NORMAL_micro_001_tm_regions.tsv
└── rep2_summary.tsv  (aggregated)
```

## Output Files

### tmbed_summary.tsv
| Column | Description |
|--------|-------------|
| protein_id | UniProt accession |
| sequence_length | Protein length (AA) |
| tm_count | Number of TM domains |
| alpha_helices | Count of alpha-helical TMs |
| beta_barrels | Count of beta-barrel TMs |
| has_tm | Has ≥1 TM domain |
| is_multipass | Has ≥2 TM domains |
| topology_string | Full topology prediction |

### tmbed_tm_regions.tsv
| Column | Description |
|--------|-------------|
| protein_id | UniProt accession |
| tm_index | TM domain number (1-indexed) |
| start | Start position (1-indexed) |
| end | End position (inclusive) |
| length | Length in amino acids |
| type | TM type (TMH/TMB/TMH_re/TMB_re) |
| topology_code | Original code (H/h/B/b) |
| sequence | Amino acid sequence |

## Example Data

Example input and output files are provided in the `examples/` directory:

- `examples/chunk_000.fasta` - Sample input proteins from UniProt
- `examples/chunk_000/tmbed_summary.tsv` - Protein-level summary statistics
- `examples/chunk_000/tmbed_tm_regions.tsv` - Detailed TM region coordinates

### How to Read the Example Files

**View FASTA file (first 3 proteins):**
```bash
head -9 examples/chunk_000.fasta
# Format: >header line, sequence line, empty line (repeated)
```

**View summary statistics:**
```bash
column -t examples/chunk_000/tmbed_summary.tsv | head -10
# Key columns:
#   - protein_id: UniProt accession (e.g., P00533)
#   - tm_count: Number of TM domains detected
#   - has_tm: True if protein has ≥1 TM domain
#   - is_multipass: True if protein has ≥2 TM domains
```

**View TM region coordinates:**
```bash
column -t examples/chunk_000/tmbed_tm_regions.tsv | head -10
# Key columns:
#   - start/end: 1-indexed amino acid positions (UniProt-style)
#   - type: TMH (alpha-helix), TMB (beta-barrel), or re-entrant variants
#   - sequence: Amino acids in the transmembrane region
```

**Quick stats from example:**
```bash
# Count proteins with TM domains
grep -c "True" examples/chunk_000/tmbed_summary.tsv

# Count total TM regions
tail -n +2 examples/chunk_000/tmbed_tm_regions.tsv | wc -l

# Find multi-pass proteins (≥2 TM domains)
awk -F'\t' '$8=="True"' examples/chunk_000/tmbed_summary.tsv | head
```

## TMbed Topology Codes

| Code | Type | Description |
|------|------|-------------|
| H | TMH | Alpha-helix transmembrane |
| h | TMH_re | Re-entrant alpha-helix |
| B | TMB | Beta-strand (outer membrane barrel) |
| b | TMB_re | Re-entrant beta-strand |
| S | SP | Signal peptide |
| I | Inside | Cytoplasmic |
| O | Outside | Extracellular |

## System Requirements

- **Python:** 3.8+
- **GPU:** CUDA, 8GB+ VRAM (tested on RTX 5070)
- **RAM:** 16GB+ for large datasets
- **Storage:** ~1GB per 1,000 proteins

**Performance:** ~2-3 min per 150-protein chunk on RTX 5070

## Troubleshooting

### Out of Memory
```bash
# Reduce batch size
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 1

# Use micro-chunking (smaller chunks)
python micro_chunker.py --input proteins.fasta --output-dir micro/ --chunk-size 25
```

### Giant Proteins Cause Failures
Filter out proteins >2000aa before processing:
```bash
python filter_long_proteins.py --input big.fasta --output-normal normal.fasta --output-giants giants.fasta
```

### Duplicate Protein IDs in Output
If all proteins show the same ID (e.g., O00762):
- **Cause:** Using `fasta-blast` format instead of `fasta` format in parsing
- **Fix:** Ensure scripts use `SeqIO.parse(path, "fasta")` not `SeqIO.parse(path, "fasta-blast")`

### File Not Found
- Use absolute paths
- Ensure output directories exist (auto-created)
- Check file permissions

### Missing TMbed
- Verify `TMBED_REPO` path in `tmbed_coords.py`
- Check CUDA: `python -c "import torch; print(torch.cuda.is_available())"`

## Pipeline Scripts Reference

| Script | Purpose | Key Features |
|--------|---------|--------------|
| `tsv_to_fasta.py` | Convert TSV to FASTA | UniProt ID extraction, contaminant filtering |
| `filter_long_proteins.py` | Length-based filtering | Separates normal (≤2000aa) from giants (>2000aa) |
| `split_fasta.py` | Standard chunking | 150 proteins/chunk default |
| `micro_chunker.py` | Fine-grained chunking | 50 proteins/chunk for OOM-prone data |
| `tmbed_coords.py` | Main TMbed runner | Embedding gen + coordinate extraction |
| `aggregate_results.py` | Combine outputs | Handles flat & recursive structures |
| `tsv_to_csv.py` | Excel conversion | Disposable CSVs for viewing/filtering |
| `batch_process.sh` | Parallel batching | Bash wrapper for multiple chunks |

## Citation

- **TMbed:** Remmert et al. (2022) *Bioinformatics* 38(22):5110-5113
- **ProtT5:** Elnaggar et al. (2021) *IEEE TPAMI*

## License

MIT License
