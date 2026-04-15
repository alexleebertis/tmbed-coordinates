# TMbed Coordinate Extractor

Extract exact transmembrane domain coordinates from protein sequences using [TMbed](https://github.com/HUBioDataLab/TMbed) and ProtT5 embeddings.

## Overview

- Input: TSV with UniProt IDs → FASTA → TMbed predictions
- Output: TSV files with exact 1-indexed TM region coordinates
- Handles large datasets via chunked processing
- No API keys required (local processing)

## Installation

**Prerequisites:**
- Python 3.8+
- CUDA GPU (8GB+ VRAM recommended)
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

### 2. Split Large FASTA (Optional)
For large datasets, split into chunks:
```bash
python split_fasta.py -i proteins.fasta -o chunks/ -n 150
```

**Arguments:**
- `--input`: Input FASTA file
- `--output`: Output directory for chunks
- `--chunk-size`: Sequences per chunk (default: 150)

### 3. Run TMbed
Single file:
```bash
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 2
```

Batch processing (multiple chunks):
```bash
./batch_process.sh chunks/ results/
```

**Arguments:**
- `--fasta`: Input FASTA
- `--output-dir`: Output directory (default: tmbed_coordinates)
- `--batch-size`: Batch size (default: 2, for RTX 5070)
- `--skip-tmbed`: Use existing predictions.txt
- `--predictions`: Path to existing predictions.txt

### 4. Aggregate Results
```bash
python aggregate_results.py results/
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

**Out of Memory:**
```bash
# Reduce batch size
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 1
```

**File Not Found:**
- Use absolute paths
- Ensure output directories exist (auto-created)
- Check file permissions

**Missing TMbed:**
- Verify `TMBED_REPO` path in `tmbed_coords.py`
- Check CUDA: `python -c "import torch; print(torch.cuda.is_available())"`

## Citation

- **TMbed:** Remmert et al. (2022) *Bioinformatics* 38(22):5110-5113
- **ProtT5:** Elnaggar et al. (2021) *IEEE TPAMI*

## License

MIT License
