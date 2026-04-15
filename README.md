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

**Configure path:** Edit `tmbed_coords.py` line 17:
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

### 2. Run TMbed
```bash
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 2
```

**Arguments:**
- `--fasta`: Input FASTA
- `--output-dir`: Output directory (default: tmbed_coordinates)
- `--batch-size`: Batch size (default: 2, for RTX 5070)
- `--skip-tmbed`: Use existing predictions.txt
- `--predictions`: Path to existing predictions.txt

### 3. Batch Processing (Large Datasets)

Split into chunks:
```bash
python split_fasta.py -i proteins.fasta -o chunks/ -n 150
mkdir -p chunks
python -c "
from Bio import SeqIO
import math
seqs = list(SeqIO.parse('proteins.fasta', 'fasta'))
chunk_size = 150
for i in range(math.ceil(len(seqs)/chunk_size)):
    chunk = seqs[i*chunk_size:(i+1)*chunk_size]
    SeqIO.write(chunk, f'chunks/chunk_{i:03d}.fasta', 'fasta')
"
```

Process all chunks:
```bash
./batch_process.sh chunks/ results/
```

Aggregate results:
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