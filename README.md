# TMbed Coordinate Extractor

Extract exact transmembrane domain coordinates from protein sequences using [TMbed](https://github.com/HUBioDataLab/TMbed) and ProtT5 embeddings.

## Overview

- **Input:** TSV with UniProt IDs → FASTA → TMbed predictions  
- **Output:** TSV files with exact 1-indexed TM region coordinates
- **Handles:** Large datasets via chunked processing, giant protein filtering (>2000aa)
- **Features:** Automatic isoform handling (e.g., `A0FGR8-2`), headerless TSV support
- **No API keys required** (local processing)

---

## Quick Start

```bash
# 1. Install dependencies
git clone https://github.com/HUBioDataLab/TMbed.git
cd TMbed && pip install -e .
pip install pandas numpy biopython requests

# 2. Configure
# Edit `tmbed_coords.py` line 17: TMBED_REPO = Path("/path/to/TMbed")

# 3. Run pipeline (see Complete Pipeline Example below)
```

---

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

---

## Complete Pipeline Example

Copy-paste this entire block and edit the `CONFIG` section for your replicate:

```bash
#!/bin/bash
# ============================================
# CONFIG - Change these for each replicate
# ============================================
REP_NAME="mem_glyco_rep3"  # Change to rep1, rep2, etc.
TSV_FILE="data/${REP_NAME}/${REP_NAME}_protein.tsv"
BASE_DIR="/home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates"

# ============================================
# 1. TSV → FASTA
# ============================================
echo "Step 1: Converting TSV to FASTA..."
mkdir -p ${BASE_DIR}/data/${REP_NAME}
python ${BASE_DIR}/tsv_to_fasta.py \
    --tsv ${BASE_DIR}/${TSV_FILE} \
    --header \
    --col 1 \
    -o ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_protein.fasta

# ============================================
# 2. Filter giants (>2000aa)
# ============================================
echo "Step 2: Filtering giant proteins..."
python ${BASE_DIR}/filter_long_proteins.py \
    --input ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_protein.fasta \
    --max-length 2000 \
    --output-normal ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_NORMAL.fasta \
    --output-giants ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_GIANTS.fasta

# ============================================
# 3. Chunk into 50-protein batches
# ============================================
echo "Step 3: Chunking FASTA..."
python ${BASE_DIR}/micro_chunker.py \
    --input ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_NORMAL.fasta \
    --output-dir ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_chunks \
    --chunk-size 50

# Count chunks
N_CHUNKS=$(ls ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_chunks/ | wc -l)
echo "Created ${N_CHUNKS} chunks"

# ============================================
# 4. TMbed Processing (Sequential for stability)
# ============================================
echo "Step 4: Running TMbed on ${N_CHUNKS} chunks..."
mkdir -p ${BASE_DIR}/results/${REP_NAME}

for i in $(seq -w 000 $((${N_CHUNKS}-1))); do
    echo "Processing chunk ${i}..."
    python ${BASE_DIR}/tmbed_coords.py \
        --fasta ${BASE_DIR}/data/${REP_NAME}/${REP_NAME}_chunks/chunk_${i}.fasta \
        --output-dir ${BASE_DIR}/results/${REP_NAME}/
done

# ============================================
# 5. Aggregate Results
# ============================================
echo "Step 5: Aggregating results..."
python ${BASE_DIR}/aggregate_results.py \
    ${BASE_DIR}/results/${REP_NAME}/ \
    --mode flat \
    --summary-name "*_summary.tsv" \
    --regions-name "*_tm_regions.tsv" \
    --output-prefix ${REP_NAME}

echo "✓ Pipeline complete for ${REP_NAME}"
echo "Output: ${BASE_DIR}/results/${REP_NAME}/${REP_NAME}_all_summary.tsv"
```

---

## Pipeline Steps (Detailed)

### Step 1: Convert TSV to FASTA

Extracts protein sequences from UniProt using IDs in your TSV file.

```bash
python tsv_to_fasta.py \
    --tsv data/your_file.tsv \
    --col 1 \
    --header \
    -o output.fasta
```

**Arguments:**
- `--tsv`: Input TSV file
- `--col`: Column **index** (0-based) containing UniProt IDs. Default: `1`
- `--header`: Use if TSV has a header row
- `--output`: Output FASTA filename

**Features:**
- ✅ Handles isoform IDs automatically (`A0FGR8-2` → fetches isoform 2)
- ✅ Extracts clean IDs from `contam_sp|O77727|K1C15_SHEEP` format
- ✅ Processes headerless proteomics TSVs

---

### Step 2: Filter Giant Proteins

TMbed crashes on proteins >2000 amino acids. Filter them out:

```bash
python filter_long_proteins.py \
    --input proteins.fasta \
    --max-length 2000 \
    --output-normal proteins_NORMAL.fasta \
    --output-giants proteins_GIANTS.fasta
```

**Arguments:**
- `--input`: Input FASTA file
- `--output-normal`: Proteins ≤ max-length (process with TMbed)
- `--output-giants`: Proteins > max-length (exclude from TMbed)
- `--max-length`: Length threshold (default: 2000)

---

### Step 3: Chunk Large FASTA

Split into manageable chunks to prevent OOM errors:

**Micro-chunks (recommended for 1000-2000aa proteins):**
```bash
python micro_chunker.py \
    --input proteins_NORMAL.fasta \
    --output-dir chunks/ \
    --chunk-size 50
```

**Standard chunks (for <1000aa proteins):**
```bash
python split_fasta.py \
    -i proteins_NORMAL.fasta \
    -o chunks/ \
    -n 150
```

**Recommendations:**
- Proteins 1500-2000aa: `--chunk-size 50`
- Proteins 500-1500aa: `--chunk-size 100`
- Proteins <500aa: `--chunk-size 150`

---

### Step 4: Run TMbed

**Single file:**
```bash
python tmbed_coords.py \
    --fasta proteins.fasta \
    --output-dir results/ \
    --batch-size 2
```

**Batch processing (sequential loop):**
```bash
N_CHUNKS=$(ls chunks/ | wc -l)
for i in $(seq -w 000 $((${N_CHUNKS}-1))); do
    python tmbed_coords.py \
        --fasta chunks/chunk_$i.fasta \
        --output-dir results/ \
        --batch-size 2
done
```

**Arguments:**
- `--fasta`: Input FASTA file
- `--output-dir`: Output directory (default: `tmbed_coordinates`)
- `--batch-size`: Batch size for embedding generation (default: 2)
- `--skip-tmbed`: Use existing `predictions.txt`

---

### Step 5: Aggregate Results

Combine all chunk outputs:

```bash
python aggregate_results.py \
    results/ \
    --mode flat \
    --summary-name "*_summary.tsv" \
    --regions-name "*_tm_regions.tsv" \
    --output-prefix my_experiment
```

**Arguments:**
- `--mode`: `flat` (all files in one dir) or `recursive` (subdirectories)
- `--summary-name`: Pattern for summary files (default: `tmbed_summary.tsv`)
- `--regions-name`: Pattern for regions files (default: `tmbed_tm_regions.tsv`)
- `--output-prefix`: Prefix for aggregated output files

**Output files:**
- `{prefix}_all_summary.tsv` - Per-protein statistics
- `{prefix}_all_tm_regions.tsv` - Per-TM-region coordinates

---

## Output File Formats

### `{prefix}_all_summary.tsv`

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt accession |
| `sequence_length` | Protein length (AA) |
| `tm_count` | Number of TM domains |
| `alpha_helices` | Count of alpha-helical TMs |
| `beta_barrels` | Count of beta-barrel TMs |
| `reentrant_helices` | Count of re-entrant helices |
| `has_tm` | Boolean: has ≥1 TM domain |
| `is_multipass` | Boolean: has ≥2 TM domains |
| `topology_string` | Full topology prediction |

### `{prefix}_all_tm_regions.tsv`

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt accession |
| `tm_index` | TM domain number (1-indexed) |
| `start` | Start position (1-indexed) |
| `end` | End position (inclusive) |
| `length` | Length in AA |
| `type` | `TMH`, `TMB`, `TMH_re`, `TMB_re` |
| `topology_code` | `H`, `h`, `B`, `b` |
| `sequence` | Amino acid sequence |

**Note:** Proteins with `has_tm=False` appear in summary but NOT in regions file (expected behavior).

---

## Troubleshooting

### Out of Memory
```bash
# Reduce batch size
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 1

# Smaller chunks
python micro_chunker.py --input proteins.fasta --output-dir micro/ --chunk-size 25
```

### Giant Proteins Cause Failures
```bash
python filter_long_proteins.py --input big.fasta --output-normal normal.fasta --output-giants giants.fasta
```

### Missing Proteins in TM Regions File
**Expected:** Proteins with `tm_count=0` (no transmembrane domains) only appear in summary file. Only `has_tm=True` proteins have coordinate entries.

### Duplicate Protein IDs in Output
Check that chunk files don't overlap:
```bash
cut -f1 results/*_summary.tsv | sort | uniq -d | wc -l
# Should show 0
```

---

## Repository Structure

```
tmbed-coordinates/
├── tsv_to_fasta.py              # Convert TSV to FASTA (handles isoforms)
├── filter_long_proteins.py      # Filter giants (>2000aa)
├── micro_chunker.py            # Create small chunks (50 proteins)
├── split_fasta.py              # Standard chunking (150 proteins)
├── tmbed_coords.py             # Main TMbed runner
├── aggregate_results.py        # Combine chunk outputs
├── data/                       # Input data (gitignored)
│   └── rep_name/
│       ├── rep_name_protein.fasta
│       ├── rep_name_NORMAL.fasta
│       ├── rep_name_GIANTS.fasta
│       └── rep_name_chunks/
│           ├── chunk_000.fasta
│           ├── chunk_001.fasta
│           └── ...
└── results/                    # Output directory (gitignored)
    └── rep_name/
        ├── rep_name_all_summary.tsv
        ├── rep_name_all_tm_regions.tsv
        └── chunk_000_summary.tsv
        └── chunk_000_tm_regions.tsv
        └── ...
```

---

## System Requirements

- **Python:** 3.8+
- **GPU:** CUDA, 8GB+ VRAM
- **RAM:** 16GB+ for large datasets
- **Storage:** ~1GB per 1,000 proteins
- **OS:** Linux/WSL

**Performance:** ~2-3 min per 150-protein chunk on RTX 5070

---

## Citation

- **TMbed:** Remmert et al. (2022) *Bioinformatics* 38(22):5110-5113
- **ProtT5:** Elnaggar et al. (2021) *IEEE TPAMI*

---

## License

MIT License
