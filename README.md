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

# 2. Configure path
# Edit `tmbed_coords.py` and set TMBED_REPO to your local TMbed installation
# Line 17: TMBED_REPO = Path("PATH_TO_YOUR_TMBED_REPO")

# 3. Run full pipeline (example: set1)
python tsv_to_fasta.py --tsv data/set1/set1_protein.tsv --header -o data/set1/set1.fasta

python filter_long_proteins.py \
    --input data/set1/set1.fasta \
    --max-length 2000 \
    --output-normal data/set1/set1_NORMAL.fasta \
    --output-giants data/set1/set1_GIANTS.fasta

python micro_chunker.py \
    --input data/set1/set1_NORMAL.fasta \
    --output-dir data/set1/set1_chunks \
    --chunk-size 50

# Process chunks (manual sequential loop for stability)
for i in $(seq -w 000 158); do
    python tmbed_coords.py \
        --fasta data/set1/set1_chunks/chunk_$i.fasta \
        --output-dir results/set1/
done

# Aggregate results
python aggregate_results.py results/set1/ --mode flat
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

**Configure path:** Edit `tmbed_coords.py` and set `TMBED_REPO`:
```python
TMBED_REPO = Path("/path/to/your/TMbed")
```

---

## Pipeline Steps

### Step 1: Convert TSV to FASTA

Extracts protein sequences from UniProt using IDs in your TSV file.

```bash
python tsv_to_fasta.py \
    --tsv data/your_file.tsv \
    --col 1 \
    --header \
    -o data/output.fasta
```

**Arguments:**
- `--tsv`: Input TSV file (no header by default, use `--header` if first row is column names)
- `--col`: Column **index** (0-based) containing UniProt IDs. Default: `1` (second column)
- `--header`: Use if TSV has a header row
- `--output`: Output FASTA filename

**Features:**
- ✅ Handles isoform IDs automatically (`A0FGR8-2` → fetches isoform 2 from UniProt)
- ✅ Extracts clean IDs from `contam_sp|O77727|K1C15_SHEEP` format
- ✅ Processes headerless proteomics TSVs (common from MaxQuant/Proteome Discoverer)

---

### Step 2: Filter Giant Proteins (Optional but Recommended)

TMbed crashes on proteins >2000 amino acids (OOM errors). Filter them out:

```bash
python filter_long_proteins.py \
    --input proteins.fasta \
    --max-length 2000 \
    --output-normal proteins_NORMAL.fasta \
    --output-giants proteins_GIANTS.fasta
```

**Arguments:**
- `--input`: Input FASTA file
- `--output-normal`: Output file for proteins ≤ max-length (process these with TMbed)
- `--output-giants`: Output file for proteins > max-length (exclude from TMbed)
- `--max-length`: Length threshold (default: 2000)

---

### Step 3: Chunk Large FASTA

For large datasets, split into manageable chunks to prevent OOM errors:

**Option A: Micro-chunks (recommended for 1000-2000aa proteins)**
```bash
python micro_chunker.py \
    --input proteins_NORMAL.fasta \
    --output-dir chunks/ \
    --chunk-size 50
```

**Option B: Standard chunks (for smaller proteins <1000aa)**
```bash
python split_fasta.py \
    -i proteins_NORMAL.fasta \
    -o chunks/ \
    -n 150
```

**Recommendations:**
- Proteins 1500-2000aa: Use `--chunk-size 50` (micro-chunks)
- Proteins 500-1500aa: Use `--chunk-size 100`
- Proteins <500aa: Use `--chunk-size 150`

---

### Step 4: Run TMbed

**Single file:**
```bash
python tmbed_coords.py \
    --fasta proteins.fasta \
    --output-dir results/ \
    --batch-size 2
```

**Batch processing (sequential loop - recommended for stability):**
```bash
# Count how many chunks you have
ls chunks/ | wc -l

# Process each chunk sequentially (change 158 to your max chunk number)
for i in $(seq -w 000 158); do
    echo "Processing chunk $i..."
    python tmbed_coords.py \
        --fasta chunks/chunk_$i.fasta \
        --output-dir results/ \
        --batch-size 2
done
```

**Arguments:**
- `--fasta`: Input FASTA file
- `--output-dir`: Output directory for predictions (default: `tmbed_coordinates`)
- `--batch-size`: Batch size for embedding generation (default: 2, for RTX 5070)
- `--skip-tmbed`: Use existing `predictions.txt` (skip re-running TMbed)
- `--predictions`: Path to existing `predictions.txt` file

---

### Step 5: Aggregate Results

Combine all chunk outputs into summary files:

**Auto-detect structure:**
```bash
python aggregate_results.py results/ --mode auto
```

**Force flat file structure:**
```bash
python aggregate_results.py results/ --mode flat --output-prefix rep1
```

**Force subdirectory structure:**
```bash
python aggregate_results.py results/ --mode recursive
```

**Output files:**
- `{prefix}_summary.tsv` - Per-protein statistics (has TM, TM count, etc.)
- `{prefix}_tm_regions.tsv` - Per-TM-region coordinates (start, end, type, sequence)

---

## Complete Pipeline Example

```bash
# Configuration
REP="set1"
TSV="data/${REP}/${REP}_protein.tsv"
RESULTS="results/${REP}"

# 1. TSV → FASTA
python tsv_to_fasta.py --tsv ${TSV} --header --col 1 -o ${REP}.fasta

# 2. Filter giants
python filter_long_proteins.py \
    --input ${REP}.fasta \
    --max-length 2000 \
    --output-normal data/${REP}/${REP}_NORMAL.fasta \
    --output-giants data/${REP}/${REP}_GIANTS.fasta

# 3. Chunk
python micro_chunker.py \
    --input data/${REP}/${REP}_NORMAL.fasta \
    --output-dir data/${REP}/${REP}_chunks \
    --chunk-size 50

# 4. Process chunks (adjust 158 to your actual chunk count)
mkdir -p ${RESULTS}
for i in $(seq -w 000 158); do
    python tmbed_coords.py \
        --fasta data/${REP}/${REP}_chunks/chunk_$i.fasta \
        --output-dir ${RESULTS}/
done

# 5. Aggregate
python aggregate_results.py ${RESULTS}/ --mode flat --output-prefix ${REP}

# 6. Verify (optional)
python diagnose_missing.py  # Check for data integrity
```

---

## Output File Formats

### `{prefix}_summary.tsv`

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt accession |
| `sequence_length` | Protein length (AA) |
| `tm_count` | Number of TM domains detected |
| `alpha_helices` | Count of alpha-helical TMs |
| `beta_barrels` | Count of beta-barrel TMs |
| `reentrant_helices` | Count of re-entrant helices |
| `has_tm` | Boolean: has ≥1 TM domain |
| `is_multipass` | Boolean: has ≥2 TM domains |
| `topology_string` | Full topology prediction code |

### `{prefix}_tm_regions.tsv`

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt accession |
| `tm_index` | TM domain number (1-indexed) |
| `start` | Start position (1-indexed, UniProt-style) |
| `end` | End position (inclusive) |
| `length` | Length in amino acids |
| `type` | TM type: `TMH` (alpha-helix), `TMB` (beta-barrel), `TMH_re`/`TMB_re` (re-entrant) |
| `topology_code` | Original code: `H`, `h`, `B`, `b` |
| `sequence` | Amino acid sequence of TM region |

**Note:** Proteins with `has_tm=False` (no transmembrane regions) appear in the summary file but NOT in the TM regions file. This is expected behavior.

---

## TMbed Topology Codes

| Code | Type | Description |
|------|------|-------------|
| `H` | TMH | Alpha-helix transmembrane |
| `h` | TMH_re | Re-entrant alpha-helix |
| `B` | TMB | Beta-strand (outer membrane barrel) |
| `b` | TMB_re | Re-entrant beta-strand |
| `S` | SP | Signal peptide |
| `I` | Inside | Cytoplasmic |
| `O` | Outside | Extracellular |

---

## Troubleshooting

### Out of Memory (OOM)
```bash
# Reduce batch size
python tmbed_coords.py -f proteins.fasta -o results/ --batch-size 1

# Use smaller chunks
python micro_chunker.py --input proteins.fasta --output-dir micro/ --chunk-size 25
```

### Giant Proteins Cause Failures
```bash
# Filter out proteins >2000aa before processing
python filter_long_proteins.py --input big.fasta --output-normal normal.fasta --output-giants giants.fasta
```

### Missing Proteins in TM Regions File
Run the diagnostic script:
```bash
python diagnose_missing.py
```

**Expected behavior:** Proteins with `tm_count=0` (no transmembrane domains) appear in the summary file but NOT in the TM regions file. Only proteins with `has_tm=True` have coordinate entries.

### File Not Found / Path Issues
- Use absolute paths when possible
- Ensure output directories exist (auto-created by most scripts)
- Check file permissions: `chmod +x *.py`

### Duplicate Protein IDs
If all proteins show the same ID (e.g., `O00762`):
- **Cause:** FASTA parsing format mismatch
- **Fix:** Scripts now use proper `SeqIO.parse(path, "fasta")` format

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
├── diagnose_missing.py         # Data integrity checker
├── examples/                   # Example inputs/outputs
├── data/                       # Input data (gitignored)
└── results/                    # Output directory (gitignored)
```

**Note:** This repo is cleaned of temporary/development scripts. Core pipeline only.

---

## System Requirements

- **Python:** 3.8+
- **GPU:** CUDA, 8GB+ VRAM (tested on RTX 5070)
- **RAM:** 16GB+ for large datasets
- **Storage:** ~1GB per 1,000 proteins
- **OS:** Linux/WSL (Windows subsystem for Linux)

**Performance:** ~2-3 min per 150-protein chunk on RTX 5070

---

## Citation

- **TMbed:** Remmert et al. (2022) *Bioinformatics* 38(22):5110-5113
- **ProtT5:** Elnaggar et al. (2021) *IEEE TPAMI*

---

## License

MIT License
