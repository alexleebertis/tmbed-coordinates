TMbed Coordinate Extractor
Extract exact transmembrane domain coordinates (start/end positions) from protein sequences using TMbed and ProtT5 embeddings.
________________________________________
Overview
Converts TSV protein lists → FASTA → TMbed predictions → TSV coordinates. Outputs exact 1-indexed positions for transmembrane regions (UniProt-style).
Key Features
•	Automatic UniProt sequence fetching from TSV input
•	Chunked processing for large datasets
•	Exact TM region coordinate extraction (1-indexed)
•	Support for multiple TM types (alpha-helices, beta-barrels, re-entrant regions)
•	TSV output format
•	Local processing with TMbed + ProtT5 (no API keys required)
________________________________________
Installation
Prerequisites
•	Python 3.8+
•	CUDA-capable GPU recommended (8GB+ VRAM)
•	TMbed repository cloned locally
Step 1: Clone TMbed
git clone https://github.com/HUBioDataLab/TMbed.git
cd TMbed
pip install -e .
Step 2: Install Pipeline Dependencies
pip install pandas numpy biopython requests
Step 3: Configure Path
Edit tmbed_coords.py line 17 to set your TMbed installation path:
TMBED_REPO = Path("/path/to/your/TMbed")
________________________________________
Usage
1. TSV to FASTA Conversion
Converts a TSV file containing UniProt IDs to FASTA format by fetching sequences from UniProt:
python tsv_to_fasta.py \
    --tsv input.tsv \
    --col "Protein ID" \
    --output proteins.fasta \
    --skip-contam
Arguments: | Argument | Description | Default | |———-|————-|———| | --tsv | Input TSV file path | Required | | --col | Column name with UniProt IDs | “Protein ID” | | --output | Output FASTA file | “output.fasta” | | --skip-contam | Skip contaminant rows | False |
2. TMbed Coordinate Extraction
Runs TMbed embed + predict and extracts exact TM region coordinates:
python tmbed_coords.py \
    --fasta proteins.fasta \
    --output-dir results/ \
    --batch-size 2
Arguments: | Argument | Description | Default |
           | --fasta | Input FASTA file | Required |
           | --output-dir | Output directory | “tmbed_coordinates” |
           | --batch-size | TMbed batch size | 2 |
           | --skip-tmbed | Skip TMbed run | False |
           | --predictions | Path to existing predictions.txt | None |
3. Batch Processing (Large Datasets)
For datasets with 1000+ proteins, process in chunks:
Split FASTA into chunks:
mkdir -p chunks
python -c "
from Bio import SeqIO
import math

seqs = list(SeqIO.parse('proteins.fasta', 'fasta'))
chunk_size = 150
num_chunks = math.ceil(len(seqs) / chunk_size)

for i in range(num_chunks):
    chunk = seqs[i*chunk_size:(i+1)*chunk_size]
    SeqIO.write(chunk, f'chunks/chunk_{i:03d}.fasta', 'fasta')
    print(f'Created chunk_{i:03d}.fasta ({len(chunk)} proteins)')
"
Process all chunks:
for f in chunks/*.fasta; do
    name=$(basename $f .fasta)
    python tmbed_coords.py -f "$f" -o "results/$name" -b 2
done
Aggregate results:
python aggregate_results.py results/
________________________________________
Output Files
Each run generates two TSV files:
tmbed_summary.tsv
Per-protein summary statistics:
Column	Description
protein_id	UniProt accession
sequence_length	Protein length in AA
tm_count	Total number of TM domains
alpha_helices	Count of alpha-helical TMs
beta_barrels	Count of beta-barrel TMs
reentrant_helices	Count of re-entrant helices
has_tm	Boolean: has ≥1 TM domain
is_multipass	Boolean: has ≥2 TM domains
topology_string	Full TMbed topology prediction
tmbed_tm_regions.tsv
Detailed TM region coordinates:
Column	Description
protein_id	UniProt accession
tm_index	TM domain number (1-indexed)
start	Start position (1-indexed, inclusive)
end	End position (inclusive)
length	Length in amino acids
type	TM type (TMH, TMB, TMH_re, TMB_re)
topology_code	Original TMbed code (H, h, B, b)
sequence	Amino acid sequence of TM region
Example Output
Summary:
protein_id    sequence_length    tm_count    has_tm    is_multipass
P00533        1210               1           True      False
P04626        1255               1           True      False
P21860        1356               7           True      True
Regions:
protein_id    tm_index    start    end    length    type    sequence
P21860        1           23       45     23        TMH     GLLLLPLLLLLLLPVLLLLGLL
P21860        2           89       111    23        TMH     VLLLLLLLLLLLLLLLLLLLLLL
________________________________________
TMbed Topology Codes
Code	Type	Description
H	TMH	Alpha-helix transmembrane domain
h	TMH_re	Re-entrant alpha-helix
B	TMB	Beta-strand (outer membrane beta barrel)
b	TMB_re	Re-entrant beta-strand
S	SP	Signal peptide
I	Inside	Cytoplasmic/interior side
O	Outside	Extracellular/exterior side
________________________________________
System Requirements
Component	Requirement
Python	3.8 or higher
GPU	CUDA-capable, 8GB+ VRAM (recommended)
RAM	16GB+ for large datasets
Storage	~1GB per 1,000 proteins for embeddings
OS	Linux (tested), macOS, Windows WSL
Performance Benchmarks
•	RTX 5070 (12GB): ~2-3 min per 150-protein chunk (batch_size=2)
•	RTX 4090 (24GB): ~1-2 min per 150-protein chunk (batch_size=4)
•	CPU-only: ~10-15 min per 150-protein chunk (not recommended)
________________________________________
File Structure
tmbed-coordinates/
├── tsv_to_fasta.py          # TSV → FASTA converter
├── tmbed_coords.py          # TMbed runner + coordinate parser
├── aggregate_results.py     # Result aggregation tool
├── batch_process.sh         # Batch processing script
├── README.md                # This file
├── requirements.txt         # Python dependencies
├── .gitignore              # Git ignore rules
└── LICENSE                 # MIT License
________________________________________
Troubleshooting
Out of Memory (OOM) Errors
•	Reduce --batch-size (try 1 instead of 2)
•	Process smaller chunks (<100 proteins)
•	Close other GPU applications
Missing Embeddings File
•	Check TMbed installation path in TMBED_REPO
•	Verify CUDA is available: python -c "import torch; print(torch.cuda.is_available())"
UniProt Fetch Failures
•	Check internet connection
•	Try smaller batch sizes in tsv_to_fasta.py (edit batch_size parameter)
•	Some accessions may be deprecated; check UniProt directly
File Not Found Errors
•	Use absolute paths for input files
•	Ensure output directories exist (scripts auto-create them)
•	Check file permissions
________________________________________
Citation
If you use this pipeline, please cite:
1.	TMbed: Remmert, M., et al. (2022). TMbed: Transmembrane protein topology prediction using ProtT5 embeddings. Bioinformatics, 38(22), 5110-5113. doi:10.1093/bioinformatics/btac669
2.	ProtT5: Elnaggar, A., et al. (2021). ProtTrans: Toward understanding the language of life through self-supervised learning. IEEE Transactions on Pattern Analysis and Machine Intelligence. doi:10.1109/TPAMI.2021.3095381
________________________________________
Contributing
Contributions welcome! Please: 1. Fork the repository 2. Create a feature branch 3. Submit a pull request
For bugs or feature requests, open an issue on GitHub.
________________________________________
License
MIT License - see LICENSE file for details.
________________________________________
Contact
For questions about this pipeline, open an issue on GitHub.
For TMbed-specific issues, visit: https://github.com/HUBioDataLab/TMbed/issues
