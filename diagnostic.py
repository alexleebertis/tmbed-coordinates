import pandas as pd

tsv_path = "/home/alexlee/work/onboarding/surfaceome_topology/tmbed-coordinates/data/mem_glyco_rep1/mem_glyco_rep1_protein.tsv"
df = pd.read_csv(tsv_path, sep='\t', header=None)

ids = df.iloc[:, 1].dropna().unique()
print(f"Total unique IDs: {len(ids)}")

# Check patterns
import re
uniprot_like = [i for i in ids if re.match(r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$', str(i))]
print(f"Standard UniProt format: {len(uniprot_like)}")

# Sample of non-matching IDs
non_standard = [i for i in ids if not re.match(r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$', str(i))]
print(f"\nFirst 20 non-standard IDs:")
for i in non_standard[:20]:
    print(f"  {i}")