# %% ============================================================================
import os
import json
from collections import defaultdict

import pandas as pd

# %% ============================================================================
project = "SC"
# Load data
print("Loading data...")
expression_data = pd.read_csv(f"{project}/normalized_expression_sparse.csv")
meta_data = pd.read_csv(f"{project}/metadata.csv", index_col=0)  # Ensure index is Cell ID

# %% ============================================================================
print("Converting data...")
# Convert metadata into a dictionary for fast lookup
cell_to_sample = meta_data["sample_id"].to_dict()

# Create a reverse mapping: sample → list of cells
sample_to_cell = defaultdict(list)
for cell, sample in cell_to_sample.items():
    sample_to_cell[sample].append(cell)

# %% ============================================================================
print("Saving samples...")
grouped_by_cell = expression_data.groupby("Cell")

# Create directory for samples
sample_dir = f"{project}/sample_jsons"
os.makedirs(sample_dir, exist_ok=True)

i = 0
for sample in sample_to_cell:
    i += 1
    print(f"Saving sample {sample} ({i}/{len(sample_to_cell)})...")
    sample_file = f"{sample_dir}/{sample}.json"
    if os.path.exists(sample_file):
        continue  # Skip if file already exists

    sample_data = defaultdict(dict)
    for cell in sample_to_cell[sample]:
        cell_genes_df = grouped_by_cell.get_group(cell)
        genes = cell_genes_df["Gene"]
        expression = cell_genes_df["Expression"]

        sample_data[cell] = dict(zip(genes, expression))

    with open(sample_file, "w") as f:
        json.dump(sample_data, f, indent=4)

print("Sample JSON files saved successfully!")
