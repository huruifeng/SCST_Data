import os
import json
import pandas as pd
from collections import defaultdict

# Load data
print("Loading data...")
expression_data = pd.read_csv("SC/normalized_expression_sparse_with_sample.csv")
meta_data = pd.read_csv("SC/metadata.csv", index_col=0)  # Ensure index is Cell ID

# Create directory for JSON output
output_dir = "SC/sample_jsons_gene_cell"
os.makedirs(output_dir, exist_ok=True)

# Group by sample
print("Grouping by sample and saving JSON files...")
grouped_by_sample = expression_data.groupby("sample_id")

# Iterate over each sample and save as JSON
for i, (sample, df) in enumerate(grouped_by_sample, 1):
    print(f"{i}. Processing sample: {sample}")

    sample_dict = defaultdict(dict)  # {Gene: {Cell: Expression}}

    for _, row in df.iterrows():
        gene = row["Gene"]
        cell = row["Cell"]
        expression = row["Expression"]

        # Store gene expression per cell
        sample_dict[gene][cell] = expression

    # Save JSON file
    file_path = os.path.join(output_dir, f"{sample}.json")
    with open(file_path, "w") as f:
        json.dump(sample_dict, f, indent=4)

print("All samples saved successfully!")
