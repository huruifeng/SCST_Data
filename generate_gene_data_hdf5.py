import pandas as pd
import json
import os
import h5py

# %% =============================================
project = "SC"
# Load normalized expression data (sparse format)
print("Reading data...")
expression_data = pd.read_csv(f"{project}/normalized_expression_sparse.csv",index_col=None, header=0)

# Group data by gene
grouped_by_gene = expression_data.groupby("Gene")

# %% =============================================
# Create directory for genes
os.makedirs(f"{project}/gene_h5", exist_ok=True)

# Save each gene as a separate JSON file
print("Saving gene...")
i=0
for gene, df in grouped_by_gene:
    i += 1
    print(i,gene)
    safe_gene_name = gene.replace("/", "_")

    # Create HDF5 file
    file_name = f"{project}/gene_h5/{safe_gene_name}"
    if not os.path.exists(file_name):
        df.to_hdf(file_name + ".h5", key='cell_expression', mode='w')
        df.to_csv(file_name + ".csv", index=False)
    else:
        print(f"File {file_name} already exists. Skipping...")

