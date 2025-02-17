# %% ============================================================================
import os
import h5py
import json
from collections import defaultdict

import pandas as pd

# %% ============================================================================
project = "SC"
#Load data
print("Loading data...HDF5")
expression_data = pd.read_csv(f"{project}/normalized_expression_sparse_with_sample.csv")

# %% ============================================================================
print("Saving samples...")
grouped_by_sample = expression_data.groupby("sample_id")

# Create directory for samples
sample_dir = f"{project}/sample_h5"
os.makedirs(sample_dir, exist_ok=True)

i = 0
for sample_id, df in grouped_by_sample:
    print(i, sample_id)
    df.drop("sample_id", axis=1, inplace=True)
    df.to_hdf(os.path.join(sample_dir, f"{sample_id}.h5"), key="cell_gene_expression", mode="w")
    df.to_csv(os.path.join(sample_dir, f"{sample_id}.csv"), index=False)

    i += 1

print("Sample HDF5 files saved successfully!")
