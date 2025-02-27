# %% ==============================
import pandas as pd
import numpy as np
import json
from collections import defaultdict
import os

# %%
##===============================
## read metadata
project = "ST"
metadata = pd.read_csv( project + "/raw_metadata.csv", index_col=0, header=0)
metadata["barcode"] = metadata.index.tolist()

spot_ids = []
subject_spot_n = {}
for index, row in metadata.iterrows():
    subject_id = row["sample_name"]
    if subject_id not in subject_spot_n:
        subject_spot_n[subject_id] = 0
    subject_spot_n[subject_id] += 1
    spot_id = subject_id + "_s" + str(subject_spot_n[subject_id])
    spot_ids.append(spot_id)
metadata.index = spot_ids

barcode_to_sid =  dict(zip(metadata["barcode"].tolist(), metadata.index.tolist()))
with open(project+"/barcode_to_spotid.txt", "w") as f:
    json.dump(barcode_to_sid, f)


# %%===============================
## change "sample_id" to "subject_id"

metadata["subject_id"] = metadata["sample_name"]
metadata.drop("sample_name", axis=1, inplace=True)

metadata.loc[metadata["diagnosis"] == "Case", "diagnosis"] = "PD"
metadata.loc[metadata["diagnosis"] == "Control", "diagnosis"] = "HC"
metadata.loc[metadata["diagnosis"] == "ILBD", "diagnosis"] = "ILB"

metadata.loc[metadata["sex"] == 1, "sex"] = "M"
metadata.loc[metadata["sex"] == 2, "sex"] = "F"

metadata["selected_spot"].fillna(True, inplace=True)

## construct new sample_id
metadata["replicate"] = "rep1"
metadata.loc[(metadata["subject_id"] == "BN0347"), "replicate"] = "rep5"
metadata.loc[(metadata["subject_id"] == "BN1554"), "replicate"] = "rep2"
metadata.loc[(metadata["subject_id"] == "BN1957"), "replicate"] = "rep2"

metadata["batch"] = metadata["batch"].apply(lambda x: "batch"+str(x))
metadata["sample_id"] = metadata["diagnosis"] + "_" + metadata["subject_id"] + "_MTG_VisiumST_" + metadata["batch"] + "_" + metadata["replicate"]

## save the new metadata
metadata.to_csv(project + "/metadata.csv")
metadata.info()
# %% ==============================================
# "spot_id","sex", "diagnosis", "selected_spot","seurat_clusters","Spatial_snn_res.0.5","layer_label_v2","smoothed_label_s5",
# "Astrocytes","Endothelial.Cells", "Fibroblast.Like.Cells","Microglia","Oligodendrocytes","OPCs"Pericytes.1","Pericytes.2","T.Cells","cell2loc_sum"

metadata["MajorCellTypes"] = metadata[["Astrocytes","Endothelial.Cells", "Fibroblast.Like.Cells","Microglia","Oligodendrocytes","OPCs","Pericytes.1","Pericytes.2","T.Cells"]].idxmax(axis=1)
# metadata["MajorType"] = metadata["MajorType"].apply(lambda x: x.replace(".",""))

## if the column data is float, keep 4 digits after the decimal point
# Round only float columns to 4 decimal places
metadata[metadata.select_dtypes(include=['float']).columns] = metadata.select_dtypes(include=['float']).round(2)


c_ls = ["sample_id","sex", "diagnosis", "selected_spot","seurat_clusters","Spatial_snn_res.0.5","layer_label_v2","smoothed_label_s5","MajorCellTypes",
        "Astrocytes","Endothelial.Cells", "Fibroblast.Like.Cells","Microglia","Oligodendrocytes","OPCs","Pericytes.1","Pericytes.2","T.Cells","cell2loc_sum"]

metadata_lite = metadata.loc[:,c_ls]
metadata_lite.to_csv(project + "/metadata_lite.csv",index_label="cs_id")

# %% ============================================================================
print("Converting data...")
# Convert metadata into a dictionary for fast lookup
spot_to_sample = metadata["sample_id"].to_dict()
# Save cell_to_sample mapping
with open(f"{project}/spot_to_sample.json", "w") as f:
    json.dump(spot_to_sample, f, indent=2)
    
# Create a reverse mapping: sample → list of spots
sample_to_spot = defaultdict(list)
for cell, sample in spot_to_sample.items():
    sample_to_spot[sample].append(cell)
# Save sample_to_cell mapping
with open(f"{project}/sample_to_spot.json", "w") as f:
    json.dump(sample_to_spot, f)


# %%
print("Loading embedding ....")
embeddings_data = pd.read_csv(project + "/raw_umap_embeddings.csv",index_col=0, header=0)

embeddings_data["UMAP_1"] = embeddings_data["UMAP_1"].apply(lambda x: round(x, 2))
embeddings_data["UMAP_2"] = embeddings_data["UMAP_2"].apply(lambda x: round(x, 2))

## reset index use barcode_cid map
# Reset index and rename using the mapping
print("Renaming....")
embeddings_data = embeddings_data.reset_index()  # Move index to a column
embeddings_data["index"] = embeddings_data["index"].map(barcode_to_sid) # Rename using mapping
embeddings_data = embeddings_data.set_index("index")  # Set the renamed column as index
embeddings_data.to_csv(project + "/umap_embeddings.csv",index_label="cs_id")

# sampling data
data_df = pd.merge(embeddings_data, metadata_lite, left_index=True, right_index=True)
data_df.to_csv(f'{project}umap_embeddings_with_meta.csv', index_label="cs_id")


## ===========================================================
# sampling data, each sample has have same number of spots, total 100k
num_samples = data_df['sample_id'].nunique()
base_rows = 100000 // num_samples  # 1063
extra_rows = 100000 % num_samples   # 78

# Get unique sample IDs and randomly select some for an extra row
sample_ids = data_df['sample_id'].unique()
extra_sample_ids = np.random.choice(sample_ids, extra_rows, replace=False)

# Define a function to sample the required number of rows for each group
def sample_group(group):
    n_to_sample = base_rows + (1 if group.name in extra_sample_ids else 0)
    return group.sample(n=n_to_sample, random_state=12)

# Apply sampling by group
data_df_100k = data_df.groupby('sample_id', group_keys=False).apply(sample_group)
print(data_df_100k.shape)  # Should be (100000, ...)

# data_df_100k = data_df.sample(n=100000, random_state=1)
data_df_100k.to_csv(f'{project}/umap_embeddings_with_meta_100k.csv', index_label="cs_id")

## add sample_id to umap_embedding data
embeddings_data["sample_id"] = embeddings_data.index.map(spot_to_sample)
embeddings_data.to_csv(f"{project}/umap_embeddings_with_sample_id.csv", index_label="cs_id")

## rename imgage file name
subject_to_sample = dict(zip(metadata["subject_id"].tolist(),metadata["sample_id"].tolist()))
files = os.listdir(project + "/images")
for file in files:
    if file.endswith(".png"):
        subject_id = file[:-4].split("_")[-1]
        if subject_id not in subject_to_sample:
            print(subject_id)
            continue
        new_name = subject_to_sample[subject_id] + ".png"
        os.rename(project + "/images/" + file, project + "/images/" + new_name)
    if file.endswith(".tiff"):
        subject_id = file[:-5].split("_")[-1]
        if subject_id not in subject_to_sample:
            print(subject_id)
            continue
        new_name = subject_to_sample[subject_id] + ".tiff"
        os.rename(project + "/images/" + file, project + "/images/" + new_name)

files = os.listdir(project + "/coordinates")
for file in files:
    if file.endswith(".csv"):
        subject_id = file[:-4].split("_")[-1]
        if subject_id not in subject_to_sample:
            print(subject_id)
            continue
        new_name = subject_to_sample[subject_id] + ".csv"
        os.rename(project + "/coordinates/" + file, project + "/coordinates/" + new_name)

files = os.listdir(project + "/coordinates")
for file in files:
    if file.endswith(".csv"):
        df = pd.read_csv(project + "/coordinates/" + file, index_col=0, header=0)
        df.rename(index=barcode_to_sid, inplace=True)
        df.to_csv(project + "/coordinates/" + file, index_label="cs_id")


# %% ============================================================================
# %%
print("Loading expression data...")
## rename_expression_data
expression_data = pd.read_csv(project + "/raw_normalized_expression_sparse.csv",index_col=None, header=0)
## rename "Cell" column use barcode_cid map
print("Renaming....")
expression_data["cs_id"] = expression_data["Spot"].map(barcode_to_sid)

## "Expression" column keep 4 digits after the decimal point
expression_data["Expression"] = expression_data["Expression"].apply(lambda x: round(x, 2))

## save the new expression data
expression_data.to_csv(project + "/normalized_expression_sparse.csv", index=False)


## add sample_id to expression data
expression_data["sample_id"] = expression_data["cs_id"].map(spot_to_sample)

## save expression data
expression_data.to_csv(f"{project}/normalized_expression_sparse_with_sample.csv", index=False)

