# %% ==============================
import pandas as pd
import json
from collections import defaultdict
import os
import numpy as np

# %%
##===============================
## read metadata
project = "SC"
metadata = pd.read_csv( project + "/raw_metadata.csv", index_col=0, header=0)
# cell_ids = ["c_"+str(i) for i in range(1,metadata.shape[0]+1)]
# metadata["cs_id"] = cell_ids

cell_ids = []
subject_cell_n = {}
for index, row in metadata.iterrows():
    subject_id = row["sample_id"]
    if subject_id not in subject_cell_n:
        subject_cell_n[subject_id] = 0
    subject_cell_n[subject_id] += 1
    c_id = subject_id + "_c" + str(subject_cell_n[subject_id])
    cell_ids.append(c_id)
metadata["cs_id"] = cell_ids


barcode_to_cid =  metadata["cs_id"].to_dict()
with open(project+"/barcode_to_cid.txt", "w") as f:
    json.dump(barcode_to_cid, f)

metadata["barcode"] = metadata.index.tolist()
metadata = metadata.set_index("cs_id")

# %%===============================
## change "sample_id" to "subject_id"
if "subject_id" not in metadata.columns:
    metadata["subject_id"] = metadata["sample_id"]
    metadata.drop("sample_id", axis=1, inplace=True)

## construct new sample_id
metadata["replicate"] = "rep1"
metadata.loc[(metadata["subject_id"] == "BN0464") & (metadata["batch"] == "batch26"), "replicate"] = "rep2"
metadata.loc[(metadata["subject_id"] == "BN0934") & (metadata["batch"] == "batch26"), "replicate"] = "rep2"
metadata.loc[(metadata["subject_id"] == "BN1076") & (metadata["batch"] == "batch26"), "replicate"] = "rep2"


metadata["sample_id"] = metadata["case"] + "_" + metadata["subject_id"] + "_MTG_snRNAseq_" +metadata["batch"] + "_" + metadata["replicate"]

all_samples = metadata["sample_id"].unique().tolist()
with open(project + "/sample_list.json", "w") as f:
    json.dump(sorted(all_samples), f)

## save the new metadata
metadata.to_csv(project + "/metadata.csv")
metadata.info()
# %% ==============================================
## "orig.ident","nCount_RNA","nFeature_RNA","sample_id",
# "case","batch","sex","RIN","PMI","age","age_bracket","PMI_bracket",
# "RIN_bracket","percent.mt","RNA_snn_res.1.5",
# "seurat_clusters","class","dblscore","MajorCellTypes",
# "CellSubtypes","lbscore","mmse","updrs","sumlbd","ncxtlbd","plaqt","tanglt",
# "Complex_Assignment","majormarker",barcode, subject_id, replicate ,sample_id
meta_list = ["sample_id","case","sex","age","seurat_clusters","MajorCellTypes","CellSubtypes","mmse","updrs"]
metadata_lite = metadata.loc[:,meta_list]
metadata_lite["mmse"].fillna(0, inplace=True)
metadata_lite["updrs"].fillna(0, inplace=True)

metadata_lite.to_csv(project + "/metadata_lite.csv")
with open(project + "/meta_list.json", "w") as f:
    json.dump(sorted(meta_list), f)


sample_meta_list = ["sample_id","case","sex","age","mmse","updrs"]
sample_meta = metadata_lite.loc[:,sample_meta_list]
sample_meta = sample_meta.drop_duplicates()
sample_meta = sample_meta.set_index("sample_id")
sample_meta.to_csv(project + "/sample_metadata.csv")

# %% ============================================================================
print("Converting data...")
# Convert metadata into a dictionary for fast lookup
cell_to_sample = metadata["sample_id"].to_dict()
# Save cell_to_sample mapping
with open(f"{project}/cell_to_sample.json", "w") as f:
    json.dump(cell_to_sample, f, indent=2)
    
# Create a reverse mapping: sample → list of cells
sample_to_cell = defaultdict(list)
for cell, sample in cell_to_sample.items():
    sample_to_cell[sample].append(cell)
# Save sample_to_cell mapping
with open(f"{project}/sample_to_cell.json", "w") as f:
    json.dump(sample_to_cell, f, indent=2)

# %%
print("Loading embedding ....")
embeddings_data = pd.read_csv(project + "/raw_umap_embeddings.csv",index_col=0, header=0)

embeddings_data["UMAP_1"] = embeddings_data["UMAP_1"].apply(lambda x: round(x, 2))
embeddings_data["UMAP_2"] = embeddings_data["UMAP_2"].apply(lambda x: round(x, 2))

## reset index use barcode_cid map
# Reset index and rename using the mapping
print("Renaming....")
embeddings_data = embeddings_data.reset_index()  # Move index to a column
embeddings_data["index"] = embeddings_data["index"].map(barcode_to_cid) # Rename using mapping
embeddings_data = embeddings_data.set_index("index")  # Set the renamed column as index
embeddings_data.to_csv(project + "/umap_embeddings.csv",index_label="cs_id")

# sampling data
data_df = pd.merge(embeddings_data, metadata_lite, left_index=True, right_index=True)
data_df.to_csv(f'{project}/umap_embeddings_with_meta.csv', index_label="cs_id")

## save each column data into a separate json file
os.makedirs(project+"/metas", exist_ok=True)
data_df_by_sample = data_df.groupby("sample_id")
for sample_id, df in data_df_by_sample:
    os.makedirs(project+"/metas/"+sample_id, exist_ok=True)
    for col in df.columns:
        if col == "sample_id":
            continue
        with open(project + f"/metas/{sample_id}/{col}.json", "w") as f:
            json.dump(df[col].to_dict(), f, indent=2)    

## ===========================================================
# sampling data, each sample has have same number of spots, total 100k
n_rows = 100000
num_samples = data_df['sample_id'].nunique()
base_rows = n_rows // num_samples 
extra_rows = n_rows % num_samples

# Get unique sample IDs and randomly select some for an extra row
sample_ids = data_df['sample_id'].unique()
extra_sample_ids = np.random.choice(sample_ids, extra_rows, replace=False)

# Define a function to sample the required number of rows for each group
def sample_group(group):
    n_to_sample = base_rows + (1 if group.name in extra_sample_ids else 0)
    if n_to_sample > group.shape[0]:
        n_to_sample = group.shape[0]
    return group.sample(n=n_to_sample, random_state=12)

# Apply sampling by group
data_df_100k = data_df.groupby('sample_id', group_keys=False).apply(sample_group)
n_x = n_rows - data_df_100k.shape[0]
if n_x > 0:
    data_df_rm_100k = data_df[~data_df.index.isin(data_df_100k.index)]
    data_df_100k = pd.concat([data_df_100k, data_df_rm_100k.sample(n=n_x, random_state=12)])
print(data_df_100k.shape)  # Should be (100000, ...)

# data_df_100k = data_df.sample(n=100000, random_state=1)
data_df_100k.to_csv(f'{project}/umap_embeddings_with_meta_100k.csv', index_label="cs_id")

meta_100k = data_df_100k.loc[:,meta_list]
meta_100k.to_csv(f"{project}/metadata_lite_100k.csv",index_label="cs_id")

## for each sample, save each column data into a separate json file
os.makedirs(project+"/metas_100k", exist_ok=True)
meta_100k_by_sample = meta_100k.groupby('sample_id')
for sample_id, df in meta_100k_by_sample:
    os.makedirs(project + f"/metas_100k/{sample_id}", exist_ok=True)
    for col in df.columns:
        if col == "sample_id":
            continue
        with open(project + f"/metas_100k/{sample_id}/{col}.json", "w") as f:
            json.dump(df[col].to_dict(), f, indent=2)


# ===========================================================
# save the 100k data
embeddings_data_100k = embeddings_data.loc[data_df_100k.index]
embeddings_data_100k.to_csv(f"{project}/umap_embeddings_100k.csv", index_label="cs_id")

## add sample_id to umap_embedding data
embeddings_data["sample_id"] = embeddings_data.index.map(cell_to_sample)
embeddings_data.to_csv(f"{project}/umap_embeddings_with_sample_id.csv", index_label="cs_id")

embeddings_data_100k = embeddings_data.loc[data_df_100k.index]
embeddings_data_100k.to_csv(f"{project}/umap_embeddings_with_sample_id_100k.csv", index_label="cs_id")

stop
# %% ============================================================================
print("Loading expression data...")
## rename_expression_data
expression_data = pd.read_csv(project + "/raw_normalized_expression_sparse.csv",index_col=None, header=0)
## rename "Cell" column use barcode_cid map
print("Renaming....")
expression_data["cs_id"] = expression_data["Cell"].map(barcode_to_cid)
expression_data.drop("Cell", axis=1, inplace=True)

## "Expression" column keep 4 digits after the decimal point
expression_data["Expression"] = expression_data["Expression"].apply(lambda x: round(x, 2))

## save the new expression data
expression_data.to_csv(project + "/normalized_expression_sparse.csv", index=False)


## add sample_id to expression data
expression_data["sample_id"] = expression_data["cs_id"].map(cell_to_sample)

## save expression data
expression_data.to_csv(f"{project}/normalized_expression_sparse_with_sample_id.csv", index=False)

