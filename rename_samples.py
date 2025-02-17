# %% ==============================
import pandas as pd
import json
from collections import defaultdict

# %%
##===============================
## read metadata
project = "SC"
metadata = pd.read_csv( project + "/raw_metadata.csv", index_col=0, header=0)
cell_ids = ["c_"+str(i) for i in range(1,metadata.shape[0]+1)]
metadata["cell_id"] = cell_ids

barcode_to_cid =  metadata["cell_id"].to_dict()
with open(project+"/barcode_to_cid.txt", "w") as f:
    json.dump(barcode_to_cid, f)

metadata["barcode"] = metadata.index.tolist()
metadata = metadata.set_index("cell_id")

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
metadata_lite = metadata.loc[:,["sample_id","case","sex","age","seurat_clusters","MajorCellTypes","CellSubtypes"]]
metadata_lite.to_csv(project + "/metadata_lite.csv")

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
embeddings_data.to_csv(project + "/umap_embeddings.csv",index_label="Cell")

# sampling data
data_df = pd.merge(embeddings_data, metadata_lite, left_index=True, right_index=True)
data_df.to_csv(f'{project}umap_embeddings_with_meta.csv', index_label="Cell")

# sampling data
data_df_20 = data_df.sample(frac=0.2, random_state=1)
data_df_20.to_csv(f'{project}/umap_embeddings_with_meta_20.csv', index_label="Cell")
data_df_50 = data_df.sample(frac=0.5, random_state=1)
data_df_50.to_csv(f'{project}/umap_embeddings_with_meta_50.csv', index_label="Cell")
data_df_100k = data_df.sample(n=100000, random_state=1)
data_df_100k.to_csv(f'{project}/umap_embeddings_with_meta_100k.csv', index_label="Cell")

## add sample_id to umap_embedding data
embeddings_data["sample_id"] = embeddings_data.index.map(cell_to_sample)
embeddings_data.to_csv(f"{project}/umap_embeddings_with_sample_id.csv", index_label="Cell")


stop
# %% ============================================================================
# %%
print("Loading expression data...")
## rename_expression_data
expression_data = pd.read_csv(project + "/raw_normalized_expression_sparse.csv",index_col=None, header=0)
## rename "Cell" column use barcode_cid map
print("Renaming....")
expression_data["Cell"] = expression_data["Cell"].map(barcode_to_cid)

## "Expression" column keep 4 digits after the decimal point
expression_data["Expression"] = expression_data["Expression"].apply(lambda x: round(x, 4))

## save the new expression data
expression_data.to_csv(project + "/normalized_expression_sparse.csv", index=False)


## add sample_id to expression data
expression_data["sample_id"] = expression_data["Cell"].map(cell_to_sample)

## save expression data
expression_data.to_csv(f"{project}/normalized_expression_sparse_with_sample.csv", index=False)

