# %%
import pandas as pd
import json
import os
import matplotlib.pyplot as plt
import seaborn as sns

# %%

project="SC"
genes=["RORB","LAMP5"]

umap_embeddings_file = os.path.join(project, 'umap_embeddings_with_meta.csv')
data_df = pd.read_csv(umap_embeddings_file, index_col=0, header=0)
## Cell,UMAP_1,UMAP_2,sample_id,case,sex,age,seurat_clusters,MajorCellTypes,CellSubtypes

# %%
for gene in genes:
    gene_expr_file = os.path.join(project, "gene_jsons",gene+".json")
    with open(gene_expr_file, 'r') as f:
        cell_expr = json.load(f)
        data_df[gene] = data_df.index.map(cell_expr).fillna(0)
data_df = data_df.loc[:, ["UMAP_1", "UMAP_2","MajorCellTypes"]+genes]

# %%
sns.violinplot(data=data_df, x="MajorCellTypes", y="RORB", palette="Set1")
plt.xticks(rotation=90)
plt.title("RORB - All (500k)")
plt.savefig("RORB_violin.png")

sns.violinplot(data=data_df, x="MajorCellTypes", y="LAMP5", palette="Set1")
plt.xticks(rotation=90)
plt.title("LAMP5 - All (500k)")
plt.savefig("LAMP5_violin.png")


# %%
# Create scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=data_df, x="UMAP_1", y="UMAP_2", hue="MajorCellTypes", palette="Set1", s=3)
plt.title("500k")
plt.legend(title="MajorCellTypes", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("MajorCellTypes_Scatter.png")

# Create scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=data_df, x="UMAP_1", y="UMAP_2", hue="RORB", s=3)
plt.title("500k")
plt.legend(title="RORB", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("RORB_Scatter.png")

# Create scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=data_df, x="UMAP_1", y="UMAP_2", hue="LAMP5", s=3)
plt.title("500k")
plt.legend(title="LAMP5", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("LAMP5_Scatter.png")
# %%
