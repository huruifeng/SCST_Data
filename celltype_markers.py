# %% ==============================
import pandas as pd
import json
from collections import defaultdict
import os
import numpy as np


# %% ============================================================================
## top 10 marker genes in each cell type
project = "SC"
marker_gene_file = project + "/celltypes/celltype_FindAllMarkers.csv"
marker_genes = pd.read_csv(marker_gene_file, index_col=None, header=0)

# Filter for significant genes
filtered_df = marker_genes[marker_genes['p_val_adj'] < 0.05]

# Rank by absolute log2FC and select top 10 per celltype
top_genes = (
    filtered_df
    .assign(abs_log2FC=filtered_df['avg_log2FC'])
    .sort_values(['cluster', 'abs_log2FC'], ascending=[True, False])
    .groupby('cluster')
    .head(10)
    .drop(columns='abs_log2FC')  # Optional: remove helper column
)

# Save or inspect the result
top_genes.loc[:,["cluster","gene","avg_log2FC","p_val_adj"]].to_csv(project+'/top10_genes_per_celltype.csv', index=False)
print(top_genes.head())

## save marker genes of each cell type into a dictionary
marker_genes_dict = defaultdict(list)
for index, row in top_genes.iterrows():
    marker_genes_dict[row["cluster"]].append(row["gene"])
print("marker genes dictionary:")
print(marker_genes_dict)
# Save the dictionary to a JSON file
with open(project + "/celltypes/celltype_markergenes.json", "w") as f:
    json.dump(marker_genes_dict, f, indent=4)

# %%============================================================================
pool_genes = top_genes["gene"].unique().tolist()
celltype_list = top_genes["cluster"].unique().tolist()
metadata = pd.read_csv(project + "/metadata_lite.csv", index_col=0, header=0)

#%% ============================================================================
## calculate the percentage of cells where the gene is detected in that cell type

pct_detected = {}

marker_genes_df = pd.DataFrame()
for celltype in celltype_list:
    print("==========================")
    print("Processing cell type: ", celltype)
    pct_detected[celltype] = {}
    cells_in_celltype = metadata[metadata["MajorCellTypes"] == celltype]
    num_cells = len(cells_in_celltype)
    pct_detected[celltype]["tatal_cells"] = num_cells

    diagnosis_sex_grouped = cells_in_celltype.groupby(["case", "sex"])
    for key, df in diagnosis_sex_grouped:
        print("_".join(key))
        key_str = "_".join(key) + "_cellcounts"
        pct_detected[celltype][key_str] = len(df)
        print(f"Number of cells in {key_str}: {len(df)}")

    
    marker_genes = {}
    for gene in pool_genes:
        marker_genes[gene] = {}
        gene_expr = json.load(open(project + "/gene_jsons/"+gene+".json"))
        gene_expr_in_cells =  list(gene_expr.keys())

        cells_with_gene_expr = [cell for cell in gene_expr_in_cells if cell in cells_in_celltype.index]
        num_cells_with_gene_expr = len(cells_with_gene_expr)

        if num_cells_with_gene_expr == 0:
            avg_expr = 0.00
        else:
            # Calculate the average expression value for the cells with gene expression and convert to float
            # avg_expr = np.mean([float(gene_expr[cell]) for cell in cells_with_gene_expr])

            # Use np.nanmean to ignore NaN values in the calculation
            # This will return NaN if all values are NaN, so we need to handle that case
            avg_expr = np.nanmean([float(gene_expr[cell]) for cell in cells_with_gene_expr])
            # If all values are NaN, set avg_expr to 0.00
            if np.isnan(avg_expr):
                avg_expr = 0.00
            else:
                avg_expr = round(avg_expr, 2)
        marker_genes[gene]["avg_expr"] = avg_expr
        marker_genes[gene]["is_marker"] = gene in marker_genes_dict[celltype]
        marker_genes[gene]["n_expr_cells"] = num_cells_with_gene_expr


        for key, df in diagnosis_sex_grouped:
            n_cells = len([cell for cell in gene_expr_in_cells if cell in df.index])
            key_str = "n_expr_cells_" + "_".join(key)
            marker_genes[gene][key_str] = n_cells
    
    marker_df = pd.DataFrame.from_dict(marker_genes, orient='index')
    marker_df["MajorCellTypes"] = celltype

    marker_genes_df = pd.concat([marker_genes_df, marker_df], axis=0)
    
marker_genes_df.to_csv(project + "/celltypes/celltype_markergenes_cellcounts.csv", index=False)

# Save the dictionary to a JSON file
with open(project + "/celltypes/celltype_cellcounts.json", "w") as f:
    json.dump(pct_detected, f, indent=4)    


# %%

for key, df in diagnosis_sex_grouped:
    print("_".join(key),df)
    # Do something with each group
    # For example, you can calculate the mean of a specific column
    # mean_value = df['your_column_name'].mean()
    # print(f"Mean value for {name}: {mean_value}")
# %%
