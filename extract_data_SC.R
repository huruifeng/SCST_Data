library(jsonlite)
library(Matrix)
library(data.table)
library(Seurat)
library(presto)


print("load RDS data...")
## Read the rds onject
seurat_obj <- readRDS("data_Jacob.rds") 
capture.output(str(seurat_obj), file = "seurat_structure_Jacob.txt")


print("Save metadata...")
# Save to CSV with index as the first column
metadata = seurat_obj@meta.data
write.csv(metadata, "SC/raw_metadata.csv", row.names = TRUE)

# save the umap embedding
print("Save umap embedding...")
umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
write.csv(umap_embeddings, "SC/raw_umap_embeddings.csv", row.names = TRUE)


# Extract normalized counts (log-normalized values)
print("Saving normalized data...")
# raw_counts <- seurat_obj@assays$RNA@counts

# Extract normalized expression data
normalized_counts <- seurat_obj@assays$RNA@data  # This is a sparse matrix

# Convert sparse matrix to triplet format (long format)
long_data <- summary(normalized_counts)

# Get row (gene) and column (cell) names
long_data$Gene <- rownames(normalized_counts)[long_data$i]
long_data$Cell <- colnames(normalized_counts)[long_data$j]
long_data$Expression <- long_data$x

# Keep only necessary columns
long_data <- long_data[, c("Gene", "Cell", "Expression")]

# Filter out zero values
nonzero_data <- long_data[long_data$Expression != 0, ]

# Convert to data.table for efficiency
nonzero_data <- as.data.table(nonzero_data)

# Save to CSV
fwrite(nonzero_data, "SC/raw_normalized_expression_sparse.csv", row.names = FALSE)


# cell type specific markers
print("Saving cell type specific markers...")
# Extract cell type specific markers
cell_type_markers <- FindAllMarkers(seurat_obj, group.by =  "MajorCellTypes")
# Convert to data.table
cell_type_markers_dt <- as.data.table(cell_type_markers)
# Save to CSV
fwrite(cell_type_markers_dt, "SC/celltypes/celltype_FindAllMarkers.csv", row.names = FALSE)
