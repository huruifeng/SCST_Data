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


## ===========================================================
# cell type specific markers
print("Saving cell type specific markers...")
# Extract cell type specific markers
cell_type_markers <- FindAllMarkers(seurat_obj, group.by =  "MajorCellTypes")
# Convert to data.table
cell_type_markers_dt <- as.data.table(cell_type_markers)
# Save to CSV
fwrite(cell_type_markers_dt, "SC/celltypes/celltype_FindAllMarkers.csv", row.names = FALSE)


## ============================================================
# Ecalcaulate differential expression within each cell type between the conditions # nolint
print("Calculating differential expression...")
# Define the cell types
cell_types <- unique(seurat_obj$MajorCellTypes)
# Initialize an empty list to store results
de_results_list <- list()
# Loop through each cell type
for (cell_type in cell_types) {
	# Print the current cell type # nolint: whitespace_linter, indentation_linter.
	print(paste("Processing cell type:", cell_type))
	# Subset the Seurat object to the current cell type # nolint
	subset_obj <- subset(seurat_obj, subset = MajorCellTypes == cell_type)
	# Calculate differential expression between conditions
	# Here, we assume that the conditions are stored in the "Condition" metadata column 
	# You may need to adjust this based on your actual metadata structure
	condition_ls <- unique(subset_obj$case)
	# if there are more than 2 conditions, compare all possible combinations
	if (length(condition_ls) > 2) {
		combinations <- combn(condition_ls, 2)
		for (i in 1:ncol(combinations)) {
			ident.1 <- combinations[1, i]
			ident.2 <- combinations[2, i]
			de_results <- FindMarkers(subset_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = "case", test.use = "DESeq2")
			# Store the results in the list
			de_results_list[[paste(cell_type, ident.1, ident.2, sep = "_")]] <- de_results
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = "case", test.use = "DESeq2")
		de_results_list[[paste(cell_type, condition_ls[1], condition_ls[2], sep = "_")]] <- de_results
	}
}
# Convert the list to a data.table
de_results_dt <- rbindlist(de_results_list, idcol = "CellType")
# Save to CSV
fwrite(de_results_dt, "SC/celltypes/celltype_DEGs.csv", row.names = FALSE)


## ============================================================
# pseudo-bulk DE analysis in each cell type
print("Calculating pseudo-bulk analysis...")
# Create a new Seurat object for pseudo-bulk analysis
pb_obj <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("sample_id", "case", "MajorCellTypes"))
tail(Cells(pseudo_ifnb))

# Define the cell types
cell_types <- unique(seurat_obj$MajorCellTypes)
# Initialize an empty list to store results
pseudo_bulk_list <- list()
# Loop through each cell type
for (cell_type in cell_types) {
	# Print the current cell type # nolint: whitespace_linter, indentation_linter.
	print(paste("Processing cell type:", cell_type))
	# Subset the Seurat object to the current cell type # nolint
	subset_obj <- subset(pb_obj, idents = cell_type)
	# Calculate differential expression between conditions
	# Here, we assume that the conditions are stored in the "Condition" metadata column 
	# You may need to adjust this based on your actual metadata structure
	condition_ls <- unique(subset_obj$case)
	# if there are more than 2 conditions, compare all possible combinations
	if (length(condition_ls) > 2) {
		combinations <- combn(condition_ls, 2)
		for (i in 1:ncol(combinations)) {
			ident.1 <- combinations[1, i]
			ident.2 <- combinations[2, i]
			de_results <- FindMarkers(subset_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = "case", test.use = "DESeq2")
			# Store the results in the list
			pseudo_bulk_list[[paste(cell_type, ident.1, ident.2, sep = "_")]] <- de_results
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = "case", test.use = "DESeq2")
		pseudo_bulk_list[[paste(cell_type, condition_ls[1], condition_ls[2], sep = "_")]] <- de_results
	}
}
# Convert the list to a data.table
pseudo_bulk_dt <- rbindlist(pseudo_bulk_list, idcol = "CellType")
# Save to CSV
fwrite(pseudo_bulk_dt, "SC/celltypes/celltype_pseudobulk_DEGs.csv", row.names = FALSE)
