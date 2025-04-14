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
de_results_topN_list <- list()

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
			de_results <- FindMarkers(subset_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = "case")
			## add a column for the gene names
			de_results$gene <- rownames(de_results)
			
			## filter out genes with padj > 0.05
			# de_results <- de_results[de_results$p_val_adj < 0.05, ]
			
			## get top 10 upregulated DE genes and downregulated, base on logFC
			de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
			# Store the results in the list
			de_results_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results
			de_results_topN_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results_topN
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = "case")
		## add a column for the gene names
		de_results$gene <- rownames(de_results)

		## filter out genes with padj > 0.05
		# de_results <- de_results[de_results$p_val_adj < 0.05, ]
		
		## get top 10 upregulated DE genes and downregulated, base on logFC
		de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
		de_results_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results
		de_results_topN_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results_topN

	}
}
# Convert the list to a data.table
de_results_dt <- rbindlist(de_results_list, idcol = "CellType_DE")
de_results_topN_dt <- rbindlist(de_results_topN_list, idcol = "CellType_DE")
# Save to CSV
fwrite(de_results_dt, "SC/celltypes/celltype_DEGs.csv", row.names = FALSE)
fwrite(de_results_topN_dt, "SC/celltypes/celltype_DEGs_top10.csv", row.names = FALSE)


## ============================================================
# pseudo-bulk DE analysis in each cell type
print("Calculating pseudo-bulk analysis...")
# Create a new Seurat object for pseudo-bulk analysis
pb_obj <- PseudobulkExpression(seurat_obj, assays = "RNA", method= "aggregate", return.seurat = T, group.by = c("sample_id", "MajorCellTypes", "case"))
tail(Cells(pb_obj))
capture.output(str(pb_obj), file = "seurat_structure_pb_Jacob.txt")
## replace "-" wuth "_" in MajorCellTypes
pb_obj$MajorCellTypes <- gsub("-", "_", pb_obj$MajorCellTypes)
pb_obj$orig.ident <- gsub("-", "_", pb_obj$orig.ident)

metadata = pb_obj@meta.data
## rename sample_id to sampleId
colnames(metadata)[colnames(metadata) == "sample_id"] <- "sampleId"
colnames(metadata)[colnames(metadata) == "case"] <- "condition"

write.csv(metadata, "SC/celltypes/metadata_sample_celltype_condition.csv", row.names = FALSE)

expr_matrix <- GetAssayData(pb_obj, assay = "RNA", slot = "data")
colnames(expr_matrix) <- gsub("-", "_", colnames(expr_matrix))
write.csv(expr_matrix, "SC/celltypes/pb_expr_matrix.csv", row.names = TRUE)


# Define the cell types
cell_types <- unique(pb_obj$MajorCellTypes)
# Initialize an empty list to store results
pseudo_bulk_list <- list()
pseudo_bulk_topN_list <- list()
# Loop through each cell type
for (cell_type in cell_types) {
	# Print the current cell type # nolint: whitespace_linter, indentation_linter.
	print(paste("Processing cell type:", cell_type))
	# Subset the Seurat object to the current cell type # nolint
	subset_obj <- subset(pb_obj, subset = MajorCellTypes == cell_type)

	# subset_obj$celltype.case <- paste(subset_obj$MajorCellTypes, subset_obj$case, sep = "_")
	# Idents(subset_obj) <- "celltype.case"
	# mono.de <- FindMarkers(subset_obj, ident.1 = "GLU_Neurons_PD", ident.2 = "GLU_Neurons_HC", verbose = FALSE)
	# fwrite(mono.de, "GLU_Neurons_PD.HC_pseudobulk_DEGs.csv", row.names = FALSE)

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
			de_results <- FindMarkers(subset_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = "case", use.method = "DEseq2")
			## add a column for the gene names
			de_results$gene <- rownames(de_results)
			
			## filter out genes with padj > 0.05
			# de_results <- de_results[de_results$p_val_adj < 0.1, ]

			## get top 10 upregulated DE genes and downregulated, base on logFC
			de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
			
			# Store the results in the list
			pseudo_bulk_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results
			pseudo_bulk_topN_list[[paste(cell_type, paste(ident.1, ident.2, sep = "vs"), sep = ".")]] <- de_results_topN
		}
	} else {
		de_results <- FindMarkers(subset_obj, ident.1 = condition_ls[1], ident.2 = condition_ls[2], group.by = "case",use.method = "DEseq2")
		## add a column for the gene names
		de_results$gene <- rownames(de_results)
		
		## filter out genes with padj > 0.05
		# de_results <- de_results[de_results$p_val_adj < 0.05, ]
		
		## get top 10 upregulated DE genes and downregulated, base on logFC
		de_results_topN <- rbind(de_results[order(de_results$avg_log2FC, decreasing = TRUE), ][1:10, ],de_results[order(de_results$avg_log2FC, decreasing = FALSE), ][1:10, ])
		
		# Store the results in the list
		pseudo_bulk_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results
		pseudo_bulk_topN_list[[paste(cell_type, paste(condition_ls[1], condition_ls[2], sep = "vs"), sep = ".")]] <- de_results_topN
	}
}
# Convert the list to a data.table
pseudo_bulk_dt <- rbindlist(pseudo_bulk_list, idcol = "CellType_DE")
pseudo_bulk_topN_dt <- rbindlist(pseudo_bulk_topN_list, idcol = "CellType_DE")
# Save to CSV
fwrite(pseudo_bulk_dt, "SC/celltypes/celltype_pseudobulk_DEGs.csv", row.names = FALSE)
fwrite(pseudo_bulk_topN_dt, "SC/celltypes/celltype_pseudobulk_DEGs_top10.csv", row.names = FALSE)

pooled_topN_DEGs = pseudo_bulk_topN_dt$gene
## remove duplicates
pooled_topN_DEGs <- unique(pooled_topN_DEGs)
## subset expr_matrix to only include pooled_topN_DEGs
expr_matrix_pooled_topN_DGEs <- expr_matrix[pooled_topN_DEGs, ]
## save the pooled_topN_DEGs expression matrix
write.csv(expr_matrix_pooled_topN_DGEs, "SC/celltypes/pb_expr_matrix_topN_DEGs.csv", row.names = TRUE)




# GLU_Neurons.PD.HC,0,2.28155130213196,0.347,0.128,0,AP001977.1
# GLU_Neurons.PD.HC,0,1.09698217771197,0.346,0.199,0,SNX31
# GLU_Neurons.PD.HC,0,3.11789946773066,0.236,0.1,0,AC105402.3
subset_obj <- subset(seurat_obj, subset = MajorCellTypes == "GLU_Neurons")
VlnPlot(subset_obj, features = c("AP001977.1", "SNX31", "AC105402.3"), idents = c("PD", "HC"), group.by = "case") 

pb_subset_obj <- subset(pb_obj, subset = MajorCellTypes == "GLU_Neurons")
VlnPlot(pb_subset_obj, features = c("AP001977.1", "SNX31", "AC105402.3"), idents = c("PD", "HC"), group.by = "case") 