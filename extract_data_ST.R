library(reticulate)

library(jsonlite)
library(Matrix)
library(data.table)

print("load RDS data...")
## Read the rds onject
seurat_obj <- readRDS("data_Jie.rds") 
capture.output(str(seurat_obj), file = "seurat_structure_Jie.txt")


print("Save metadata...")
# Save to CSV with index as the first column
metadata = seurat_obj@meta.data
write.csv(metadata, "ST/raw_metadata.csv", row.names = FALSE)

# save the umap embedding
x = seurat_obj@reductions
print("Save umap embedding...")
umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
write.csv(umap_embeddings, "ST/raw_umap_embeddings.csv", row.names = TRUE)


library(Seurat)
library(EBImage)
np <- import("numpy")
images = seurat_obj@images
all_names = names(images)

## Create the directory
dir.create("ST/images", showWarnings = FALSE)
dir.create("ST/coordinates", showWarnings = FALSE)

#Loop and extract the image data
for (i in 1:length(all_names)) {
    print(all_names[i])

    image_name = all_names[i]
    spatial_data = images[[image_name]]
    image_array = spatial_data@image
    coordinates = spatial_data@coordinates

    # Extract the scale.factors from the Seurat object
    scale_factors <- spatial_data@scale.factors

    # Remove the custom class attributes
    scale_factors_plain <- unclass(scale_factors)
    scale_factors_plain$spot.radius = spatial_data@spot.radius

    # Convert to JSON with pretty formatting
    json_output <- toJSON(scale_factors_plain, pretty = TRUE, auto_unbox = TRUE)

    # Save to a JSON file
    write(json_output, paste0("ST/coordinates/raw_scalefactors_", image_name, ".json"))

    # Convert to EBImage format
    image_eb <- Image(image_array, colormode = "Color")

    # Convert to NumPy array and save
    
    np$save(paste0("ST/images/raw_array_", image_name, ".npy"), image_array)

    # Save as PNG (best for analysis)
    writeImage(image_eb, paste0("ST/images/raw_image_", image_name, ".png"), type = "png")
    writeImage(image_eb, paste0("ST/images/raw_image_", image_name, ".tiff"), type = "tiff")

    # Save coordinates

    write.csv(coordinates, paste0("ST/coordinates/raw_coordinates_", image_name, ".csv"), row.names = TRUE)
}

exit(0)
# ===================================================
# Extract normalized counts (log-normalized values)
print("Saving normalized data...")
# raw_counts <- seurat_obj@assays$Spatial@counts

# Extract normalized expression data
normalized_counts <- seurat_obj@assays$Spatial@data  # This is a sparse matrix


# Convert sparse matrix to triplet format (long format)
long_data <- summary(normalized_counts)

# Get row (gene) and column (spot) names
long_data$Gene <- rownames(normalized_counts)[long_data$i]
long_data$Spot <- colnames(normalized_counts)[long_data$j]
long_data$Expression <- long_data$x

# Keep only necessary columns
long_data <- long_data[, c("Gene", "Spot", "Expression")]

# Filter out zero values
nonzero_data <- long_data[long_data$Expression != 0, ]

# Convert to data.table for efficiency
nonzero_data <- as.data.table(nonzero_data)

# Save to CSV
fwrite(nonzero_data, "ST/raw_normalized_expression_sparse.csv", row.names = FALSE)
