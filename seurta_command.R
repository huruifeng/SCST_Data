# Load required libraries
library(Seurat)
library(ggplot2)

# 1. Load the RDS file
seurat_object <- readRDS("data_Jie.rds")

# 2. Verify the structure of the object
# Check available samples (assuming samples are stored in a list)
print(names(seurat_object@images))

# 3. Select a specific sample
selected_sample <- "BN2003"  # Replace with your actual sample name
subset_obj <- subset(seurat_object, subset = orig.ident == selected_sample)

# 4. Plot using SpatialFeaturePlot
SpatialFeaturePlot(
  subset_obj,
  features = "nCount_Spatial",      # Replace with your target gene
  alpha = c(0.8, 1),          # Adjust transparency
  image.alpha = 0.8,            # Image transparency
  crop = F,                  # Crop the image
) + 
theme(legend.position = "right")  # Customize legend
