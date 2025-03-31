library(ggplot2)
library(jsonlite)
library(data.table)

project <- "ST"
genes <- c("RORB", "LAMP5")

# Read UMAP embeddings
umap_file <- file.path(project, "umap_embeddings_with_meta.csv")
data_df <- fread(umap_file)
setnames(data_df, names(data_df)[1], "Spot")  # Set first column as 'Spot'

# Add gene expression data
for(gene in genes) {
  gene_file <- file.path(project, "gene_jsons", paste0(gene, ".json"))
  gene_data <- fromJSON(gene_file)
  
  # Convert JSON to data.frame
  gene_df <- stack(gene_data)
  colnames(gene_df) <- c(gene, "Spot")
  
  # Merge with main data
  data_df <- merge(data_df, gene_df, by = "Spot", all.x = TRUE)
  
  # Replace NA with 0
  data_df[is.na(get(gene)), (gene) := 0]
}

# Create violin plots
create_violin <- function(gene) {
  ggplot(data_df, aes(x = smoothed_label_s5, y = .data[[gene]], fill=smoothed_label_s5)) +
	geom_violin(trim = FALSE, scale = "width") +
	# geom_boxplot(width=0.1) +
	# geom_jitter(width = 0.1, alpha = 0.5) +
	labs(title = paste(gene, "Expression - 50k"), x = "smoothed_label_s5", y = "Expression Level") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate and save plots
for(gene in genes) {
  ggsave(
	filename = paste0(gene, "_violin_R.pdf"),
	plot = create_violin(gene),
	device = "pdf",
	width = 8,
	height = 6
  )
}

# Create scatter plots
create_scatter_continuous <- function(feature) {
  ggplot(data_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
	geom_point(size = 0.1) +
	scale_color_gradient(low = "gray", high = "red") +  # Color gradient from gray to red
	labs(title = "50k", x = "UMAP_1", y = "UMAP_2") +
	theme_minimal()
}

create_scatter<- function(feature) {
  ggplot(data_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
	geom_point(size = 0.1) +
	labs(title = "50k", x = "UMAP_1", y = "UMAP_1") +
	theme_minimal()
}

# Generate and save plots
for(gene in genes){
  print(gene)
  ggsave(
	filename = paste0(gene, "_scatter_R.pdf"),
	plot = create_scatter_continuous(gene),
	device = "pdf",
	width = 8,
	height = 6
  )
}


ggsave(
	filename = paste0("smoothed_label_s5", "_scatter_R.pdf"),
	plot = create_scatter("smoothed_label_s5"),
	device = "pdf",
	width = 8,
	height = 6
	)
