library(ggplot2)
library(jsonlite)
library(data.table)

project <- "SC"
genes <- c("RORB", "LAMP5")

# Read UMAP embeddings
umap_file <- file.path(project, "umap_embeddings_with_meta_100k.csv")
data_df <- fread(umap_file)
setnames(data_df, names(data_df)[1], "Cell")  # Set first column as 'Cell'

# Add gene expression data
for(gene in genes) {
  gene_file <- file.path(project, "gene_jsons", paste0(gene, ".json"))
  gene_data <- fromJSON(gene_file)
  
  # Convert JSON to data.frame
  gene_df <- stack(gene_data)
  colnames(gene_df) <- c(gene, "Cell")
  
  # Merge with main data
  data_df <- merge(data_df, gene_df, by = "Cell", all.x = TRUE)
  
  # Replace NA with 0
  data_df[is.na(get(gene)), (gene) := 0]
}

# Create violin plots
create_violin <- function(gene) {
  ggplot(data_df, aes(x = MajorCellTypes, y = .data[[gene]])) +
    geom_violin(trim = FALSE, scale = "width") +
    labs(title = paste(gene, "Expression"), x = "Major Cell Types", y = "Expression Level") +
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