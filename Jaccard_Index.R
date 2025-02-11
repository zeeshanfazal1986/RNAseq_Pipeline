library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

setwd("/Users/fazalz2/Desktop/2024-11-22-scRNAseq_Viral/")

seurat_standard<-readRDS("integrated_seurat_LabeledTransfered_standardworkflow_3000.rds")
seurat_sct<-readRDS("integrated_seurat_LabeledTransfered_SCT_3000.rds")

# Extract cell barcodes and their annotations
meta_standard <- seurat_standard@meta.data %>%
  select(mouselung) %>%
  rename(annotation_standard = mouselung)

meta_sct <- seurat_sct@meta.data %>%
  select(mouselung) %>%
  rename(annotation_sct = mouselung)

# Merge metadata based on cell barcodes (rownames should be identical cell barcodes)
merged_meta <- merge(meta_standard, meta_sct, by = "row.names", all = FALSE)
colnames(merged_meta)[1] <- "cell_barcode"

# Identify matching and mismatching annotations
merged_meta <- merged_meta %>%
  mutate(match = ifelse(annotation_standard == annotation_sct, "Match", "Mismatch"))

# Compute Jaccard Index
jaccard_index <- table(merged_meta$annotation_standard, merged_meta$annotation_sct)
jaccard_index <- as.matrix(jaccard_index)
write.table(jaccard_index, file="Jaccard_matrix.txt", quote=FALSE, sep="\t")

# Convert to Jaccard similarity (intersection / union)
jaccard_sim <- outer(
  rownames(jaccard_index), colnames(jaccard_index),
  Vectorize(function(x, y) {
    intersection <- jaccard_index[x, y]
    union <- sum(jaccard_index[x, ]) + sum(jaccard_index[, y]) - intersection
    if (union == 0) return(0) else return(intersection / union)
    
  })
)

# Convert to data frame for plotting
jaccard_df <- melt(jaccard_sim)
colnames(jaccard_df) <- c("Annotation_Standard", "Annotation_SCT", "Jaccard_Index")


# Plot heatmap using ggplot2
plot1<-ggplot(jaccard_df, aes(x = Annotation_SCT, y = Annotation_Standard, fill = Jaccard_Index)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Standard vs SCT 3000 Genes",
       x = "SCT Normalization Annotation",
       y = "Standard Normalization Annotation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("standard_vs_sct_3000.pdf", plot=plot1, height = 10, width = 10, dpi=600)

