library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(Seurat)

rm(list=ls())

setwd("/home/lin/c_group/cellchat/")

HSCs <- readRDS("/home/lin/c_group/subpopulation_analysis_filtered_complete/HSCs/HSCs_processed_filtered_complete.rds")
CD4T <- readRDS("/home/lin/c_group/CD4+T_ident.rds")
Hepatocytes <- readRDS("/home/lin/c_group/subpopulation_analysis/Hepatocytes/Hepatocytes_processed.rds")
CD8T <- readRDS("/home/lin/c_group/CD8+T_ident.rds")

# Check the features in each dataset
HSCs_features <- rownames(HSCs)
CD4T_features <- rownames(CD4T)
Hepatocytes_features <- rownames(Hepatocytes)
CD8T_features <- rownames(CD8T)
# Find common features across all datasets
common_features <- Reduce(intersect, list(HSCs_features, CD4T_features, Hepatocytes_features, CD8T_features))
cat("Number of common features:", length(common_features), "\n")

# Add cell type information to metadata before merging
HSCs$cell_type <- "HSC"
CD4T$cell_type <- "CD4+T"
Hepatocytes$cell_type <- "Hepatocyte"
CD8T$cell_type <- "CD8+T"

# Optionally, add dataset information to cell IDs to avoid duplicate cell names
HSCs <- RenameCells(HSCs, add.cell.id = "HSC")
CD4T <- RenameCells(CD4T, add.cell.id = "CD4+T")
Hepatocytes <- RenameCells(Hepatocytes, add.cell.id = "Hep")
CD8T <- RenameCells(CD8T, add.cell.id = "CD8+T")
# Merge the datasets
merged_seurat <- merge(HSCs, y = c(CD4T, Hepatocytes, CD8T), 
                       add.cell.ids = NULL,  # Already added cell IDs
                       project = "Liver_Cells")
# Save the merged dataset
saveRDS(merged_seurat, file = "/home/lin/c_group/merged_liver_cells.rds")
