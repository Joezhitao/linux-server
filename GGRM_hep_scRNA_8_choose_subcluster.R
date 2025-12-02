library(Seurat)

rm(list = ls())

pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_ident_11m27d.rds")

levels(pbmc@meta.data$celltype)

celltype <- c("Epithelial_cells")
pbmc <- pbmc[, pbmc@meta.data$celltype %in% celltype]
pbmc@meta.data$celltype <- droplevels(pbmc@meta.data$celltype)

saveRDS(pbmc, "/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")
