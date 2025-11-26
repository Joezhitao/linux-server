library(Seurat)

rm(list = ls())

setwd("/home/lin/GGRM_hep/result/")

pbmc <- readRDS("/home/lin/GGRM_hep/result/cluster_optimization/filtered_both_threshold0.2.rds")

p1 <- DimPlot(pbmc,group.by = "celltype")
p2 <- DimPlot(pbmc,group.by = "orig.ident", reduction = "tsne")
p1+p2

#添加标签

levels(pbmc@meta.data$orig.ident)
pbmc@meta.data$group <- pbmc@meta.data$orig.ident
levels(pbmc@meta.data$group) <- c("resistance","Non-resistant","resistance","Non-resistant")

DimPlot(pbmc, group.by = "group", reduction = "tsne", pt.size = 0.8)
DimPlot(pbmc, group.by = "seurat_clusters", reduction = "tsne", pt.size = 0.8, label = TRUE)
DimPlot(pbmc, group.by = "celltype", reduction = "tsne", pt.size = 0.8)

saveRDS(pbmc, "/home/lin/GGRM_hep/result/GGRM_hep_ident_11m27d.rds")
