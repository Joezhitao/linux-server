library(Seurat)
library(readxl)
setwd("/home/lin/c_group/")
pbmc <- readRDS("/home/lin/c_group/optimization_results/umap/filtered_all_clusters_threshold0.2.rds")

#手工标注
ident <- read_excel("ident.xlsx", sheet = "Sheet1")
new.cluster.ids <- c(ident$ident)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)

#删除分簇
pbmc <- pbmc[, pbmc@meta.data$seurat_clusters != "删除"]
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

saveRDS(pbmc,"/home/lin/c_group/ident.rds")
