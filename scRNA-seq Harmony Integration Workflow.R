library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(cowplot)
library(clustree)
library(kableExtra)
library(readr)
library(here) # 如果用here包管理路径

rm(list = ls())
################
## 1. 加载数据
################

pbmc <- readRDS("/home/lin/CRC/CRC_scRNA.rds")
levels(pbmc@meta.data$orig.ident)
pbmc@meta.data$orig.ident <- pbmc@meta.data$samples
cat("Loaded data with", ncol(pbmc), "cells and", length(unique(pbmc$orig.ident)), "samples.")

################
## 2. NormalizeData标准化 & Harmony批次校正
################

pc.num <- 1:30

pbmc <- pbmc %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.3)

cat("Harmony integration completed. Number of clusters:", length(unique(pbmc$seurat_clusters)))

################
## 3. 可视化批次校正效果
################

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", pt.size=0.4) + ggtitle("By Sample")
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", pt.size=0.4, label=T, repel=T) + ggtitle("By Cluster")

plot_grid(p1, p2, ncol = 2)

################
## 4. 不同分辨率聚类树可视化
################

res_values <- seq(0.1, 1.2, by=0.1)
for(res in res_values){
  pbmc <- FindClusters(pbmc, resolution = res)
}
clustree(pbmc)
pbmc <- FindClusters(pbmc, resolution = 0.7)
################
## 5. Marker基因分析并导出
################

## seurat V5版本运行前JoinLayers
pbmc[["RNA"]] <- JoinLayers(pbmc[["RNA"]])

markers <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25)

top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

top_markers %>%
  kbl() %>%
  kable_styling()

write_csv(markers, here("/home/lin/CRC/combined_all_markers_res0.7.csv"))

################
## 6. 保存结果
################

saveRDS(pbmc, here("/home/lin/CRC/CRC_scRNA_harmony.rds"))
cat("Harmony-integrated object saved.")

################
## 7. 关键基因表达可视化
################

VlnPlot(pbmc, features = "ZG16", layer = "counts", log = TRUE) +
  scale_fill_viridis_d() +
  theme_minimal()

################
## 8. UMAP图导出
################

ggsave(here("umap_harmony_res0.7.tiff"), 
       DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", pt.size=0.5, label=TRUE, repel=TRUE), 
       width=12, height=10, dpi=600, compression="lzw")
