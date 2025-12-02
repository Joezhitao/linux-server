library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(clustree)
library(kableExtra)
library(readr)

rm(list = ls())

pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")

pc.num <- 1:20

pbmc <- pbmc %>%
  ScaleData(vars.to.regress = c("mt_percent")) %>%
  FindVariableFeatures(nfeatures = 2000) %>% 
  RunPCA(npcs = 50, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.3) %>% 
  RunTSNE(reduction = "pca", dims = pc.num)
#微调
pbmc <- FindClusters(pbmc, resolution = 0.3)

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE,pt.size=0.4) + ggtitle("By Sample")

saveRDS(pbmc, "/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")
