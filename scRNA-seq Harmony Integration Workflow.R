library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

rm(list = ls())
################
## 1. 加载数据
################

pbmc <- readRDS("/home/lin/CRC/CRC_scRNA.rds")
levels(pbmc@meta.data$orig.ident)
pbmc@meta.data$orig.ident <- pbmc@meta.data$samples
cat("Loaded data with", ncol(pbmc), "cells and", length(unique(pbmc$orig.ident)), "samples.")

################
## 2. SCTransform标准化 & Harmony批次校正
################

pc.num <- 1:30

pbmc <- pbmc %>%
  SCTransform(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.3)

cat("Harmony integration completed. Number of clusters:", length(unique(pbmc$seurat_clusters)))
