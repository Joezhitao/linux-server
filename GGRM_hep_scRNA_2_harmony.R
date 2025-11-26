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

setwd("/home/lin/GGRM_hep/result/")
################
## 1. 加载数据
################

pbmc <- readRDS("/home/lin/GGRM_hep/result/GGRM_hep_scRNAseq_11m26d.rds")
pbmc@meta.data$orig.ident <- as.factor(pbmc@meta.data$orig.ident)
levels(pbmc@meta.data$orig.ident)

################
## 2. NormalizeData标准化 & Harmony批次校正
################

pc.num <- 1:20

pbmc <- pbmc %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("mt_percent")) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = pc.num) %>%
  RunTSNE(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.3)


ElbowPlot(pbmc, ndims = 60)

cat("Harmony integration completed. Number of clusters:", length(unique(pbmc$seurat_clusters)))

################
## 3. 可视化批次校正效果
################

p1 <- DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident", pt.size=0.4) + ggtitle("By Sample")
p2 <- DimPlot(pbmc, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.4, label=T, repel=T) + ggtitle("By Cluster")

plot_grid(p1, p2, ncol = 2)

################
## 4. 保存结果
################

saveRDS(pbmc, here("/home/lin/GGRM_hep/result/CRC_scRNA_harmony.rds"))
cat("Harmony-integrated object saved.")
