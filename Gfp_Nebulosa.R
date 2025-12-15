# 加载必要的包
library(Seurat)
library(schard)
library("Nebulosa")

rm(list = ls())
setwd("/home/lin/gfpdata")

# 读取h5ad文件并转换为Seurat对象
seu_obj <- schard::h5ad2seurat('/home/lin/gfpdata/CD8 T cells_processed.h5ad')
seu_obj@reductions$umap <- seu_obj@reductions$Xumap_
# 绘制CD4和CD8a的联合密度图
m3 <- plot_density(
  seu_obj, 
  features = c("Pdcd1"),  # 两个基因
  joint = FALSE,                  # 显示联合密度（共表达）
  reduction = "umap",          # 使用的降维方法
  size = 0.8                     # 点的大小
)
m3
