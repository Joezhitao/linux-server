rm(list = ls())
# 加载包
library("Seurat")
library("Nebulosa")

# 设置工作目录和读取数据
setwd("/home/lin/c_group/subpopulation_analysis")

pbmc <- readRDS("/home/lin/c_group/hep.rds")

# 绘制CD4和CD8a的联合密度图
m3 <- plot_density(
  pbmc, 
  features = c("Sod1", "Atp5mc1", "Npas2", "Zbtb16", "Nrg1", "Cpt1a"),  # 两个基因
  joint = FALSE,                  # 显示联合密度（共表达）
  reduction = "umap",          # 使用的降维方法
  size = 0.8                     # 点的大小
)
m3
# 保存联合密度图
ggsave(
  plot = m3, 
  filename = "m3.png", 
  width = 16,          # 宽度设置较大以容纳三个子图
  height = 4.5, 
  dpi = 600,
  bg = "white"
)
