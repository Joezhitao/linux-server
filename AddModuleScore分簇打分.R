library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(xlsx)

setwd("/home/lin/B_group")

rm(list = ls())

# 读取已处理的Seurat对象
pbmc <- readRDS('/home/lin/B_group/Hepatocytes.RDS')

# 读取衰老基因集
geneset <- read.xlsx("/home/lin/B_group/gene.xlsx", sheetIndex = 1, header = T, encoding = "UTF-8")
geneset <- geneset$Gene
# 选择seurat对象中有的基因
geneset <- geneset[geneset %in% rownames(pbmc)]

# 计算衰老基因集的模块评分
pbmc <- AddModuleScore(
  object = pbmc,
  features = list(geneset),
  ctrl.size = 100,
  name = "Senescence"  # 更改为"Senescence"
)

# 定义簇的颜色
my_colors <- c('#5470c6', '#91cc75')

# 创建小提琴图比较不同簇之间的衰老评分
p <- ggviolin(pbmc@meta.data, 
              x = "seurat_clusters", 
              y = "Senescence1",  # 使用新的评分名称
              color = "seurat_clusters", 
              fill = "seurat_clusters",
              add = "mean_sd", 
              add.params = list(color = "black")) +
  stat_compare_means(comparisons = list(c("S1", "S2")), 
                     label = "p.signif") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "", 
       y = "Senescence Geneset Score", 
       title = "Senescence Geneset Score") +  # 更改标题
  theme(legend.position = "none")

# 保存小提琴图
ggsave("Senescence_geneset_score_violin.png", p, width = 6, height = 4.5)

# 在UMAP上展示衰老基因评分
p2 <- FeaturePlot(pbmc, features = "Senescence1") +  # 使用新的评分名称
  labs(title = "Senescence Geneset Score")  # 更改标题
ggsave("UMAP_Senescence_geneset_score.png", p2, width = 6, height = 4.5)
