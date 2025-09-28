library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(xlsx)

setwd("/home/lin/B_group")

rm(list = ls())

# ============ 参数设置区域 ============
# 设置基因集参数
geneset_file <- "gene1.xlsx"  # 基因集文件名
geneset_sheet <- 1            # Excel中的sheet索引
geneset_name <- "Proliferation" # 基因集名称 (例如: "Senescence", "Regeneration", etc.)
output_prefix <- "Proliferation_geneset" # 输出文件前缀
# ====================================

# 读取已处理的Seurat对象
pbmc <- readRDS('/home/lin/B_group/Hepatocytes.RDS')

# 读取基因集
geneset <- read.xlsx(paste0("/home/lin/B_group/", geneset_file), 
                     sheetIndex = geneset_sheet, 
                     header = T, 
                     encoding = "UTF-8")
geneset <- geneset$gene
# 选择seurat对象中有的基因
geneset <- geneset[geneset %in% rownames(pbmc)]

# 显示找到的基因数量
cat(paste0("在Seurat对象中找到", length(geneset), "个", geneset_name, "基因\n"))

# 计算基因集的模块评分
pbmc <- AddModuleScore(
  object = pbmc,
  features = list(geneset),
  ctrl.size = 100,
  name = geneset_name
)

# 评分结果列名
score_column <- paste0(geneset_name, "1")

# 定义簇的颜色
my_colors <- c('#5470c6', '#91cc75')

# 创建小提琴图比较不同簇之间的评分
p <- ggviolin(pbmc@meta.data, 
              x = "seurat_clusters", 
              y = score_column,
              color = "seurat_clusters", 
              fill = "seurat_clusters",
              add = "mean_sd", #mean
              add.params = list(color = "black")) +
  stat_compare_means(comparisons = list(c("S1", "S2")), 
                     label = "p.signif") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "", 
       y = paste0(geneset_name, " Geneset Score"), 
       title = paste0(geneset_name, " Geneset Score")) +
  theme(legend.position = "none")

# 保存小提琴图
violin_filename <- paste0(output_prefix, "_score_violin.png")
ggsave(violin_filename, p, width = 6, height = 4.5)
cat(paste0("小提琴图已保存为: ", violin_filename, "\n"))

# 在UMAP上展示基因评分
p2 <- FeaturePlot(pbmc, features = score_column) +
  labs(title = paste0(geneset_name, " Geneset Score"))
umap_filename <- paste0("UMAP_", output_prefix, "_score.png")
ggsave(umap_filename, p2, width = 6, height = 4.5)
cat(paste0("UMAP图已保存为: ", umap_filename, "\n"))

# 创建小提琴图比较不同簇之间的评分 - 去掉对比线和p值
p <- ggviolin(pbmc@meta.data, 
              x = "seurat_clusters", #选择分簇进行分析
              y = score_column,
              color = "seurat_clusters", 
              fill = "seurat_clusters",
              add = "mean_sd", 
              add.params = list(color = "black")) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "", 
       y = paste0(geneset_name, " Geneset Score"), 
       title = paste0(geneset_name, " Geneset Score")) +
  theme(legend.position = "none")

# 保存小提琴图
violin_filename <- paste0(output_prefix, "_score_violin.png")
ggsave(violin_filename, p, width = 6, height = 4.5)
cat(paste0("小提琴图已保存为: ", violin_filename, "\n"))
