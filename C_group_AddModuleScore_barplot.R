library(Seurat)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
rm(list = ls())
# 读取数据
seu_obj <- readRDS("/home/lin/c_group/CD8+T_ident.rds")
gene_list <- read_excel("/home/lin/c_group/gene.xlsx")

# 设置样本顺序和重命名
sample_order <- c("Con", "Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R")
seu_obj$orig.ident <- factor(seu_obj$orig.ident, 
                             levels = c("Ad1R", "Cd1L", "Cd1R", "Cd15L", "Cd15R", "Cd30L", "Cd30R"))
levels(seu_obj$orig.ident) <- sample_order

# 计算耗竭评分
exhaustion_genes <- na.omit(gene_list$gene)
seu_obj <- AddModuleScore(seu_obj,
                          features = list(exhaustion_genes),
                          name = "Exhaustion_Score")

# 准备绘图数据
plot_data <- data.frame(
  Sample = factor(seu_obj$orig.ident, levels = sample_order),  # 确保顺序
  Score = seu_obj$Exhaustion_Score1
)

# 绘制箱线图
p <- ggplot(plot_data, aes(x = Sample, y = Score)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", 
       y = "Exhaustion Score", 
       title = "CD8+ T Cell Exhaustion Score by Sample")

# 显示和保存
print(p)
ggsave("/home/lin/c_group/exhaustion_score_boxplot.pdf", p, width = 8, height = 6)
