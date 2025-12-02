#查看基因表达情况
library(scRNAtoolVis)
library(Seurat)
library(readxl)  # 用于读取Excel文件

rm(list = ls())

# 设置工作目录
setwd("/home/lin/c_group/")
# 读取之前处理好的单细胞数据
pbmc <- readRDS("/home/lin/c_group/ident.rds")

# 读取基因列表
gene_df <- read_excel("/home/lin/c_group/gene.xlsx")
genes <- gene_df$gene

# 计算需要多少组图（每4个基因一组）
n_groups <- ceiling(length(genes) / 4)

for (i in 1:n_groups) {
  # 计算当前组的基因索引
  start_idx <- (i-1)*4 + 1
  end_idx <- min(i*4, length(genes))
  
  # 获取当前组的基因
  current_genes <- genes[start_idx:end_idx]
  
  # 直接绘制图形
  p <- featureCornerAxes(object = pbmc, reduction = 'umap',
                         groupFacet = NULL,
                         features = current_genes,
                         cornerTextSize = 3,
                         themebg = 'bwCorner')
  print(p)  # 显示图形
  
  # 输出信息
  cat("Displayed plot for genes:", paste(current_genes, collapse=", "), "\n")
  
  # 可选：等待用户按Enter键继续查看下一组
  cat("Press Enter to continue to next group...\n")
  readline()
}
featureCornerAxes(object = pbmc, reduction = 'umap',
                  groupFacet = NULL,
                  features = c("Pdcd1"),
                  cornerTextSize = 3,
                  themebg = 'bwCorner')
