# ======================================================
# C_group_拟时许分析(monocle2)
# 作者: Joe
# 日期: 2025-11-23
# ======================================================

# 清空环境
rm(list = ls())

# 设置工作目录
setwd("/home/lin/c_group/monocle/")

# 加载必要的包
library(monocle)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(ggplot2)
library(schard)
library(showtext)
library(igraph)
library(readxl)
# 加载数据
seu_obj <- readRDS("/home/lin/c_group/CD8+T_ident.rds")

# 转换Seurat对象为CellDataSet
sample_ann <- seu_obj@meta.data
expr_matrix <- as(as.matrix(GetAssayData(seu_obj, layer = "counts")), "sparseMatrix")

# 创建特征注释
feature_ann <- data.frame(
  gene_id = rownames(expr_matrix),
  gene_short_name = rownames(expr_matrix)
)
rownames(feature_ann) <- rownames(expr_matrix)

# 创建AnnotatedDataFrame对象
fd <- new("AnnotatedDataFrame", data = feature_ann)
pd <- new("AnnotatedDataFrame", data = sample_ann)

# 创建CellDataSet
cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

# 预处理CellDataSet
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, cores = 4)
cds <- detectGenes(cds, min_expr = 0.1)

# 选择用于排序的基因
# 只使用在至少5%细胞中表达的基因
expressed_genes <- row.names(subset(fData(cds), 
                                    num_cells_expressed >= 0.05 * ncol(cds)))

# 设置排序基因
cds <- setOrderingFilter(cds, expressed_genes)

# 降维和排序
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

# 绘制轨迹图
p1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
p3 <- plot_cell_trajectory(cds, color_by = "orig.ident")
p4 <- plot_cell_trajectory(cds, markers = c("Foxp3"), use_color_gradient = TRUE)
plot_cell_trajectory(cds, color_by = "TREM2_EFFEROCYTOSIS_score1")
# 组合图像
combined_plot <- p1 | p2
dev.off()
# 保存图像
#pdf("Bcells_pseudotime.pdf", width = 30, height = 6)
print(combined_plot)
#dev.off()

png("Myeloid_cells_pseudotime.png", width = 18, height = 6, res = 300, units = 'in')
print(combined_plot)
dev.off()

# 特定基因的表达图
gene_of_interest <- read_excel("/home/lin/c_group/gene.xlsx", sheet = "Sheet1")
gene_of_interest <- gene_of_interest$gene
gene_of_interest <- gene_of_interest[gene_of_interest %in% rownames(seu_obj)]
for (i in gene_of_interest) {
  gene_plot <- plot_cell_trajectory(cds, markers = i, use_color_gradient = TRUE)
  
  png(paste0("CD8T_", i, "_pseudotime.png"), 
      width = 6, height = 6, units = 'in', res = 300)
  print(gene_plot)
  dev.off()
}
cat("拟时序分析完成!\n")
