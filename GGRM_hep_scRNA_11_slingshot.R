# 加载必要的包
library(Seurat)
library(slingshot)
library(RColorBrewer)
library(tradeSeq)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grDevices)

rm(list = ls())

#----------------------------------------------------------------
# 1. 数据准备
#----------------------------------------------------------------

# 读取预处理好的Seurat对象
scRNAsub <- readRDS("/home/lin/GGRM_hep/result/GGRM_Ep_11m27d.rds")
Idents(scRNAsub) <- scRNAsub$seurat_clusters
levels(scRNAsub)
# 设置输出路径
output_dir <- "/home/lin/GGRM_hep/result/"
dir.create(output_dir, showWarnings = FALSE)

#----------------------------------------------------------------
# 2. 轨迹分析
#----------------------------------------------------------------

# 转为SingleCellExperiment对象
scRNAsub_sce <- as.SingleCellExperiment(scRNAsub, assay = "RNA") 

# 运行Slingshot轨迹分析
sce_slingshot <- slingshot(scRNAsub_sce,     
                           reducedDim = 'umap',  
                           start.clus = 'S2',
                           clusterLabels = scRNAsub_sce$seurat_clusters)

# 查看Slingshot结果
print(SlingshotDataSet(sce_slingshot))

#----------------------------------------------------------------
# 3. 轨迹图 (与教程相同的样式)
#----------------------------------------------------------------
output_file1 <- paste0(output_dir, "slingshot_trajectory.png")
png(output_file1, width = 2000, height = 2000, res = 300)

# 设置伪时间颜色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_slingshot$slingPseudotime_1, breaks=100)]

# 绘制基础图
plot(reducedDims(sce_slingshot)$UMAP, col = plotcol, pch=16, asp = 1,
     main = "Trajectory Analysis with Pseudotime",
     xlab = "UMAP 1", ylab = "UMAP 2")

# 绘制轨迹线
lines(SlingshotDataSet(sce_slingshot), lwd=2, col='black')

dev.off()

#----------------------------------------------------------------
# 4. 基因动态表达分析
#----------------------------------------------------------------

# 准备伪时间数据
pseudotimeED <- slingPseudotime(sce_slingshot, na = FALSE)
cellWeightsED <- slingCurveWeights(sce_slingshot)

# 获取计数矩阵
counts <- assay(scRNAsub_sce, "counts")

# 找差异基因
diff.genes_results <- FindAllMarkers(scRNAsub, only.pos = TRUE, 
                                     min.pct = 0.1, logfc.threshold = 0.25)

# 选择显著的差异基因
sig_genes <- diff.genes_results %>% 
  filter(p_val_adj < 0.05) %>%
  top_n(50, avg_log2FC) %>%
  pull(gene)

# 验证基因存在于数据中
valid_genes <- sig_genes[sig_genes %in% rownames(counts)]
cat("使用", length(valid_genes), "个差异基因进行轨迹分析\n")

# 提取这些基因的计数
counts_subset <- counts[valid_genes, ]

# 拟合GAM模型
sce_gam <- fitGAM(counts = as.matrix(counts_subset), 
                  pseudotime = pseudotimeED, 
                  cellWeights = cellWeightsED, 
                  nknots = 5, 
                  verbose = TRUE)

# 评估模型收敛情况
convergence_rate <- mean(rowData(sce_gam)$tradeSeq$converged, na.rm = TRUE)
cat("模型收敛率:", convergence_rate, "\n")

# 关联测试找出随伪时间变化最显著的基因
ATres <- associationTest(sce_gam)

#----------------------------------------------------------------
# 5. 简单热图 (与教程相同的样式)
#----------------------------------------------------------------

# 选择热图显示的基因
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:50]

# 按伪时间排序细胞
pst.ord <- order(sce_slingshot$slingPseudotime_1, na.last = NA)

# 提取热图数据
heatdata <- assay(scRNAsub_sce, "counts")[topgenes, pst.ord]
heatclus <- scRNAsub_sce$seurat_clusters[pst.ord]
heatdata <- as.matrix(heatdata)

# 绘制热图
output_file2 <- paste0(output_dir, "gene_expression_heatmap.png")
png(output_file2, width = 2500, height = 3000, res = 300)

# 确保heatclus是因子并且有正确的水平
heatclus <- factor(heatclus)
cluster_colors <- brewer.pal(min(9, length(levels(heatclus))), "Set1")

heatmap(log1p(heatdata), 
        Colv = NA,  # 不对列进行聚类
        ColSideColors = cluster_colors[as.numeric(heatclus)],
        labCol = "",  # 不显示列标签
        main = "Gene Expression along Pseudotime")

dev.off()

#----------------------------------------------------------------
# 6. 基因动态曲线图 (与教程相同的样式)
#----------------------------------------------------------------

# 选择要绘制的基因
# 根据您的数据选择适当的基因
genes_plot <- rownames(ATres)[order(ATres$pvalue)][1:2]  # 选择最显著的两个基因

# 创建基因动态图
gene_dynamic <- list()
for(i in 1:length(genes_plot)){
  p = plotSmoothers(sce_gam, assays(sce_slingshot)$counts,
                    gene = genes_plot[i], alpha = 0.6, border = TRUE, lwd = 2) +
    ggtitle(genes_plot[i])
  gene_dynamic[[i]] <- p
}

# 保存基因动态图
output_file3 <- paste0(output_dir, "gene_dynamics.png")
png(output_file3, width = 2500, height = 1200, res = 300)
cowplot::plot_grid(plotlist = gene_dynamic, ncol = 2)
dev.off()

cat("轨迹分析完成，所有图像已保存到:", output_dir, "\n")
